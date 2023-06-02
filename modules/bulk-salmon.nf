#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastp{
    container params.FASTP_CONTAINER
    label 'cpus_8'
    label 'mem_8'
    tag "${meta.library_id}-bulk"
    input:
        tuple val(meta), path(read1), path(read2)
    output:
        tuple val(meta), path(trimmed_reads)
    script:
        trimmed_reads = "${meta.library_id}_trimmed"
        fastp_report = "${meta.library_id}_fastp.html"
        """
        mkdir -p ${trimmed_reads}
        fastp --in1 <(gunzip -c ${read1}) --out1 ${trimmed_reads}/${meta.library_id}_R1_trimmed.fastq.gz \
        ${meta.technology == 'paired_end' ? "--in2 <(gunzip -c ${read2}) --out2 ${trimmed_reads}/${meta.library_id}_R2_trimmed.fastq.gz" : ""} \
        --length_required 20 \
        --thread ${task.cpus}
        """
    stub:
        trimmed_reads = "${meta.library_id}_trimmed"
        fastp_report = "${meta.library_id}_fastp.html"
        """
        mkdir -p ${trimmed_reads}
        """
}

process salmon{
    container params.SALMON_CONTAINER
    label 'cpus_12'
    label 'mem_24'
    tag "${meta.library_id}-bulk"
    publishDir "${meta.salmon_publish_dir}", mode: 'copy'
    input:
        tuple val(meta), path(read_dir), path(index)
    output:
        tuple val(meta), path(salmon_dir)
    script:
        salmon_dir = file(meta.salmon_results_dir).name
        // get meta to write as file
        meta_json = Utils.makeJson(meta)
        """
        salmon quant -i ${index} \
        -l A \
        ${meta.technology == 'paired_end' ? "-1": "-r"} ${read_dir}/*_R1_*.fastq.gz \
        ${meta.technology == 'paired_end' ? "-2 ${read_dir}/*_R2_*.fastq.gz" : "" } \
        -o ${salmon_dir} \
        --validateMappings \
        --rangeFactorizationBins 4 \
        --gcBias \
        --seqBias \
        --threads ${task.cpus}

        echo '${meta_json}' > ${salmon_dir}/scpca-meta.json
        """
    stub:
        salmon_dir = file(meta.salmon_results_dir).name
        meta_json = Utils.makeJson(meta)
        """
        mkdir -p ${salmon_dir}
        echo '${meta_json}' > ${salmon_dir}/scpca-meta.json
        """

}

process merge_bulk_quants {
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    publishDir "${params.results_dir}/${meta.project_id}", mode: 'copy'
    input:
        tuple val(meta), path(salmon_directories), path(t2g_bulk)
        path(library_metadata)
    output:
        path(tximport_file), emit: bulk_counts
        path(bulk_metadata_file), emit: bulk_metadata
    script:
        tximport_file = "${meta.project_id}_bulk_quant.tsv"
        bulk_metadata_file = "${meta.project_id}_bulk_metadata.tsv"
        workflow_url = workflow.repository ?: workflow.manifest.homePage
        """
        ls -d ${salmon_directories} > salmon_directories.txt

        merge_counts_tximport.R \
          --project_id ${meta.project_id} \
          --salmon_dirs salmon_directories.txt \
          --output_file ${tximport_file} \
          --tx2gene ${t2g_bulk}

        generate_bulk_metadata.R \
         --project_id ${meta.project_id} \
         --salmon_dirs salmon_directories.txt \
         --library_metadata_file ${library_metadata} \
         --metadata_output ${bulk_metadata_file} \
         --genome_assembly ${meta.ref_assembly} \
         --workflow_url "${workflow_url}" \
         --workflow_version "${workflow.revision}" \
         --workflow_commit "${workflow.commitId}"
        """
    stub:
        tximport_file = "${meta.project_id}_bulk_quant.tsv"
        bulk_metadata_file = "${meta.project_id}_bulk_metadata.tsv"
        """
        touch ${tximport_file}
        touch ${bulk_metadata_file}
        """

}

workflow bulk_quant_rna {
    take: bulk_channel
    // a channel with a map of metadata for each rna library to process
    main:
        bulk_channel = bulk_channel
          // add salmon directory and salmon file location to meta
          .map{it.salmon_publish_dir = "${params.checkpoints_dir}/salmon";
               it.salmon_results_dir = "${it.salmon_publish_dir}/${it.library_id}";
               it}
          // split based on whether repeat_mapping is false and the salmon quant.sf file exists
          // and whether the assembly matches the current assembly
          .branch{
              has_quants: (!params.repeat_mapping
                           && file(it.salmon_results_dir).exists()
                           && Utils.getMetaVal(file("${it.salmon_results_dir}/scpca-meta.json"), "ref_assembly") == "${it.ref_assembly}"
                           && Utils.getMetaVal(file("${it.salmon_results_dir}/scpca-meta.json"), "t2g_bulk_path") == "${it.t2g_bulk_path}"
                          )
              make_quants: true
          }

        // If the quants are current and repeat_mapping is false
        // create tuple of metadata map (read from output), salmon output directory to use as input to merge_bulk_quants
        quants_ch = bulk_channel.has_quants
          .map{meta -> tuple(Utils.readMeta(file("${meta.salmon_results_dir}/scpca-meta.json")),
                             file(meta.salmon_results_dir)
                             )}

        // If we need to run salmon, create tuple of (metadata map, [Read 1 files], [Read 2 files])
        bulk_reads_ch = bulk_channel.make_quants
          .map{meta -> tuple(meta,
                             file("${meta.files_directory}/*_R1_*.fastq.gz"),
                             file("${meta.files_directory}/*_R2_*.fastq.gz")
                             )}

        // run fastp and salmon for libraries that are not skipping salmon
        fastp(bulk_reads_ch)
        salmon_ch = fastp.out
          .map{it.toList() + [file(it[0].salmon_bulk_index)]}
        salmon(salmon_ch)

        // group libraries together by project
        grouped_salmon_ch = salmon.out.mix(quants_ch)
          .map{[it[0].project_id,
                it[0],
                it[1]]} // salmon directories
          .groupTuple(by: 0)
          .map{[it[1][0], // meta; relevant data should all be the same by project, so take the first
                it[2], // salmon directories
                file(it[1][0].t2g_bulk_path)]}

        // create tsv file and combined metadata for each project containing all libraries
        merge_bulk_quants(grouped_salmon_ch, file(params.run_metafile))

    emit:
        bulk_counts = merge_bulk_quants.out.bulk_counts
        bulk_metadata = merge_bulk_quants.out.bulk_metadata
}
