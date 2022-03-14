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

}

process salmon{
    container params.SALMON_CONTAINER
    label 'cpus_12'
    label 'mem_24'
    tag "${meta.library_id}-bulk"
    publishDir "${meta.salmon_publish_dir}"
    input: 
        tuple val(meta), path(read_dir)
        path (index)
    output: 
        tuple val(meta), path(salmon_results_dir)
    script:
        salmon_results_dir = "${meta.library_id}"
        """
        salmon quant -i ${index} \
        -l A \
        ${meta.technology == 'paired_end' ? "-1": "-r"} ${read_dir}/*_R1_*.fastq.gz \
        ${meta.technology == 'paired_end' ? "-2 ${read_dir}/*_R2_*.fastq.gz" : "" } \
        -o ${salmon_results_dir} \
        --validateMappings \
        --rangeFactorizationBins 4 \
        --gcBias \
        --seqBias \
        --threads ${task.cpus}
        """

}

process merge_bulk_quants {
    container params.SCPCATOOLS_CONTAINER
    publishDir "${params.outdir}/publish/${project_id}"
    input:
        tuple val(project_id), path(salmon_directories)
        path(tx2gene)
        path(library_metadata)
    output:
        path(tximport_file), emit: bulk_counts
        path(bulk_metadata_file), emit: bulk_metadata
    script:
        tximport_file = "${project_id}_bulk_quant.tsv"
        bulk_metadata_file = "${project_id}_bulk_metadata.tsv"
        workflow_url = workflow.repository ?: workflow.manifest.homePage
        """
        ls -d ${salmon_directories} > salmon_directories.txt

        merge_counts_tximport.R \
          --project_id ${project_id} \
          --salmon_dirs salmon_directories.txt \
          --output_file ${tximport_file} \
          --tx2gene ${tx2gene}

        generate_bulk_metadata.R \
         --project_id ${project_id} \
         --salmon_dirs salmon_directories.txt \
         --library_metadata_file ${library_metadata} \
         --metadata_output ${bulk_metadata_file} \
         --genome_assembly ${params.assembly} \
         --workflow_url "${workflow_url}" \
         --workflow_version "${workflow.revision}" \
         --workflow_commit "${workflow.commitId}"
        """
}

workflow bulk_quant_rna {
    take: bulk_channel 
    // a channel with a map of metadata for each rna library to process
    main: 
        bulk_channel = bulk_channel
          // add salmon directory and salmon file location to meta 
          .map{it.salmon_publish_dir = "${params.outdir}/internal/salmon";
               it.salmon_results_dir = "${it.salmon_publish_dir}/${it.library_id}";
               it}
          // split based on whether repeat_mapping is false and the salmon quant.sf file exists 
          .branch{
              has_quants: !params.repeat_mapping && file(it.salmon_results_dir).exists()
              make_quants: true
          }
        
        // If we need to run salmon, create tuple of (metadata map, [Read 1 files], [Read 2 files])
        bulk_reads_ch = bulk_channel.make_quants
          .map{meta -> tuple(meta,
                             file("${meta.files_directory}/*_R1_*.fastq.gz"),
                             file("${meta.files_directory}/*_R2_*.fastq.gz")
                             )}

        // If the quant.sf file from salmon exits and repeat_mapping is false 
        // create tuple of metadata map, salmon output directory to use as input to merge_bulk_quants
        quants_ch = bulk_channel.has_quants
          .map{meta -> tuple(meta,
                             file(meta.salmon_results_dir)
                             )}

        
        // run fastp and salmon for libraries that are not skipping salmon
        fastp(bulk_reads_ch)
        salmon(fastp.out, params.bulk_index)

        // group libraries together by project
        grouped_salmon_ch = salmon.out.mix(quants_ch)
            .map{[it[0]["project_id"], it[1]]}
            .groupTuple(by: 0)

        // create tsv file and combined metadata for each project containing all libraries
        merge_bulk_quants(grouped_salmon_ch, params.t2g_bulk_path, params.run_metafile)

    emit: 
        bulk_counts = merge_bulk_quants.out.bulk_counts
        bulk_metadata = merge_bulk_quants.out.bulk_metadata
}
