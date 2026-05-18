// formatting checks for sce files

include { pullthroughContainer } from '../lib/utils.nf'

process check_sce {
  container "${pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)}"
  label 'mem_16'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), 
          path(sce_file), 
          val(object_type)
    path(reference_sce_file)
  output:
    tuple val(meta), 
          path("${meta.unique_id}_formatting_errors.txt"),
          val(object_type)
  script:
    """
    sce_formatting_checks.R \
        --sce_file ${sce_file} \
        --object_type ${object_type} \
        --reference_file ${reference_sce_file} \
        --output_file ${meta.unique_id}_formatting_errors.txt
    """
  stub:
    """
    touch ${meta.unique_id}_formatting_errors.txt
    """
}

workflow format_checks {
  take: 
    sce_ch // [ meta, unfiltered sce, filtered sce, processed sce, metadata json]
  main: 

    sce_format_ch = sce_ch
      // drop the metadata json
      .map{ meta, unfiltered_sce, filtered_sce, processed_sce, metadata_json -> 
      [meta, [unfiltered_sce, filtered_sce, processed_sce]]
      } 
      .transpose() // [ [meta, unfiltered sce], [meta, filtered sce], [meta, processed sce] ]
      .map{ meta, sce -> 
        def object_type = sce.baseName.replace("${meta.library_id}_", "")
        [meta, sce, object_type]
      }

    check_sce(sce_format_ch, file(params.sce_format_reference_file))

    // check the contents of each error file and print any errors if present
    check_sce.out
      .map{ meta, error_file, object_type ->
        if (error_file.size() > 0) {
          log.error "Formatting errors for library ${meta.library_id} ${object_type}:\n${error_file.text}"
        }
      }

}
