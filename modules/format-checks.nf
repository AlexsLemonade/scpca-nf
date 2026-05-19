// formatting checks for sce files

include { pullthroughContainer } from '../lib/utils.nf'

process check_sce {
  container "${pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)}"
  label 'mem_16'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), 
          path(sce_file)
    path(reference_sce_file)
  output:
    tuple val(meta), 
          path("${meta.unique_id}_formatting_errors.txt")
  script:
    def object_type = sce_file.baseName.replace("${meta.library_id}_", "")
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

process compile_errors {
  container "${pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)}"
  publishDir "${params.outdir}", mode: 'copy'
  label 'mem_8'
  input:
   // more than one error file per library (unfiltered, filtered, processed), so use stageAs
   path error_files, stageAs: "error_file_*.txt"
  output:
    path "format_check_results.txt"
  script:
    """
    # check if all files are empty, if so print a success message
    # otherwise concatenate all error files into one output file
    errors="\$(cat ${error_files})"
    if [ -z \${errors} ]; then
      echo "No formatting errors found." > format_check_results.txt
    else
      echo \${errors} > format_check_results.txt
    fi
    """
  stub: 
    """
    touch format_check_results.txt
    """
}

workflow format_checks {
  take: 
    sce_ch // [ meta, unfiltered sce, filtered sce, processed sce, metadata json]
    sce_format_reference_file
  main: 

    sce_format_ch = sce_ch
      // drop the metadata json
      .map{ meta, unfiltered_sce, filtered_sce, processed_sce, metadata_json -> 
        [meta, [unfiltered_sce, filtered_sce, processed_sce]]
      } 
      .transpose() // [ [meta, unfiltered sce], [meta, filtered sce], [meta, processed sce] ]

    check_sce(sce_format_ch, sce_format_reference_file)

    // collect all error files and concatenate to print to a formatting errors output file
    error_input_ch = check_sce.out
      .collect{ meta, error_file -> error_file } // collect into a list of just the error files
    
    compile_errors(error_input_ch)

    // collect all error files and print out an error to the log file
    check_sce.out
      .filter{ meta, error_file -> error_file.size() > 0 } // only warn about libraries with errors
      .subscribe{ meta, error_file -> 
        log.error "Formatting errors for ${error_file.text}"
      }
    

}
