// Utility functions for scpca-nf

/**
  * Read a metadata file in JSON format
  *
  * @param file A Nextflow file object with metadata in JSON format
  * @return Metadata as a groovy Map object
  */
def readMeta(file) {
  def meta = new groovy.json.JsonSlurper().parse(file)
  meta = meta.each { key, value -> meta[key] = parseNA(value) }

  if (meta.sample_id && meta.library_id && !meta.unique_id) {
    meta.unique_id = meta.technology.contains("_multi") ? "${meta.library_id}-${meta.sample_id}" : meta.library_id
  }

  return (meta)
}

/**
  * Write a metadata map to a JSON string
  *
  * @param meta A metadata map
  * @param updates A map of updates to the metadata
  * @return A JSON string
  */
def makeJson(meta, updates = [:]) {
  def meta_out = meta + updates
  // merge/overwrite the metadata with any updates
  def meta_json = groovy.json.JsonOutput.toJson(meta_out)
  meta_json = groovy.json.JsonOutput.prettyPrint(meta_json)
  return (meta_json)
}

/**
  * Read an element from a metadata file in JSON format
  *
  * @param file A Nextflow file object with metadata in JSON format
  * @return A value from the metadata
  */
def getMetaVal(file, key) {
  if (!file.exists()) {
    return (null)
  }

  def obj = new groovy.json.JsonSlurper().parse(file)
  def value = obj[key]

  if (value instanceof String) {
    value = parseNA(value)
  }

  return (value)
}


/**
  * Replace a string with an NA value with ""
  * (which evaluates as false in boolean contexts)
  *
  * @param str A string
  * @return The input string unless it was NA or a variant thereof, in which case returns ""
  */
def parseNA(str) {
  if (str) {
    if (str instanceof String) {
      // has to be a string to have NA vals replaced
      str.toLowerCase() in ['na', 'n/a', 'nan'] ? '' : str
    }
    else {
      // not a string, so just return the unmodified value
      str
    }
  }
  else {
    // all falsey values get turned into empty strings
    ''
  }
}


/**
  * Make a map of versions for the workflow and nextflow, for adding to metadata
  *
  * @param workflow A nextflow WorkflowMetadata object containing the revision and manifest
  * @param nextflow A NextflowMeta object (created by the nextflow runtime)
  * @return A map with keys "scpca_version" and "nextflow_version"
  */
def getVersions(WorkflowMetadata workflow, NextflowMeta nextflow) {
  [
    scpca_version: workflow.revision ?: workflow.manifest.version,
    nextflow_version: nextflow.version.toString(),
  ]
}
