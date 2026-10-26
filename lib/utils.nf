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

  if (meta.technology) {
    meta.technology = meta.technology.toLowerCase()
  }

  return (meta)
}

/**
  * Update a metadata map with values from a stored metadata file in JSON format
  * Stored values are used unless they are blank or NA, in which case the current value is kept
  * Keys that are only present in the stored file are kept, even if blank
  *
  * @param meta The current metadata map
  * @param meta_file A Nextflow file object with metadata in JSON format
  * @return A new metadata map
  */
def updateMetaFromFile(meta, meta_file) {
  // readMeta converts null and NA values to ""
  def stored_meta = readMeta(meta_file)
  // current values are only used where the stored value is blank or missing
  return (stored_meta + meta.findAll { key, _value -> stored_meta[key] == null || stored_meta[key] == '' })
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
  * Replace a null or a string with an NA value with ""
  * (which evaluates as false in boolean contexts)
  *
  * @param str A string, or any other value
  * @return The input unless it was null or an NA string or variant thereof, in which case returns ""
  */
def parseNA(str) {
  if (str == null) {
    ''
  }
  else if (str instanceof String) {
    // has to be a string to have NA vals replaced
    str.toLowerCase() in ['na', 'n/a', 'nan'] ? '' : str
  }
  else {
    // not a string, so just return the unmodified value (including false and 0)
    str
  }
}


/**
  * Make a map of versions for the workflow and nextflow, for adding to metadata
  *
  * @return A map with keys "scpca_version" and "nextflow_version"
  */
def getVersions() {
  [
    scpca_version: workflow.revision ?: workflow.manifest.version,
    nextflow_version: nextflow.version.toString(),
  ]
}


/**
  * Make a pullthrough container URL if a pullthrough URL is provided
  *
  * @param image_url The original image URL
  * @param pullthrough_url The pullthrough URL to use
  * @return The modified image URL with the pullthrough prefix
  */
def pullthroughContainer(image_url, pullthrough_url = ""){
  def container = image_url
  def pullthrough_prefixes = [
    "public.ecr.aws": "public_ecr_aws",
    "quay.io": "quay_io",
    "ghcr.io": "ghcr_io",
  ]
  if (pullthrough_url) {
    def registry = container.tokenize('/')[0]
    if (registry in pullthrough_prefixes.keySet()) {
      container = container.replaceFirst(registry, "${pullthrough_url}/${pullthrough_prefixes[registry]}")
    }
  }
  return container
}


/**
  * For a given directory, return its space ranger image files
  *
  * @param dir Directory to search image files(s) in
  * @param cytaimage Boolean if the image file(s) to grab are the cytaimage
  * @return Array of file paths
  */
def getImageFiles(dir, cytaimage = false) {

  // https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/image-image-recommendation
  def image_extensions = cytaimage
    ? ["tif", "tiff"]
    : ["tif", "tiff", "btf", "tf2", "tf8", "jpg", "jpeg", "qptiff"]

  def f = files("${dir}/*")
    .findAll { it.name.tokenize('.').last().toLowerCase() in image_extensions }

  return (f)
}
