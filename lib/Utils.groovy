
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

/**
 * Utility functions for scpca-nf
 */
class Utils {

   /**
     * Read a metadata file in JSON format
     *
     * @param file A Nextflow file object with metadata in JSON format
     * @return Metadata as a groovy Map object
     */
  static def readMeta(file) {
    def meta = new JsonSlurper().parse(file)
    return(meta)
  }

  /**
   * Write a metadata map to a JSON string
   *
   * @param meta A metadata map
   * @return A JSON string
   */
  static def makeJson(meta) {
    def meta_json = JsonOutput.toJson(meta)
    meta_json = JsonOutput.prettyPrint(meta_json)
    return(meta_json)
  }

  /**
   * Read an element from a metadata file in JSON format
   *
   * @param file A Nextflow file object with metadata in JSON format
   * @return A value from the metadata
   */
  static def getMetaVal(file, key){
    def obj = new JsonSlurper().parse(file)

    return(obj[key] ?: "")
  }


  /**
   * Replace a string with an NA value with ""
   * (which evaluates as false in boolean contexts)
   *
   * @param str A string
   * @return The input string unless it was NA or a variant thereof, in which case returns ""
   */
  static def parseNA(str) {
    if (str){
      str.toLowerCase() in ["na","n/a","nan"]? "" : str
    } else {
      ""
    }
   }
}
