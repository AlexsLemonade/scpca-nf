// Groovy functions used in multiple modules/processes

import groovy.json.JsonSlurper
import groovy.json.JsonOutput


class Utils {
  static def readMeta(file) {
    def meta = new JsonSlurper().parse(file)
    return(meta)
  }

  static def makeJson(meta) {
    def meta_json = JsonOutput.toJson(meta)
    meta_json = JsonOutput.prettyPrint(meta_json)
    return(meta_json)
  }
}
