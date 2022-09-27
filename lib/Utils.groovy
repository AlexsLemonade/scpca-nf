// Groovy functions used in multiple modules/processes

import groovy.json.JsonSlurper
import groovy.json.JsonOutput

class Utils {
  def readMeta(path) {
    meta = new JsonSlurper().parse(file(path))
    return(meta)
  }

  def makeJson(meta) {
    meta_json = JsonOutput.toJson(meta)
    meta_json = JsonOutput.prettyPrint(meta_json)
    return(meta_json)
  }
}
