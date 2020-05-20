#' \code{SUNGEO} package
#'
#' Sub-National Geospatial Data Archive System: Geoprocessing Toolkit
#'
#' See the README on
#' \href{https://github.com/zhukovyuri/SUNGEO#readme}{GitHub}
#'
#' @docType package
#' @name SUNGEO
#' @importFrom dplyr %>%
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(".",":=","AREA_W","POP_INT","POP_TOTAL","POP_W","TERM","alternatenames","asciiname","feature_code","geometry","geonameid","gn_type","latitude","longitude","match_type","name","sgeoz_length","strdist","area_code","col.id","overlap_1","overlap_2"))
  }
