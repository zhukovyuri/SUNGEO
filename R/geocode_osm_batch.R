#' Batch geocode addresses with OpenStreetMap
#'
#' Function to find geographic coordinates of multiple addresses and place names, using OpenStreetMap's Nominatum API.
#'
#' @param query Addresses or place names to be geocoded. Character string.
#' @param delay Delay between requests. Default is 1 second. Numeric.
#' @param match_num If query matches multiple locations, which match to return? Default is 1 (highest-ranking match, by relevance). Numeric.
#' @param return_all Should all matches be returned? Overrides \code{match_num} if \code{TRUE}. Default is \code{FALSE}. Logical.
#' @param details Should detailed results be returned? Default is \code{FALSE}. Logical.
#' @param user_agent Valid User-Agent identifying the application for OSM-Nominatum. If none supplied, function will attempt to auto-detect. Character string.
#' @param verbose Print status messages and progress? Default is \code{FALSE}. Logical.
#' @return A \code{data.table} object. If \code{details=FALSE}, contains fields
#' \itemize{
##'  \item{"query". }{User-supplied address query(ies). Character string.}
##'  \item{"osm_id". }{OpenStreetMap ID. Character string.}
##'  \item{"address". }{OpenStreetMap address. Character string.}
##'  \item{"longitude". }{Horizontal coordinate. Numeric.}
##'  \item{"latitude". }{Vertical coordinate. Numeric.}
##'  }
##' If \code{details=TRUE}, contains additional fields
#' \itemize{
##'  \item{"osm_type". }{OpenStreetMap ID. Character string.}
##'  \item{"importance". }{Relevance of Nominatum match to query, from 0 (worst) to 1 (best). Numeric.}
##'  \item{"bbox_ymin". }{Minimum vertical coordinate of bounding box. Numeric.}
##'  \item{"bbox_ymax". }{Maximum vertical coordinate of bounding box. Numeric.}
##'  \item{"bbox_xmin". }{Minimum horizontal coordinate of bounding box. Numeric.}
##'  \item{"bbox_xmax". }{Maximum horizontal coordinate of bounding box. Numeric.}
##'  }
#' @details Wrapper function for \code{\link[SUNGEO]{geocode_osm}}. Because Nominatim Usage Policy stipulates an absolute maximum of 1 request per second, this function facilitates batch geocoding by adding a small delay between queries (\url{https://operations.osmfoundation.org/policies/nominatim/}).
#' @import data.table tidyverse RCurl jsonlite
#' @importFrom data.table last first between
#' @importFrom rvest html_session
#' @importFrom dplyr bind_rows
#' @examples
#' # Geocode multiple addresses (top matches only)
#' \dontrun{
#' geocode_osm_batch(c("Ann Arbor","East Lansing","Columbus"))
#' }
#' # With progress reports
#' \dontrun{
#' geocode_osm_batch(c("Ann Arbor","East Lansing","Columbus"), verbose = TRUE)
#' }
#' # Return detailed results for all matches
#' \dontrun{
#' geocode_osm_batch(c("Ann Arbor","East Lansing","Columbus"),
#'                   details = TRUE, return_all = TRUE)
#' }
#' @export

geocode_osm_batch <- function(
  query,
  delay = 1,
  return_all = FALSE,
  match_num = 1,
  details = FALSE,
  user_agent = NULL,
  verbose = FALSE){

  # Enforce minimum delay time
  if(delay<1){delay <- 1; warning("OSM Nominatum Usage Policy requires maximum of 1 request per second.")}

  # Loop over queries
  osm_mat <- lapply(1:length(query),function(i){
    geo_i <- data.frame(query = query[i], osm_id=NA_character_,osm_type=NA_character_,importance=NA_real_,address=NA_character_,longitude=NA_real_,latitude=NA_real_,bbox_ymin=NA_real_,bbox_ymax=NA_real_,bbox_xmin=NA_real_,bbox_xmax=NA_real_,stringsAsFactors = FALSE)
    tryCatch({
      geo_i <- geocode_osm(query[i], match_num = match_num, return_all = return_all, details = details, user_agent = user_agent)
    }, error=function(e){})
    Sys.sleep(delay)
    if(verbose){print(paste0(i,"/",length(query)," ",geo_i$address %>% paste0(collapse="; ")))}
    geo_i
  }) %>% dplyr::bind_rows() %>% as.data.table()

  # Output
  return(osm_mat)
}
