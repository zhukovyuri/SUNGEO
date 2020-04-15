#' Geocode addresses with OpenStreetMap
#'
#' Function to find geographic coordinates of addresses and place names, using OpenStreetMap's Nominatum API.
#'
#' @param query Address or place name to be geocoded. Character string.
#' @param match_num If query matches multiple locations, which match to return? Default is 1 (highest-ranking match, by relevance). Numeric.
#' @param return_all Should all matches be returned? Overrides \code{match_num} if \code{TRUE}. Default is \code{FALSE}. Logical.
#' @param details Should detailed results be returned? Default is \code{FALSE}. Logical.
#' @param user_agent Valid User-Agent identifying the application for OSM-Nominatum. If none supplied, function will attempt to auto-detect. Character string.
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
#' @details Note that Nominatim Usage Policy stipulates an absolute maximum of 1 request per second (\url{https://operations.osmfoundation.org/policies/nominatim/}). For batch geocoding of multiple addresses, please use \code{\link[SUNGEO]{geocode_osm_batch}} or \code{\link[SUNGEO]{geocode_gn}}.
#' @import data.table tidyverse RCurl jsonlite
#' @importFrom data.table last first between
#' @importFrom rvest html_session
#' @examples
#' # Geocode an address (top match only)
#' \dontrun{
#' geocode_osm("Michigan Stadium")
#' }
#' # Return detailed results for top match
#' \dontrun{
#' geocode_osm("Michigan Stadium", details=TRUE)
#' }
#' # Return detailed results for all matches
#' \dontrun{
#' geocode_osm("Michigan Stadium", details=TRUE, return_all = TRUE)
#' }
#' @export

geocode_osm <- function(
  query,
  match_num = 1,
  return_all = FALSE,
  details = FALSE,
  user_agent = NULL
){

  # Batch geocoding warning
  if(length(query)>1){query <- query[1]; warning("Returning first result only. Please use geocode_osm_batch() or geocode_gn() to geocode multiple addresses.")}

  # Error handling
  downloadFail <- FALSE
  tryCatch({

  # Create empty objects
  osm_id=NA_character_; osm_type = NA_character_; address=NA_character_; longitude=NA_real_; latitude=NA_real_; importance=NA_real_; bbx=data.frame(bbox_ymin=NA_real_,bbox_ymax=NA_real_,bbox_xmin=NA_real_,bbox_xmax=NA_real_,stringsAsFactors = FALSE)

  # User-Agent
  if(length(user_agent)==0){user_agent <- rvest::html_session("https://httpbin.org/user-agent") %>%
    (function(.){.["response"][[1]]["request"][[1]]["options"][[1]]["useragent"][[1]]})}

  # Send query to OSM Nominatum API
  root_url <- "https://nominatim.openstreetmap.org/search?q="
  sufx_url <- "&format=json&polygon=1&addressdetails=1"
  doc <- paste0(root_url, query, sufx_url) %>% URLencode() %>% RCurl::getURL(httpheader = c('User-Agent' = user_agent))

  }, warning = function(w) {
    downloadFail <<- TRUE
  }, error = function(e) {
    message("Cannot access OSM server. Please check your internet connection and try again.")
    downloadFail <<- TRUE
  }, finally = {
  })
  if(downloadFail){
    cat("Cannot access OSM server. Please check your internet connection and try again.")
    return()
  } else {

  # Parse results
  if(nchar(doc)>2){
    dat <- jsonlite::fromJSON(doc)
    if(return_all){match_num <- 1:nrow(dat)}
    if(nrow(dat)>0){
      osm_id <- dat$osm_id[match_num] %>% as.character()
      osm_type <- dat$osm_type[match_num] %>% as.character()
      address <- dat$display_name[match_num] %>% as.character()
      longitude <- dat$lon[match_num] %>% as.character() %>% as.numeric()
      latitude <- dat$lat[match_num] %>% as.character() %>% as.numeric()
      importance <- dat$importance[match_num] %>% as.character() %>% as.numeric()
      bbx <- dat$boundingbox %>% do.call(rbind,.) %>% as.data.frame(stringsAsFactors = FALSE) %>% data.table::setnames(names(.),c("bbox_ymin","bbox_ymax","bbox_xmin","bbox_xmax")) %>% dplyr::mutate_all(function(.){as.numeric(.)}) %>% dplyr::slice(match_num)
    }
  }

  # Output
  out <- data.frame(query=query,
                    osm_id=osm_id,
                    osm_type=osm_type,
                    importance=importance,
                    address=address,
                    longitude=longitude,
                    latitude=latitude,
                    stringsAsFactors = FALSE) %>% cbind(bbx)
  if(!details){out <- out %>% dplyr::mutate(osm_type=NULL,importance=NULL) %>% dplyr::select(-dplyr::starts_with("bbox"))}
  return(out %>% as.data.table())

  }
}
