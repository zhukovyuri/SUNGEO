#' Geocode addresses with OpenStreetMap
#'
#' Function to find geographic coordinates of addresses and place names, using OpenStreetMap's Nominatum API.
#'
#' @param query Address or place name to be geocoded. Character string.
#' @param match_num If query matches multiple locations, which match to return? Default is 1 (highest-ranking match, by relevance). Numeric.
#' @param return_all Should all matches be returned? Overrides \code{match_num} if \code{TRUE}. Default is \code{FALSE}. Logical.
#' @param details Should detailed results be returned? Default is \code{FALSE}. Logical.
#' @param user_agent Valid User-Agent identifying the application for OSM-Nominatum. If none supplied, function will attempt to auto-detect. Character string.
#' @return A \code{data.frame} object. If \code{details=FALSE}, contains fields
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
#' @details Note that Nominatim Usage Policy stipulates an absolute maximum of 1 request per second (\url{https://operations.osmfoundation.org/policies/nominatim/}). For batch geocoding of multiple addresses, please use \code{\link[SUNGEO]{geocode_osm_batch}}.
#' @import RCurl
#' @importFrom dplyr mutate_all slice mutate select starts_with
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom httr handle config status_code headers GET warn_for_status
#' @importFrom rlang abort
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

  # Internal functions
  request_GET <- function(x, url, ...) {
    x$response <- httr::GET(url, x$config, ..., handle = x$handle)
    x$html <- new.env(parent = emptyenv(), hash = FALSE)
    x$url <- x$response$url
    httr::warn_for_status(x$response)
    x
  }

  html_sessionSunGeo <- function(url, ...) {
    session <- structure(
      list(
        handle   = httr::handle(url),
        config   = c(..., httr::config(autoreferer = 1L)),
        url      = NULL,
        back     = character(),
        forward  = character(),
        response = NULL,
        html     = new.env(parent = emptyenv(), hash = FALSE)
      ),
      class = "session"
    )
    request_GET(session, url)
  }

  print.session <- function(x, ...) {
    cat("<session> ", x$url, "\n", sep = "")
    cat("  Status: ", httr::status_code(x$response), "\n", sep = "")
    cat("  Type:   ", httr::headers(x)$`Content-Type`, "\n", sep = "")
    cat("  Size:   ", length(x$response$content), "\n", sep = "")
  }

  # Batch geocoding warning
  if(length(query)>1){query <- query[1]; warning("Returning first result only. Please use geocode_osm_batch() to geocode multiple addresses.")}

  # Error handling
  downloadFail <- FALSE
  tryCatch({

  # Create empty objects
  osm_id=NA_character_; osm_type = NA_character_; address=NA_character_; longitude=NA_real_; latitude=NA_real_; importance=NA_real_; bbx=data.frame(bbox_ymin=NA_real_,bbox_ymax=NA_real_,bbox_xmin=NA_real_,bbox_xmax=NA_real_,stringsAsFactors = FALSE)

  # User-Agent
  if(length(user_agent)==0){user_agent <- html_sessionSunGeo("https://httpbin.org/user-agent")["response"][[1]]["request"][[1]]["options"][[1]]["useragent"][[1]]}

  # Send query to OSM Nominatum API
  root_url <- "https://nominatim.openstreetmap.org/search?q="
  sufx_url <- "&format=json&polygon=1&addressdetails=1"
  doc <- RCurl::getURL(utils::URLencode(paste0(root_url, query, sufx_url)),httpheader = c('User-Agent' = user_agent))

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
      osm_id <- as.character(dat$osm_id[match_num])
      osm_type <- as.character(dat$osm_type[match_num])
      address <- as.character(dat$display_name[match_num])
      longitude <- as.numeric(as.character(dat$lon[match_num]))
      latitude <- as.numeric(as.character(dat$lat[match_num]))
      importance <- as.numeric(as.character(dat$importance[match_num]))
      bbx <- dplyr::slice(dplyr::mutate_all(as.data.frame(do.call(rbind,dat$boundingbox),stringsAsFactors = FALSE),function(.){as.numeric(.)}),match_num)
      colnames(bbx) <- c("bbox_ymin","bbox_ymax","bbox_xmin","bbox_xmax")
    }
  }

  # Output
  out <- cbind(data.frame(query=query,
                    osm_id=osm_id,
                    osm_type=osm_type,
                    importance=importance,
                    address=address,
                    longitude=longitude,
                    latitude=latitude,
                    stringsAsFactors = FALSE), bbx)
  if(!details){out <- dplyr::select(dplyr::mutate(out,osm_type=NULL,importance=NULL),-dplyr::starts_with("bbox"))}
  return(as.data.frame(out))
  }
}
