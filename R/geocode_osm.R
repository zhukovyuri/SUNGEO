#' Geocode addresses and coordinates with OpenStreetMap
#'
#' Finds geographic coordinates of addresses and place names (forward geocoding),
#' or converts longitude/latitude coordinates to place names and administrative
#' units (reverse geocoding), using OpenStreetMap's Nominatim API.
#'
#' @param query Address or place name to be geocoded. Character string. Ignored
#'   if \code{reverse=TRUE}.
#' @param match_num If query matches multiple locations, which match to return?
#'   Default is 1 (highest-ranking match, by relevance). Numeric.
#' @param return_all Should all matches be returned? Overrides \code{match_num}
#'   if \code{TRUE}. Default is \code{FALSE}. Logical.
#' @param details Should detailed results be returned? Default is \code{FALSE}.
#'   Logical. Ignored if \code{reverse=TRUE}.
#' @param reverse Should reverse geocoding be performed (coordinates to address)?
#'   Default is \code{FALSE}. Logical.
#' @param lon Longitude of the point to reverse geocode. Required if
#'   \code{reverse=TRUE}. Numeric.
#' @param lat Latitude of the point to reverse geocode. Required if
#'   \code{reverse=TRUE}. Numeric.
#' @param zoom Zoom level for reverse geocoding, controlling the level of detail
#'   returned (0 = country, 10 = city, 18 = building). Default is 10. Numeric.
#' @param user_agent Valid User-Agent identifying the application for
#'   OSM Nominatim. If none supplied, function will attempt to auto-detect.
#'   Character string.
#' @return A \code{data.frame} object.
#'
#'   For forward geocoding (\code{reverse=FALSE}):
#'   \itemize{
#'   \item \code{query}. User-supplied address query. Character string.
#'   \item \code{osm_id}. OpenStreetMap ID. Character string.
#'   \item \code{address}. OpenStreetMap address. Character string.
#'   \item \code{longitude}. Horizontal coordinate. Numeric.
#'   \item \code{latitude}. Vertical coordinate. Numeric.
#'   }
#'   If \code{details=TRUE}, also contains:
#'   \itemize{
#'   \item \code{osm_type}. OpenStreetMap feature type. Character string.
#'   \item \code{importance}. Relevance of Nominatim match to query, from 0
#'     (worst) to 1 (best). Numeric.
#'   \item \code{bbox_ymin}. Minimum vertical coordinate of bounding box. Numeric.
#'   \item \code{bbox_ymax}. Maximum vertical coordinate of bounding box. Numeric.
#'   \item \code{bbox_xmin}. Minimum horizontal coordinate of bounding box. Numeric.
#'   \item \code{bbox_xmax}. Maximum horizontal coordinate of bounding box. Numeric.
#'   }
#'
#'   For reverse geocoding (\code{reverse=TRUE}):
#'   \itemize{
#'   \item \code{osm_name}. Full display name from OpenStreetMap. Character string.
#'   \item \code{osm_adm0_code}. ISO country code. Character string.
#'   \item \code{osm_adm0}. Country name. Character string.
#'   \item \code{osm_adm1}. State or first-level administrative division. Character string.
#'   \item \code{osm_adm2}. District or county. Character string.
#'   \item \code{osm_adm3}. Municipality. Character string.
#'   \item \code{osm_adm4}. Borough, village or street. Character string.
#'   }
#' @details Note that the Nominatim Usage Policy stipulates an absolute maximum
#'   of 1 request per second
#'   (\url{https://operations.osmfoundation.org/policies/nominatim/}).
#'   For batch geocoding of multiple addresses, please use
#'   \code{\link[SUNGEO]{geocode_osm_batch}}.
#' @importFrom curl curl_version
#' @importFrom data.table data.table :=
#' @importFrom dplyr mutate_all slice mutate select starts_with
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom httr handle config status_code headers GET warn_for_status add_headers content
#' @importFrom rlang abort
#' @examples
#' # Geocode an address (top match only)
#' geocode_osm("Michigan Stadium")
#' # Return detailed results for top match
#' geocode_osm("Michigan Stadium", details = TRUE)
#' # Return detailed results for all matches
#' geocode_osm("Michigan Stadium", details = TRUE, return_all = TRUE)
#' # Reverse geocode a coordinate pair
#' geocode_osm(reverse = TRUE, lon = -83.74868, lat = 42.26587, zoom = 18)
#' @export

geocode_osm <- function(
    query,
    match_num   = 1,
    return_all  = FALSE,
    details     = FALSE,
    reverse     = FALSE,
    lon         = NULL,
    lat         = NULL,
    zoom        = 10,
    user_agent  = NULL
) {

  # Internal helpers

  request_GET <- function(x, url, ...) {
    x$response <- httr::GET(url, x$config, ..., handle = x$handle)
    x$html     <- new.env(parent = emptyenv(), hash = FALSE)
    x$url      <- x$response$url
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

  #' @exportS3Method print session
  print.session <- function(x, ...) {
    cat("<session> ", x$url, "\n", sep = "")
    cat("  Status: ", httr::status_code(x$response), "\n", sep = "")
    cat("  Type:   ", httr::headers(x)$`Content-Type`, "\n", sep = "")
    cat("  Size:   ", length(x$response$content), "\n", sep = "")
  }

  # Reverse geocoder

  if (isTRUE(reverse)) {
    if (is.null(lon) || is.null(lat)) {
      rlang::abort("When reverse=TRUE, both `lon` and `lat` must be supplied.")
    }

    reverse_osm <- function(lon, lat, zoom = 10, user_agent = NULL) {
      url <- sprintf(
        "https://nominatim.openstreetmap.org/reverse?lat=%f&lon=%f&format=json&zoom=%d&addressdetails=1",
        lat, lon, zoom
      )
      res <- tryCatch(
        httr::GET(url, httr::add_headers("User-Agent" = user_agent %||% "r-geocoder")),
        error = function(e) NULL
      )
      if (is.null(res)) return(NULL)

      doc <- jsonlite::fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))
      if (is.null(doc$address)) return(NULL)

      out <- data.table::data.table(
        osm_name      = doc$display_name,
        osm_adm0_code = NA_character_,
        osm_adm0      = NA_character_,
        osm_adm1      = NA_character_,
        osm_adm2      = NA_character_,
        osm_adm3      = NA_character_,
        osm_adm4      = NA_character_
      )
      if ("country_code" %in% names(doc$address)) {
        out <- out[, osm_adm0_code := doc$address$country_code]
      }
      if ("country" %in% names(doc$address)) {
        out <- out[, osm_adm0 := doc$address$country]
      }
      if ("state" %in% names(doc$address)) {
        out <- out[, osm_adm1 := doc$address$state]
      }
      if ("district" %in% names(doc$address)) {
        out <- out[, osm_adm2 := doc$address$district]
      }
      if ("county" %in% names(doc$address) & !"district" %in% names(doc$address)) {
        out <- out[, osm_adm2 := doc$address$county]
      }
      if ("municipality" %in% names(doc$address)) {
        out <- out[, osm_adm3 := doc$address$municipality]
      }
      if ("city" %in% names(doc$address)) {
        out <- out[, osm_adm3 := doc$address$city]
      }
      if ("borough" %in% names(doc$address)) {
        out <- out[, osm_adm4 := doc$address$borough]
      }
      if ("village" %in% names(doc$address)) {
        out <- out[, osm_adm4 := doc$address$village]
      }
      if ("road" %in% names(doc$address) & !"village" %in% names(doc$address)) {
        out <- out[, osm_adm4 := doc$address$road]
      }
      rm(lon, lat, user_agent, doc, res, url, zoom)
      return(out)
    }

    # User-Agent auto-detect
    if (length(user_agent) == 0) {
      user_agent <- tryCatch({
        jsonlite::fromJSON(httr::content(httr::GET("https://httpbin.org/user-agent"), as = "text", encoding = "UTF-8"))$`user-agent`
      }, error = function(e) "r-geocoder")
    }

    result <- reverse_osm(lon = lon, lat = lat, zoom = zoom, user_agent = user_agent)
    if (is.null(result)) {
      message("Cannot access OSM server or no result returned. Please check your internet connection and coordinates.")
      return(invisible(NULL))
    }
    return(as.data.frame(result))
  }

  # Forward geocoder

  # Batch geocoding warning
  if (length(query) > 1) {
    query <- query[1]
    warning("Returning first result only. Please use geocode_osm_batch() to geocode multiple addresses.")
  }

  # Error handling
  downloadFail <- FALSE
  tryCatch({

    # Create empty objects
    osm_id     <- NA_character_
    osm_type   <- NA_character_
    address    <- NA_character_
    longitude  <- NA_real_
    latitude   <- NA_real_
    importance <- NA_real_
    bbx <- data.frame(
      bbox_ymin = NA_real_, bbox_ymax = NA_real_,
      bbox_xmin = NA_real_, bbox_xmax = NA_real_,
      stringsAsFactors = FALSE
    )

    # User-Agent
    if (length(user_agent) == 0) {
      user_agent <- tryCatch({
        jsonlite::fromJSON(httr::content(httr::GET("https://httpbin.org/user-agent"), as = "text", encoding = "UTF-8"))$`user-agent`
      }, error = function(e) "r-geocoder")
    }


    # Send query to OSM Nominatim API
    root_url <- "https://nominatim.openstreetmap.org/search?q="
    sufx_url <- "&format=json&polygon=1&addressdetails=1"
    response <- httr::GET(
      utils::URLencode(paste0(root_url, query, sufx_url), repeated = TRUE),
      httr::add_headers("User-Agent" = user_agent)
    )
    doc <- httr::content(response, as = "text", encoding = "UTF-8")

  }, warning = function(w) {
    downloadFail <<- TRUE
  }, error = function(e) {
    message("Cannot access OSM server. Please check your internet connection and try again.")
    downloadFail <<- TRUE
  }, finally = {
  })

  if (downloadFail) {
    cat("Cannot access OSM server. Please check your internet connection and try again.")
    return()
  } else {

    # Parse results
    if (nchar(doc) > 2) {
      dat <- jsonlite::fromJSON(doc)
      if (return_all) { match_num <- 1:nrow(dat) }
      if (nrow(dat) > 0) {
        osm_id     <- as.character(dat$osm_id[match_num])
        osm_type   <- as.character(dat$osm_type[match_num])
        address    <- as.character(dat$display_name[match_num])
        longitude  <- as.numeric(as.character(dat$lon[match_num]))
        latitude   <- as.numeric(as.character(dat$lat[match_num]))
        importance <- as.numeric(as.character(dat$importance[match_num]))
        bbx <- dplyr::slice(
          dplyr::mutate_all(
            as.data.frame(do.call(rbind, dat$boundingbox), stringsAsFactors = FALSE),
            function(.) { as.numeric(.) }
          ),
          match_num
        )
        colnames(bbx) <- c("bbox_ymin", "bbox_ymax", "bbox_xmin", "bbox_xmax")
      }
    }

    # Output
    out <- cbind(
      data.frame(
        query      = query,
        osm_id     = osm_id,
        osm_type   = osm_type,
        importance = importance,
        address    = address,
        longitude  = longitude,
        latitude   = latitude,
        stringsAsFactors = FALSE
      ),
      bbx
    )
    if (!details) {
      out <- dplyr::select(
        dplyr::mutate(out, osm_type = NULL, importance = NULL),
        -dplyr::starts_with("bbox")
      )
    }
    return(as.data.frame(out))
  }
}
