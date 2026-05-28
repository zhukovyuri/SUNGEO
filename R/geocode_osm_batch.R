#' Batch geocode addresses and coordinates with OpenStreetMap
#'
#' Finds geographic coordinates of multiple addresses and place names (forward
#' geocoding), or converts multiple longitude/latitude coordinate pairs to place
#' names and administrative units (reverse geocoding), using OpenStreetMap's
#' Nominatim API.
#'
#' @param query Addresses or place names to be geocoded. Character string.
#'   Ignored if \code{reverse=TRUE}.
#' @param delay Delay between requests, in seconds. Default is 1. Numeric.
#' @param match_num If query matches multiple locations, which match to return?
#'   Default is 1 (highest-ranking match, by relevance). Numeric.
#' @param return_all Should all matches be returned? Overrides \code{match_num}
#'   if \code{TRUE}. Default is \code{FALSE}. Logical.
#' @param details Should detailed results be returned? Default is \code{FALSE}.
#'   Logical. Ignored if \code{reverse=TRUE}.
#' @param reverse Should reverse geocoding be performed (coordinates to
#'   address)? Default is \code{FALSE}. Logical.
#' @param lon Longitudes of points to reverse geocode. Required if
#'   \code{reverse=TRUE}. Must be the same length as \code{lat}. Numeric vector.
#' @param lat Latitudes of points to reverse geocode. Required if
#'   \code{reverse=TRUE}. Must be the same length as \code{lon}. Numeric vector.
#' @param zoom Zoom level for reverse geocoding, controlling the level of detail
#'   returned (0 = country, 10 = city, 18 = building). Default is 10. Numeric.
#' @param user_agent Valid User-Agent identifying the application for
#'   OSM Nominatim. If none supplied, function will attempt to auto-detect.
#'   Character string.
#' @param verbose Print status messages and progress? Default is \code{FALSE}.
#'   Logical.
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
#'   \item \code{osm_adm4}. Borough or village. Character string.
#'   }
#' @details Wrapper function for \code{\link[SUNGEO]{geocode_osm}}. Because the
#'   Nominatim Usage Policy stipulates an absolute maximum of 1 request per
#'   second, this function facilitates batch geocoding by adding a small delay
#'   between queries
#'   (\url{https://operations.osmfoundation.org/policies/nominatim/}).
#' @importFrom dplyr bind_rows
#' @examples
#' \donttest{
#' # Geocode multiple addresses (top matches only)
#' geocode_osm_batch(c("Ann Arbor", "East Lansing", "Columbus"))
#' }
#' \donttest{
#' # With progress reports
#' geocode_osm_batch(c("Ann Arbor", "East Lansing", "Columbus"), verbose = TRUE)
#' }
#' \donttest{
#' # Return detailed results for all matches
#' geocode_osm_batch(c("Ann Arbor", "East Lansing", "Columbus"),
#'   details = TRUE, return_all = TRUE)
#' }
#' \donttest{
#' # Reverse geocode multiple coordinates (city level)
#' geocode_osm_batch(
#'   reverse = TRUE,
#'   lon = c(-83.743, -84.556, -82.999),
#'   lat = c(42.278,  42.701,  39.961),
#'   zoom = 10
#' )
#' }
#' @export

geocode_osm_batch <- function(
    query      = NULL,
    delay      = 1,
    return_all = FALSE,
    match_num  = 1,
    details    = FALSE,
    reverse    = FALSE,
    lon        = NULL,
    lat        = NULL,
    zoom       = 10,
    user_agent = NULL,
    verbose    = FALSE
) {

  # Enforce minimum delay time
  if (delay < 1) {
    delay <- 1
    warning("OSM Nominatim Usage Policy requires a maximum of 1 request per second.")
  }

  # Reverse geocoder

  if (isTRUE(reverse)) {
    if (is.null(lon) || is.null(lat)) {
      stop("When reverse=TRUE, both `lon` and `lat` must be supplied.")
    }
    if (length(lon) != length(lat)) {
      stop("`lon` and `lat` must be the same length.")
    }

    osm_mat <- as.data.frame(
      dplyr::bind_rows(
        lapply(seq_along(lon), function(i) {
          geo_i <- data.frame(
            osm_name      = NA_character_,
            osm_adm0_code = NA_character_,
            osm_adm0      = NA_character_,
            osm_adm1      = NA_character_,
            osm_adm2      = NA_character_,
            osm_adm3      = NA_character_,
            osm_adm4      = NA_character_,
            stringsAsFactors = FALSE
          )
          tryCatch({
            geo_i <- geocode_osm(
              reverse    = TRUE,
              lon        = lon[i],
              lat        = lat[i],
              zoom       = zoom,
              user_agent = user_agent
            )
          }, error = function(e) {})
          Sys.sleep(delay)
          if (verbose) {
            print(paste0(i, "/", length(lon), " ", paste0(geo_i$osm_name, collapse = "; ")))
          }
          geo_i
        })
      )
    )
    return(osm_mat)
  }

  # Forward geocoder

  osm_mat <- as.data.frame(
    dplyr::bind_rows(
      lapply(seq_along(query), function(i) {
        geo_i <- data.frame(
          query      = query[i],
          osm_id     = NA_character_,
          osm_type   = NA_character_,
          importance = NA_real_,
          address    = NA_character_,
          longitude  = NA_real_,
          latitude   = NA_real_,
          bbox_ymin  = NA_real_,
          bbox_ymax  = NA_real_,
          bbox_xmin  = NA_real_,
          bbox_xmax  = NA_real_,
          stringsAsFactors = FALSE
        )
        tryCatch({
          geo_i <- geocode_osm(
            query[i],
            match_num  = match_num,
            return_all = return_all,
            details    = details,
            user_agent = user_agent
          )
        }, error = function(e) {})
        Sys.sleep(delay)
        if (verbose) {
          print(paste0(i, "/", length(query), " ", paste0(geo_i$address, collapse = "; ")))
        }
        geo_i
      })
    )
  )

  return(osm_mat)
}
