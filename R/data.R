#' Constituency level results for lower chamber legislative elections, Germany 2009.
#'
#' A simple feature collection containing the spatial geometries of electoral contituency
#' borders, and data on turnout levels, votes shares and other attributes of lower chamber
#' legislative elections.
#'
#' @format Simple feature collection with 16 features and 10 fields.
#' geometry type:  MULTIPOLYGON.
#' dimension:      XY.
#' bbox:           xmin: 5.867281 ymin: 47.27096 xmax: 15.04388 ymax: 55.05902.
#' epsg (SRID):    4326.
#' proj4string:    +proj=longlat +datum=WGS84 +no_defs.
#' \describe{
#'   \item{cst }{Constituency number. Numeric.}
#'   \item{cst_n }{Constituency name. Character.}
#'   \item{ctr }{Country number. Numeric.}
#'   \item{ctr_n }{Country name. Character.}
#'   \item{yrmo }{Year and month of election (YYYYMM). Character.}
#'   \item{to1 }{Turnout in first round. Numeric.}
#'   \item{vv1 }{Number of valid votes in first round. Numeric.}
#'   \item{pvs1_margin }{Popular vote share margin in first round. Numeric.}
#'   \item{incumb_pty_n }{Incumbent party name.}
#'   \item{win1_pty_n }{Party name of popular vote share winner in first round. Character.}
#' }
#' @source Constituency-Level Elections Archive (CLEA) \url{http://www.electiondataarchive.org/}
"clea_deu2009"

#' Constituency level results for lower chamber legislative elections, Germany 2009.
#'
#' A simple feature collection containing the geographic centroids of electoral
#' contituencies, and data on turnout levels, votes shares and other attributes of
#' lower chamber legislative elections.
#'
#' @format Simple feature collection with 16 features and 10 fields.
#' geometry type:  POINT.
#' dimension:      XY.
#' bbox:           xmin: 6.953882 ymin: 48.54535 xmax: 13.40315 ymax: 54.18635.
#' epsg (SRID):    4326.
#' proj4string:    +proj=longlat +datum=WGS84 +no_defs.
#' \describe{
#'   \item{cst }{Constituency number. Numeric.}
#'   \item{cst_n }{Constituency name. Character.}
#'   \item{ctr }{Country number. Numeric.}
#'   \item{ctr_n }{Country name. Character.}
#'   \item{yrmo }{Year and month of election (YYYYMM). Character.}
#'   \item{to1 }{Turnout in first round. Numeric.}
#'   \item{vv1 }{Number of valid votes in first round. Numeric.}
#'   \item{pvs1_margin }{Popular vote share margin in first round. Numeric.}
#'   \item{incumb_pty_n }{Incumbent party name.}
#'   \item{win1_pty_n }{Party name of popular vote share winner in first round. Character.}
#' }
#' @source Constituency-Level Elections Archive (CLEA) \url{http://www.electiondataarchive.org/}
"clea_deu2009_pt"

#' Constituency level results for lower chamber legislative elections, Germany 2009.
#'
#' A data.frame object containing the geographic centroids of electoral
#' contituencies, and data on turnout levels, votes shares and other attributes of
#' lower chamber legislative elections.
#'
#' @format data.frame with 16 observations and 12 variables.
#' \describe{
#'   \item{cst }{Constituency number. Numeric.}
#'   \item{cst_n }{Constituency name. Character.}
#'   \item{ctr }{Country number. Numeric.}
#'   \item{ctr_n }{Country name. Character.}
#'   \item{yrmo }{Year and month of election (YYYYMM). Character.}
#'   \item{to1 }{Turnout in first round. Numeric.}
#'   \item{vv1 }{Number of valid votes in first round. Numeric.}
#'   \item{pvs1_margin }{Popular vote share margin in first round. Numeric.}
#'   \item{incumb_pty_n }{Incumbent party name.}
#'   \item{win1_pty_n }{Party name of popular vote share winner in first round. Character.}
#'   \item{longitude }{Longitude of constituency centroid. Numeric.}
#'   \item{latitude }{Latitude of constituency centroid. Numeric.}
#' }
#' @source Constituency-Level Elections Archive (CLEA) \url{http://www.electiondataarchive.org/}
"clea_deu2009_df"

#' Population count raster for Germany, 2010.
#'
#' 2.5 arc-minute resolution raster of estimates of human population (number of persons per pixel),
#' consistent with national censuses and population registers, for the year 2010.
#'
#' @format
#' class       : RasterLayer.
#' dimensions  : 186, 220, 40920  (nrow, ncol, ncell).
#' resolution  : 0.04166667, 0.04166667  (x, y).
#' extent      : 5.875, 15.04167, 47.29167, 55.04167  (xmin, xmax, ymin, ymax).
#' coord. ref. : +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0.
#' data source : in memory.
#' names       : gpw_v4_population_count_rev11_2010_2pt5_min.
#' values      : 0, 92915.66  (min, max).
#' @source Gridded Population of the World (GPW) v4: Population Count, v4.11 \url{https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-rev11}
"gpw4_deu2010"

#' Hexagonal grid for Germany.
#'
#' Regular hexagonal grid of 0.5 degree diameter cells, covering territory of Germany (2020 borders).
#'
#' @format Simple feature collection with 257 features and 3 fields.
#' geometry type:  POLYGON.
#' dimension:      XY.
#' bbox:           xmin: 5.375001 ymin: 46.76568 xmax: 15.375 ymax: 55.13726.
#' epsg (SRID):    4326.
#' proj4string:    +proj=longlat +datum=WGS84 +no_defs.
#' \describe{
#'   \item{HEX_ID }{Unique cell identifier. Character.}
#'   \item{HEX_X }{Longitude of cell centroid. Numeric.}
#'   \item{HEX_Y }{Latitude of cell centroid. Numeric.}
#' }
#' @source SUNGEO
"hex_05_deu"

#' Roads polylines for Germany, 1992
#'
#' Roads thematic layer from Digital Chart of the World. Subset: divided multi-lane highways.
#'
#' @format
#' Simple feature collection with 1741 features and 5 fields.
#' geometry type:  MULTILINESTRING.
#' dimension:      XY.
#' bbox:           xmin: 5.750933 ymin: 47.58799 xmax: 14.75109 ymax: 54.80712
#' epsg (SRID):    4326.
#' proj4string:    +proj=longlat +datum=WGS84 +no_defs.
#' \describe{
#'   \item{MED_DESCRI }{Is the road a divided multi-lane highway with a median? Character string.}
#'   \item{RTT_DESCRI }{Primary or secondary route? Character string.}
#'   \item{F_CODE_DES }{Feature code description (road or trail). Character string.}
#'   \item{ISO }{ISO 3166-1 alpha-3 country code. Character string.}
#'   \item{ISOCOUNTRY }{Country name. Character string.}
#' }
#' @source Defense Mapping Agency (DMA), 1992. Digital Chart of the World. Defense Mapping Agency, Fairfax, Virginia. (Four CD-ROMs). Available through DIVA-GIS: \url{http://www.diva-gis.org/gData} (accessed April 15, 2020).
"highways_deu1992"
