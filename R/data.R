#' Constituency level results for lower chamber legislative elections, Germany 2009.
#'
#' A simple feature collection containing the spatial geometries of electoral constituency
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
#' @source Constituency-Level Elections Archive (CLEA) \url{https://electiondataarchive.org/}
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
#' @source Constituency-Level Elections Archive (CLEA) \url{https://electiondataarchive.org/}
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
#' @source Constituency-Level Elections Archive (CLEA) \url{https://electiondataarchive.org/}
"clea_deu2009_df"

#' Population count raster for Germany, 2010.
#'
#' 2.5 arc-minute resolution raster of estimates of human population (number of persons per pixel),
#' consistent with national censuses and population registers, for the year 2010.
#'
#' @format
#' class       : SpatRaster
#' dimensions  : 186, 220, 1  (nrow, ncol, nlyr)
#' resolution  : 0.04166667, 0.04166667  (x, y)
#' extent      : 5.875, 15.04167, 47.29167, 55.04167  (xmin, xmax, ymin, ymax)
#' coord. ref. : lon/lat WGS 84 (EPSG:4326)
#' source(s)   : memory
#' name        : gpw_v4_population_count_rev11_2010_2pt5_min
#' min value   :                                        0.00
#' max value   :                                    92915.66
#' @source Gridded Population of the World (GPW) v4: Population Count, v4.11 <doi:10.7927/H4JW8BX5>.
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
#' @format Simple feature collection with 1741 features and 5 fields.
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
#' @source Defense Mapping Agency (DMA), 1992. Digital Chart of the World. Defense Mapping Agency, Fairfax, Virginia. (Four CD-ROMs). Available through DIVA-GIS: \url{http://www.diva-gis.org/gData} (accessed August 12, 2021).
"highways_deu1992"

#' Data availability through SUNGEO API
#'
#' Census of geospatial and processed data files available to download using \code{SUNGEO::get_data()}.
#'
#' @format List of 42 data.table objects
#' Geoset:GADM :Classes ‘data.table’ and 'data.frame': 249 obs. of  4 variables
#' Geoset:GAUL :Classes ‘data.table’ and 'data.frame': 242 obs. of  4 variables
#' Geoset:geoBoundaries :Classes ‘data.table’ and 'data.frame': 197 obs. of  4 variables
#' Geoset:GRED :Classes ‘data.table’ and 'data.frame': 74 obs. of  4 variables
#' Geoset:HEXGRID :Classes ‘data.table’ and 'data.frame': 199 obs. of  4 variables
#' Geoset:MPIDR :Classes ‘data.table’ and 'data.frame': 52 obs. of  4 variables
#' Geoset:NHGIS :Classes ‘data.table’ and 'data.frame': 1 obs. of  4 variables
#' Geoset:PRIOGRID :Classes ‘data.table’ and 'data.frame': 199 obs. of  4 variables
#' Geoset:SHGIS :Classes ‘data.table’ and 'data.frame': 68 obs. of  4 variables
#' \describe{
#'   \item{country_iso3 }{Codes for available countries (ISO 3166-1 alpha-3). Character string.}
#'   \item{country_name }{Names of available countries. Character string.}
#'   \item{geoset_years }{Years of available historical boundary files. Character string.}
#'   \item{space_units }{Available spatial units of analysis. Character string.}
#' }
#' Elections:LowerHouse:CLEA :Classes ‘data.table’ and 'data.frame': 168 obs. of  6 variables
#' Demographics:Ethnicity:EPR :Classes ‘data.table’ and 'data.frame': 180 obs. of  6 variables
#' Demographics:Ethnicity:GREG :Classes ‘data.table’ and 'data.frame': 234 obs. of  6 variables
#' Demographics:Population:GHS :Classes ‘data.table’ and 'data.frame': 257 obs. of  6 variables
#' Events:PoliticalViolence:ABADarfur :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:ACLED :Classes ‘data.table’ and 'data.frame': 100 obs. of  6 variables
#' Events:PoliticalViolence:BeissingerProtest :Classes ‘data.table’ and 'data.frame': 15 obs. of  6 variables
#' Events:PoliticalViolence:BeissingerRiot :Classes ‘data.table’ and 'data.frame': 15 obs. of  6 variables
#' Events:PoliticalViolence:BeissingerUkraine :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:COCACW :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:ESOCAfghanistanWITS :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:ESOCIraqSIGACT :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:ESOCIraqWITS :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:ESOCMexicoDrugRelatedMurders :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:ESOCMexicoHomicide :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:ESOCPakistanBFRS :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:ESOCPakistanWITS :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:GED :Classes ‘data.table’ and 'data.frame': 121 obs. of  6 variables
#' Events:PoliticalViolence:Lankina :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:NIRI :Classes ‘data.table’ and 'data.frame': 12 obs. of  6 variables
#' Events:PoliticalViolence:NVMS :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:PITF :Classes ‘data.table’ and 'data.frame': 133 obs. of  6 variables
#' Events:PoliticalViolence:SCAD :Classes ‘data.table’ and 'data.frame': 60 obs. of  6 variables
#' Events:PoliticalViolence:yzCaucasus2000 :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:yzChechnya :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:yzLibya :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Events:PoliticalViolence:yzUkraine2014 :Classes ‘data.table’ and 'data.frame': 1 obs. of  6 variables
#' Infrastructure:Roads:gRoads :Classes ‘data.table’ and 'data.frame': 240 obs. of  6 variables
#' Infrastructure:NightLights:DMSP :Classes ‘data.table’ and 'data.frame': 257 obs. of  6 variables
#' PublicHealth:Covid19:JHUCSSEC19 :Classes ‘data.table’ and 'data.frame': 207 obs. of  6 variables
#' Terrain:Elevation:ETOPO1 :Classes ‘data.table’ and 'data.frame': 256 obs. of  6 variables
#' Terrain:LandCover:GLCC :Classes ‘data.table’ and 'data.frame': 257 obs. of  6 variables
#' Weather:AirTemperatureAndPrecipitation:NOAA :Classes ‘data.table’ and 'data.frame': 209 obs. of  6 variables
#' \describe{
#'   \item{country_iso3 }{Codes for available countries (ISO 3166-1 alpha-3). Character string.}
#'   \item{country_name }{Names of available countries. Character string.}
#'   \item{year_range }{Range of available years for data topic. Character string.}
#'   \item{time_units }{Available time units. Character string.}
#'   \item{space_units }{Available spatial units. Character string.}
#'   \item{geosets }{Names of available geographic boundary data sources. Character string.}
#' }
#' @source Sub-National Geospatial Data Archive System: Geoprocessing Toolkit (updated March 17, 2023).
"available_data"


#' Country code dictionary
#'
#' Reference table of country names and ISO-3166 codes, adapted from \code{countrycode} package.
#'
#' @format data.table object, with	8626 obs. of  3 variables:
#' \describe{
#'   \item{country_name }{Country names. Character string.}
#'   \item{country_name_alt }{Alternative spellings of country names, ASCII characters only. Character string.}
#'   \item{country_iso3 }{Country codes (ISO 3166-1 alpha-3). Character string.}
#' }
#' @source Vincent Arel-Bundock. Package countrycode: Convert Country Names and Country Code, version 1.40. CRAN (October 12, 2022).
"cc_dict"
