#' Line-in-polygon analysis
#'
#' Function for basic geometry calculations on polyline features, within an overlapping destination polygon layer.
#'
#' @param linez Source polyline layer. \code{sf} object.
#' @param polyz Destination polygon layer. Must have identical CRS to \code{linez}. \code{sf} object.
#' @param poly_id Name of unique ID column for destination polygon layer. Character string.
#' @param measurez Desired measurements. Could be any of "length" (sum of line lengths by polygon), "density" (sum of line lengths divided by area of polygon) and/or "distance" (distance from each polygon to nearest line feature). Default is to report all three. Character string or vector of character strings.
#' @param outvar_name Name (root) to be given to output variable. Default is \code{"line"}. Character string.
#' @param unitz Units of measurement (linear). Defaul is \code{"km"}. Character string.
#' @param reproject Temporarily reproject layers to planar projection for geometric operations? Defaul is \code{TRUE}. Logical.
#' @param na_val Value to be assigned to missing values (line lengths and densities only). Defaul is \code{NA}. Logical or list.
#' @param verbose Print status messages and progress? Default is \code{TRUE}. Logical.
#' @return An \code{sf} polygon object, with summary statisics of \code{linez} features aggregated to the geometries of \code{polyz}.
#'
#' If \code{measurez = "lengths"}, contains fields with suffixes
#' \itemize{
##'  \item{"\code{_length}". }{Sum of line lengths within each polygon, in km or other units supplied in \code{unitz}.}
##'  }
#' If \code{measurez = "density"}, contains fields with suffixes
#' \itemize{
##'  \item{"\code{_length}". }{Sum of line lengths within each polygon, in km or other units supplied in \code{unitz}.}
##'  \item{"\code{_area}". }{Area of each polygon, in km^2 or the square of linear units supplied in \code{unitz}.}
##'  \item{"\code{_density}". }{Sum of line lengths divided by area of each polygon, in km/km^2 or other units supplied in \code{unitz}.}
##'  }
#' If \code{measurez = "distance"}, contains fields with suffixes
#' \itemize{
##'  \item{"\code{_distance}". }{Distance from each polygon to nearest line feature, in km or other units supplied in \code{unitz}.}
##'  }
#' If \code{measurez = c("length","density","distance")} (default), contains all of the above.
#' @import sf data.table tidyverse
#' @importFrom stats as.dist
#' @importFrom raster extract pointDistance raster projectRaster
#' @importFrom methods as
#' @importFrom udunits2 ud.convert
#' @importFrom dplyr select everything
#' @examples
#' # Road lengths, densities and distance from polygon to nearest highway
#' \dontrun{
#' data(hex_05_deu)
#' data(highways_deu1992)
#' out_1 <- line2poly(linez = highways_deu1992,
#'                    polyz = hex_05_deu,
#'                    poly_id = "HEX_ID")
#' plot(out_1["line_length"])
#' plot(out_1["line_density"])
#' plot(out_1["line_distance"])
#' }
#'
#' # Replace missing road lengths and densities with 0's, rename variables
#' \dontrun{
#' out_2 <- line2poly(linez = highways_deu1992,
#'                    polyz = hex_05_deu,
#'                    poly_id = "HEX_ID",
#'                    outvar_name = "road",
#'                    na_val = 0)
#' plot(out_2["road_length"])
#' plot(out_2["road_density"])
#' plot(out_2["road_distance"])
#' }
#' @export

line2poly <- function(linez,
                      polyz,
                      poly_id,
                      measurez=c("length","density","distance"),
                      outvar_name="line",
                      unitz="km",
                      reproject=TRUE,
                      na_val=NA,
                      verbose=TRUE){

  # Find optimal planar projection for map
  if(reproject){
    if(verbose){print("Finding optimal planar projection...")}
    suppressWarnings({polyz_tr <- crs_select(polyz = polyz %>% fix_geom())})
    polyz_tr_ <- polyz_tr[["sf"]]
    epsg <- polyz_tr[["epsg_best"]]
    linez_tr_ <- linez %>% st_transform(st_crs(paste0("EPSG:",epsg)))
  }
  if(!reproject){
    polyz_tr_ <- polyz
    linez_tr_ <- linez
  }

  # Set units
  unitz_sq <- paste0(unitz,"^2")

  # Overlay
  if(verbose){print("Conducting overlay operations...")}
  lp_o <- suppressWarnings({linez_tr_ %>% st_intersection(polyz_tr_)})
  lp_o$sgeoz_length <- st_length(lp_o) %>% udunits2::ud.convert("m", unitz) %>% as.numeric()
  lp_o <- lp_o %>% as.data.table()
  lp_o <- lp_o[,sum(sgeoz_length),by=poly_id] %>% data.table::setnames("V1","sgeoz_length")

  # Distance
  if("distance"%in%measurez){
    polyz_tr_$sgeoz_distance <- polyz_tr_ %>% st_distance(linez_tr_ %>% st_union()) %>% udunits2::ud.convert("m", unitz) %>% as.numeric()
  }

  # Length
  if("length"%in%measurez|"density"%in%measurez){
    polyz_tr_ <- polyz_tr_ %>% merge(lp_o %>% dplyr::select(poly_id %>% tidyselect::all_of(),"sgeoz_length"),by=poly_id,all.x=T,all.y=F)
    polyz_tr_$sgeoz_length <- polyz_tr_$sgeoz_length %>% replace(is.na(.), na_val)
  }

  # Density
  if("density"%in%measurez){
    polyz_tr_$sgeoz_area <- polyz_tr_ %>% st_area() %>% udunits2::ud.convert("m^2", unitz_sq) %>% as.numeric()
    polyz_tr_$sgeoz_density <- polyz_tr_$sgeoz_length/polyz_tr_$sgeoz_area %>% replace(is.na(.), na_val)
  }

  # Rename variables
  polyz_tr_ <- polyz_tr_ %>% dplyr::select(poly_id %>% tidyselect::all_of() ,names(.)[grep("^sgeoz",names(.))],dplyr::everything()) %>% data.table::setnames(names(.),gsub("^sgeoz",outvar_name,names(.)))

  # Reproject
  if(reproject){
    if(verbose){print("Restoring original projection...")}
    polyz_ <- polyz_tr_ %>% st_transform(st_crs(polyz))
  }
  if(!reproject){
    polyz_ <- polyz_tr_
  }

  # Output
  return(polyz_)
}
