#' Point-to-polygon interpolation, tessellation method
#'
#' Function for interpolating values from a source point layer to a destination polygon layer, using Voronoi tessellation and area/population weights.
#'
#' @param pointz Source points layer. \code{sf} object.
#' @param polyz Destination polygon layer. Must have identical CRS to \code{pointz}. \code{sf} object.
#' @param poly_id Name of unique ID column for destination polygon layer. Character string.
#' @param methodz Area interpolation method(s). Could be either of "aw" (areal weighting, default) and/or "pw" (population weighting). See "details". Character string or vector of character strings.
#' @param pop_raster Population raster to be used for population weighting, Must be supplied if \code{methodz="pw"}. Must have identical CRS to \code{poly_from}. \code{raster} object.
#' @param varz Names of numeric variable(s) to be interpolated from source polygon layer to destination polygons. Character string or list of character strings.
#' @param funz Aggregation function to be applied to variables specified in \code{varz}. Must take as an input a numeric vector \code{x} and vector of weights \code{w}. Function or list of functions.
#' @param char_varz  Names of character string variables to be interpolated from source polygon layer to destination polygons. Character string or vector of character strings.
#' @param char_assign Assignment rule to be used for variables specified in \code{char_varz}. Could be either "biggest_overlap" (default) or "all_overlap". See "details". Character string or vector of character strings.
#' @param return_tess Return Voronoi polygons, in addition to destinaton polygon layer? Defaul is \code{FALSE}. Logical.
#' @return If \code{return_tess=FALSE}, returns a \code{sf} polygon object, with variables from \code{pointz} interpolated to the geometries of \code{polyz}.
#'
#' If \code{return_tess=TRUE}, returns a list, containing
#' \itemize{
##'  \item{"result". }{The destination polygon layer. \code{sf} object.}
##'  \item{"tess". }{The (intermediate) Voronoi tessellation polygon layer. \code{sf} object.}
##'  }
#' @details
#' This function interpolates point data to polygons with a two-step process. In the first step (tessellation), each point is assigned a Voronoi cell, drawn such that (a) the distance from its borders to the focal point is less than or equal to the distance to any other point, and (b) no gaps between cells remain. The second step (interpolation) performs a polygon-in-polygon interpolation, using the Voronoi cells as source polygons.
#'
#' Currently supported integration methods in the second step (\code{methodz}) include:
#' \itemize{
##'  \item{Areal weighting ("aw"). }{Values from \code{poly_from} weighted in proportion to relative area of spatial overlap between source features and geometries of \code{poly_to}.}
##'  \item{Population weighting ("pw"). }{Values from \code{poly_from} weighted in proportion to relative population sizes in areas of spatial overlap between source features and geometries of \code{poly_to}. This routine uses a third layer (supplied in \code{pop_raster}) to calculate the weights.}
##' }
#' It is possible to pass multiple arguments to \code{methodz} (e.g. \code{methodz=c("aw","pw")}), in which case the function will calculate both sets of weights, and append the resulting columns to the output.
#'
#' Interpolation procedures are handled somewhat differently for numeric and character string variables. For numeric variables supplied in \code{varz}, "aw" and/or "pw" weights are passed to the function specified in \code{funz}. If different sets of numeric variables are to be aggregated with different functions, both \code{varz} and \code{funz} should be specified as lists (see examples below).
#'
#' For character string (and any other) variables supplied in \code{char_varz}, "aw" and/or "pw" weights are passed to the assignment rule(s) specified in \code{char_assign}. Note that the \code{char_varz} argument may include numerical variables, but \code{varz} cannot include character string variables.
#'
#' Currently supported assignment rules for character strings (\code{char_assign}) include:
#' \itemize{
##'  \item{"biggest_overlap". }{For each variable in \code{char_varz}, the features in \code{poly_to} are assigned a single value from overlapping \code{poly_from} features, corresponding to the intersection with largest area and/or population weight.}
##'  \item{"all_overlap". }{For each variable in \code{char_varz}, the features in \code{poly_to} are assigned all values from overlapping \code{poly_from} features, ranked by area and/or population weights (largest-to-smallest) of intersections.}
##' }
#' It is possible to pass multiple arguments to \code{char_assign} (e.g. \code{char_assign=c("biggest_overlap","all_overlap")}), in which case the function will calculate both, and append the resulting columns to the output.
#' @import sf maptools data.table tidyverse
#' @importFrom stats as.dist
#' @importFrom raster extract pointDistance raster projectRaster
#' @importFrom methods as
#' @importFrom rmapshaper ms_dissolve
#' @examples
#' # Interpolation of a single variable, with area weights
#' \dontrun{
#' data(hex_05_deu)
#' data(clea_deu2009_pt)
#' out_1 <- point2poly_tess(pointz = clea_deu2009_pt,
#'                              polyz = hex_05_deu,
#'                              poly_id = "HEX_ID",
#'                              varz = "to1")
#' plot(out_1["to1_aw"])
#' }
#'
#' # Extract and inspect tessellation polygons
#' \dontrun{
#' out_2 <- point2poly_tess(pointz = clea_deu2009_pt,
#'                              polyz = hex_05_deu,
#'                              poly_id = "HEX_ID",
#'                              varz = "to1",
#'                              return_tess = TRUE)
#' plot(out_2$tess["to1"])
#' plot(out_2$result["to1_aw"])
#' }
#'
#' # Interpolation of multiple variables, with area and population weights
#' \dontrun{
#' data(gpw4_deu2010)
#' out_3 <- point2poly_tess(pointz = clea_deu2009_pt,
#'                          polyz = hex_05_deu,
#'                          poly_id = "HEX_ID",
#'                          methodz = c("aw","pw"),
#'                          varz = list(
#'                            c("to1","pvs1_margin"),
#'                            c("vv1")
#'                          ),
#'                          funz = list(
#'                            function(x,w){weighted.mean(x,w)},
#'                            function(x,w){sum(x*w)}
#'                            ),
#'                          char_varz = c("incumb_pty_n","win1_pty_n"),
#'                          pop_raster = gpw4_deu2010)
#' plot(out_3["pvs1_margin_pw"])
#' }
#' @export

point2poly_tess <- function(
  pointz,
  polyz,
  poly_id,
  methodz="aw",
  pop_raster=NULL,
  varz=NULL,
  char_varz=NULL,
  char_assign="biggest_overlap",
  funz=function(x,w){weighted.mean(x,w,na.rm=T)},
  return_tess=FALSE){

  # Put variables and functions into list
  if(class(varz)=="character"){varz <- list(varz)}
  if(class(funz)=="function"){funz <- list(funz)}

  # Stop if no population raster
  if("pw"%in%methodz & length(pop_raster)==0){stop("No population raster provided.")}

  # Create union layer
  suppressWarnings({
    polyz_u <- polyz %>% rmapshaper::ms_dissolve() %>% fix_geom(self_int = FALSE) %>% dplyr::select(-1)
  })

  # Jitter and crop by polygon
  suppressWarnings({
    pointz_crop <-  suppressMessages(
      pointz %>% st_jitter() %>% st_crop(st_bbox(polyz_u))
    )
  })

  # Convert to multipoint
  pointz_geom <- pointz_crop  %>% st_geometry() %>% st_union()

  # # Buffer layer
  # buffsize <- 1
  # polyz_buff <- suppressMessages(
  #   polyz_u %>% st_buffer(dist=buffsize)
  # )

  # Create Voronoi polygons
  suppressWarnings({
    geo_vor <- suppressMessages(
      st_voronoi(pointz_geom) %>% st_cast() %>% st_as_sf() %>% st_intersection(.,polyz_u) %>% as.data.table() %>% setnames(old="x",new="geometry") %>% st_as_sf()
    )
  })
  # Exception for single-point geometries
  if(nrow(pointz_crop)==1){
    geo_vor <- suppressMessages(
      polyz_u %>% st_geometry() %>% st_as_sf() %>% st_intersection(polyz_u$geometry) %>% setnames("x","geometry")
    )
  }
  st_geometry(geo_vor) <- "geometry"

  # Combine with point feature attributes
  int <- suppressMessages(
    st_intersects(geo_vor,pointz_crop) %>% as.data.table()
  )
  geo_vor <- dplyr::bind_cols(geo_vor  %>% dplyr::slice(int$row.id),pointz_crop %>% as.data.table() %>% dplyr::select(-geometry) %>% dplyr::slice(int$col.id))

  # Crop
  suppressWarnings({
    geo_vor <-  suppressMessages(
      st_intersection(geo_vor,polyz_u) %>% fix_geom(self_int=FALSE)
    )
  })

  # Poly-in-poly
  polyz_ <- poly2poly_ap(
    poly_from = geo_vor,
    poly_to = polyz,
    poly_to_id = poly_id,
    methodz = methodz,
    pop_raster = pop_raster,
    varz = varz,
    char_varz = char_varz,
    char_assign = char_assign,
    funz = funz
  )

  # Output
  if(!return_tess){return(polyz_)}
  if(return_tess){return(list(result=polyz_,tess=geo_vor))}

}
