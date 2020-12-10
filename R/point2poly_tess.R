#' Point-to-polygon interpolation, tessellation method
#'
#' Function for interpolating values from a source point layer to a destination polygon layer, using Voronoi tessellation and area/population weights.
#'
#' @param pointz Source points layer. \code{sf} object.
#' @param polyz Destination polygon layer. Must have identical CRS to \code{pointz}. \code{sf} object.
#' @param poly_id Name of unique ID column for destination polygon layer. Character string.
#' @param methodz Interpolation method(s) for numeric covariates. Could be either of "aw" (areal weighting, default) and/or "pw" (population weighting). See "details". Character string or vector of character strings.
#' @param char_methodz Interpolation method(s) for character strings. Could be either of "aw" (areal weighting, default) or "pw" (population weighting). See "details". Character string.
#' @param pop_raster Population raster to be used for population weighting, Must be supplied if \code{methodz="pw"}. Must have identical CRS to \code{poly_from}. \code{raster} object.
#' @param varz Names of numeric variable(s) to be interpolated from source polygon layer to destination polygons. Character string or list of character strings.
#' @param funz Aggregation function to be applied to variables specified in \code{varz}. Must take as an input a numeric vector \code{x} and vector of weights \code{w}. Function or list of functions.
#' @param char_varz  Names of character string variables to be interpolated from source polygon layer to destination polygons. Character string or vector of character strings.
#' @param char_assign Assignment rule to be used for variables specified in \code{char_varz}. Could be either "biggest_overlap" (default) or "all_overlap". See "details". Character string or vector of character strings.
#' @param return_tess Return Voronoi polygons, in addition to destinaton polygon layer? Default is \code{FALSE}. Logical.
#' @param seed Seed for generation of random numbers. Default is 1. Numeric.
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
#' When a list of variables are supplied and one methods argument specified, then the chosen method will be applied to all variables.
#'
#' When a list of of variables are supplied and multiple methods arguments specified, then weighting methods will be applied in a pairwise order. For example, specifying \code{varz = list(c("to1","pvs1_margin"), c("vv1"))} and \code{methodz = c('aw', 'pw')} will apply areal weighting to the the first set of variables (to1 and pvs1_margin) and population weighing to the second set (vv1).
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
#' @import sf
#' @importFrom stats as.dist weighted.mean
#' @importFrom data.table data.table rbindlist as.data.table setnames
#' @importFrom raster extract pointDistance raster projectRaster
#' @importFrom methods as
#' @importFrom rmapshaper ms_dissolve
#' @importFrom dplyr select bind_cols left_join
#' @importFrom purrr reduce
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
#'                            function(x,w){stats::weighted.mean(x,w)},
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
  char_methodz = "aw",
  methodz="aw",
  pop_raster=NULL,
  varz=NULL,
  char_varz=NULL,
  char_assign="biggest_overlap",
  funz=function(x,w){stats::weighted.mean(x,w,na.rm=T)},
  return_tess=FALSE,
  seed = 1){

  set.seed(seed)

  ###########################################################
  #Section A - Preparing the Gemetries for Tessalation Prcess
  ###########################################################
  # Put variables and functions into list
  if(class(varz)=="character"){varz <- list(varz)}
  if(class(funz)=="function"){funz <- list(funz)}

  # Stop if no population raster
  if("pw"%in%methodz & length(pop_raster)==0){stop("No population raster provided.")}

  #
  #
  #
  #

  ###########################################################
  #Section B - Preparing the Gemetries for Tessalation Prcess
  ###########################################################
  ############################
  #Part 1 - Create union layer
  ###########################
  suppressWarnings({
    # polyz_u <- dplyr::select(fix_geom(rmapshaper::ms_dissolve(polyz)),-1)
    polyz_u <- fix_geom(sf::st_union(polyz))
  })


  #####################################
  #Part 2 -  Jitter and crop by polygon
  #####################################
  #Part i -
  suppressWarnings({
    suppressMessages(
      pointz_crop <-  sf::st_crop(sf::st_jitter(pointz),sf::st_bbox(polyz_u))
    )
  })

  #Part ii -
  pointz_crop$Unique_ID <- 1:nrow(pointz_crop)

  #
  #

  ###############################
  #Part 3 - Convert to multipoint
  ###############################
  suppressWarnings({
    suppressMessages(
      pointz_geom <- sf::st_union(sf::st_geometry(pointz_crop))
    )
  })

  #
  #

  #################################
  #Part 4 - Create Voronoi polygons
  #################################
  suppressWarnings({
    geo_vor <- suppressMessages(
      sf::st_as_sf(data.table::setnames(data.table::as.data.table(sf::st_intersection(sf::st_as_sf(sf::st_cast(sf::st_voronoi(pointz_geom))),polyz_u)), old="x",new="geometry"))
    )
  })

  #
  #

  ###############################################
  #Part 5 - Exception for single-point geometries
  ###############################################
  if(nrow(pointz_crop)==1){
    suppressWarnings({
      geo_vor <- suppressMessages(
         data.table::setnames(sf::st_intersection(sf::st_as_sf(sf::st_geometry(polyz_u)),polyz_u$geometry),"x","geometry")
      )
    })
  }
  sf::st_geometry(geo_vor) <- "geometry"

  #
  #

  ###############################################
  #Part 6 - Combine with point feature attributes
  ###############################################
  #Part i -
  suppressWarnings({
    suppressMessages(
      int <- sf::st_intersects(geo_vor,pointz_crop)
    )
  })
  Errors <- sapply(int, function(x) length(x))

  #Part ii -
  geo_vor$Unique_ID <- NA
  geo_vor$Unique_ID[Errors == 1] <- pointz_crop$Unique_ID[unlist(int[Errors == 1])]

  #
  #

  #############################
  #Part 7 - Crop the Geometries
  #############################
  suppressWarnings({
    geo_vor <-  suppressMessages(
      fix_geom(sf::st_intersection(geo_vor,polyz_u))
    )
  })

  #
  #

  ################################
  #Part 8 - Merge with Point Layer
  ################################
  pointz_crop_dt <- pointz_crop

  sf::st_geometry(pointz_crop_dt) <- NULL

  suppressWarnings({
    geo_vor <- suppressMessages(
      dplyr::left_join(geo_vor, pointz_crop_dt, by = 'Unique_ID')
    )
  })

  #
  #
  #
  #

  ###########################################################
  #Section C - Preparing the Gemetries for Tessalation Prcess
  ###########################################################
  # Poly-in-poly
  polyz_ <- poly2poly_ap(
    poly_from = geo_vor,
    poly_to = polyz,
    poly_to_id = poly_id,
    geo_vor = geo_vor,
    methodz = methodz,
    pop_raster = pop_raster,
    varz = varz,
    char_methodz = char_methodz,
    char_varz = char_varz,
    char_assign = char_assign,
    funz = funz,
    seed = seed
  )


  # Output
  polyz_$Unique_ID <- NULL; geo_vor$Unique_ID <- NULL
  if(!return_tess){return(polyz_)}
  if(return_tess){return(list(result=polyz_,tess=geo_vor))}

  #Equal Area Projection as a Default

}


