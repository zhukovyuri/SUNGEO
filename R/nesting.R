#' Relative scale and nesting coefficients
#'
#' Function to calculate relative scale and nesting metrics for changes of support from a source polygon layer to an overlapping (but spatially misaligned) destination polygon layer.
#'
#' @param poly_from Source polygon layer. \code{sf} object (polygon or multipolygon).
#' @param poly_to Destination polygon layer. Must have identical CRS to \code{poly_from}. \code{sf} object (polygon or multipolygon).
#' @param metrix Requested scaling and nesting metrics. See "details". Default is "all". Character string or vector of character strings.
#' @param tol_ Minimum area of polygon intersection, in square meters. Default is 0.001. Numeric.
#' @param by_unit Include a by-unit decomposition of requested nesting metrics (if available)? Default is FALSE. Logical.
#' @return Named list, with numeric values for each requested metric in \code{metrix}. If \code{by_unit==TRUE}, last element of list is a data.table, with nesting metrics disaggregated by source unit, where the first column is a row index for the source polygon layer.
#' @details Currently supported metrics (\code{metrix}) include:
#' \itemize{
##'  \item Relative scale ("rs"). Measures whether a change-of-support (CoS) task is one of aggregation or disaggregation, by calculating the share of source units that are smaller than destination units. Its range is from 0 to 1, where values of 1 indicate pure aggregation (all source units are smaller than destination units) and values of 0 indicate no aggregation (all source units are at least as large as destination units). Values between 0 and 1 indicate a hybrid (i.e. some source units are smaller, others are larger than target units).
##'  \item Relative nesting ("rn"). Measures how closely source and destination boundaries align, by calculating the share of source units that cannot be split across multiple destination units. Its range is from 0 to 1, where values of 0 indicate no nesting (every source unit can be split across multiple destination units) and values of 1 indicate full nesting (no source unit can be split across multiple destination units).
##'  \item Relative scale, symmetric ("rs_sym"). Alternative measure of "rs", which ranges from -1 to 1. It calculates a difference between two proportions: the share of source units that is smaller than destination units (i.e. "rs" from standpoint of source units), and the share that is larger (i.e. "rs" from standpoint of destination units). Values of -1 indicate pure disaggregation (all source units are larger than destination units), 1 indicates pure aggregation (all source units are smaller than destination units). Values of 0 indicate that all source units are the same size as target units.
##'  \item Relative nesting, symmetric ("rn_sym"). Alternative measure of "rn", which ranges from -1 to 1. It calculates a difference between two components: the nesting of source units within destination units (i.e. "rn" from standpoint of source units), and the nesting of destination units within source units (i.e. "rn" from standpoint of destination units. Values of 1 indicate that source units are perfectly nested within destination units; -1 indicates that destination units are perfectly nested within source units.
##'  \item Relative scale, alternative ("rs_alt"). Alternative measure of "rs", rescaled as a proportion of destination unit area. This measure can take any value on the real line, with positive values indicating aggregation and negative values indicating disaggregation.
##'  \item Relative nesting, alternative ("rn_alt"). Alternative measure of "rn", which places more weight on areas of maximum overlap. The main difference between this measure and "rn" is its use of the maximum intersection area for each source polygon instead of averaging over the quadratic term. Two sets of polygons are considered nested if one set is completely contained within another, with as few splits as possible. If none or only a sliver of a source polygon area falls outside a single destination polygon, those polygons are "more nested" than a case where half of a source polygon falls in destination polygon A and half falls into another polygon B.
##'  \item Relative scale, conditional ("rs_nn"). Alternative measure of "rs", calculated for the subset of source units that are not fully nested within destination units.
##'  \item Relative nesting, conditional ("rn_nn"). Alternative measure of "rn", calculated for the subset of source units that are not fully nested within destination units.
##'  \item Proportion intact ("p_intact"). A nesting metric that requires no area calculations at all. This measure ranges from 0 to 1, where 1 indicates full nesting (i.e. every source unit is intact/no splits), and 0 indicates no nesting (i.e. no source unit is intact/all are split).
##'  \item Proportion fully nested ("full_nest"). A stricter version of "p_intact". This measure ranges from 0 to 1, where 1 indicates full nesting (i.e. every source unit is intact/no splits AND falls completely inside the destination layer), and 0 indicates no nesting (i.e. no source unit is both intact and falls inside destination layer).
##'  \item Relative overlap ("ro"). Assesses extent of spatial overlap between source and destination polygons. This measure is scaled between -1 and 1. Values of 0 indicate perfect overlap (there is no part of source units that fall outside of destination units, and vice versa). Values between 0 and 1 indicate a "source underlap" (some parts of source polygons fall outside of destination polygons; more precisely, a larger part of source polygon area falls outside destination polygons than the other way around). Values between -1 and 0 indicate a "destination underlap" (some parts of destination polygons fall outside of source polygons; a larger part of destination polygon area falls outside source polygons than the other way around). Values of -1 and 1 indicate no overlap (all source units fall outside destination units, and vice versa). This is a theoretical limit only; the function returns an error if there is no overlap.
##'  \item Gibbs-Martin index of diversification ("gmi"). Inverse of "rn", where values of 1 indicate that every source unit is evenly split across multiple destination units, and 0 indicates that no source unit is split across any destination units.
#' }
#' It is possible to pass multiple arguments to \code{metrix} (e.g. \code{metrix=c("rn","rs")}). The default (\code{metrix="all"}) returns all of the above metrics.
#'
#' The function automatically reprojects source and destination geometries to Lambert Equal Area prior to calculation, with map units in meters.
#'
#' Values of \code{tol_} can be adjusted to increase or decrease the sensitivity of these metrics to small border misalignments. The default value discards polygon intersections smaller than 0.001 square meters in area.
#' @importFrom sf sf_use_s2 st_geometry st_crs st_transform st_coordinates st_as_sfc st_intersection st_bbox st_area st_buffer st_difference st_union st_combine
#' @importFrom data.table data.table as.data.table
#' @importFrom stats median
#' @examples
#' # Calculate all scale and nesting metrics for two sets of polygons
#' \dontrun{
#' data(clea_deu2009)
#' data(hex_05_deu)
#' nest_1 <- nesting(
#'               poly_from = clea_deu2009,
#'               poly_to = hex_05_deu
#'               )
#' nest_1
#' }
#'
#' # Calculate just Relative Nesting, in the opposite direction
#' \dontrun{
#' nest_2 <- nesting(
#'               poly_from = hex_05_deu,
#'               poly_to = clea_deu2009,
#'               metrix = "rn"
#'               )
#' nest_2
#' }
#' @export

nesting <- function(
  poly_from = NULL,
  poly_to = NULL,
  metrix = "all",
  tol_ = 0.001,
  by_unit=FALSE){

  # Null out undefined variables
  .N <- ID_1_ <- area_in_ <- V1 <- index <- area_1_ <- area_ix_ <-NULL

  # Turn off s2 processing
  suppressMessages({
    sf::sf_use_s2(FALSE)
  })

  # Check input geometries
  if(is.null(poly_from)|is.null(poly_to)){
    base::stop("poly_from and poly_to geographies must both be sfc_MULTIPOLYGON or sfc_POLYGON.")
  }
  if(!("sfc_MULTIPOLYGON"%in%base::class(sf::st_geometry(poly_from))|"sfc_POLYGON"%in%base::class(sf::st_geometry(poly_from)))|!("sfc_MULTIPOLYGON"%in%base::class(sf::st_geometry(poly_to))|"sfc_POLYGON"%in%base::class(sf::st_geometry(poly_to)))){
    base::stop("poly_from and poly_to geographies must both be sfc_MULTIPOLYGON or sfc_POLYGON.")
  }

  # Select scale and nesting metrics
  if("all"%in%metrix){
    metrix <- c("rs","rn","rs_sym","rn_sym","rs_alt","rn_alt","rs_nn","rn_nn","p_intact","full_nest","ro","gmi")
  }
  metrix <- metrix[metrix%in%c("rs","rn","rs_sym","rn_sym","rs_alt","rn_alt","rs_nn","rn_nn","p_intact","full_nest","ro","gmi")]

  # Ensure layers have same CRS
  if(sf::st_crs(poly_from)!=sf::st_crs(poly_to)){
    poly_from <- sf::st_transform(poly_from,crs=sf::st_crs(poly_to))
  }

  # Reproject to equal area
  if(!grepl("4326",sf::st_crs(poly_from)["input"]$input)){
    poly_from <- sf::st_transform(poly_from,crs=4326)
    poly_to <- sf::st_transform(poly_to,crs=4326)
  }
  xy0 <- data.table::as.data.table(sf::st_coordinates(poly_to))[,lapply(.(get("X"),get("Y")),stats::median)]
  xy0 <- data.table::setnames(xy0, c("X","Y"))
  poly_from <- sf::st_transform(poly_from, crs = paste0("+proj=laea +lon_0=",xy0[,get("X")]," +lat_0=",xy0[,get("Y")]))
  poly_to <- sf::st_transform(poly_to, crs = paste0("+proj=laea +lon_0=",xy0[,get("X")]," +lat_0=",xy0[,get("Y")]))

  # Stop if no overlap
  if(length(sf::st_intersection(sf::st_as_sfc(sf::st_bbox(poly_from)),sf::st_as_sfc(sf::st_bbox(poly_to))))==0){
    base::stop("Zero overlap between polygons.")
  } else {

    # Assign temporary IDs
    poly_from$ID_1_ <- 1:nrow(poly_from)
    poly_to$ID_2_ <- 1:nrow(poly_to)

    # Calculate area of original units
    poly_from$area_1_ <- sf::st_area(poly_from)
    poly_to$area_2_ <- sf::st_area(poly_to)

    # Intersect
    suppressWarnings({
      o_sf <- sf::st_intersection(poly_from,poly_to)[,c("ID_1_","ID_2_","area_1_","area_2_")]
    })
    o_sf$area_ix_ <- sf::st_area(o_sf)
    o_sf$area_in_ <- o_sf$area_ix_/o_sf$area_1_

    # Create unit matrix
    if(by_unit==TRUE){
      nesting_dt <- data.table::data.table(index = 1:nrow(poly_from))
    }

    # Relative scale
    if("rs"%in%metrix){
      rs <- data.table::data.table(o_sf)[as.numeric(get("area_in_"))>tol_,base::mean(1*(get("area_1_")<get("area_2_")),na.rm=TRUE)]
      if(by_unit==TRUE){
        rs_byun <- data.table::data.table(o_sf)[as.numeric(get("area_in_"))>tol_,base::mean(1*(get("area_1_")<get("area_2_")),na.rm=TRUE),by="ID_1_"]
        nesting_dt[,rs := rs_byun[match(nesting_dt[,index],rs_byun[,ID_1_]),V1]]
      }
    }

    # Relative nesting
    if("rn"%in%metrix|"rn_sym"%in%metrix){
      rn <- data.table::as.data.table(o_sf)[as.numeric(get("area_in_"))>tol_,base::sum(as.numeric(get("area_ix_")/get("area_1_"))^2),by="ID_1_"][,base::mean(get("V1"),na.rm=TRUE)]
      if(by_unit==TRUE){
        rn_byun <- data.table::as.data.table(o_sf)[as.numeric(get("area_in_"))>tol_,base::sum(as.numeric(get("area_ix_")/get("area_1_"))^2),by="ID_1_"]
        nesting_dt[,rn := rn_byun[match(nesting_dt[,index],rn_byun[,ID_1_]),V1]]
      }
    }

    # Relative scale, symmetric
    if("rs_sym"%in%metrix){
      rs_sym <- data.table::data.table(o_sf)[as.numeric(get("area_in_"))>tol_,base::mean(1*(get("area_1_")<get("area_2_"))-1*(get("area_1_")>get("area_2_")),na.rm=TRUE)]
    }

    # Relative nesting, symmetric
    if("rn_sym"%in%metrix){
      suppressWarnings({
        o_sf_2 <- sf::st_intersection(poly_to,poly_from)[,c("ID_1_","ID_2_","area_1_","area_2_")]
      })
      o_sf_2$area_ix_ <- sf::st_area(o_sf_2)
      o_sf_2$area_in_ <- o_sf_2$area_ix_/o_sf_2$area_2_
      rn_sym <- data.table::as.data.table(o_sf_2)[as.numeric(get("area_in_"))>tol_,rn - base::sum(as.numeric(get("area_ix_")/get("area_2_"))^2),by="ID_2_"][,base::mean(get("V1"))]
    }

    # Relative scale, alternative
    if("rs_alt"%in%metrix){
      rs_alt <- data.table::data.table(o_sf)[as.numeric(get("area_in_"))>tol_,as.numeric(base::mean((get("area_2_")-get("area_1_"))/base::mean(get("area_2_")),na.rm=TRUE))]
      if(by_unit==TRUE){
        rs_alt_byun <- data.table::data.table(o_sf)[as.numeric(get("area_in_"))>tol_,as.numeric(base::mean((get("area_2_")-get("area_1_"))/base::mean(get("area_2_")),na.rm=TRUE)),by="ID_1_"]
        nesting_dt[,rs_alt := rs_alt_byun[match(nesting_dt[,index],rs_alt_byun[,ID_1_]),V1]]
      }
    }

    # Relative nesting, alternative
    if("rn_alt"%in%metrix){
      rn_alt <- data.table::as.data.table(o_sf)[as.numeric(get("area_in_"))>tol_,as.numeric(base::max(area_in_)),by="ID_1_"][,mean(V1)]
      if(by_unit==TRUE){
        rn_alt_byun <- data.table::as.data.table(o_sf)[as.numeric(get("area_in_"))>tol_,as.numeric(base::max(area_in_)),by="ID_1_"]
        nesting_dt[,rn_alt := rn_alt_byun[match(nesting_dt[,index],rn_alt_byun[,ID_1_]),V1]]
      }
    }

    # Relative scale, conditional
    if("rs_nn"%in%metrix){
      full_nest_byun <- data.table::data.table(o_sf)[as.numeric(get("area_in_"))>tol_,1*(as.numeric(sum(area_ix_)/mean(area_1_))==1&.N==1),by=ID_1_]
      nn <- full_nest_byun[data.table::data.table(o_sf)[,match(ID_1_,full_nest_byun[,ID_1_])],V1]!=1
      rs_nn <- data.table::data.table(o_sf)[nn&as.numeric(get("area_in_"))>tol_,base::mean(1*(get("area_1_")<get("area_2_")),na.rm=TRUE)]
      if(by_unit==TRUE){
        rs_byun <- data.table::data.table(o_sf)[as.numeric(get("area_in_"))>tol_,base::mean(1*(get("area_1_")<get("area_2_")),na.rm=TRUE),by="ID_1_"]
        nesting_dt[,rs_nn := rs_byun[match(nesting_dt[,index],rs_byun[,ID_1_]),V1]]
        nesting_dt[index%in%full_nest_byun[V1==1,ID_1_], rs_nn := NA_real_]
      }
    }

    # Relative nesting, conditional
    if("rn_nn"%in%metrix){
      full_nest_byun <- data.table::data.table(o_sf)[as.numeric(get("area_in_"))>tol_,1*(as.numeric(sum(area_ix_)/mean(area_1_))==1&.N==1),by=ID_1_]
      nn <- full_nest_byun[data.table::data.table(o_sf)[,match(ID_1_,full_nest_byun[,ID_1_])],V1]!=1
      rn_nn <- data.table::as.data.table(o_sf)[nn&as.numeric(get("area_in_"))>tol_,base::sum(as.numeric(get("area_ix_")/get("area_1_"))^2),by="ID_1_"][,base::mean(get("V1"),na.rm=TRUE)]
      if(by_unit==TRUE){
        rn_byun <- data.table::as.data.table(o_sf)[as.numeric(get("area_in_"))>tol_,base::sum(as.numeric(get("area_ix_")/get("area_1_"))^2),by="ID_1_"]
        nesting_dt[,rn_nn := rn_byun[match(nesting_dt[,index],rn_byun[,ID_1_]),V1]]
        nesting_dt[index%in%full_nest_byun[V1==1,ID_1_], rn_nn := NA_real_]
      }
    }

    # Proportion intact
    if("p_intact"%in%metrix){
      p_intact <- data.table::as.data.table(o_sf)[as.numeric(get("area_in_"))>tol_,base::mean(base::table(get("ID_1_"))==1)]
      if(by_unit==TRUE){
        p_intact_byun <- data.table::as.data.table(o_sf)[as.numeric(get("area_in_"))>tol_,1*(.N==1),by="ID_1_"]
        nesting_dt[,p_intact := p_intact_byun[match(nesting_dt[,index],p_intact_byun[,ID_1_]),V1]]
      }
    }

    # Proportion fully nested
    if("full_nest"%in%metrix){
        full_nest_byun <- data.table::data.table(o_sf)[as.numeric(get("area_in_"))>tol_,1*(as.numeric(sum(area_ix_)/mean(area_1_))==1&.N==1),by=ID_1_]
        full_nest <- full_nest_byun[,mean(V1)]
      if(by_unit==TRUE){
        nesting_dt[,full_nest := full_nest_byun[match(nesting_dt[,index],full_nest_byun[,ID_1_]),V1]]
      }
    }

    # Relative overlap
    if("ro"%in%metrix){
      suppressWarnings({
        suppressMessages({
          area_all_1 <- sum(sf::st_area(sf::st_buffer(poly_from,0)))
          area_under_1 <- sum(sf::st_area(sf::st_difference(sf::st_buffer(poly_from,0), sf::st_buffer(sf::st_union(sf::st_combine(sf::st_buffer(poly_to,0))),0))))
          area_all_2 <- sum(sf::st_area(sf::st_buffer(poly_to,0)))
          area_under_2 <- sum(sf::st_area(sf::st_difference(sf::st_buffer(poly_to,0), sf::st_buffer(sf::st_union(sf::st_combine(sf::st_buffer(poly_from,0))),0))))
          ro <- (1-as.numeric(area_under_2/area_all_2))-(1-as.numeric(area_under_1/area_all_1))
        })
      })
    }

    # Gibbs-Martin index
    if("gmi"%in%metrix){
      gmi <- data.table::as.data.table(o_sf)[as.numeric(get("area_in_"))>tol_,1-sum(as.numeric(get("area_ix_")/get("area_1_"))^2),by="ID_1_"][,base::mean(get("V1"))]
      if(by_unit==TRUE){
        gmi_byun <- data.table::as.data.table(o_sf)[as.numeric(get("area_in_"))>tol_,1-sum(as.numeric(get("area_ix_")/get("area_1_"))^2),by="ID_1_"][,base::mean(get("V1")),by="ID_1_"]
        nesting_dt[,gmi := gmi_byun[match(nesting_dt[,index],gmi_byun[,ID_1_]),V1]]
      }
    }

    # Output
    out_list <- lapply(metrix,function(m0){
      get(m0)
    })
    names(out_list) <- metrix
    if(by_unit==TRUE){
      out_list <- append(out_list,list(nesting_dt))
      names(out_list)[length(out_list)] <- "by_unit"
    }
    return(out_list)

  }
}
