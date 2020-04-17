#' Area/population weighted polygon-to-polygon interpolation
#'
#' Function for interpolating values from a source polygon layer to an overlapping (but spatially misaligned) destination polygon layer, using area and/or population weights.
#'
#' @param poly_from Source polygon layer. \code{sf} object.
#' @param poly_to Destination polygon layer. Must have identical CRS to \code{poly_from}. \code{sf} object.
#' @param poly_to_id Name of unique ID column for destination polygon layer. Character string.
#' @param methodz Area interpolation method(s). Could be either of "aw" (areal weighting, default) and/or "pw" (population weighting). See "details". Character string or vector of character strings.
#' @param pop_raster Population raster to be used for population weighting, Must be supplied if \code{methodz="pw"}. Must have identical CRS to \code{poly_from}. \code{raster} object.
#' @param varz Names of numeric variable(s) to be interpolated from source polygon layer to destination polygons. Character string or vector of character strings.
#' @param funz Aggregation function to be applied to variables specified in \code{varz}. Must take as an input a numeric vector \code{x} and vector of weights \code{w}. Function or list of functions.
#' @param char_varz  Names of character string variables to be interpolated from source polygon layer to destination polygons. Character string or vector of character strings.
#' @param char_assign Assignment rule to be used for variables specified in \code{char_varz}. Could be either "biggest_overlap" (default) or "all_overlap". See "details". Character string or vector of character strings.
#' @return \code{sf} polygon object, with variables from \code{poly_from} interpolated to the geometries of \code{poly_to}.
#' @details Currently supported integration methods (\code{methodz}) include:
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
#' @examples
#' # Interpolation of a single variable, with area weights
#' \dontrun{
#' data(clea_deu2009)
#' data(hex_05_deu)
#' out_1 <- poly2poly_ap(poly_from = clea_deu2009,
#'               poly_to = hex_05_deu,
#'               poly_to_id = "HEX_ID",
#'               varz = "to1"
#'              )
#' }
#'
#' # Interpolation of multiple variables, with area weights
#' \dontrun{
#' out_2 <- poly2poly_ap(
#'               poly_from = clea_deu2009,
#'               poly_to = hex_05_deu,
#'               poly_to_id = "HEX_ID",
#'               varz = list(
#'                 c("to1","pvs1_margin"),
#'                 c("vv1") ),
#'               funz = list(
#'                 function(x,w){weighted.mean(x,w)},
#'                 function(x,w){sum(x*w)} ),
#'               char_varz = c("incumb_pty_n","win1_pty_n")
#'              )
#' }
#'
#' # Interpolation of a single variable, with population weights
#' \dontrun{
#' data(gpw4_deu2010)
#' out_3 <- poly2poly_ap(poly_from = clea_deu2009,
#'                          poly_to = hex_05_deu,
#'                          poly_to_id = "HEX_ID",
#'                          varz = "to1",
#'                          methodz = "pw",
#'                          pop_raster = gpw4_deu2010)
#' }
#'
#' # Interpolation of a single variable, with area and population weights
#' \dontrun{
#' out_4 <- poly2poly_ap(poly_from = clea_deu2009,
#'                          poly_to = hex_05_deu,
#'                          poly_to_id = "HEX_ID",
#'                          varz = "to1",
#'                          methodz = c("aw","pw"),
#'                          pop_raster = gpw4_deu2010)
#' }
#' @export


poly2poly_ap <- function(
  poly_from,
  poly_to,
  poly_to_id,
  methodz="aw",
  pop_raster=NULL,
  varz=NULL,
  char_varz=NULL,
  char_assign="biggest_overlap",
  funz=function(x,w){weighted.mean(x,w,na.rm=T)}
){

  # Put variables and functions into list
  if(length(varz)>0){
    if(class(varz)=="character"){varz <- list(varz)}
  }
  if(class(funz)=="function"){funz <- list(funz)}

  # Stop if no population raster
  if("pw"%in%methodz & length(pop_raster)==0){stop("No population raster provided.")}

  # Area-weights (part 1)
  if("aw"%in%methodz){
    # Calculate polygon areas
    poly_from$AREA_TOTAL <- st_area(poly_from) %>% as.numeric()
  }

  # Population-weights (part 1)
  if("pw"%in%methodz){
    # Calculate polygon total
    poly_from$POP_TOTAL <- pop_raster %>% raster::extract(poly_from %>% as("Spatial"),fun=sum,na.rm=T)
  }

  # Intersection
  suppressWarnings({
    int_1 <- suppressMessages(
      st_intersection(poly_from,poly_to) %>% st_buffer(dist=0)
    )
  })
  # # Fix geometry types
  # if(int_1 %>% st_geometry_type() %>% grepl("GEOMETRY",.) %>% sum() > 0){
  #   int_1 <- int_1[!(int_1 %>% st_geometry_type() %>% grepl("GEOMETRY",.)),]
  #   format(object.size(int_1),units="Mb")
  # }

  # Area-weights (part 2)
  if("aw"%in%methodz){
    # Calculate weights
    int_1$AREA_INT <- st_area(int_1) %>% as.numeric()
    int_1$AREA_W <- int_1$AREA_INT/int_1$AREA_TOTAL
  }

  # Population-weights (part 2)
  if("pw"%in%methodz){
    # Calculate weights
    int_1$POP_INT <- pop_raster %>% raster::extract(int_1 %>% as("Spatial"),fun=sum,na.rm=T)
    int_1$POP_W <- int_1$POP_INT/int_1$POP_TOTAL
  }

  # Convert to dt
  int_1_dt <- data.table(int_1)

  # Interpolate missing population values
  if("pw"%in%methodz){
    if(sum(is.na(int_1_dt$POP_INT))==1){
      w <- raster::pointDistance(int_1 %>% st_centroid() %>% st_geometry() %>% unlist() %>% matrix(ncol=2,byrow=T) %>% as.data.frame() %>% data.table::setnames(c("lon","lat")),lonlat=T) %>% as.dist() %>% as.matrix()
      diag(w) <- NA
      int_1_dt[is.na(POP_INT),POP_INT := int_1_dt$POP_INT[t(apply(w, 1, order)[ 1:min(10,nrow(int_1)), is.na(int_1$POP_INT)])] %>% mean(na.rm=T)]
      int_1_dt[is.na(POP_W),POP_W := POP_INT/POP_TOTAL]
    }
    if(sum(is.na(int_1_dt$POP_INT))>1){
      w <- raster::pointDistance(int_1 %>% st_centroid() %>% st_geometry() %>% unlist() %>% matrix(ncol=2,byrow=T) %>% as.data.frame() %>% data.table::setnames(c("lon","lat")),lonlat=T) %>% as.dist() %>% as.matrix()
      diag(w) <- NA
      int_1_dt[is.na(POP_INT),POP_INT := t(apply(w, 1, order)[ 1:min(10,nrow(int_1)), is.na(int_1$POP_INT)]) %>% apply(1,function(.){mean(int_1$POP_INT[.],na.rm=T)})]
      int_1_dt[is.na(POP_W),POP_W := POP_INT/POP_TOTAL]
    }
  }

  # Aggregate
  int_1_w <- int_1_w0 <- int_1_dt %>% dplyr::select(poly_to_id) %>% unique()

  # Character string variables
  if(length(char_varz)>0){
    if("biggest_overlap"%in%char_assign){
      int_1_w <- lapply(seq_along(char_varz),function(j0){
        if("aw"%in%methodz){
          int_1_w <- int_1_w0 %>% merge(int_1_dt[,list(w = get(char_varz[j0])[which.max(AREA_W)]),by=poly_to_id] %>% data.table::setnames("w",paste0(char_varz[j0],"_aw")),by=poly_to_id,suffixes=c("","_sngz"))
        }
        if("pw"%in%methodz){
          int_1_w <- int_1_w %>% merge(int_1_dt[,list(w = get(char_varz[j0])[which.max(POP_W)]),by=poly_to_id] %>% data.table::setnames("w",paste0(char_varz[j0],"_pw")),by=poly_to_id,suffixes=c("","_sngz"))
        }
        if(j0>1){int_1_w <- int_1_w %>% dplyr::select(-poly_to_id)}
        int_1_w
      }) %>% dplyr::bind_cols()
    }
    int_1_w1 <- int_1_w
    if("all_overlap"%in%char_assign){
      int_1_w <- lapply(seq_along(char_varz),function(j0){
        if("aw"%in%methodz){
          int_1_w <- int_1_w0 %>% merge(int_1_dt[order(AREA_W) %>% rev(),list(w = get(char_varz[j0]) %>% unlist() %>% unique() %>% paste0(collapse=" | ")),by=poly_to_id] %>% data.table::setnames("w",paste0(char_varz[j0],"_all_aw")),by=poly_to_id,suffixes=c("","_sngz"))
        }
        if("pw"%in%methodz){
          int_1_w <- int_1_w %>% merge(int_1_dt[order(POP_W) %>% rev(),list(w = get(char_varz[j0]) %>% unlist() %>% unique() %>% paste0(collapse=" | ")),by=poly_to_id] %>% data.table::setnames("w",paste0(char_varz[j0],"_all_pw")),by=poly_to_id,suffixes=c("","_sngz"))
        }
        if(j0>1){int_1_w <- int_1_w %>% dplyr::select(-poly_to_id)}
        int_1_w
      }) %>% dplyr::bind_cols()
    }
    int_1_w <- int_1_w1 %>% merge(int_1_w,by=poly_to_id,suffixes=c("","_sngz"))
  }

  # Numeric variables
  if(length(varz)>0){
    int_1_w <- int_1_w %>% merge(
      {
        lapply(seq_along(varz),function(v0){
          valz_agg_2 <- lapply(seq_along(varz[[v0]]),function(j0){
            if("aw"%in%methodz){
              int_1_w0 <- int_1_w0 %>% merge(int_1_dt[,list(w = funz[[v0]](x=get(varz[[v0]][j0]),w=AREA_W)),by=poly_to_id] %>% data.table::setnames("w",paste0(varz[[v0]][j0],"_aw")),by=poly_to_id,suffixes=c("","_sngz"))
            }
            if("pw"%in%methodz){
              int_1_w0 <- int_1_w0 %>% merge(int_1_dt[,list(w = funz[[v0]](get(varz[[v0]][j0]),w=POP_W)),by=poly_to_id] %>% data.table::setnames("w",paste0(varz[[v0]][j0],"_pw")),by=poly_to_id,suffixes=c("","_sngz"))
            }
            if(j0>1){int_1_w0 <- int_1_w0 %>% dplyr::select(-poly_to_id)}
            int_1_w0
          }) %>% dplyr::bind_cols()
          if(v0>1){valz_agg_2 <- valz_agg_2 %>% dplyr::select(-poly_to_id)}
          valz_agg_2
        }) %>% dplyr::bind_cols()
      },by=poly_to_id,suffixes=c("","_sngz"))
  }
  int_1_w

  # Merge with sf
  polyz_ <- merge(poly_to %>% as.data.table(),int_1_w %>% as.data.table(),by=poly_to_id,suffixes=c("","_sngz")) %>% dplyr::select(-grep("_sngz$",names(.))) %>% st_as_sf()


  # Output
  return(polyz_)

}
