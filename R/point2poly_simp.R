#' Point-to-polygon interpolation, simple overlay method
#'
#' Function for assigning values from a source point layer to a destination polygon layer, using simple point-in-polygon overlays
#'
#' @param pointz Source points layer. \code{sf} object.
#' @param polyz Destination polygon layer. Must have identical CRS to \code{pointz}. \code{sf} object.
#' @param varz Names of variable(s) to be assigned from source polygon layer to destination polygons. Character string or vector of character strings.
#' @param funz Aggregation function to be applied to variables specified in \code{varz}. Must take as an input a vector \code{x}. Function or list of functions.
#' @param na_val Value to be assigned to missing values. Defaul is \code{NA}. Logical or list.
#' @return Returns a \code{sf} polygon object, with variables from \code{pointz} assigned to the geometries of \code{polyz}.
#' @details Assignment procedures are the same for numeric and character string variables. All variables supplied in \code{varz} are passed directly to the function specified in \code{funz}. If different sets of variables are to be aggregated with different functions, both \code{varz} and \code{funz} should be specified as lists (see examples below).
#' @import sf maptools data.table tidyverse
#' @importFrom stats as.dist
#' @importFrom raster extract pointDistance raster
#' @importFrom methods as
#' @importFrom rmapshaper ms_dissolve
#' @examples
#' # Assignment of a single variable (sums)
#' \dontrun{
#' data(hex_05_deu)
#' data(clea_deu2009_pt)
#' out_1 <- point2poly_simp(pointz=clea_deu2009_pt,
#'                          polyz=hex_05_deu,
#'                          varz="vv1")
#' plot(out_1["vv1"])
#' }
#'
#' # Replace NA's with 0's
#' \dontrun{
#' out_2 <- point2poly_simp(pointz = clea_deu2009_pt,
#'                          polyz = hex_05_deu,
#'                          varz = "vv1",
#'                          na_val = 0)
#' plot(out_2["vv1"])
#' }
#'
#' # Multiple variables, with different assignment functions
#' \dontrun{
#' out_3 <- point2poly_simp(pointz = clea_deu2009_pt,
#'                          polyz = hex_05_deu,
#'                          varz = list(
#'                            c("to1","pvs1_margin"),
#'                            c("vv1"),
#'                            c("incumb_pty_n","win1_pty_n")),
#'                          funz = list(
#'                            function(x){mean(x,na.rm=T)},
#'                            function(x){sum(x,na.rm=T)},
#'                            function(x){paste0(unique(na.omit(x)),collapse=" | ") }),
#'                          na_val = list(NA_real_,0,NA_character_))
#' }
#' @export
#'
point2poly_simp <- function(pointz,polyz,varz,funz=function(x){sum(x,na.rm=T)},na_val=NA){

  # Put variables and functions into list
  if(class(varz)=="character"){varz <- list(varz)}
  if(class(funz)=="function"){funz <- list(funz)}
  if(class(na_val)=="logical"){na_val <- list(na_val)}

  # Empty points layer
  pointz_dt <- pointz[1,] %>% as.data.table()
  for(v0 in seq_along(varz)){
    pointz_dt[,c(varz[[v0]]) := lapply(1:length(varz[[v0]]),function(.){na_val[[v0]]})] %>% (function(.){.[,o0 := 1]})
  }
  # Match points to polygons
  o0 <- suppressMessages(
    pointz %>% st_within(polyz) %>% as.data.table()
  )
  # Add polygon index to point layer
  if(nrow(o0)>0){
    pointz_dt <- pointz %>% as.data.table() %>% (function(.){.[o0$row.id,o0 := o0$col.id]}) %>% dplyr::select(-geometry)
  }
  # Aggregate (numeric variables)
  pointz_agg <- lapply(seq_along(varz),function(v0){
    pointz_agg_ <- lapply(seq_along(varz[[v0]]),function(j0){
      int_1_ <- pointz_dt[,list(w = funz[[v0]](get(varz[[v0]][j0]))),by=o0] %>% data.table::setnames("w",paste0(varz[[v0]][j0]))
      if(j0>1){int_1_ <- int_1_ %>% dplyr::select(-o0)}
      int_1_
    }) %>% dplyr::bind_cols()
    pointz_agg_
  }) %>% dplyr::bind_cols()


  # Merge with polygons
  polyz$o0 <- 1:nrow(polyz)
  polyz_ <- merge(polyz,pointz_agg,by="o0",all.x=T,all.y=F, suffixes=c("_","")) %>% (function(.){.[order(.$o0 %>% as.numeric()),]}) %>% dplyr::select(-o0)
  for(v0 in seq_along(varz)){
    polyz_[,varz[[v0]]] <- polyz_ %>% as.data.table() %>% dplyr::select(varz[[v0]]) %>% replace(., is.na(.), na_val[[v0]])
  }
  if(length(grep("^o0",names(polyz_)))>0){polyz_[,grep("^o0",names(polyz_))] <- NULL}

  # Output
  return(polyz_)
}
