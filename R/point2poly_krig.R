#' Point-to-polygon interpolation, Ordinary Kriging
#'
#' Function for interpolating values from a source points layer to an overlapping destination polygon layer, using Ordinary Kriging with automatic variogram fitting
#'
#' @param pointz Source points layer. \code{sf} object.
#' @param polyz Destination polygon layer. Must have identical CRS to \code{pointz}. \code{sf} object.
#' @param polyz2 Optional polygon layer defining study region, used for identification of optimal planar projection. Must have identical CRS to \code{pointz}. \code{sf} object.
#' @param varz Names of numeric variable(s) to be interpolated from source points layer to destination polygons. Character string or vector of character strings.
#' @param cellsize Target cell size to be used in Kriging estimation. Default is 25 km. Numeric.
#' @param messagez Optional message to be printed during Kriging estimation. Character string.
#' @return \code{sf} polygon object, with variables from \code{pointz} interpolated to the geometries of \code{polyz}. For each variable "x" in \code{varz}, two estimates are generated, with suffixes
#' \itemize{
##'  \item{"_kr". }{Mean of predicted values within each destination polygon.}
##'  \item{"_kr_sd". }{Average standard deviation of predicted values within each destination polygon.}
##' }
#' @details This function performs Ordinary Kriging, autmatically selecting a variogram model with the smallest residual sum of squares with the sample variogram.
#'
#' Unlike other available point-to-polygon interpolation techniques, this function currently only accepts numeric variables in \code{varz} and does not support interpolation of character strings.
#' @export
#' @import sf maptools data.table tidyverse automap
#' @importFrom stats as.dist
#' @importFrom raster extract pointDistance raster projectRaster
#' @importFrom methods as
#' @importFrom rmapshaper ms_dissolve
#' @importFrom sp SpatialPixelsDataFrame coordinates proj4string
#' @examples
#' # Ordinary Kriging with one variable
#' \dontrun{
#' data(hex_05_deu)
#' data(clea_deu2009_pt)
#' out_1 <- point2poly_krig(pointz = clea_deu2009_pt,
#'                          polyz = hex_05_deu,
#'                          varz = "to1")
#' plot(out_1["to1_kr"])
#' plot(out_1["to1_kr_sd"])
#' }
#'
#' # Ordinary Kriging with multiple variable
#' \dontrun{
#' out_2 <- point2poly_krig(pointz = clea_deu2009_pt,
#'                          polyz = hex_05_deu,
#'                          varz = c("to1","pvs1_margin"))
#' }

point2poly_krig <- function(pointz,polyz,polyz2=NULL,varz=NULL,cellsize=25000,messagez=""){

  if(length(varz)==0){stop("Please supply at least one variable in varz.")}

  # Union layer
  polyz_u <- polyz %>% ms_dissolve()
  if(length(polyz2)==0){polyz2 <- polyz_u}

  # Find optimal planar projection for map
  pointz_tr <- geo2planar(
    polyz = polyz2 %>% fix_geom(),
    sf_layer = pointz %>% st_transform(st_crs(polyz)) %>% fix_geom()
  )
  pointz_tr_ <- pointz_tr[["sf"]]
  epsg <- pointz_tr[["epsg_best"]]

  # Create regular grid
  k_grid <- st_make_grid(polyz_u %>% st_transform(crs=epsg),cellsize=cellsize,what="centers")
  if(length(k_grid[polyz_u %>% st_transform(crs=epsg)])==1){
    k_grid <- st_make_grid(polyz_u %>% st_transform(crs=epsg),n=25,what="centers")
  }

  # Ordinary kriging
  krig_mat <- lapply(seq_along(varz),function(v0){
    yvar <- varz[v0]
    a1 <- data.frame(x = rep(NA,nrow(polyz)), x_sd = rep(NA,nrow(polyz)), stringsAsFactors = F) %>% as.data.table() %>% data.table::setnames(old=c("x","x_sd"),new=c(varz[v0],paste0(varz[v0],"_sd")))
    if(pointz_tr_[yvar] %>% na.omit() %>% nrow() > 2 & pointz_tr_ %>% as.data.table() %>% dplyr::select(yvar) %>% unlist() %>% unique() %>% length() > 2){
      print(paste0("Krige ",messagez," ",v0,"/",length(varz)))
      krig_form <- as.formula(paste0(yvar,"~1"))
      krig_out <- suppressMessages(
        automap::autoKrige(krig_form, as(pointz_tr_[yvar] %>% na.omit(),"Spatial"), new_data=as(k_grid,"Spatial"))
      )
      krig_pix <- sp::SpatialPixelsDataFrame(krig_out$krige_output %>% sp::coordinates(),data = krig_out$krige_output %>% as.data.frame())
      sp::proj4string(krig_pix) <- st_crs(epsg)$proj4string
      krig_fit <- krig_pix[,"var1.pred"] %>% raster::raster() %>% raster::projectRaster( crs = sp::proj4string(as(polyz_u,"Spatial")))
      krig_sd <- krig_pix[,"var1.stdev"] %>% raster::raster() %>% raster::projectRaster( crs = sp::proj4string(as(polyz_u,"Spatial")))
      a1 <- data.frame(
        x = raster::extract(krig_fit,polyz %>% as("Spatial"),fun=mean,na.rm=T),
        x_sd = raster::extract(krig_sd,polyz %>% as("Spatial"),fun=mean,na.rm=T),
        stringsAsFactors = FALSE
      ) %>% as.data.table() %>% data.table::setnames(old=c("x","x_sd"),new=c(varz[v0],paste0(varz[v0],"_sd")))
    }
    a1
  }) %>% dplyr::bind_cols() %>% as.data.table() %>% data.table::setnames(1:(2*length(varz)),paste0(rep(varz,each=2),c("_kr","_kr_sd")))

  # Merge with polygons
  krig_mat <- polyz %>% dplyr::bind_cols(krig_mat)

  # Output
  return(krig_mat)
}
