#' Point-to-polygon interpolation, ordinary and universal Kriging method
#'
#' Function for interpolating values from a source points layer to an overlapping destination polygon layer, using ordinary and universal kriging with automatic variogram fitting
#'
#' @param pointz Source points layer. \code{sf}, \code{sp}, or data frame object.
#' @param polyz Destination polygon layer. Must have identical CRS to \code{pointz}. \code{sf}, \code{sp}, or data frame object.
#' @param rasterz Source raster layer (or list of raster), with covariate(s) used for universal kriging. Must have identical CRS to \code{polyz}.  \code{RasterLayer} object or list of \code{RasterLayer} objects.
#' @param yvarz Names of numeric variable(s) to be interpolated from source points layer to destination polygons. Character string or vector of character strings.
#' @param xvarz Names of numeric variable(s) for universal Kriging, in which yvarz is linearly dependent. Character string or vector of character strings.
#' @param pointz_x_coord Name of numeric variable corresponding to a measure of longitude (Easting) in a data frame object for \code{pointz}. Character string.
#' @param pointz_y_coord Name of numeric variable corresponding to a measure of Latitude (Northing) in a data frame object for \code{pointz}. Character string.
#' @param polyz_x_coord Name of numeric variable corresponding to a measure of longitude (Easting) in a data frame object for \code{polyz}. Character string.
#' @param polyz_y_coord Name of numeric variable corresponding to a measure of Latitude (Northing) in a data frame object for \code{polyz}. Character string.
#' @param messagez Optional message to be printed during Kriging estimation. Character string.
#' @return \code{sp} polygon object, with variables from \code{pointz} interpolated to the geometries of \code{polyz}.
#' @details This function performs Ordinary and Universal Kriging, automatically selecting a variogram model with the smallest residual sum of squares from the sample variogram. See \link[automap]{autofitVariogram}.
#'
#' Unlike other available point-to-polygon interpolation techniques, this function currently only accepts numeric variables in \code{varz} and does not support interpolation of character strings.
#' @import sf sp
#' @importFrom raster extract
#' @importFrom automap autoKrige
#' @importFrom stats as.formula
#' @examples
#' # Ordinary Kriging with one variable
#' \dontrun{
#' data(clea_deu2009)
#' data(clea_deu2009_pt)
#' out_1 <- point2poly_krige(pointz = clea_deu2009_pt,
#'                          polyz = clea_deu2009,
#'                          yvarz = "to1")
#' par(mfrow=c(1,2))
#' plot(clea_deu2009["to1"], key.pos = NULL, reset = FALSE)
#' plot(out_1["to1.pred"], key.pos = NULL, reset = FALSE)
#' }
#'
#' # Ordinary Kriging with multiple variables
#' \dontrun{
#' out_2 <- point2poly_krige(pointz = clea_deu2009_pt,
#'                          polyz = clea_deu2009,
#'                          yvarz = c("to1","pvs1_margin"))
#' par(mfrow=c(1,2))
#' plot(clea_deu2009["pvs1_margin"], key.pos = NULL, reset = FALSE)
#' plot(out_2["pvs1_margin.pred"], key.pos = NULL, reset = FALSE)
#' }
#'
#' # Universal Kriging with one variable from a raster
#' \dontrun{
#' data(gpw4_deu2010)
#' data(clea_deu2009)
#' data(clea_deu2009_pt)
#' out_3 <- point2poly_krige(pointz = clea_deu2009_pt,
#'                          polyz = clea_deu2009,
#'                          yvarz = "to1",
#'                          rasterz = gpw4_deu2010)
#' par(mfrow=c(1,2))
#' plot(clea_deu2009["to1"], key.pos = NULL, reset = FALSE)
#' plot(out_3["to1.pred"], key.pos = NULL, reset = FALSE)
#' }
#' @export


point2poly_krige <- function(pointz,
                             polyz,
                             rasterz=NULL,
                             yvarz=NULL,
                             xvarz=NULL,
                             pointz_x_coord=NULL,
                             pointz_y_coord=NULL,
                             polyz_x_coord=NULL,
                             polyz_y_coord=NULL,
                             messagez=""){

  # Convert SF/DataFrame to SP
  if(class(pointz)[1] == 'sf'){
    krig_pointz <- as(pointz, "Spatial")
  } else if(class(pointz)[1] == 'data.frame'){
    if(is.null(pointz_x_coord) == TRUE & is.null(pointz_y_coord) == TRUE){stop("Please supply both a pointz_x_coord and pointz_y_coord.")}
    colnames(pointz)[which(colnames(pointz)== pointz_x_coord)] <- "x"
    colnames(pointz)[which(colnames(pointz)== pointz_y_coord)] <- "y"
    sf.convert_pointz <- df2sf(pointz$x, pointz$y, pointz)
    krig_pointz <- as(sf.convert_pointz, "Spatial")
  } else if(attr(class(pointz), 'package')[1] == 'sp'){
    krig_pointz <- pointz
  } else{stop("Please supply either a data frame, sf, or sp object for pointz.")}

  if(class(polyz)[1] == 'sf'){
    krig_polyz <- as(polyz, "Spatial")
  } else if(class(polyz)[1] == 'data.frame'){
    if(is.null(polyz_x_coord) == TRUE & is.null(polyz_y_coord) == TRUE){stop("Please supply both a pointz_x_coord and pointz_y_coord.")}
    colnames(polyz)[which(colnames(polyz)== polyz_x_coord)] <- "x"
    colnames(polyz)[which(colnames(polyz)== polyz_y_coord)] <- "y"
    sf.convert_polyz <- df2sf(polyz$x, polyz$y, polyz)
    krig_polyz <- as(sf.convert_polyz, "Spatial")
  } else if(attr(class(polyz), 'package') == 'sp'){
    krig_polyz <- polyz
  } else{stop("Please supply either a sf or sp object for polyz.")}

  # Check that xvarz entry is in both pointz and polyz
  if (is.null(xvarz) == FALSE){
    pointz_check <- ifelse(xvarz %in% colnames(krig_pointz@data), TRUE, FALSE)
    poly_check <- ifelse(xvarz %in% colnames(krig_polyz@data), TRUE, FALSE)
    if (pointz_check == FALSE | poly_check == FALSE) {stop("Please ensure that column names in xvarz appear in both pointz and polyz.")}
  }

  # Extract raster information
  if(is.null(rasterz) == FALSE){
    finalrasterz2polyz <- NULL
    if(class(rasterz)=="list"){
      for(v0 in seq_along(rasterz)){
        krig_polyz_ <- sp::spTransform(krig_polyz, sp::proj4string(rasterz[[v0]]))
        extractdata <- suppressWarnings(raster::extract(rasterz[[v0]], krig_polyz_))
        rastervarz <- unlist(lapply(extractdata, function(x) if (!is.null(x)) sum(x, na.rm=TRUE) else NA))
        finalrasterz2polyz<- cbind(finalrasterz2polyz, rastervarz)
        rm(krig_polyz_)
      }} else{
        krig_polyz_ <- sp::spTransform(krig_polyz, sp::proj4string(rasterz))
        extractdata <- suppressWarnings(raster::extract(rasterz, krig_polyz_))
        rastervarz <- unlist(lapply(extractdata, function(x) if (!is.null(x)) sum(x, na.rm=TRUE) else NA))
        finalrasterz2polyz<- cbind(finalrasterz2polyz, rastervarz)
        rm(krig_polyz_)
      }

    #c Cmbine with pointz and polyz
    col_length <- ncol(krig_polyz@data)
    krig_polyz <- cbind(krig_polyz, finalrasterz2polyz)
    krig_pointz <- cbind(krig_pointz, finalrasterz2polyz)

    #add raster variables to xvarz
    xvarz <- c(xvarz, colnames(krig_polyz@data)[(col_length+1):ncol(krig_polyz@data)])
  }

  # Find optimal planar projection for map
  suppressMessages({
    suppressWarnings({
      polyz_layer <- utm_select(krig_polyz)
      pointz_layer <- sp::spTransform(krig_pointz, sp::proj4string(polyz_layer))
    })
  })

  #autokrige
  krige_mat <- lapply(seq_along(yvarz), function(v1){
    if(is.null(xvarz) == TRUE){
      krige_form <- stats::as.formula(paste(yvarz[v1],1, sep = "~"))
      krige_result <- suppressWarnings(automap::autoKrige(krige_form, pointz_layer, polyz_layer))
    }else if (length(xvarz) >= 1){
      krige_form <- stats::as.formula(paste(yvarz[v1],paste0(xvarz, collapse = "+"), sep = "~"))
      krige_result <- suppressWarnings(automap::autoKrige(krige_form, pointz_layer, polyz_layer))
    }
    colnames(krige_result[["krige_output"]]@data)[3:5] <- c(paste0(yvarz[v1],".pred"), paste0(yvarz[v1],".var"),paste0(yvarz[v1],".stdev"))
    out <- krige_result[["krige_output"]]@data[,c(3,5)]
    return(out)
  })

  #Mbind with polygons
  if(class(polyz)[1] == "sf"){
    krige_out <- cbind(polyz_layer, krige_mat)
    krige_out <- sf::st_as_sf(krige_out)
  }else if(class(polyz)[1] == "data.frame"){
    krige_out <- cbind(polyz_layer, krige_mat)
    krige_out <- sf::st_as_sf(krige_out)
  }else if(attr(class(polyz), 'package') == 'sp'){
    krige_out <- cbind(polyz_layer, krige_mat)
  }

  #Output
  return(krige_out)
}



