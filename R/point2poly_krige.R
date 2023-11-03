#' Point-to-polygon interpolation, ordinary and universal Kriging method
#'
#' Function for interpolating values from a source points layer to an overlapping destination polygon layer, using ordinary and universal kriging with automatic variogram fitting
#'
#' @param pointz Source points layer. \code{sf}, \code{sp}, or data frame object.
#' @param polyz Destination polygon layer. Must have identical CRS to \code{pointz}. \code{sf}, \code{sp}, or data frame object.
#' @param rasterz Source raster layer (or list of raster), with covariate(s) used for universal kriging. Must have identical CRS to \code{polyz}.  \code{RasterLayer} object or list of \code{RasterLayer} objects.
#' @param yvarz Names of numeric variable(s) to be interpolated from source points layer to destination polygons. Character string or vector of character strings.
#' @param xvarz Names of numeric variable(s) for universal Kriging, in which yvarz is linearly dependent. Character string or vector of character strings.
#' @param pycno_yvarz  Names of spatially extensive numeric variables for which the pycnophylactic (mass-preserving) property should be preserved. Must be a subset of \code{yvarz}. Character string or vector of character strings.
#' @param funz Aggregation function to be applied to values in \code{rasterz} and to interpolated values. Must take as an input a vector \code{x}. Default is mean.  Function.
#' @param use_grid Use regular grid as destination layer for interpolation, before aggregating to polygons? Default is FALSE.
#' @param nz_grid Number of grid cells in x and y direction (columns, rows). Integer of length 1 or 2. Default is 25. Ignored if use_grid=FALSE.
#' @param blockz Size of blocks used for Block Kriging, in meters. Integer of length 1 or 2. Default is 0.
#' @param pointz_x_coord Name of numeric variable corresponding to a measure of longitude (Easting) in a data frame object for \code{pointz}. Character string.
#' @param pointz_y_coord Name of numeric variable corresponding to a measure of Latitude (Northing) in a data frame object for \code{pointz}. Character string.
#' @param polyz_x_coord Name of numeric variable corresponding to a measure of longitude (Easting) in a data frame object for \code{polyz}. Character string.
#' @param polyz_y_coord Name of numeric variable corresponding to a measure of Latitude (Northing) in a data frame object for \code{polyz}. Character string.
#' @param messagez Optional message to be printed during Kriging estimation. Character string.
#' @return \code{sf} polygon object, with variables from \code{pointz} interpolated to the geometries of \code{polyz}.
#' @details This function performs Ordinary and Universal Kriging, automatically selecting a variogram model with the smallest residual sum of squares from the sample variogram. See \link[automap]{autofitVariogram}.
#'
#' Unlike other available point-to-polygon interpolation techniques, this function currently only accepts numeric variables in \code{varz} and does not support interpolation of character strings.
#' @importFrom sf sf_use_s2 st_crs st_as_sf st_transform st_make_grid st_intersects st_nearest_feature st_make_valid
#' @importFrom terra extract crs values
#' @importFrom automap autoKrige
#' @importFrom stats as.formula aggregate sd
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
#'
#' # Block Kriging with block size of 100 km
#' \dontrun{
#' data(clea_deu2009)
#' data(clea_deu2009_pt)
#' out_4 <- point2poly_krige(pointz = clea_deu2009_pt,
#'                          polyz = clea_deu2009,
#'                          yvarz = "to1",
#'                          blockz = 100000)
#' par(mfrow=c(1,2))
#' plot(clea_deu2009["to1"], key.pos = NULL, reset = FALSE)
#' plot(out_4["to1.pred"], key.pos = NULL, reset = FALSE)
#' }
#' @export


point2poly_krige <- function(pointz,
                             polyz,
                             rasterz=NULL,
                             yvarz=NULL,
                             xvarz=NULL,
                             pycno_yvarz=NULL,
                             funz=base::mean,
                             use_grid=FALSE,
                             nz_grid=25,
                             blockz = 0,
                             pointz_x_coord=NULL,
                             pointz_y_coord=NULL,
                             polyz_x_coord=NULL,
                             polyz_y_coord=NULL,
                             messagez=""){


  # Turn off s2 processing
  suppressMessages({
    sf::sf_use_s2(FALSE)
  })


  # Prepare point layer
  if(class(pointz)[1] %in% 'sf'){
    krig_pointz <- pointz
    if(is.na(sf::st_crs(krig_pointz)$input)){sf::st_crs(krig_pointz)<-"EPSG:4326"}
  } else if(class(pointz)[1] %in% 'data.frame'){
    if(is.null(pointz_x_coord) == TRUE & is.null(pointz_y_coord) == TRUE){stop("Please supply both a pointz_x_coord and pointz_y_coord.")}
    colnames(pointz)[which(colnames(pointz)== pointz_x_coord)] <- "x"
    colnames(pointz)[which(colnames(pointz)== pointz_y_coord)] <- "y"
    krig_pointz <- df2sf(x_coord="x", y_coord="y", input_data=pointz)
    if(is.na(sf::st_crs(krig_pointz)$input)){sf::st_crs(krig_pointz)<-"EPSG:4326"}
  } else if(attr(class(pointz), 'package')[1] == 'sp'){
    krig_pointz <- sf::st_as_sf(pointz)
    if(is.na(sf::st_crs(krig_pointz)$input)){sf::st_crs(krig_pointz)<-"EPSG:4326"}
  } else{stop("Please supply either a data frame, sf, or sp object for pointz.")}
  if(class(polyz)[1] %in% 'sf'){
    krig_polyz <- polyz
    if(is.na(sf::st_crs(krig_polyz)$input)){sf::st_crs(krig_polyz)<-"EPSG:4326"}
  } else if(class(polyz)[1] %in% 'data.frame'){
    if(is.null(polyz_x_coord) == TRUE & is.null(polyz_y_coord) == TRUE){stop("Please supply both a pointz_x_coord and pointz_y_coord.")}
    colnames(polyz)[which(colnames(polyz)== polyz_x_coord)] <- "x"
    colnames(polyz)[which(colnames(polyz)== polyz_y_coord)] <- "y"
    krig_polyz <- df2sf(x_coord="x", y_coord="y", input_data=polyz)
    if(is.na(sf::st_crs(krig_polyz)$input)){sf::st_crs(krig_polyz)<-"EPSG:4326"}
  } else if(attr(class(polyz), 'package') == 'sp'){
    krig_polyz <- sf::st_as_sf(polyz)
    if(is.na(sf::st_crs(krig_polyz)$input)){sf::st_crs(krig_polyz)<-"EPSG:4326"}
  } else{stop("Please supply either a sf or sp object for polyz.")}

  # Check that xvarz entry is in both pointz and polyz
  if(is.null(xvarz) == FALSE){
    pointz_check <- ifelse(xvarz %in% colnames(krig_pointz), TRUE, FALSE)
    poly_check <- ifelse(xvarz %in% colnames(krig_polyz), TRUE, FALSE)
    if (pointz_check == FALSE | poly_check == FALSE) {stop("Please ensure that column names in xvarz appear in both pointz and polyz.")}
  }

  # Extract raster information
  if(is.null(rasterz) == FALSE){
    finalrasterz2polyz <- NULL
    finalrasterz2pointz <- NULL
    if(inherits(rasterz,"list")){
      for(v0 in seq_along(rasterz)){
        if(stats::sd(terra::values(rasterz[[v0]]),na.rm=TRUE)==0){
          terra::values(rasterz[[v0]]) <- jitter(terra::values(rasterz[[v0]]))
        }
        krig_polyz_ <- sf::st_transform(krig_polyz, sf::st_crs(rasterz[[v0]]))
        extractdata <- suppressWarnings(terra::extract(rasterz[[v0]], krig_polyz_))
        rastervarz <- unlist(lapply(extractdata, function(x) if (!is.null(x)) funz(x, na.rm=TRUE) else NA))
        finalrasterz2polyz <- cbind(finalrasterz2polyz, rastervarz)
        rm(krig_polyz_)
        krig_pointz_ <- sf::st_transform(krig_pointz, sf::st_crs(rasterz[[v0]]))
        extractdata <- suppressWarnings(terra::extract(rasterz, krig_pointz_))
        rastervarz <- unlist(lapply(extractdata, function(x) if (!is.null(x)) funz(x, na.rm=TRUE) else NA))
        finalrasterz2pointz <- cbind(finalrasterz2pointz, rastervarz)
        rm(krig_pointz_)
      }} else{
        if(stats::sd(terra::values(rasterz),na.rm=TRUE)==0){
          terra::values(rasterz) <- jitter(terra::values(rasterz))
        }

        krig_polyz_ <- sf::st_transform(krig_polyz, sf::st_crs(rasterz))
        extractdata <- suppressWarnings(terra::extract(rasterz, krig_polyz_))
        rastervarz <- unlist(lapply(extractdata, function(x) if (!is.null(x)) funz(x, na.rm=TRUE) else NA))
        finalrasterz2polyz<- cbind(finalrasterz2polyz, rastervarz)
        rm(krig_polyz_)
        krig_pointz_ <- sf::st_transform(krig_pointz, sf::st_crs(rasterz))
        extractdata <- suppressWarnings(terra::extract(rasterz, krig_pointz_))
        rastervarz <- unlist(lapply(extractdata, function(x) if (!is.null(x)) funz(x, na.rm=TRUE) else NA))
        finalrasterz2pointz <- cbind(finalrasterz2pointz, rastervarz)
        rm(krig_pointz_)
      }

    # Combine with pointz and polyz
    col_length <- ncol(krig_polyz)
    krig_polyz <- cbind(krig_polyz, finalrasterz2polyz)
    krig_pointz <- cbind(krig_pointz, finalrasterz2pointz)

    #add raster variables to xvarz
    xvarz <- c(xvarz, colnames(finalrasterz2polyz))
  }

  # Create empty prediction grid
  if(use_grid==TRUE){
    suppressMessages({
      suppressWarnings({
        # k_grid <- sf::st_make_grid(sf::st_as_sf(krig_polyz),n=nz_grid,what="centers")
        k_grid <- sf::st_make_grid(krig_polyz,n=nz_grid,what="centers")
        krig_polyz$ID_kriggrid <- 1:nrow(krig_polyz)
        k_ix <- sf::st_intersects(k_grid,krig_polyz)
        if(any(lengths(k_ix)==0)&any(lengths(k_ix)>1)){
          suppressMessages({sf::sf_use_s2(TRUE)})
          k_mat <- rbind(as.data.frame(k_ix),data.frame(row.id=which(lengths(k_ix)==0),col.id=sf::st_nearest_feature(k_grid[lengths(k_ix)==0],sf::st_make_valid( sf::st_as_sf(krig_polyz) ))))
          k_mat <- k_mat[order(k_mat$row.id),]
          k_mat <- k_mat[!duplicated(k_mat$row.id),]
          suppressMessages({sf::sf_use_s2(FALSE)})
        } else if(any(lengths(k_ix)==0)){
          suppressMessages({sf::sf_use_s2(TRUE)})
          k_mat <- rbind(as.data.frame(k_ix),data.frame(row.id=which(lengths(k_ix)==0),col.id=sf::st_nearest_feature(k_grid[lengths(k_ix)==0],sf::st_make_valid( sf::st_as_sf(krig_polyz) ))))
          k_mat <- k_mat[order(k_mat$row.id),]
          suppressMessages({sf::sf_use_s2(FALSE)})
        }else if(any(lengths(k_ix)>1)){
          k_mat <- as.data.frame(k_ix)
          k_mat <- k_mat[!duplicated(k_mat$row.id),]
        }else{
          k_mat <- as.data.frame(k_ix)
        }
        k_grid <- sf::st_as_sf(cbind(k_grid,as.data.frame(krig_polyz)[k_mat$col.id,]))
      })
    })
  }

  # Find optimal planar projection for map
  suppressMessages({
    suppressWarnings({
      sf::st_crs(krig_polyz)$input
      if(grepl("4326",sf::st_crs(krig_polyz)$input)){
        polyz_layer <- utm_select(krig_polyz)
        pointz_layer <- sf::st_transform(krig_pointz, sf::st_crs(polyz_layer))
        if(use_grid==TRUE){
          gridz_layer <- sf::st_transform(k_grid, sf::st_crs(polyz_layer))
        }
      }else{
        polyz_layer <- krig_polyz
        pointz_layer <- sf::st_transform(krig_pointz, sf::st_crs(polyz_layer))
        if(use_grid==TRUE){
          gridz_layer <- sf::st_transform(k_grid, sf::st_crs(polyz_layer))
        }
      }
    })
  })

  if(use_grid==TRUE){
    # autokrige (to grid)
    krige_mat <- lapply(seq_along(yvarz), function(v1){
      if(is.null(xvarz) == TRUE){
        krige_form <- stats::as.formula(paste(yvarz[v1],1, sep = "~"))
        krige_result <- suppressWarnings({automap::autoKrige(krige_form, pointz_layer[!is.na(as.data.frame(pointz_layer)[,yvarz[v1]]),], gridz_layer,block=blockz)})
      }else if (length(xvarz) >= 1){
        krige_form <- stats::as.formula(paste(yvarz[v1],paste0(xvarz, collapse = "+"), sep = "~"))
        krige_result <- suppressWarnings(automap::autoKrige(krige_form, pointz_layer[!is.na(as.data.frame(pointz_layer)[,xvarz])&!is.na(as.data.frame(pointz_layer)[,yvarz[v1]]),], gridz_layer,block=blockz))
      }
      colnames(krige_result[["krige_output"]])[grep("pred$|var$|stdev$",colnames(krige_result[["krige_output"]]))] <- c(paste0(yvarz[v1],".pred"), paste0(yvarz[v1],".var"),paste0(yvarz[v1],".stdev"))
      out <- as.data.frame(krige_result[["krige_output"]][,grep("pred$|var$|stdev$",colnames(krige_result[["krige_output"]]))])
      out$geometry <- NULL
      return(out)
    })

    # bind with polygons
    krige_agg <- stats::aggregate(krige_mat,by=list(ID_kriggrid=gridz_layer$ID_kriggrid),FUN=funz, na.rm=TRUE, na.action=NULL)

    # Pycno
    if(!is.null(pycno_yvarz)){
      # Loop over pycno_varz
      for(p0 in 1:length(pycno_yvarz)){
        # Find sum of original variable
        sum_from <- data.table::as.data.table(krig_pointz)[,sum(get(pycno_yvarz[p0]),na.rm=TRUE)]
        # Find matching processed variables in target geometry
        pycno_yvarz_to <- grep(paste0("^",pycno_yvarz[p0]),names(krige_agg),value=TRUE)
        # Rescale variables in target geometry
        krige_agg_dt <- data.table::as.data.table(krige_agg)
        for(p00 in 1:length(pycno_yvarz_to)){
          krige_agg_dt[,eval(pycno_yvarz_to[p00]) := get(pycno_yvarz_to[p00])*sum_from/sum(get(pycno_yvarz_to[p00]),na.rm = TRUE)]
        }
        krige_agg <- as.data.frame(krige_agg_dt)
      }
    }
    if(class(polyz)[1] %in% "sf"){
      krige_out <- merge(polyz_layer, krige_agg)
      krige_out$ID_kriggrid <- NULL
      if(!"geometry"%in%names(krige_out)){
        krige_out$geometry <- polyz$geometry
      }
      # krige_out <- suppressWarnings(sf::st_as_sf(krige_out))
    }else if(class(polyz)[1] %in% "data.frame"){
      krige_out <- merge(polyz_layer, krige_agg)
      krige_out$ID_kriggrid <- NULL
      krige_out <- sf::st_as_sf(krige_out)
    }else if(attr(class(polyz), 'package') == 'sp'){
      krige_out <- merge(polyz_layer, krige_agg)
      krige_out$ID_kriggrid <- NULL
    }
  } else {
    # autokrige (to polyz)
    krige_mat <- lapply(seq_along(yvarz), function(v1){
      if(is.null(xvarz) == TRUE){
        krige_form <- stats::as.formula(paste(yvarz[v1],1, sep = "~"))
        krige_result <- suppressWarnings(automap::autoKrige(krige_form, pointz_layer[!is.na(as.data.frame(pointz_layer)[,yvarz[v1]]),], polyz_layer,block=blockz))
      }else if (length(xvarz) >= 1){
        krige_form <- stats::as.formula(paste(yvarz[v1],paste0(xvarz, collapse = "+"), sep = "~"))
        krige_result <- suppressWarnings(automap::autoKrige(krige_form, pointz_layer[!is.na(as.data.frame(pointz_layer)[,xvarz])&!is.na(as.data.frame(pointz_layer)[,yvarz[v1]]),], polyz_layer,block=blockz))
      }
      colnames(krige_result[["krige_output"]])[grep("pred$|var$|stdev$",colnames(krige_result[["krige_output"]]))] <- c(paste0(yvarz[v1],".pred"), paste0(yvarz[v1],".var"),paste0(yvarz[v1],".stdev"))
      out <- as.data.frame(krige_result[["krige_output"]][,grep("pred$|var$|stdev$",colnames(krige_result[["krige_output"]]))])
      out$geometry <- NULL
      return(out)
    })
    # bind with polygons
    if(class(polyz)[1] %in% "sf"){
      krige_out <- cbind(polyz_layer, krige_mat)
      krige_out <- sf::st_as_sf(krige_out)
    }else if(class(polyz)[1] %in% "data.frame"){
      krige_out <- cbind(polyz_layer, krige_mat)
      krige_out <- sf::st_as_sf(krige_out)
    }else if(attr(class(polyz), 'package') == 'sp'){
      krige_out <- cbind(polyz_layer, krige_mat)
    }
  }

  #Output
  return(krige_out)
}
