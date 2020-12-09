#' Automatically calculate Local G hot spot intensity
#'
#' Function automatically calculates the Local G hot spot intensity measure for spatial points, spatial polygons, and single raster layers. Uses RANN for efficient nearest neighbor calculation (spatial points and single raster layers only); users can specify the number of neighbors (k). Users can specify the neighborhood style (see spdep::nb2listw) with default being standardized weight matrix (W).
#'
#' @param insert Spatial point, spatial polygon, or single raster layer object. Acceptable formats include \code{sf}, \code{SpatialPolygonsDataFrame}, \code{SpatialPointsDataFrame}, and \code{RasterLayer}.
#' @param variable Column name or numeric vector containing the variable from which the local G statistic will be calculated. Must possess a natural scale that orders small and large observations (i.e. number, percentage, ratio and not model residuals).
#' @param style Style can take values \code{'W'}, \code{'B'}, \code{'C'}, \code{'U'}, \code{'mimax'}, \code{'S'} (see  \code{\link[spdep]{nb2listw}}). Character string.
#' @param k Number of neighbors. Default is 9. Numeric.
#' @param remove_missing Whether to calculate statistic without missing values. If \code{FALSE}, substitute value must be supplied to \code{NA_Value}.
#' @param NA_Value Substitute for missing values. Default value is 0. Numeric.
#' @param include_Moran Calculate local Moran's I statistics. Default is \code{FALSE}. Logical.
#' @return If \code{input} is \code{sf}, \code{SpatialPolygonsDataFrame} or \code{SpatialPointsDataFrame} object, returns \code{sf} object with same geometries and columns as \code{input}, appended with additional column containing Local G estimates (\code{LocalG}). If \code{input} is \code{RasterLayer} object, returns \code{RasterBrick} object containing original values (\code{Original}) and Local G estimates (\code{LocalG}).
#' @import sp
#' @importFrom sf st_as_sf st_is st_coordinates
#' @importFrom raster raster brick values as.data.frame
#' @importFrom spdep poly2nb nb2listw localG localmoran
#' @importFrom RANN nn2
#' @examples
#' # Calculate Local G for sf point layer
#'
#' \dontrun{
#' data(clea_deu2009_pt)
#' out_1 <- hot_spot(insert=clea_deu2009_pt, variable = clea_deu2009_pt$to1)
#' class(out_1)
#' plot(out_1["LocalG"])
#' }
#'
#' # Calculate Local G for sf polygon layer (variable as numeric vector)
#'
#' \dontrun{
#' data(clea_deu2009)
#' out_2 <- hot_spot(insert=clea_deu2009, variable = clea_deu2009$to1)
#' summary(out_2$LocalG)
#' plot(out_2["LocalG"])
#' }
#'
#' # Calculate Local G for sf polygon layer (variable as column name)
#'
#' \dontrun{
#' out_3 <- hot_spot(insert=clea_deu2009, variable = "to1")
#' summary(out_3$LocalG)
#' plot(out_3["LocalG"])
#' }
#'
#' # Calculate Local G for sf polygon SpatialPolygonsDataFrame (variable as column name)
#'
#' \dontrun{
#' out_4 <- hot_spot(insert=as(clea_deu2009,"Spatial"), variable = "to1")
#' summary(out_4$LocalG)
#' plot(out_4["LocalG"])
#' }
#'
#' # Calculate Local G for RasterLayer
#' \dontrun{
#' data(gpw4_deu2010)
#' out_5 <- hot_spot(insert=gpw4_deu2010)
#' class(out_5)
#' raster::plot(out_5)
#' }
#' @export
hot_spot <- function(insert,
                     variable = NULL,
                     style = 'W',
                     k = 9,
                     remove_missing = TRUE,
                     NA_Value = 0,
                     include_Moran = FALSE){
  #Section 1 -
  if(any(class(insert)%in%c("SpatialPolygonsDataFrame", "SpatialPointsDataFrame"))){
    insert <- sf::st_as_sf(insert)
  }

  #Section 2 -
  if(is.character(variable) && 'sf'%in%class(insert)){
    insert_df <- raster::as.data.frame(insert)

    variable <- insert_df[,names(insert_df)%in%variable]
    rm(insert_df)

    variable <- as.vector(variable)
  }

  #Section 3 -
  if('sf'%in%class(insert)){
    if(all(sf::st_is(insert, 'POLYGON')) | all(sf::st_is(insert, 'MULTIPOLYGON'))){
      #Part A -
      NNobj <- spdep::poly2nb(insert)

      #Part B -
      NNobj <- spdep::nb2listw(NNobj, style = style)

      #Part C -
      LocalG_Output <- spdep::localG(variable, NNobj)

      if(include_Moran%in%T){
        LocalM_Output <- spdep::localmoran(variable, NNobj)[,1]
      }

      #Part D -
      insert$LocalG <- LocalG_Output

      if(include_Moran%in%T){
        insert$LocalM <- LocalM_Output
      }

      #Part E -
      RANN_Act <- 0

      #RETURN
      return(insert)
    }
  }

  #Section 4 -
  if('RasterLayer'%in%class(insert)){
    #Part A -
    variable <- raster::values(insert)

    #Part B -
    IDxy_Matrix <- raster::as.data.frame(insert, xy = T)

    #Part C -
    IDxy_Matrix <- IDxy_Matrix[,c(1:2)]

    #Part D -
    names(IDxy_Matrix) <- c('x', 'y')

    #Part E -
    Location_Missing <- which(is.na(variable))

    #Part F -
    if(T%in%remove_missing && length(Location_Missing) > 0){
      IDxy_Matrix <- IDxy_Matrix[-Location_Missing,]

      variable <- variable[-Location_Missing]
    }

    #Part G -
    if(F%in%remove_missing && length(Location_Missing) > 0){
      variable[Location_Missing] <- NA_Value
    }

    #Part H -
    RANN_Act <- 1

  }

  #Section 5 -
  if('sf'%in%class(insert)){
    if(all(sf::st_is(insert, 'POINT'))){
      #Part A -
      IDxy_Matrix <- raster::as.data.frame(sf::st_coordinates(insert), xy = T)

      #Part B -
      Location_Missing <- which(is.na(variable))

      #Part C -
      if(T%in%remove_missing && length(Location_Missing) > 0){
        IDxy_Matrix <- IDxy_Matrix[-Location_Missing,]
      }

      #Part D -
      if(F%in%remove_missing && length(Location_Missing) > 0){
        variable[Location_Missing] <- NA_Value
      }

      #Part E -
      RANN_Act <- 1

    }
  }

  #Section 6 -
  if(RANN_Act == 1){
    #Part A -
    KNN <- RANN::nn2(IDxy_Matrix[,1:2], IDxy_Matrix[,1:2], k = k)

    #Part B -
    KNN <- lapply(1:nrow(KNN$nn.idx[,1:k]), function(x){KNN$nn.idx[x,1:k]})

    #Part C -
    class(KNN) <- 'nb'

    #Part D -
    NNobj <- spdep::nb2listw(KNN, style = style)

    #Part E -
    LocalG_Output <- spdep::localG(variable, NNobj)

    if(include_Moran%in%T){
      LocalM_Output <- spdep::localmoran(variable, NNobj)[,1]
    }
  }

  #Section 7 -
  if('sf'%in%class(insert)){
    if(all(sf::st_is(insert, 'POINT'))){
      #Part A -
      insert$LocalG <- NA

      #Part B -
      if(length(Location_Missing) > 0){
        insert$LocalG[-Location_Missing] <- as.vector(LocalG_Output)
      } else {
        insert$LocalG <- as.vector(LocalG_Output)
      }


      #RETURN
      return(insert)
    }
  }

  #Section 8 -
  if('RasterLayer'%in%class(insert)){
    #Part A -
    IDxy_Matrix <- raster::as.data.frame(insert, xy = T)

    ##########################
    #Part B -
    ##########################
    IDxy_Matrix$LocalG <- NA

    #Part C -
    if(length(Location_Missing) > 0){
      IDxy_Matrix$LocalG[-Location_Missing] <- as.vector(LocalG_Output)
    } else {
      IDxy_Matrix$LocalG <- as.vector(LocalG_Output)
    }

    #Part D -
    TemporaryRaster <- insert

    #Part E -
    raster::values(TemporaryRaster) <- IDxy_Matrix$LocalG

    ##########################
    #Part C -
    ##########################
    if(include_Moran%in%T){

      IDxy_Matrix$LocalM <- NA

      #Part C -
      if(length(Location_Missing) > 0){
        IDxy_Matrix$LocalM[-Location_Missing] <- as.vector(LocalM_Output)
      } else {
        IDxy_Matrix$LocalM <- as.vector(LocalM_Output)
      }

      #Part D -
      TemporaryRaster_Moran <- insert

      #Part E -
      raster::values(TemporaryRaster_Moran) <- IDxy_Matrix$LocalM


    }

    if(include_Moran%in%F){
      #Part F -
      returnBrick <- raster::brick(insert, TemporaryRaster)

      #Part G -
      names(returnBrick) <- c('Original', 'LocalG')
    } else {
      #Part F -
      returnBrick <- raster::brick(insert, TemporaryRaster, TemporaryRaster_Moran)

      #Part G -
      names(returnBrick) <- c('Original', 'LocalG', 'Moran')
    }


    #RETURN
    return(returnBrick)
  }

  gc()
}

