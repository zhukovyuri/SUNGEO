#' Convert data.frame object into simple features object
#'
#' Function takes in x-, y-coordinates, and a data.frame of variables (optional) and returns an SFC object
#'
#' @param x_coord Numeric vector with longitude or easting projected coordinates. When \code{input_data} or \code{file} is supplied, can be either column name or numeric vector of the same length as \code{nrow(input_data)}.
#' @param y_coord Numeric vector with latitude or northing projected coordinates. Must be equal to the vector length of \code{x_coord}. When \code{input_data} or \code{file} is supplied, can be either column name or numeric vector of the same length as \code{nrow(input_data)}.
#' @param input_data Optional data frame object, containing \code{x_coord} and \code{y_coord}. \code{nrow(input_data)} must be equal to the vector length of \code{x_coord}. NOTE: Rows corresponding to non-usable coordinates are removed from the final output.
#' @param file Optional path to csv file. Overrides \code{input_data}.
#' @param n_max Maximum number of rows to read in \code{file}. Default is \code{Inf}.
#' @param start Number of rows to skip in \code{file}. Default is 0 (start on first row).
#' @param projection_input Projection string associated with \code{x_coord} and \code{y_coord}. Default is \code{'+proj=longlat'}.
#' @param zero.policy If \code{TRUE}, removes rows where corresponding coordinates equals (0,0). Default is \code{FALSE}.
#' @param show_removed If \code{TRUE}, returns a vector of indices corresponding to non-usable coordinates. Default is \code{FALSE}.
#' @return If \code{show_removed==FALSE}, returns an \code{sf} object, with rows corresponding to non-usable coordinates removed. If \code{show_removed==TRUE}, returns a list, with an \code{sf} object (\code{Spatial_Coordinates}), and a vector of indices corresponding to non-usable coordinates removed (\code{Removed_Rows}).
#' @importFrom data.table fread
#' @importFrom sf st_as_sf st_crs st_geometry
#' @examples
#' # Coordinates supplied as vectors
#' \dontrun{
#' data(clea_deu2009_df)
#' out_1 <- df2sf(x_coord=clea_deu2009_df$longitude,y_coord = clea_deu2009_df$latitude)
#' class(out_1)
#' plot(out_1$geometry)
#' }
#' # Coordinates supplied as column mames
#' \dontrun{
#' out_2 <- df2sf(x_coord="longitude",y_coord ="latitude", input_data = clea_deu2009_df)
#' plot(out_2["geometry"])
#' }
#' # Load from external file
#' \dontrun{
#' tmp <- tempfile()
#' write.csv(clea_deu2009_df,file=tmp)
#' out_3 <- df2sf(x_coord="longitude",y_coord ="latitude", file=tmp)
#' plot(out_3["geometry"])
#' }
#' @export
#'

df2sf <- function(x_coord, y_coord, input_data = NULL, file = NULL, n_max = Inf, start = 0, projection_input = 'EPSG:4326', zero.policy = FALSE, show_removed = FALSE){
  #Part 1 -
  if(is.null(file) == F){
    input_data <- as.data.frame(data.table::fread(file, nrows = n_max, skip = start))

  }

  #Part 2 -
  if(!is.null(input_data)&(length(x_coord)>1|length(y_coord)>1)){
    message('Warning: When input_data is supplied, x_coord and y_coord should be either column names or numeric vectors of the same length as nrow(input_data)')
  }

  if(is.null(input_data) == F & (length(x_coord)==1&length(y_coord)==1)){
    locationCoords_x <- which(names(input_data)%in%x_coord)

    locationCoords_y <- which(names(input_data)%in%y_coord)

    Coordinate_Matrix <- cbind(input_data[locationCoords_x], input_data[locationCoords_y])

  } else {
    Coordinate_Matrix <- cbind(x_coord, y_coord)

  }

  #Part 3 - Identify NA Results
  if(zero.policy == T){
    Location_Missing_Rows <- which(
      is.na(Coordinate_Matrix[,1]) |
        is.na(Coordinate_Matrix[,2]) |

        is.nan(Coordinate_Matrix[,1]) |
        is.nan(Coordinate_Matrix[,2]) |

        Coordinate_Matrix[,1] == 0 |
        Coordinate_Matrix[,2] == 0)
  } else {
    Location_Missing_Rows <- which(
      is.na(Coordinate_Matrix[,1]) |
        is.na(Coordinate_Matrix[,2]) |

        is.nan(Coordinate_Matrix[,1]) |
        is.nan(Coordinate_Matrix[,2])
    )
  }

  #Part 4 -
  if(length(Location_Missing_Rows) > 0){
    Coordinate_Matrix <- Coordinate_Matrix[-Location_Missing_Rows,]
  }

  #Part 5 -
  Coordinate_Matrix <- as.data.frame(Coordinate_Matrix)

  #Part 6 -
  ST_SFC <- sf::st_as_sf(Coordinate_Matrix, coords = c(1,2))

  #Part 7 -
  sf::st_crs(ST_SFC) <- projection_input

  #Part 8 -
  if(is.null(input_data) == F){
    if(length(Location_Missing_Rows) > 0){
      input_data <- input_data[-Location_Missing_Rows,]

      sf::st_geometry(input_data) <- ST_SFC$geometry
    } else {
      sf::st_geometry(input_data) <- ST_SFC$geometry
    }
  }

  #RETURN
  if(show_removed){
    if(is.null(input_data)){
      return(list('Spatial_Coordinates' = ST_SFC, 'Removed_Rows' = Location_Missing_Rows))
    } else {
      return(list('Spatial_Coordinates' = input_data, 'Removed_Rows' = Location_Missing_Rows))
    }
  } else {
    if(is.null(input_data)){
      return(ST_SFC)
    } else {
      return(input_data)
    }
  }

  rm(Coordinate_Matrix, Location_Missing_Rows, ST_MULTIPOINT, ST_SFC)
  gc()
}


