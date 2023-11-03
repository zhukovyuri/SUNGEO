#' Convert simple features object into regularly spaced raster
#'
#' This function takes in an \code{sf} spatial object (polygon or point) and returns a regularly spaced RasterLayer. Reverse translation option allows users to create an \code{sf} polygon object from the regularly spaced RasterLayer. This function can also conver the \code{sf} object into a cartogram with a user-specified variable name.
#'
#' @param polyz_from Source polygon layer. \code{sf} object.
#' @param pointz_from Source point layer. \code{sf} object.
#' @param input_variable Name of input variable from source layer. Character string.
#' @param aggregate_function Aggregation function to be applied to variables specified in \code{input_variable}. Must take as an input a numeric vector \code{x}. Function or list of functions. Default is \code{mean}.
#' @param reverse Reverse translation from raster layer to \code{sf} polygon object (polygon features only). Default is \code{FALSE}.
#' @param poly_to Destination polygon layer for reverse conversion. Must be specified if \code{reverse=TRUE}. \code{sf} object.
#' @param return_output Return output for reverse conversion. Must be specified if \code{reverse=TRUE}.
#' @param return_field Return field for reverse conversion. Must be specified if \code{reverse=TRUE}.
#' @param reverse_function  Aggregation function for reverse conversion. Must be specified if \code{reverse=TRUE}. Function or list of functions. Default is \code{mean}.
#' @param grid_dim Dimensions of raster grid. Numerical vector of length 2 (number of rows, number of columns). Default is \code{c(1000,1000)}.
#' @param cartogram Cartogram transformation. Logical. Default is \code{FALSE}.
#' @param carto_var Input variable for cartogram transformation. Must be specified if \code{cartogram=TRUE}. Character string.
#' @param message_out Print informational messages. Logical. Default is \code{TRUE}.
#' @param return_list Return full set of results, including input polygons, centroids and field raster. Default is \code{FALSE}. Logical.
#' @return If \code{return_list=FALSE} (default) and \code{reverse=FALSE} (default), returns \code{RasterLayer} object, with cell values corresponding to \code{input_variable}.
#'
#' If \code{return_list=TRUE} and input layer is polygon, returns a list containing
#' \itemize{
##'  \item{"return_output". }{Output raster, with values corresponding to \code{input_variable}. \code{RasterLayer} object.}
##'  \item{"return_centroid". }{Raster of centroids, with values corresponding to \code{input_variable}. \code{RasterLayer} object.}
##'  \item{"poly_to". }{Source polygons, with columns corresponding to \code{input_variable} and auto-generated numerical ID \code{Field}. \code{sf} object.}
##'  \item{"return_field". }{Output raster, with values corresponding to auto-generated numerical ID \code{Field}. \code{RasterLayer} object.}
##'  }
##' If \code{return_list=TRUE} and  input layer is points, returns a list containing
#' \itemize{
##'  \item{"return_output". }{Output raster, with values corresponding to \code{input_variable}. \code{RasterLayer} object.}
##'  \item{"return_point". }{Source points, with column corresponding to \code{input_variable}.}
##'  }
##' If \code{reverse=TRUE}, returns an \code{sf} polygon layer, with columns corresponding to \code{input_variable} and auto-generated numerical ID \code{Field}.
#' @import packcircles cartogram Rcpp
#' @importFrom terra rasterize ncell res ext vect as.data.frame res extract values mask plot
#' @importFrom raster raster
#' @importFrom cartogram cartogram_cont
#' @importFrom sf st_collection_extract st_crs st_geometry st_centroid st_buffer st_transform st_convex_hull st_union
#' @importFrom dplyr group_by summarize left_join
#' @importFrom purrr reduce
#' @examples
#' # Rasterization of polygon layer.
#' \dontrun{
#' data(clea_deu2009)
#' out_1 <- sf2raster(polyz_from = utm_select(clea_deu2009),
#'                    input_variable = "to1")
#' terra::plot(out_1)
#' }
#' # Rasterization of point layer
#' \dontrun{
#' data(clea_deu2009_pt)
#' out_2 <- sf2raster(pointz_from = utm_select(clea_deu2009_pt),
#'                    input_variable = "to1",
#'                    grid_dim = c(25,25))
#' terra::plot(out_2)
#' }
#' # Cartogram (vote turnout scaled by number of valid votes)
#' \dontrun{
#' out_3 <- sf2raster(polyz_from = utm_select(clea_deu2009),
#'                    input_variable = "to1",
#'                    cartogram = TRUE,
#'                    carto_var = "vv1")
#' terra::plot(out_3)
#' }
#' # Polygonization of cartogram raster
#' \dontrun{
#' out_4a <- sf2raster(polyz_from = utm_select(clea_deu2009),
#'                     input_variable = "to1",
#'                     cartogram = TRUE,
#'                     carto_var = "vv1",
#'                     return_list = TRUE)
#' out_4 <- sf2raster(reverse = TRUE,
#'                    poly_to = out_4a$poly_to,
#'                    return_output = out_4a$return_output,
#'                    return_field = out_4a$return_field)
#' terra::plot(out_4)
#' }
#' @export


sf2raster <- function(polyz_from = NULL,
                       pointz_from = NULL,
                       input_variable = NULL,
                       reverse = FALSE,
                       poly_to = NULL,
                       return_output = NULL,
                       return_field = NULL,
                       aggregate_function = list(function(x) mean(x, na.rm = TRUE)),
                       reverse_function = list(function(x) mean(x, na.rm = TRUE)),
                       grid_dim = c(1000, 1000),
                       cartogram = FALSE,
                       carto_var = NULL,
                       message_out = TRUE,
                       return_list = FALSE){

  # Turn off s2 processing
  suppressMessages({
    sf::sf_use_s2(FALSE)
  })


  if (!inherits(aggregate_function, "list")) {
    aggregate_function <- list(aggregate_function)
  }
  if (!inherits(reverse_function, "list")) {
    reverse_function <- list(reverse_function)
  }

  if (!is.null(polyz_from) & message_out && !reverse && !any(grepl("meter|metre",sf::st_crs(polyz_from),ignore.case = TRUE))) {
    message("Note: Please project the SF object into a meter-based projection coordinate system prior to using this function. Default grid resolution is in meters.")
  }
  if (!is.null(pointz_from) & message_out && !reverse && !any(grepl("meter|metre",sf::st_crs(pointz_from),ignore.case = TRUE))) {
    message("Note: Please project the SF object into a meter-based projection coordinate system prior to using this function. Default grid resolution is in meters.")
  }
  if (!is.null(polyz_from) && !is.null(pointz_from) && !reverse) {
    stop("Error: Please enter only one type of SF object {polygon or point}.")
  }
  if (!reverse) {
    CheckObj <- if (!is.null(polyz_from)) {
      polyz_from
    }
    else {
      pointz_from
    }
  } else {
    CheckObj <- poly_to
  }
  if (all(!("sf" %in% class(CheckObj)))) {
    stop("Error: Only SF formatted objects are accepted. Please convert to an sf object and then re-run function")
  }
  if (reverse && (is.null(poly_to) || is.null(return_output) ||
                  is.null(return_field))) {
    stop("Error: If translating back into SF Polygon object, then ensure poly_to, return_output, and return_field are specified")
  }

  if (!reverse) {
    if (message_out) {
      message(paste("Converting SF object into a regularly spaced raster grid of dimensions:",
                    grid_dim[1], "by", grid_dim[2]))
    }
    if (cartogram && !is.null(polyz_from)) {
      if (is.null(carto_var)) {
        stop("Please supply a input_variable for cartogram transformation")
      }
      polyz_from <- cartogram::cartogram_cont(polyz_from,
                                              carto_var)
      polyz_from <- sf::st_collection_extract(polyz_from,
                                              "POLYGON")
    }
    if (cartogram && is.null(polyz_from)) {
      stop("Cartogram translation is only available for spatial polygon objects. Please insert the appropriate spatial feature.")
    }
    if (!is.null(polyz_from)) {
      polyz_from$Field <- 1:nrow(polyz_from)
      sel.col <- which(names(polyz_from) %in% c("Field",
                                                input_variable))
      polyz_from <- polyz_from[, sel.col]
      polyz_from_SP_keep <- polyz_from
    }
    if (!is.null(pointz_from)) {
      sel.col <- which(names(pointz_from) %in% c(input_variable))
      pointz_from <- pointz_from[, sel.col]
      names(pointz_from)[names(pointz_from) %in% input_variable] <- "Var"
      pointz_from_SP_keep <- pointz_from
      names(pointz_from_SP_keep)[names(pointz_from_SP_keep) %in%
                                   "Var"] <- input_variable
    }
    ProjectionObj <- if (!is.null(polyz_from)) {
      sf::st_crs(polyz_from)
    } else {
      sf::st_crs(pointz_from)
    }
    if (grepl("^EPSG", ProjectionObj$input) & !grepl("^\\+",
                                                     ProjectionObj$input)) {
      ProjectionObj$input <- tolower(paste0("+init=", ProjectionObj$input))
    }
    if (!is.null(polyz_from)){
      Template_Raster <- terra::rast(ncol = grid_dim[2], nrow = grid_dim[1], crs = ProjectionObj$input, extent=terra::ext(polyz_from))
    } else {
      Template_Raster <- terra::rast(ncol = grid_dim[2], nrow = grid_dim[1], crs = ProjectionObj$input, extent=terra::ext(pointz_from))
    }
    suppressWarnings({
      terra::values(Template_Raster) <- 1:terra::ncell(Template_Raster)
    })
    names(Template_Raster) <- "ID"
    impute_layer <- Template_Raster
    if (!is.null(pointz_from)) {
      pointz_from$Field <- unlist(terra::extract(x=Template_Raster,y=terra::vect(pointz_from),ID=FALSE))
      pointz_from_xy <- pointz_from
      sf::st_geometry(pointz_from_xy) <- NULL
      suppressMessages({
        pointz_from_xy <- dplyr::group_by(pointz_from_xy,
                                          Field)
        pointz_from_xy <- dplyr::summarize(pointz_from_xy,
                                           Var = aggregate_function[[1]](Var))
        pointz_from_xy <- as.data.frame(pointz_from_xy)
      })
      Template_Raster_xy <- terra::as.data.frame(Template_Raster,na.rm=FALSE)
      # head(Template_Raster_xy)
      Template_Raster_xy <- dplyr::left_join(Template_Raster_xy,
                                             pointz_from_xy, by = c(ID = "Field"))
      output_data <- Template_Raster
      terra::values(output_data) <- Template_Raster_xy$Var
      names(output_data) <- input_variable
      ReturnList <- list(return_output = output_data, return_point = pointz_from_SP_keep)
    }
    if (!is.null(polyz_from)) {
      polyz_from_layer <- terra::rasterize(polyz_from,
                                           Template_Raster, field = "Field")
      # terra::res(polyz_from_layer) <- c(grid_res[1], grid_res[2])
      # terra::plot(polyz_from_layer)
      field_layer <- polyz_from_layer
      names(field_layer) <- "Field"
      polyz_from_xy <- terra::as.data.frame(polyz_from_layer,na.rm=FALSE)
      names(polyz_from_xy)[1] <- "Field"
      polyz_from_df <- polyz_from
      sf::st_geometry(polyz_from_df) <- NULL
      dim(polyz_from_xy)
      dim(polyz_from_df)
      terra::ncell(polyz_from_layer)
      terra::ncell(Template_Raster)
      polyz_from_xy <- dplyr::left_join(polyz_from_xy,
                                        polyz_from_df, by = "Field")
      # polyz_from_xy <- base::merge(polyz_from_xy,polyz_from_df, by = "Field",all.x=TRUE,all.y=FALSE)

      dim(polyz_from_xy)
      names(polyz_from_xy)[names(polyz_from_xy) %in% c(input_variable)] <- "Var"
      terra::values(polyz_from_layer) <- polyz_from_xy$Var
      output_data <- polyz_from_layer
      names(output_data) <- input_variable
      suppressWarnings({
        centroid_polygon <- sf::st_centroid(polyz_from_SP_keep)
        centroid_polygon <- sf::st_buffer(centroid_polygon,
                                          dist = terra::res(Template_Raster)[1])
      })
      centroid_layer <- terra::mask(output_data, centroid_polygon)
      # terra::plot(centroid_layer)
      ReturnList <- list(return_output = output_data, return_centroid = centroid_layer,
                         poly_to = polyz_from_SP_keep, return_field = field_layer)
    }
    if (return_list == FALSE) {
      ReturnList <- ReturnList["return_output"][[1]]
    }
    return(ReturnList)
  }
  if (reverse) {
    if (message_out) {
      message("Converting regularly spaced raster grid cell into SF polygon object")
    }
    namesOutput <- names(return_output)
    Polygon_Frame <- c(return_field, return_output)
    names(Polygon_Frame)[1] <- "Field"
    names(Polygon_Frame)[-1] <- namesOutput
    Polygon_Frame <- terra::as.data.frame(Polygon_Frame)
    Polygon_Frame <- lapply(2:ncol(Polygon_Frame), function(col) {
      Selected_Matrix <- Polygon_Frame[, c(1, col)]
      names(Selected_Matrix)[2] <- "Var"
      suppressMessages({
        Selected_Matrix <- dplyr::group_by(Selected_Matrix,
                                           Field)
        Selected_Matrix <- dplyr::summarize(Selected_Matrix,
                                            Var = reverse_function[[1]](Var))
      })
      names(Selected_Matrix)[2] <- names(Polygon_Frame)[col]
      return(Selected_Matrix)
    })
    if (length(Polygon_Frame) > 1) {
      Polygon_Frame <- purrr::reduce(.x = Polygon_Frame,
                                     .f = full_join, by = "Field")
    }
    else {
      Polygon_Frame <- Polygon_Frame[[1]]
    }
    poly_to <- poly_to[, which(names(poly_to) %in% "Field")]
    poly_to <- dplyr::left_join(poly_to, Polygon_Frame, by = "Field")
    return(poly_to)
  }
}
