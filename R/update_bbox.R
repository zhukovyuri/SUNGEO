#' Update bounding box of sf object
#'
#' Function to update the coordinates of the bounding box of sf vector data objects (e.g. after cropping or subsetting).
#'
#' @param sfobj Layer to be updated. \code{sf} object.
#' @return \code{sf} object, with corrected bounds.
#' @importFrom sf st_coordinates st_geometry
#' @importFrom data.table data.table as.data.table
#' @examples
#' # Update bbox for subset of sf object
#' \dontrun{
#' data(clea_deu2009)
#' out_1 <- update_bbox(clea_deu2009[clea_deu2009$cst_n%in%c("Berlin"),])
#' out_1
#' 
#' # Bounding box of full dataset
#' data.table::as.data.table(clea_deu2009)[,sf::st_bbox(geometry)]
#' 
#' # Bounding box of subset (incorrect)
#' data.table::as.data.table(clea_deu2009)[cst_n%in%c("Berlin"),sf::st_bbox(geometry)]
#' 
#' # Corrected bounding box
#' data.table::as.data.table(out_1)[,sf::st_bbox(geometry)]
#' }
#' @export

update_bbox <- function(
  sfobj
  ){
  # Manually calculate bounds from coordinates
  new_bb <- data.table::as.data.table(sf::st_coordinates(sfobj))[,c(min(X),min(Y),max(X),max(Y))]
  # Rename columns
  names(new_bb) <- c("xmin", "ymin", "xmax", "ymax")
  # Change object class
  attr(new_bb, "class") <- "bbox"
  # Assign to bbox slot of sfobj
  attr(sf::st_geometry(sfobj), "bbox") <- new_bb

  return(sfobj)
}
