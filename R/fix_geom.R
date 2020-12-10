#' Fix polygon geometries
#'
#' Function to check validity and fix broken geometries in simple features polygon objects
#'
#' @param x Polygon layer to be checked and fixed. \code{sf} object.
#' @param n_it Number of iterations. Default is 10. Numeric..
#' @return Returns a \code{sf} polygon object, with self-intersections and other geometry problems fixed.
#' @importFrom sf st_make_valid st_is_valid st_buffer st_cast
#' @examples
#' # Assignment of a single variable (sums)
#' \dontrun{
#' data(clea_deu2009)
#' out_1 <- fix_geom(clea_deu2009)
#' }
#' @export

fix_geom <- function(x, n_it = 10){
  #Part 1 -
  ValidityVector <- sf::st_is_valid(x)

  #Part 2 -
  countIter <- 0

  while(all(ValidityVector) == F || countIter < n_it){
    suppressMessages(
      x <- sf::st_buffer(x,dist=0)
    )

    ValidityVector <- sf::st_is_valid(x)

    countIter <- 1 + countIter

    if(all(ValidityVector)){
      break
    }

    if(all(ValidityVector) == F && countIter == n_it){
      stop('Non-Convergence for Polygon; Please increase number of iterations or check to see whether polygon is corrupted')
    }

  }

  message(paste0('Completed ', countIter, ' Round of Buffering'))

  # Simplify geometry type
  x <- sf::st_cast(x)

  return(x)
}
