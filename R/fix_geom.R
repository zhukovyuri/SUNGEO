#' Polygon geometry fix
#'
#' Function to check and fix broken geometries in simple features polygon objects
#'
#' @param x Polygon layer to be checked and fixed. \code{sf} object.
#' @param self_int Look only for self-intersections? Logical.
#' @return Returns a \code{sf} polygon object, with self-intersections and other geometry problems fixed.
#' @import sf tidyverse
#' @examples
#' # Assignment of a single variable (sums)
#' \dontrun{
#' data(clea_deu2009)
#' out_1 <- fix_geom(clea_deu2009)
#' }
#' @export

fix_geom <- function(x,self_int=TRUE){
  if(self_int){
    if(sum(grepl("Self-inter",st_is_valid(x,reason=T)))>0){
      suppressMessages(
        x <- x %>% st_buffer(dist=0)
      )
    }}
  if(!self_int){
    if(sum(!st_is_valid(x))>0){
      suppressMessages(
        x <- x %>% st_buffer(dist=0)
      )
    }}
  x
}
