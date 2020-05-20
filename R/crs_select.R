#' Automatic planar coordinate reference system (CRS) selection
#'
#' Function to automatically convert simple features with geodetic coordinates (longitude, latitude / WGS 1984, EPSG:4326) to planar coordinates.
#'
#' @param polyz Polygon layer definiting the study region. \code{sf} object.
#' @param sf_layer Layer to be reprojected, if different from \code{polyz}. \code{sf} object.
#' @return List object, containing
#' \itemize{
##'  \item{"sf". }{The re-projected layer. \code{sf} object.}
##'  \item{"epsg_best". }{EPSG code of the projection. Character string.}
##'  }
#' @details Optimal map projection for the object \code{sf_layer} is defined as one that maximizes areal overlap between the study area and the spatial extent of planar coordinate reference systems and transformations in the EPSG Geodetic Parameter Dataset (polygons version 9.8, 2019-09-19).
#' @import sf data.table tidyverse
#' @importFrom stats dist
#' @importFrom utils read.delim
#' @importFrom dplyr select bind_rows
#' @importFrom rlang .data
#' @examples
#' # Find a planar projection for an unprojected (WSG 1984) hexagonal grid of Germany
#' \dontrun{
#' data(hex_05_deu)
#' hex_tr <- crs_select(hex_05_deu)
#' }
#' @export

crs_select <- function(polyz,sf_layer=polyz){

  # Convert to single polygon
  polyz_u <- polyz %>% st_union()

  # Find EPSG polygons that overlap with sf_layer
  suppressWarnings({
    epsg_candidates <- suppressMessages(
      list(sf_layer %>% st_centroid(),sf_layer)[c(grepl("POLYGON",sf_layer %>% st_geometry_type() %>% (function(.){.[1]})),!grepl("POLYGON",sf_layer %>% st_geometry_type() %>% (function(.){.[1]})))][[1]] %>% st_within(.,epsg_poly %>% st_transform(st_crs(sf_layer))) %>% as.data.table()
    )
  })
  epsg_candidates <- epsg_candidates[,col.id] %>% table() %>% sort(decreasing=F) %>% names() %>% as.numeric()

  # Area overlap between polyz and EPSG region
  epsg_mat <- lapply(seq_along(epsg_candidates),function(e0){
    area_1 <- polyz_u %>% st_area() %>% as.numeric()
    area_2 <- epsg_poly[epsg_candidates[e0],] %>% st_area() %>% as.numeric()
    area_ix <- suppressMessages(
      polyz_u %>% st_intersection(epsg_poly[epsg_candidates[e0],]) %>% st_area() %>% as.numeric()
    )
    data.frame(
      overlap_1 = area_ix/area_1,
      overlap_2 = area_ix/area_2,
      area_code = epsg_poly[epsg_candidates[e0],] %>% as.data.table() %>% dplyr::select(.data$AREA_CODE) %>% unlist(),
      area_name = epsg_poly[epsg_candidates[e0],] %>% as.data.table() %>% dplyr::select(.data$AREA_NAME) %>% unlist(),
      region = epsg_poly[epsg_candidates[e0],] %>% as.data.table() %>% dplyr::select(.data$REGION) %>% unlist(),
      stringsAsFactors = FALSE
    ) %>% as.data.table()
  }) %>% dplyr::bind_rows()
  epsg_mat <- epsg_mat[overlap_1<=1&overlap_2<=1,]
  epsg_dist <- epsg_mat[,c("overlap_1","overlap_2"),with=FALSE] %>% as.matrix() %>% rbind(c(1,1),.) %>% dist() %>% as.matrix()
  epsg_best <- epsg_mat[epsg_dist[-1,1] %>% order(),area_code]
  epsg_best <- lapply(epsg_best,function(x){st_crs(x)$epsg}) %>% unlist() %>% na.omit()

  # Find and assign optimal projection
  crs_bad <- TRUE; x0 <- 1
  while(crs_bad&x0<=length(epsg_best)){
    epsg_best_ <- epsg_best[x0]
    suppressWarnings({
      sf_layer_tr <- sf_layer %>% st_transform(crs=st_crs(paste0("EPSG:",epsg_best[x0])))
    })
    suppressWarnings({
      crs_bad <- (sf_layer %>% st_transform(crs=st_crs(paste0("EPSG:",epsg_best[x0]))) %>% st_bbox() %>% as.numeric() %>% is.na() %>% mean())==1 | st_is_longlat(sf_layer_tr)
    })
    x0 <- x0+1
  }

  # Output
  return(list(sf=sf_layer_tr,epsg_best=epsg_best_))
}
