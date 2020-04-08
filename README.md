# `SUNGEO` / Sub-National Geospatial Data Archive System: Geoprocessing Toolkit
R package for integrating spatially-misaligned GIS datasets.

Version 0.1.0 (April 8, 2020)

* `fix_geom` / Function to check and fix broken geometries in simple features polygon objects
* `geo2planar` / Automatic planar coordinate reference system (CRS) selection
* `point2poly_krig` / Point-to-polygon interpolation, ordinary kriging method
* `point2poly_simp` / Point-to-polygon interpolation, simple overlay method
* `point2poly_tess` / Point-to-polygon interpolation, tessellation method
* `poly2poly_ap` / Area/population weighted polygon-to-polygon interpolation

Dependencies: `sf`, `data.table`, `tidyverse`

Very preliminary. Feedback, bug reports welcome: zhukov-at-umich-dot-edu

To install in R:

> library(devtools)

> devtools::install_github("zhukovyuri/SUNGEO")
