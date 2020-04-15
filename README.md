# `SUNGEO` / Sub-National Geospatial Data Archive System: Geoprocessing Toolkit
R package for integrating spatially-misaligned GIS datasets.

Version 0.1.0 (April 8, 2020)

Dependencies: `sf`, `data.table`, `tidyverse`

* `clean_geonames` / Function to clean GeoNames gazetteer files
* `crs_select` / Automatic planar coordinate reference system (CRS) selection
* `fix_geom` / Function to check and fix broken geometries in simple features polygon objects
* `geocode_osm` / Geocode addresses with OpenStreetMap
* `geocode_gn` / Geocode addresses with GeoNames
* `point2poly_krig` / Point-to-polygon interpolation, ordinary kriging method
* `point2poly_simp` / Point-to-polygon interpolation, simple overlay method
* `point2poly_tess` / Point-to-polygon interpolation, tessellation method
* `poly2poly_ap` / Area/population weighted polygon-to-polygon interpolation

Very preliminary. Feedback, bug reports welcome: zhukov-at-umich-dot-edu

To install in R:

```
library(devtools)
devtools::install_github("zhukovyuri/SUNGEO")
```

Load package:

```
library(SUNGEO)
```

Read help files:

```
?poly2poly_ap
?geo2planar
```

Example: geocode an address

```
# Get geographic coordinates for the Big House (top match only)
geocode_osm("Michigan Stadium")

# Return detailed results for top match
geocode_osm("Michigan Stadium", details=TRUE)

# Return detailed results for all matches
geocode_osm("Michigan Stadium", details=TRUE, return_all = TRUE)

```

Example: area-weighted polygon-to-polygon interpolation

```
# Load legislative election results (from CLEA)
data(clea_deu2009)

# Visualize voter turnout at constituency level
plot(clea_deu2009["to1"])

# Load 0.5 degree hexagonal grid
data(hex_05_deu)

# Interpolate
out_1 <- poly2poly_ap(poly_from = clea_deu2009,
                      poly_to = hex_05_deu,
                      poly_to_id = "HEX_ID",
                      varz = "to1"
                      )

# Visualize voter turnout at grid cell level
plot(out_1["to1_aw"])
```

Example: population-weighted polygon-to-polygon interpolation

```
# Load population raster (from GPW v4)
data(gpw4_deu2010)

# Interpolate
out_2 <- poly2poly_ap(poly_from = clea_deu2009,
                      poly_to = hex_05_deu,
                      poly_to_id = "HEX_ID",
                      varz = "to1",
                      methodz = "pw",
                      pop_raster = gpw4_deu2010)

# Visualize voter turnout at grid cell level
plot(out_2["to1_pw"])
```

Example: point-to-polygon interpolation using ordinary kriging (with automatic variogram model selection)

```
# Load point-level legislative election results (from CLEA)
data(clea_deu2009_pt)

# Interpolate
out_3 <- point2poly_krig(pointz = clea_deu2009_pt,
                         polyz = hex_05_deu,
                         varz = "to1")

# Visualize predicted voter turnout values
plot(out_3["to1_kr"])

# Visualize standard deviation of predicted values
plot(out_3["to1_kr_sd"])
```

Example: point-to-polygon interpolation using tessellation method and area weights

```
# Interpolate
out_4 <- point2poly_tess(pointz = clea_deu2009_pt,
                         polyz = hex_05_deu,
                         poly_id = "HEX_ID",
                         varz = "to1",
                         return_tess = TRUE)

# Visualize voter turnout at grid cell level 
plot(out_4$result["to1_aw"]))

# Visualize Voronoi polygons used in estimation
plot(out_4$tess["to1"])
```

Example: Automatically find a planar CRS for a GIS dataset

```
# Find a suitable CRS and re-project
clea_tr <- geo2planar(clea_deu2009)

# EPSG code of transformed data
clea_tr$epsg_best

# Visualize transformed geometries
plot(clea_tr$sf$geometry)
```

Additional examples in help files of individual functions.
