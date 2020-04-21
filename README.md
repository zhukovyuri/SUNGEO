# `SUNGEO` / Sub-National Geospatial Data Archive System: Geoprocessing Toolkit
R package for integrating spatially-misaligned GIS datasets.

Version 0.1.0 (April 15, 2020)

Very preliminary. Feedback, bug reports welcome: zhukov-at-umich-dot-edu

* `clean_geonames` / Function to clean GeoNames gazetteer files
* `crs_select` / Automatic planar coordinate reference system (CRS) selection
* `fix_geom` / Function to check and fix broken geometries in simple features polygon objects
* `geocode_gn` / Batch geocode addresses with GeoNames
* `geocode_osm` / Geocode addresses with OpenStreetMap
* `geocode_osm_batch` / Batch geocode addresses with OpenStreetMap
* `line2poly` / Line-in-polygon analysis
* `point2poly_krig` / Point-to-polygon interpolation, ordinary kriging method
* `point2poly_simp` / Point-to-polygon interpolation, simple overlay method
* `point2poly_tess` / Point-to-polygon interpolation, tessellation method
* `poly2poly_ap` / Area/population weighted polygon-to-polygon interpolation

To install in R:

```
library(devtools)
devtools::install_github("zhukovyuri/SUNGEO", dependencies = TRUE)
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

Example: geocode an address with OpenStreetMap

```
# Get geographic coordinates for the Big House (top match only)
geocode_osm("Michigan Stadium")

# Return detailed results for top match
geocode_osm("Michigan Stadium", details=TRUE)

# Return detailed results for all matches
geocode_osm("Michigan Stadium", details=TRUE, return_all = TRUE)

```

Example: geocode multiple addresses

```
# Geocode multiple addresses (top matches only)
geocode_osm_batch(c("Ann Arbor","East Lansing","Columbus"))

# ... with progress reports
geocode_osm_batch(c("Ann Arbor","East Lansing","Columbus"), 
                  verbose = TRUE)

# Return detailed results for all matches
geocode_osm_batch(c("Ann Arbor","East Lansing","Columbus"),
                  details = TRUE, return_all = TRUE)

```

Example: geocode addresses with GeoNames

```
# Geocode an address
geocode_gn("Chisinau", country_name = "Moldova")

# Return detailed results
geocode_gn("Chisinau", country_name = "Moldova", details = TRUE)

# Return detailed results for multiple addresses, with progress reports
geocode_gn(query = c("Chisinau","Buiucani, Chisinau","Chisinau centru"),
           country_name = "Moldova", details = TRUE, verbose = TRUE)
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

Example: line-in-polygon analysis

```
# Load road data (from Digital Chart of the World) and extract highways
data(dcwroad_deu1992)
highways <- dcwroad_deu1992[dcwroad_deu1992$MED_DESCRI%in%"With Median",]

# Basic map overlay
plot(hex_05_deu$geometry)
plot(highways$geometry, add=TRUE, col = "blue", lwd=2)

# Calculate road lengths, densities and distances from each polygon to nearest highway
out_5 <- line2poly(linez = highways,
                   polyz = hex_05_deu,
                   poly_id = "HEX_ID")
                   
# Visualize results
plot(out_5["line_length"])
plot(out_5["line_density"])
plot(out_5["line_distance"])

# Replace missing road lengths and densities with 0's, rename variables
out_6 <- line2poly(linez = highways,
                   polyz = hex_05_deu,
                   poly_id = "HEX_ID",
                   outvar_name = "road",
                   na_val = 0)

# Visualize results
plot(out_6["road_length"])
plot(out_6["road_density"])
plot(out_6["road_distance"])
```

Example: Automatically find a planar CRS for a GIS dataset

```
# Find a suitable CRS and re-project
clea_tr <- crs_select(clea_deu2009)

# EPSG code of transformed data
clea_tr$epsg_best

# Visualize transformed geometries
plot(clea_tr$sf$geometry)
```

Additional examples in help files of individual functions.
