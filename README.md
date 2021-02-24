# `SUNGEO` / Sub-National Geospatial Data Archive: Geoprocessing Toolkit
R package for integrating spatially-misaligned GIS datasets.

Version 0.2.21 (February 23, 2021)

Jason Byers, Marty Davidson, Yuri M. Zhukov

Center for Political Studies, Institute for Social Research

University of Michigan

Feedback, bug reports welcome: zhukov-at-umich-dot-edu

* `df2sf` / Convert data.frame object into simple features object
* `fix_geom` / Fix polygon geometries
* `geocode_osm` / Geocode addresses with OpenStreetMap
* `geocode_osm_batch` / Batch geocode addresses with OpenStreetMap
* `hot_spot` / Automatically calculate Local G hot spot intensity
* `line2poly` / Line-in-polygon analysis
* `point2poly_simp` / Point-to-polygon interpolation, simple overlay method
* `point2poly_tess` / Point-to-polygon interpolation, tessellation method
* `poly2poly_ap` / Area and population weighted polygon-to-polygon interpolation
* `sf2raster` / Convert simple features object into regularly spaced raster
* `utm_select` / Automatically convert geographic (degree) to planar coordinates (meters)


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
?utm_select
```

Example: geocode an address with OpenStreetMap

```
# Get geographic coordinates for the Big House (top match only)
geocode_osm("Michigan Stadium")
geocode_osm("Big House")

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

Example: point-to-polygon interpolation using tessellation method and area weights

```
# Load point-level election results
data(clea_deu2009_pt)

# Interpolate
out_4 <- point2poly_tess(pointz = clea_deu2009_pt,
                         polyz = hex_05_deu,
                         poly_id = "HEX_ID",
                         varz = "to1",
                         return_tess = TRUE)

# Visualize voter turnout at grid cell level 
plot(out_4$result["to1_aw"])

# Visualize Voronoi polygons used in estimation
plot(out_4$tess["to1"])
```

Example: line-in-polygon analysis

```
# Load road data (from Digital Chart of the World) and extract highways
data(highways_deu1992)

# Basic map overlay
plot(hex_05_deu["geometry"])
plot(highways_deu1992$geometry, add=TRUE, col = "blue", lwd=2)

# Calculate road lengths, densities and distances from each polygon to nearest highway
out_5 <- line2poly(linez = highways_deu1992,
                   polyz = hex_05_deu,
                   poly_id = "HEX_ID")
                   
# Visualize results
plot(out_5["line_length"])
plot(out_5["line_density"])
plot(out_5["line_distance"])

# Replace missing road lengths and densities with 0's, rename variables
out_6 <- line2poly(linez = highways_deu1992,
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
# Visualize original geometries (WGS1984, degrees)
plot(clea_deu2009["geometry"], axes=TRUE)

# Find a suitable CRS and re-project
out_7 <- utm_select(clea_deu2009)

# Visualize transformed geometries (UTM 32N, meters)
plot(out_7["geometry"], axes=TRUE)

# proj4string of transformed data
utm_select(clea_deu2009, return_list=TRUE)$proj_out
```

Example: Rasterization of polygons

```
# Transform sf polygon layer into 1km-by-1km RasterLayer (requires planar CRS)
out_8 <- sf2raster(polyz_from = utm_select(clea_deu2009),
                   input_variable = "to1")
# Visualize raster
raster::plot(out_8)

# 25km-by-25km RasterLayer (requires planar CRS)
out_9 <- sf2raster(polyz_from = utm_select(clea_deu2009),
                   input_variable = "to1",
                   grid_res = c(25000, 25000))
# Visualize raster
raster::plot(out_9)
```

Example: Create cartogram

```
# Cartogram of turnout scaled by number of valid votes
out_10 <- sf2raster(polyz_from = utm_select(clea_deu2009),
                  input_variable = "to1",
                  cartogram = TRUE,
                  carto_var = "vv1")
raster::plot(out_10)
```

Example: reverse rasterization

```
# Polygonization of cartogram raster
out_11a <- sf2raster(polyz_from = utm_select(clea_deu2009),
                    input_variable = "to1",
                    cartogram = TRUE,
                    carto_var = "vv1",
                    return_list = TRUE)
out_11 <- sf2raster(reverse = TRUE,
                   poly_to = out_11a$poly_to,
                   return_output = out_11a$return_output,
                   return_field = out_11a$return_field)
plot(out_11["to1"])

```

Additional examples in help files of individual functions.
