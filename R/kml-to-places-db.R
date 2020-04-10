# import data from kml to place
require(tidyverse)
require(sf)
require(sfc)
require(SFtools)
require(tidykml)
# devtools::install_github("briatte/tidykml")

## R functions written for this project
# projection
source("./code/crs_2193.R")

#Already added `.zip` and unzipped into folder called seed-locations
#in raw_data folder of data....
map.bounds <- tidykml::kml_bounds("./data/research/doc.kml")
map.info <- tidykml::kml_info("./data/research/doc.kml")
tidy.kat <- tidykml::kml_points("./data/research/doc.kml")
tidy.poly <- tidykml::kml_polygons("./data/research/doc.kml")
tidy.read <- tidykml::kml_read("./data/research/doc.kml")

full_tidy <- dplyr::bind_rows(tidy.kat, tidy.poly) %>%
  #need to doo but doesnt work currently
  # st_transform(crs = 4326) %>%
  mutate(x = latitude,
         y = longitude)

full_tidy_sf <- dplyr::bind_rows(tidy.kat, tidy.poly) %>%
  #need to doo but doesnt work currently
  # st_transform(crs = 4326) %>%
  mutate(x = latitude,
         y = longitude) %>%
  st_as_sf(coords = c("latitude", "longitude"), crs = 4326)

str(full_tidy_sf)
lat <- full_tidy$latitude
long <- full_tidy$longitude
# #need to doo but doesnt work currently
# # st_transform(crs = 4326) %>%
#   mutate(x = latitude,
#          y = longitude) %>%

seedling_polygon_range <- tidykml::kml_polygons("./data/research/doc.kml") %>%
  #need to doo but doesnt work currently
  # st_transform(crs = 4326) %>%
  mutate(x = latitude,
         y = longitude)
