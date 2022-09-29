## code to prepare `landwater` dataset goes here
library(dplyr)
library(USAboundaries)
library(USAboundariesData)
library(terra)
library(sf)

# state of interest
st <-
  USAboundaries::us_states() %>%
  dplyr::filter(state_abbr == "CA") %>%
  terra::vect() %>%
  terra::project("epsg:5070")

# # Socal/Santa Barbara bounding box extent to crop to
bb <- terra::ext(c(-2317958.20937478, -1973298.78204874, 1309390.77983874, 1623941.40966438))

# landwater raster
landwater  <- get_landwater(
  r              = st,
  res            = 2000,
  verbose        = TRUE
) %>%
  terra::crop(bb)

usethis::use_data(landwater, overwrite = TRUE)
