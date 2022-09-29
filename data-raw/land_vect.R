## code to prepare `land_vect` dataset goes here
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

# bounding box to crop to
bb <- terra::ext(c(-2317958.20937478, -1973298.78204874, 1309390.77983874, 1623941.40966438))

# # rasterize land polygon
land_vect <- terra::rasterize(
  x = st,
  make_grid(
    r   = st,
    res = 2000
  )
) %>%
  terra::crop(bb) %>%
  terra::as.polygons()

usethis::use_data(land_vect, overwrite = TRUE)
