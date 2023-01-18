library(AOI)
library(sf)
library(terra)
library(usethis)

# state of interest
st <-
  AOI::aoi_get(state = "CA") %>%
  terra::vect() %>%
  terra::project("epsg:5070")

# bounding box to crop to
bb <- terra::ext(c(-2338958.20937478, -1973298.78204874, 1245390.77983874, 1623941.40966438))
# bb <- terra::ext(c(-2317958.20937478, -1973298.78204874, 1309390.77983874, 1623941.40966438))

# Template raster of HSI grid
r_grid <- terra::rast(
  crs         = terra::crs(st),
  extent      = terra::ext(st),
  resolution  = c(2000, 2000),
  vals        = 1
)

# # rasterize land polygon
land_vect <- terra::rasterize(
  x = st,
  r_grid
) %>%
  terra::crop(bb) %>%
  terra::as.polygons() %>%
  sf::st_as_sf()

usethis::use_data(land_vect, overwrite = TRUE)
