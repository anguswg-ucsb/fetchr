library(AOI)
library(sf)
library(terra)
library(raster)
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
land <- terra::rasterize(
  x = st,
  r_grid
) %>%
  terra::crop(bb) %>%
  raster::raster()

usethis::use_data(land, overwrite = TRUE)
