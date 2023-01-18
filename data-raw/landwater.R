## code to prepare `landwater` dataset goes here
#' @title Binary landwater raster
#' @description  Example binary landwater raster. Values of 1 indicate land and values of 0 indicate water.
#' @importFrom AOI aoi_get
#' @importFrom terra vect project ext rast crs rasterize crop setValues values
#' @importFrom raster raster
#' @export

# state of interest
st <-
  AOI::aoi_get(state = "CA") %>%
  terra::vect() %>%
  terra::project("epsg:5070")

# # Socal/Santa Barbara bounding box extent to crop to
bb <- terra::ext(c(-2338958.20937478, -1973298.78204874, 1245390.77983874, 1623941.40966438))
# bb <- terra::ext(c(-2317958.20937478, -1973298.78204874, 1309390.77983874, 1623941.40966438))

# Template raster of HSI grid
r_grid <- terra::rast(
  crs         = terra::crs(st),
  extent      = terra::ext(st),
  resolution  = c(2000, 2000),
  vals        = 1
)

# rasterize land polygon
r <- terra::rasterize(
  x = st,
  r_grid
) %>%
  terra::crop(bb)

# set all water cells to a given value and all other cells to another. Presumes raster is just land and water binary.
landwater <- terra::setValues(
    r, ifelse(is.na(terra::values(r)), 1, 0)
  ) %>%
  raster::raster()

usethis::use_data(landwater, overwrite = TRUE)

