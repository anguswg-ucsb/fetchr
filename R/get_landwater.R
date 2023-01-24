#' Prepare a spatial object for fetch calculations
#' @description Provides a method for preprocessing sf/terra polygons or raster/terra rasters for calculating fetch using get_fetch(). When provided raster data, the default behavior is to treat all NA cells as water cells and all non NA cells as land. If you want to specify a different value for water cells, use the water_value argument to specify the grid cells that represent water. Returns a binary landwater SpatRaster with land values of 0 and water values of 1.
#'  A water value should be specified to indicate which pixel value should represent water in the final raster. The rest of the pixels will indicate land cells.
#' @param r sf POLYGON/MULTIPOLYGON, terra SpatVect POLYGON, SpatRaster, RasterLayer, or a file path (character) to a vector/raster data file. Data must be on a projected CRS.
#' @param water_value numeric, value indicating the value of the water cells within the input raster. If a sf/terra polygon is given, this argument is ignored. Default is NA, such that NA cells are treated as water cells and all other cells are land cells.
#' @param res numeric resolution for template raster. If a polygon is given, the polygon must be rasterized at a specific resolution. This argument is ignored if a raster/terra raster is given. Default is 2000, which will rasterize polygons onto a grid with a grid cell resolution of 2000 x 2000.
#' @param verbose logical, whether messages should be printed. Default is TRUE, so messages will print.
#' @importFrom methods is
#' @return SpatRaster, binary landwater raster with land values of 0 and water values of 1
#' @export
#' @examples
#' \dontrun{
#' # Create a binary landwater SpatRaster from a sf polygon at a 2000m x 2000m resolution
#' lw <- fetchr::get_landwater(
#'                 r   = fetchr::land_vect,
#'                 res = 2000
#'                 )
#' }
get_landwater <- function(
    r,
    water_value    = NA,
    res            = 2000,
    verbose        = TRUE
) {

  if(verbose == TRUE) {
    message(paste0("Creating landwater raster..."))
  }

  # if none of the valid 'r' arguments are provided
  if(!any(methods::is(r, "character"), methods::is(r, "sf"),
          methods::is(r, "RasterLayer"), methods::is(r, "SpatRaster"),
          methods::is(r, "SpatVector"))) {

    stop(
      r_error_msg()
    )

  }
  # if r is a character
  if(methods::is(r, "character")) {

    # load and prep file path data
    r <- prep_filepath(
      r           = r,
      water_value = water_value,
      res         = res,
      verbose     = FALSE
    )

    return(r)

  }
  # if terra SpatVector or sf object are given
  if(any(methods::is(r, "sf"),  methods::is(r, "SpatVector"))) {


    r <- prep_poly(
          r       = r,
          res     = res,
          verbose = FALSE
          )

    return(r)

  }

  # if terra SpatRaster or raster RasterLayer are given
  if(any(methods::is(r, "RasterLayer"), methods::is(r, "SpatRaster"))) {

    # convert RasterLayer to SpatRaster
   r <- prep_raster(
         r           = r,
         water_value = water_value,
         verbose     = FALSE
         )

   return(r)

  }

}

#' Make a land raster from binary land type raster
#' @description returns a SpatRaster of only land cells. Land cells have a value of 2 by default, all other cells are NA. Fits the HSI raster specifications.
#' @param r SpatRast of binary land type raster. Land and water values only.
#' @param land_value numeric indicating what value land is represented by in raster. Default is 0 to conform to get_landwater() output
#' @param out_value numeric indicating what value land should be represented by in output raster. Default is 2
#' @param as_polygon logical, should polygon be returned instead of a raster. Default is FALSE, returns a terra rast object
#' @param as_sf logical, should polygon be returned as an SF object. Default is FALSE, returns a terra vect object
#' @importFrom terra setValues values
#' @return SpatRaster with only land cells having a value, all other cells are NA
# get_land = function(
#     r,
#     land_value = 0,
#     out_value  = 2,
#     as_polygon = FALSE,
#     as_sf      = FALSE
# ) {
#
#   # if cell value is 1 (land) make it a value of 2, otherwise make cell NA
#   land <- terra::setValues(
#     r, ifelse(terra::values(r) == land_value, out_value, NA)
#   )
#
#   # if a polygon should be returned instead of a raster
#   if(as_polygon == TRUE) {
#
#     # if polygon should be returned as an SF object
#     if(as_sf == TRUE) {
#
#       message("Generating land polygon...")
#       message("Type: SF")
#
#       # SF land polygon
#       land <-
#         land %>%
#         terra::as.polygons() %>%
#         sf::st_as_sf()  %>%
#         sf::st_make_valid() %>%
#         sf::st_cast() %>%
#         stats::setNames(c("land", "geometry")) %>%
#         sf::st_transform(26915)
#
#       # return as terra vect object
#     } else {
#
#       message("Generating land polygon...")
#       message("Type: SpatVector")
#
#       # spatvector land polygon
#       land <-
#         land %>%
#         terra::as.polygons() %>%
#         stats::setNames(c("land"))
#
#     }
#
#
#   } else {
#
#     message("Generating land raster...")
#     message("Type: SpatRaster")
#
#   }
#
#   return(land)
#
#   # view map
#   # land %>%
#   #   raster::raster()%>%
#   #   mapview::mapview()
#
# }

#' Make a water raster from binary land type raster
#' @description returns a SpatRaster of only water cells. Water cells have a value of 1 by default, all other cells are NA. Fits the HSI raster specifications.
#' @param r SpatRast of binary land type raster. Land and water values only.
#' @param water_value numeric indicating what value water is represented by in raster. Default is 1 to conform to get_landwater() output
#' @param out_value numeric indicating what value land should be represented by in output raster. Default is 1
#' @param as_polygon logical, should polygon be returned instead of a raster. Default is FALSE, returns a terra rast object
#' @param as_sf logical, should polygon be returned as an SF object. Default is FALSE, returns a terra vect object
#' @importFrom terra setValues values
#' @return SpatRaster with only land cells having a value, all other cells are NA
# get_water = function(
#     r,
#     water_value = 1,
#     out_value   = 1,
#     as_polygon  = FALSE,
#     as_sf       = FALSE
# ) {
#
#   # if cell value is 1 (land) make it a value of 2, otherwise make cell NA
#   water <- terra::setValues(
#     r, ifelse(terra::values(r) == water_value, out_value, NA)
#   )
#
#   # if a polygon should be returned instead of a raster
#   if(as_polygon == TRUE) {
#
#     # if polygon should be returned as an SF object
#     if(as_sf == TRUE) {
#
#       # make polygon from water raster
#       water <-
#         water %>%
#         terra::as.polygons() %>%
#         sf::st_as_sf()  %>%
#         sf::st_make_valid() %>%
#         sf::st_cast() %>%
#         stats::setNames(c("water", "geometry")) %>%
#         sf::st_transform(26915)
#
#       # return as terra vect object
#     } else {
#
#       water <-
#         water %>%
#         terra::as.polygons() %>%
#         stats::setNames(c("water"))
#
#     }
#   }
#
#   return(water)
#
# }
#
