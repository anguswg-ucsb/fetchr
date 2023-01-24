utils::globalVariables(c(".", "res", "cell_count", "diag_cell_count", "side_cell_count"))

#' Calculate wind fetch distances from spatial objects
#' @description This function will calculate the mean, maximum, minimum, or total fetch distance in 8 directions (horizontally, vertically, and diagonally). When provided raster data, the default behavior is to treat all NA cells as water cells and all non NA cells as land. If you want to specify a different value for water cells, use the water_value argument to specify the grid cells that represent water. When a sf/terra polygon is provided, a cell resolution must be provided to specify the desired resolution of the output raster. Returns a binary landwater SpatRaster with land values of 0 and water values of 1.Returns a raster of the same CRS, resolution, and extent as the input raster or vector.
#' @param r terra SpatRaster, raster RasterLayer, sf POLYGON/MULTIPOLYGON geometry, terra SpatVector POLYGON/MULTIPOLYGON geometry, or a character filepath pointing to a .shp, .gpkg, .geojson, or .tif file representing land and water. Provide a SpatRaster or RasterLayer should be a binary landwater raster with land values of 0 and water values of 1.
#' @param max_dist numeric, maximum distance (meters) to calculate fetch values out to. Default is 100,000m (100km)
#' @param res numeric, specifying the desired resolution of the output raster when providing a polygon. If a polygon is given, the polygon must be rasterized at a specific resolution. This argument is ignored if a raster/terra raster is given. Default is 2000, which will rasterize polygons onto a grid with a grid cell resolution of 2000 x 2000.
#' @param water_value numeric, value indicating the value of the water cells within the input raster. If a sf/terra polygon is given, this argument is ignored. Default is NA, such that NA cells are treated as water cells and all other cells are land cells.
#' @param func character R function to apply to distances from each cell. Either mean, min, max, or sum.
#' @param ncores numeric indicating how many cores to use during parallel processing. Default is 2.
#' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
#' @return SpatRaster of mean, max, min, or total wind fetch distances from land calculated for every water pixel.
#' @importFrom dplyr arrange group_by summarise left_join
#' @importFrom terra as.matrix rast nrow ncol crs ext plot
#' @export
#' @examples
#' \dontrun{
#' # mean fetch distance using a binary raster RasterLayer as the input
#' mean_fetch <- fetchr::get_fetch(
#' r      = fetchr::landwater
#' )
#' # plot mean fetch
#' plot(mean_fetch)
#'
#' # maximum fetch distance using a binary raster RasterLayer as the input
#' max_fetch <- fetchr::get_fetch(
#' r      = fetchr::landwater,
#' func   = "max"
#' )
#'
#' # plot max fetch
#' plot(max_fetch)
#'
#' # minimum fetch distance using a binary raster RasterLayer as the input
#'
#' min_fetch <- fetchr::get_fetch(
#' r      = fetchr::landwater,
#' func   = "min"
#' )
#'
#' # plot min fetch
#'
#' plot(min_fetch)
#'
#' # mean fetch distance using a sf MULTIPOLYGON object as the input
#' fetch_poly <- fetchr::get_fetch(
#' r = fetchr::land_vect
#' )
#' # plot mean fetch from a polygon object
#' plot(fetch_poly)
#'
#' # increase the maximum distance to calculate fetch lengths out to 200km
#' increase_max_dist <- fetchr::get_fetch(
#' r        = fetchr::landwater,
#' max_dist = 200000,
#' ncores   = 4                      # increase the number of cores used on computer
#' )
#'
#' # plot mean fetch going out 200km
#' plot(increase_max_dist)
#' }
get_fetch = function(
    r              = NULL,
    max_dist       = 100000,
    res            = 2000,
    water_value    = NA,
    func           = "mean",
    ncores         = 2,
    verbose        = TRUE
) {

  # error if no path or raster are given
  if(is.null(r)) {

    stop(r_error_msg())
  }

  # process input data to meet specifications for calculating fetch
  r <- get_landwater(
          r              = r,
          water_value    = water_value,
          res            = res,
          verbose        = FALSE
        )

  if(verbose == TRUE) {
    message(paste0("Calculating Fetch..."))
  }

  # cell resolution for distance calcs
  cell_res  <- terra::res(r)[1]

  # max count distance
  max_count <- floor(max_dist/(cell_res))

  # make landwater raster a wide matrix
  rmat      <- terra::as.matrix(r, wide = TRUE)

  # get indices of water points
  idx <-
    dplyr::arrange(
      as.data.frame(
        which(rmat == 1, arr.ind = T)
        ),
      row
      )
  # # get indices of water points
  # idx <-
  #   which(rmat == 1, arr.ind = T) %>%
  #   as.data.frame() %>%
  #   dplyr::arrange(row)

  # verbose = T
  if(verbose == TRUE) {
    message(paste0("Calculating north/south/east/west distances"))
  }

  # side cell totals list
  side_lst <- side_summary(
    m          = rmat,
    indices_df = idx,
    cell_res   = cell_res,
    max_dist   = max_dist,
    func       = "mean",
    ncores     = ncores,
    verbose    = T
  )

  if(verbose == TRUE) {
    message(paste0("Calculating diagonal distances"))
  }

  # diagonal totals list
  diag_lst <- diag_summary(
    m          = rmat,
    indices_df = idx,
    func       = "mean",
    max_count  = max_count,
    ncores     = ncores,
    verbose    = FALSE
  )

  # add mean diagonal cell counts
  idx$diag_cell_count <- unlist(diag_lst)

  # add mean diagonal cell counts
  idx$side_cell_count <- unlist(side_lst)

  # calculate cell distances from cell counts
  dist_mat <- cell_dist(
    indices_df = idx,
    cell_res   = cell_res,
    max_dist   = max_dist,
    func       = func,
    na_remove  = TRUE
  )

  # create fetch raster from gridded matrix
  fetch_r <-
    terra::rast(
      x = matrix(
        data = dplyr::arrange(
          dplyr::left_join(
            expand.grid(
              row = 1:terra::nrow(r),
              col = 1:terra::ncol(r)
            ),
            dist_mat[, c("row", "col", "distance")],
            by = c("row", "col")
          ),
          row
        )$distance,
        nrow  = terra::nrow(r),
        ncol  = terra::ncol(r),
        byrow = TRUE
      ),
      crs    = terra::crs(r),
      extent = terra::ext(r)
    )

  return(fetch_r)
}
