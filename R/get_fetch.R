utils::globalVariables(c(".", "res", "cell_count", "diag_cell_count", "side_cell_count"))

#' Calculate wind fetch distances from spatial objects
#' @description This function will calculate the mean, maximum, minimum, or total fetch distance in 8 directions (horizontally, vertically, and diagonally). When provided raster data, the default behavior is to treat all NA cells as water cells and all non NA cells as land. If you want to specify a different value for water cells, use the water_value argument to specify the grid cells that represent water. When a sf/terra polygon is provided, a cell resolution must be provided to specify the desired resolution of the output raster. Returns a binary landwater SpatRaster with land values of 0 and water values of 1.Returns a raster of the same CRS, resolution, and extent as the input raster or vector.
#' @param r terra SpatRaster, raster RasterLayer, sf POLYGON/MULTIPOLYGON geometry, terra SpatVector POLYGON/MULTIPOLYGON geometry, or a character filepath pointing to a .shp, .gpkg, .geojson, or .tif file representing land and water. Provide a SpatRaster or RasterLayer should be a binary landwater raster with land values of 0 and water values of 1.
#' @param max_dist numeric, maximum distance (meters) to calculate fetch values out to. Default is 100,000m (100km)
#' @param res numeric, specifying the desired resolution of the output raster when providing a polygon. If a polygon is given, the polygon must be rasterized at a specific resolution. This argument is ignored if a raster/terra raster is given. Default is 2000, which will rasterize polygons onto a grid with a grid cell resolution of 2000 x 2000.
#' @param water_value numeric, value indicating the value of the water cells within the input raster. If a sf/terra polygon is given, this argument is ignored. Default is NA, such that NA cells are treated as water cells and all other cells are land cells.
#' @param func character R function to apply to distances from each cell. Either mean, min, max, or sum.
#' @param in_parallel logical whether calculations should be done using parallel workers. If TRUE, (default), processes are run in parallel.
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
#' r           = fetchr::landwater,
#' max_dist    = 200000,
#' in_parallel = TRUE                 # Run in parallel
#' )
#'
#' # plot mean fetch going out 200km
#' plot(increase_max_dist)
#' }
get_fetch <- function(
    r              = NULL,
    max_dist       = 100000,
    res            = 2000,
    water_value    = NA,
    func           = "mean",
    in_parallel    = TRUE,
    verbose        = TRUE
) {

  # r      = fetchr::landwater %>% terra::unwrap()
  # r <- AOI::aoi_get(state = "Hawaii") %>%
  #   sf::st_transform(507)
  # func   = "mean"
  # max_dist       = 100000
  # res            = 2000
  # water_value    = NA
  # verbose        = TRUE
  # in_parallel    = TRUE

  # error if no path or raster are given
  if(is.null(r)) {

    stop(r_error_msg())
  }

  # r      = fetchr::landwater %>% terra::unwrap()
  # # func   = "max"
  # max_dist       = 100000
  # res            = 2000
  # water_value    = NA
  # func           = "mean"
  # in_parallel    = TRUE
  # verbose        = TRUE

  if(!func %in% c("mean", "max", "min", "sum")) {

    stop("Invalid input to 'func' argument: '", func, "', must be 'mean', 'max', 'min', or 'sum'")

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
  rmat      <- terra::as.matrix(
    r,
    wide = TRUE
  )

  # get indices of water points
  idx <-
    dplyr::arrange(
      as.data.frame(
        which(rmat == 1, arr.ind = T)
      ),
      row
    )

  # system.time({

  # get summarized distance for each cell in matrix
  dists <- do_cell_distances(
    grid        = rmat,
    index_df    = idx,
    func        = func,
    max_count   = max_count,
    in_parallel = in_parallel
  )

  # })

  # unlist result
  dists <- unlist(dists)

  # # max_count, func, in_parallel
  # # match function to R function
  # func <- match.fun(func)
  # # get average distance for each cell in matrix
  # dists <-
  #   lapply(1:nrow(idx), function(i) {
  #     func(
  #       get_cell_distances(
  #         grid = rmat,
  #         row  =   idx$row[i],
  #         col  =   idx$col[i],
  #         search_value = 0,
  #         max_count = max_count
  #       ),
  #       na.rm = TRUE
  #       )
  #   })
  #
  # # unlist list
  # dists <- unlist(dists)

  # add distances back to index dataframe
  idx$distance <- dists * cell_res

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
            idx[, c("row", "col", "distance")],
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

  # plot(fetch_r)

  return(fetch_r)

  }

#' Calculate wind fetch distances in 8 directions from spatial objects
#' @description This function will calculate fetch values for 8 directions (horizontally, vertically, and diagonally). The function will create an 8 layer SpatRaster with a raster layer for north, south, east, west, northeast, southeast, southwest, and northwest fetch values. When provided raster data, the default behavior is to treat all NA cells as water cells and all non NA cells as land. If you want to specify a different value for water cells, use the water_value argument to specify the grid cells that represent water. When a sf/terra polygon is provided, a cell resolution must be provided to specify the desired resolution of the output raster. Returns a binary landwater SpatRaster with land values of 0 and water values of 1.Returns a raster of the same CRS, resolution, and extent as the input raster or vector.
#' @param r terra SpatRaster, raster RasterLayer, sf POLYGON/MULTIPOLYGON geometry, terra SpatVector POLYGON/MULTIPOLYGON geometry, or a character filepath pointing to a .shp, .gpkg, .geojson, or .tif file representing land and water. Provide a SpatRaster or RasterLayer should be a binary landwater raster with land values of 0 and water values of 1.
#' @param max_dist numeric, maximum distance (meters) to calculate fetch values out to. Default is 100,000m (100km)
#' @param res numeric, specifying the desired resolution of the output raster when providing a polygon. If a polygon is given, the polygon must be rasterized at a specific resolution. This argument is ignored if a raster/terra raster is given. Default is 2000, which will rasterize polygons onto a grid with a grid cell resolution of 2000 x 2000.
#' @param water_value numeric, value indicating the value of the water cells within the input raster. If a sf/terra polygon is given, this argument is ignored. Default is NA, such that NA cells are treated as water cells and all other cells are land cells.
#' @param in_parallel logical whether calculations should be done using parallel workers. If TRUE, (default), processes are run in parallel.
#' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
#' @return SpatRaster of mean, max, min, or total wind fetch distances from land calculated for every water pixel.
#' @importFrom dplyr arrange group_by summarise left_join
#' @importFrom terra as.matrix rast nrow ncol crs ext plot
#' @export
#' @examples
#' \dontrun{
#' # calculate 8 fetch directions using a binary raster RasterLayer as the input
#' fetch_dirs <- fetchr::get_fetch_directions(
#' r      = fetchr::landwater
#' )
#' # plot 8 fetch directions fetch
#' plot(fetch_dirs)

#' # calculate 8 fetch directions using a sf MULTIPOLYGON object as the input
#' fetch_poly <- fetchr::get_fetch_directions(
#' r = fetchr::land_vect
#' )
#' # plot mean fetch from a polygon object
#' plot(fetch_poly)
#' }
get_fetch_directions <- function(
    r              = NULL,
    max_dist       = 100000,
    res            = 2000,
    water_value    = NA,
    in_parallel    = TRUE,
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
  # max_count = 100

  # make landwater raster a wide matrix
  rmat      <- terra::as.matrix(
    r,
    wide = TRUE
  )

  # get indices of water points
  idx <-
    dplyr::arrange(
      as.data.frame(
        which(rmat == 1, arr.ind = T)
      ),
      row
    )

  # get distances for each direction
  dists <- do_dir_distances(
    grid        = rmat,
    index_df    = idx,
    max_count   = max_count,
    in_parallel = in_parallel
  )

  # bind into dataframe as rows
  dists <- dplyr::bind_rows(dists)

  # system.time({
  # # get average distance for each cell in matrix
  # dists <-
  #   lapply(1:nrow(idx), function(i) {
  #     # i = 900
  #     # i = 335
  #
  #     get_direction_distances(
  #       grid = rmat,
  #       row  =   idx$row[i],
  #       col  =   idx$col[i],
  #       search_value = 0,
  #       max_count = max_count
  #     )
  #   })
  #   })
  # # bind distances into a dataframe
  # dists <- dplyr::bind_rows(dists)
  #
  # })

  # loop through each of the columns in the dists dataframe and create a raster with each directions distance values
  res <- lapply(1:ncol(dists), function(i) {
    # message(i, "/", ncol(dists))

    # use for raster layer name
    direction <- names(dists[,i])

    # calculate distance using cell counts
    idx$distance <- as.integer(dists[[direction]]) * cell_res

    # create fetch raster from gridded matrix
    fr <-
      terra::rast(
        x = matrix(
          data = dplyr::arrange(
            dplyr::left_join(
              expand.grid(
                row = 1:terra::nrow(r),
                col = 1:terra::ncol(r)
              ),
              idx[, c("row", "col", "distance")],
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

    # set direction name
    terra::set.names(fr, direction)

    fr
  })

  # put list of rasters into a single SpatRaster
  res <- terra::rast(res)

  return(res)
}

#' #' Calculate wind fetch distances from spatial objects
#' #' @description This function will calculate the mean, maximum, minimum, or total fetch distance in 8 directions (horizontally, vertically, and diagonally). When provided raster data, the default behavior is to treat all NA cells as water cells and all non NA cells as land. If you want to specify a different value for water cells, use the water_value argument to specify the grid cells that represent water. When a sf/terra polygon is provided, a cell resolution must be provided to specify the desired resolution of the output raster. Returns a binary landwater SpatRaster with land values of 0 and water values of 1.Returns a raster of the same CRS, resolution, and extent as the input raster or vector.
#' #' @param r terra SpatRaster, raster RasterLayer, sf POLYGON/MULTIPOLYGON geometry, terra SpatVector POLYGON/MULTIPOLYGON geometry, or a character filepath pointing to a .shp, .gpkg, .geojson, or .tif file representing land and water. Provide a SpatRaster or RasterLayer should be a binary landwater raster with land values of 0 and water values of 1.
#' #' @param max_dist numeric, maximum distance (meters) to calculate fetch values out to. Default is 100,000m (100km)
#' #' @param res numeric, specifying the desired resolution of the output raster when providing a polygon. If a polygon is given, the polygon must be rasterized at a specific resolution. This argument is ignored if a raster/terra raster is given. Default is 2000, which will rasterize polygons onto a grid with a grid cell resolution of 2000 x 2000.
#' #' @param water_value numeric, value indicating the value of the water cells within the input raster. If a sf/terra polygon is given, this argument is ignored. Default is NA, such that NA cells are treated as water cells and all other cells are land cells.
#' #' @param func character R function to apply to distances from each cell. Either mean, min, max, or sum.
#' #' @param ncores numeric indicating how many cores to use during parallel processing. Default is 2.
#' #' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
#' #' @return SpatRaster of mean, max, min, or total wind fetch distances from land calculated for every water pixel.
#' #' @importFrom dplyr arrange group_by summarise left_join
#' #' @importFrom terra as.matrix rast nrow ncol crs ext plot
#' #' @examples
#' #' \dontrun{
#' #' # mean fetch distance using a binary raster RasterLayer as the input
#' #' mean_fetch <- fetchr::get_fetch(
#' #' r      = fetchr::landwater
#' #' )
#' #' # plot mean fetch
#' #' plot(mean_fetch)
#' #'
#' #' # maximum fetch distance using a binary raster RasterLayer as the input
#' #' max_fetch <- fetchr::get_fetch(
#' #' r      = fetchr::landwater,
#' #' func   = "max"
#' #' )
#' #'
#' #' # plot max fetch
#' #' plot(max_fetch)
#' #'
#' #' # minimum fetch distance using a binary raster RasterLayer as the input
#' #'
#' #' min_fetch <- fetchr::get_fetch(
#' #' r      = fetchr::landwater,
#' #' func   = "min"
#' #' )
#' #'
#' #' # plot min fetch
#' #'
#' #' plot(min_fetch)
#' #'
#' #' # mean fetch distance using a sf MULTIPOLYGON object as the input
#' #' fetch_poly <- fetchr::get_fetch(
#' #' r = fetchr::land_vect
#' #' )
#' #' # plot mean fetch from a polygon object
#' #' plot(fetch_poly)
#' #'
#' #' # increase the maximum distance to calculate fetch lengths out to 200km
#' #' increase_max_dist <- fetchr::get_fetch(
#' #' r        = fetchr::landwater,
#' #' max_dist = 200000,
#' #' ncores   = 4                      # increase the number of cores used on computer
#' #' )
#' #'
#' #' # plot mean fetch going out 200km
#' #' plot(increase_max_dist)
#' #' }
#' get_fetch = function(
#'     r              = NULL,
#'     max_dist       = 100000,
#'     res            = 2000,
#'     water_value    = NA,
#'     func           = "mean",
#'     ncores         = 2,
#'     verbose        = TRUE
#' ) {
#'
#'   # r      = fetchr::landwater %>% terra::unwrap()
#'   # # func   = "max"
#'   # max_dist       = 100000
#'   # res            = 2000
#'   # water_value    = NA
#'   # func           = "mean"
#'   # ncores         = 0
#'   # verbose        = TRUE
#'   # xmin = -2339582
#'   # xmax = -2156582
#'   # ymin = 1435369
#'   # ymax = 1624369
#'   #
#'   # cext <- ext(xmin, xmax, ymin, ymax)
#'   # r <- terra::crop(r, cext)
#'
#'   # error if no path or raster are given
#'   if(is.null(r)) {
#'
#'     stop(r_error_msg())
#'   }
#'
#'   # process input data to meet specifications for calculating fetch
#'   r <- get_landwater(
#'           r              = r,
#'           water_value    = water_value,
#'           res            = res,
#'           verbose        = FALSE
#'         )
#'
#'   if(verbose == TRUE) {
#'     message(paste0("Calculating Fetch..."))
#'   }
#'
#'   # cell resolution for distance calcs
#'   cell_res  <- terra::res(r)[1]
#'
#'   # max count distance
#'   max_count <- floor(max_dist/(cell_res))
#'
#'   # make landwater raster a wide matrix
#'   rmat      <- terra::as.matrix(r, wide = TRUE)
#'
#'   # get indices of water points
#'   idx <-
#'     dplyr::arrange(
#'       as.data.frame(
#'         which(rmat == 1, arr.ind = T)
#'         ),
#'       row
#'       )
#'   # # get indices of water points
#'   # idx <-
#'   #   which(rmat == 1, arr.ind = T) %>%
#'   #   as.data.frame() %>%
#'   #   dplyr::arrange(row)
#'
#'   # verbose = T
#'   if(verbose == TRUE) {
#'     message(paste0("Calculating north/south/east/west distances"))
#'   }
#'
#'   # side cell totals list
#'   side_lst <- side_summary(
#'     m          = rmat,
#'     indices_df = idx,
#'     cell_res   = cell_res,
#'     max_dist   = max_dist,
#'     func       = "mean",
#'     ncores     = ncores,
#'     verbose    = T
#'   )
#'
#'   if(verbose == TRUE) {
#'     message(paste0("Calculating diagonal distances"))
#'   }
#'
#'   # diagonal totals list
#'   diag_lst <- diag_summary(
#'     m          = rmat,
#'     indices_df = idx,
#'     func       = "mean",
#'     max_count  = max_count,
#'     ncores     = ncores,
#'     verbose    = FALSE
#'   )
#'
#'   # add mean diagonal cell counts
#'   idx$diag_cell_count <- unlist(diag_lst)
#'
#'   # add mean diagonal cell counts
#'   idx$side_cell_count <- unlist(side_lst)
#'
#'   # calculate cell distances from cell counts
#'   dist_mat <- cell_dist(
#'     indices_df = idx,
#'     cell_res   = cell_res,
#'     max_dist   = max_dist,
#'     func       = func,
#'     na_remove  = TRUE
#'   )
#'
#'   # create fetch raster from gridded matrix
#'   fetch_r <-
#'     terra::rast(
#'       x = matrix(
#'         data = dplyr::arrange(
#'           dplyr::left_join(
#'             expand.grid(
#'               row = 1:terra::nrow(r),
#'               col = 1:terra::ncol(r)
#'             ),
#'             dist_mat[, c("row", "col", "distance")],
#'             by = c("row", "col")
#'           ),
#'           row
#'         )$distance,
#'         nrow  = terra::nrow(r),
#'         ncol  = terra::ncol(r),
#'         byrow = TRUE
#'       ),
#'       crs    = terra::crs(r),
#'       extent = terra::ext(r)
#'     )
#'
#'   return(fetch_r)
#' }
