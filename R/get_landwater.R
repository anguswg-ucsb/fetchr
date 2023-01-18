#' Create grid template
#' @description returns a template SpatRaster. All cell values are 1 so the template raster is not actually "empty
#' @param r SpatRaster, SpatVector or SF object with a Projected Coordinate System CRS.
#' @param res numeric resolution for template raster. Default is 2000, generates a template raster with a grid cell resolution of 2000x2000.
#' @importFrom terra rast crs ext vect
#' @return SpatRaster of template grid
make_grid = function(
    r,
    res = 2000
    ) {

  # if SF object given, convert to SpatVector
  if(is(r, "sf")) {

    r <- terra::vect(r)

  }

  # Template raster grid
  r_grid <- terra::rast(
    crs         = terra::crs(r),
    extent      = terra::ext(r),
    resolution  = c(res, res),
    vals        = 0
  )

  return(r_grid)

}

#' Internal r error handling message
#' @return character error message
r_error_msg = function() {

  paste0(
    "Invalid 'r' argument, 'r' must be one of the following:\n",
    "1. Filepath to a .shp, .gpkg, .geojson, or .tif file\n",
    "2. terra SpatRaster\n",
    "3. raster RasterLayer\n",
    "4. sf POLYGON/MULTIPOLYGON geometry\n",
    "5. terra SpatVector POLYGON/MULTIPOLYGON geometry\n\n",
    "Notes:\n - Spatial object must be in a PROJECTED CRS\n",
    " - For raster data, the default behavior is to classify NA values as water cells and to classify all other values as land cells. If water cells should be classified differently, see the 'water_value' argument"
  )
}

#' Internal filepath error handling message
#' @return character error message for invalid file paths
filepath_error_msg = function() {

  paste0(
    "Invalid 'r' file path argument, file path must point to a .shp, .gpkg, .geojson, or .tif file\n"
    )
}

#' Prepare filepath vector/raster data for fetch calculations
#' @description Internal function for checking and aligning raster data to specifications for calculating fetch using get_fetch
#' @param r character, file path to a shp, gpkg, geojson, or tif file in a projected CRS
#' @param water_value numeric, value indicating the value of the water cells within the provided tif file. If a shp, gpkg, or geojson is given, this argument is ignored. Default is NA, such that NA cells are treated as water cells and all other cells are land cells.
#' @param res numeric resolution of the aligned output raster. Default is 2000.
#' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
#' @importFrom terra vect rast
#' @importFrom methods is
#' @return SpatRaster, binary landwater raster with land values of 0 and water values of 1
prep_filepath <- function(
    r,
    water_value = NA,
    res         = 2000,
    verbose     = FALSE
) {

  # if r is a character
  if(methods::is(r, "character")) {

    # if r is an existing file
    if(file.exists(file.path(r))) {

      if(verbose == TRUE) {

        message(paste0("Loading data from file path: ", r))

      }
        # if shp, gpkg, or geojson is provided, read in 'r' as a SpatVector
        if(any(c("shp", "gpkg", "geojson") %in% strsplit(basename(r), split="\\.")[[1]] )) {

          # read in shapefile as SpatVector
          r <- terra::vect(r)

          # prepare vector object
          r <- prep_poly(
            r           = r,
            res         = res,
            verbose     = TRUE
          )

          return(r)
        }

        # if tif is provided, read in 'r' as a SpatRaster
        if(any(c("tif") %in% strsplit(basename(r), split="\\.")[[1]] )) {

          # read in raster as SpatRaster
          r <- terra::rast(r)

          # prepare vector object
          r <- prep_raster(
            r           = r,
            water_value = water_value ,
            verbose     = TRUE
          )
          return(r)
        }

        # if file path is not too a shapefile, geopackage, geojson, or tif file
        if(!any(c("tif", "shp", "gpkg", "geojson") %in% strsplit(basename(r), split="\\.")[[1]] )) {

          stop(
            filepath_error_msg()
          )

        }
    } else {
      stop(
        paste0("Invalid 'r' file path argument\n", r, " does not exist")
      )
    }
  # if 'r' is NOT a character, return the original data
  } else {

    return(r)

    }
}

#' Prepare polygon vector data for fetch calculations
#' @description Internal function for checking and aligning polygon data to specifications for calculating fetch using get_fetch
#' @param r terra SpatVector POLYGON/MULTIPOLYGON or sf POLYGON/MULTIPOLYGON object in a projected CRS
#' @param res numeric resolution of the aligned output raster. Default is 2000.
#' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
#' @importFrom terra vect is.lonlat rasterize setValues values
#' @importFrom methods is
#' @return SpatRaster, binary landwater raster with land values of 0 and water values of 1
prep_poly <- function(
    r,
    res     = 2000,
    verbose = FALSE
    ) {

  # if none of the required 'r' arguments are provided
  if(!any(methods::is(r, "character"), methods::is(r, "sf"), methods::is(r, "SpatVector"))) {

    return(r)

  }

  # double check r is not a RasterLayer or SpatRaster
  if(any(methods::is(r, "RasterLayer"), methods::is(r, "SpatRaster"))) {

    return(r)

  }

  # if SF object given, convert to SpatVector
  if(methods::is(r, "sf")) {

    # add land value of 0
    r$land_value <- 0

    r <- terra::vect(r)

  }

  # if given a SpatVect polygon
  if(methods::is(r, "SpatVector")) {

    # make sure  is in a projected CRS
    if(terra::is.lonlat(r)) {

      stop(paste0("Incorrect lat/long CRS, spatial object must be in a PROJECTED CRS"))

    }

    # add land value of 0
    r$land_value <- 0

    if(verbose == TRUE) {
      message(paste0("Rasterizing SpatVect to a ", res, "x", res, " grid"))
    }

    # rasterize land polygon
    r <- terra::rasterize(
      x = r,
      make_grid(
        r   = r,
        res = res
      ),
      field = "land_value"
    )

  }

  # fix any gridf
  r <- fix_grid(
    r       = r,
    res     = res,
    verbose = FALSE
  )

  # convert all cells == water_value to 1 and non water_value cells to 0
  r <- terra::setValues(
    r, ifelse(is.na(terra::values(r)), 1, 0)
  )

  return(r)

}

#' Prepare gridded raster data for fetch calculations
#' @description Internal function for checking and aligning raster data to specifications for calculating fetch using get_fetch
#' @param r terra SpatRaster or raster RasterLayer object in a projected CRS.
#' @param water_value numeric, value indicating the value of the water cells within the input raster. If a sf/terra polygon is given, this argument is ignored. Default is NA, such that NA cells are treated as water cells and all other cells are land cells.
#' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
#' @importFrom terra rast is.lonlat setValues values
#' @importFrom methods is
#' @return SpatRaster, binary landwater raster with land values of 0 and water values of 1
prep_raster <- function(
    r,
    water_value = NA,
    verbose     = FALSE
    ) {

  # if none of the required 'r' arguments are provided
  if(!any(methods::is(r, "character"), methods::is(r, "RasterLayer"), methods::is(r, "SpatRaster"))) {

    return(r)

  }

  # double check r is not a sf or SpatVector
  if(any(methods::is(r, "sf"), methods::is(r, "SpatVector"))) {

    return(r)

  }

  # if none of the required 'r' arguments are provided
  if(methods::is(r, "RasterLayer")) {

    # convert RasterLayer to SpatRaster
    r <- terra::rast(r)

  }

  # make sure  is in a projected CRS
  if(terra::is.lonlat(r)) {

    stop(paste0("Incorrect lat/long CRS, spatial object must be in a PROJECTED CRS"))

  }

  # if a raster is provided with only values of 0 and 1s is provided
  if(all(c(0, 1) %in% unique(terra::values(r)))) {

    return(r)

  }

  # fix any grid
  r <- fix_grid(
          r       = r,
          res     = res,
          verbose = FALSE
          )

  # set all water cells to a given value and all other cells to another. Presumes raster is just land and water binary.
  if(is.na(water_value)) {

    if(verbose == TRUE) {
      message(paste0("Replacing NA values with 1 and non NA with 0"))
    }

    # convert all NA cells values to 1 and non NA to 0
    r <- terra::setValues(
      r, ifelse(is.na(terra::values(r)), 1, 0)
    )

  } else {

    if(verbose == TRUE) {
      message(paste0("Replacing ", water_value,  " with 1, and non ", water_value, " with 0"))
    }

    # convert all cells == water_value to 1 and non water_value cells to 0
    r <- terra::setValues(
      r, ifelse(terra::values(r) == water_value, 1, 0)
    )

  }

  return(r)

}

#' Align irregular grid to fit a regular grid with a specified resolution
#' @description Internal utils function for generating landwater rasters
#' @param r terra SpatRaster or raster RasterLayer object in a projected CRS.
#' @param res numeric resolution of the aligned output raster. Default is 2000.
#' @param verbose logical, whether messages should be printed. Default is FALSE, so no messages will print.
#' @importFrom terra yres xres resample
#' @importFrom methods is
#' @return SpatRaster on a regular grid at a specified resolution
fix_grid <- function(
    r,
    res     = 2000,
    verbose = FALSE
    ) {

  # difference in x and y resolutions
  res_diff <- terra::xres(r) - terra::yres(r)

  # if raster resolution is irregular, resample to fit a regular grid with a specified resolution
  if(abs(res_diff) > 0.1) {

    if(verbose == TRUE) {
      message(paste0("Correcting irregular grid resolution..."))
    }

    # resample grid to regular grid
    r <-
      r %>%
      terra::resample(
        y      = make_grid(
          r   = r,
          res = res
        ),
        method = "near"
      )

    return(r)

  # if NOT irregular, return input raster
  } else {

    if(verbose == TRUE) {
      message(paste0("Raster is on a regular grid, returning input raster"))
    }

    return(r)

  }

}

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
get_land = function(
    r,
    land_value = 0,
    out_value  = 2,
    as_polygon = FALSE,
    as_sf      = FALSE
) {

  # if cell value is 1 (land) make it a value of 2, otherwise make cell NA
  land <- terra::setValues(
    r, ifelse(terra::values(r) == land_value, out_value, NA)
  )

  # if a polygon should be returned instead of a raster
  if(as_polygon == TRUE) {

    # if polygon should be returned as an SF object
    if(as_sf == TRUE) {

      message("Generating land polygon...")
      message("Type: SF")

      # SF land polygon
      land <-
        land %>%
        terra::as.polygons() %>%
        sf::st_as_sf()  %>%
        sf::st_make_valid() %>%
        sf::st_cast() %>%
        stats::setNames(c("land", "geometry")) %>%
        sf::st_transform(26915)

      # return as terra vect object
    } else {

      message("Generating land polygon...")
      message("Type: SpatVector")

      # spatvector land polygon
      land <-
        land %>%
        terra::as.polygons() %>%
        stats::setNames(c("land"))

    }


  } else {

    message("Generating land raster...")
    message("Type: SpatRaster")

  }

  return(land)

  # view map
  # land %>%
  #   raster::raster()%>%
  #   mapview::mapview()

}

#' Make a water raster from binary land type raster
#' @description returns a SpatRaster of only water cells. Water cells have a value of 1 by default, all other cells are NA. Fits the HSI raster specifications.
#' @param r SpatRast of binary land type raster. Land and water values only.
#' @param water_value numeric indicating what value water is represented by in raster. Default is 1 to conform to get_landwater() output
#' @param out_value numeric indicating what value land should be represented by in output raster. Default is 1
#' @param as_polygon logical, should polygon be returned instead of a raster. Default is FALSE, returns a terra rast object
#' @param as_sf logical, should polygon be returned as an SF object. Default is FALSE, returns a terra vect object
#' @importFrom terra setValues values
#' @return SpatRaster with only land cells having a value, all other cells are NA
get_water = function(
    r,
    water_value = 1,
    out_value   = 1,
    as_polygon  = FALSE,
    as_sf       = FALSE
) {

  # if cell value is 1 (land) make it a value of 2, otherwise make cell NA
  water <- terra::setValues(
    r, ifelse(terra::values(r) == water_value, out_value, NA)
  )

  # if a polygon should be returned instead of a raster
  if(as_polygon == TRUE) {

    # if polygon should be returned as an SF object
    if(as_sf == TRUE) {

      # make polygon from water raster
      water <-
        water %>%
        terra::as.polygons() %>%
        sf::st_as_sf()  %>%
        sf::st_make_valid() %>%
        sf::st_cast() %>%
        stats::setNames(c("water", "geometry")) %>%
        sf::st_transform(26915)

      # return as terra vect object
    } else {

      water <-
        water %>%
        terra::as.polygons() %>%
        stats::setNames(c("water"))

    }
  }

  return(water)

}

