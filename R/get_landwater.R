#' Create grid template
#' @description returns a template SpatRaster. All cell values are 1 so the template raster is not actually "empty
#' @param r SpatRaster, SpatVector or SF object with a Projected Coordinate System CRS.
#' @param res numeric resolution for template raster. Default is 2000, generates a template raster with a grid cell resolution of 2000x2000.
#' @return SpatRaster of template grid
#' @export
#' @importFrom terra rast crs ext
make_grid = function(
    r,
    res = 2000
    ) {

  # if SF object given, convert to SpatVector
  if(any(class(r) == "sf")) {

    r <- terra::vect(r)

  }

  # Template raster of HSI grid
  r_grid <- terra::rast(
    crs         = terra::crs(r),
    extent      = terra::ext(r),
    resolution  = c(res, res),
    vals        = 1
  )

  return(r_grid)

}

#' Make a binary land water raster from multiple land types raster.
#' @description returns a binary land water raster with land cells set to 1 and water cells set to 0. A water value should be specified to indicate which pixel value should represent water in the final raster. The rest of the pixels will indicate land cells.
#' @param r SpatRaster, RasterLayer, SpatVect or a filepath (character) to raster data. Data must be on a projected CRS.
#' @param water_value numeric, value of water cells in input raster. Default is NA.
#' @param res numeric resolution for template raster. Default is 2000, generates a raster with a grid cell resolution of 2000x2000.
#' @param verbose logical, whether messages should be printed. Default is TRUE, no messages print.
#' @return SpatRaster, a binary land water raster with land cells set to 1 and water cells set to 0.
#' @export
#' @importFrom terra rast is.lonlat rasterize yres xres resample setValues values
get_landwater = function(
    r,
    water_value    = NA,
    res            = 2000,
    verbose        = TRUE
){

  # res = 2000
  # r <- st
  # r

  # IF filepath given
  if(any(class(r) == "character")) {

    if(verbose == TRUE) {
      message(paste0("Loading ", r, " as SpatRaster"))
    }

    # read in raster as SpatRaster
    r <- terra::rast(r)

  }

  # if SF object given, convert to SpatVector
  if(any(class(r) == "sf")) {
    if(verbose == TRUE) {
      message(paste0("Converting SF object to SpatVect..."))
    }

    r <- terra::vect(r)

  }

  # make sure raster/shape is in a projected CRS
  if(terra::is.lonlat(r) == TRUE) {

    stop(paste0("Incorrest lat/long CRS, spatial object must be in a projected CRS"))

  }

  # if RasterLayer given
  if(class(r) == "RasterLayer") {

    if(verbose == TRUE) {
      message(paste0("Converting ", class(r), " to ", "SpatRaster"))
    }

    # convert RasterLayer to SpatRaster
    r <- terra::rast(r)

  }

  # if object isn't a spatraster yet at this point, convert to SpatRaster
  if(any(class(r) != "SpatRaster")) {

    # if given a SpatVect polygon
    if(class(r) == "SpatVector") {

      if(verbose == TRUE) {
        message(paste0("Rasterizing SpatVect..."))
        message(paste0("Resolution: ", res, " x ", res))
      }

      # rasterize land polygon
      r <- terra::rasterize(
        x = r,
        make_grid(
          r   = r,
          res = res
          )
        )
    }
  } else {

    if(verbose == TRUE) {
      message(paste0(class(r), " provided..."))
    }

  }
    # difference in x and y resolutions
    res_diff <- terra::xres(r) - terra::yres(r)

    # if raster resolution is irregular, resample to fit a regular grid with a specified resolution
    if(abs(res_diff) > 0.1) {
    # if(terra::xres(r) != terra::yres(r)) {
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
                     # res = max(c(terra::xres(r), terra::yres(r)))
                   ),
         method = "near"
       )
     }

  # set all water cells to a given value and all other cells to another. Presumes raster is just land and water binary.
  if(is.na(water_value)) {

    r <- terra::setValues(
      r, ifelse(is.na(terra::values(r)), 1, 0)
    )

  } else {

    r <- terra::setValues(
      r, ifelse(terra::values(r) == water_value, 1, 0)
    )

  }

  return(r)

}

#' Make a land raster from binary land type raster
#' @description returns a SpatRaster of only land cells. Land cells have a value of 2 by default, all other cells are NA. Fits the HSI raster specifications.
#' @param r SpatRast of binary land type raster. Land and water values only.
#' @param land_value numeric indicating what value land is represented by in raster. Default is 0 to conform to get_landwater() output
#' @param out_value numeric indicating what value land should be represented by in output raster. Default is 2
#' @param as_polygon logical, should polygon be returned instead of a raster. Default is FALSE, returns a terra rast object
#' @param as_sf logical, should polygon be returned as an SF object. Default is FALSE, returns a terra vect object
#' @return SpatRaster with only land cells having a value, all other cells are NA
#' @importFrom terra setValues values
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
#' @return SpatRaster with only land cells having a value, all other cells are NA
#' @importFrom terra setValues values
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
