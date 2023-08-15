# version 1.1.1.00

-   big improvements to speed and efficiency for `get_fetch`

-   providing `func =` inputs other than "mean", "min", "max", or "sum" raises an error telling the user to use one of those 4 functions

-   `ncores` argument has been replaced with `in_parallel` which is a logical that when TRUE, will run the calculations in parallel using the number of cores on the machine - 1. The number of cores is determined by the `parallel::detectCores()` function

-   Added `get_fetch_directions` function which will generate an 8 layer `SpatRaster` with each layer being a wind fetch raster for one of the 8 cardinal directions (north, south, east, west, northeast, northwest, southeast, and southwest)
