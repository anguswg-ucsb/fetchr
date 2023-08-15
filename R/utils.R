utils::globalVariables(c(".", "res", "cell_count", "diag_cell_count", "side_cell_count"))

#' Get the primary and secondary diagonals for a cell in a matrix
#'
#' @param grid matrix
#' @param row numeric, row number
#' @param col numeric, column number
#' @param NROW numeric, number of rows in `grid`
#' @param NCOL numeric, number of columns in `grid`
#' @noRd
#' @keywords internal
#' @return list of length 2, containing the numeric vector of values for the primary and secondary diagonals and the index of the given cell within that numeric vector
get_diag_vectors <- function(grid,
                             row,
                             col,
                             NROW = NULL,
                             NCOL = NULL) {
  # grid = rmat
  # row  =   idx$row[i]
  # col  =   idx$col[i]
  # search_value = 0
  # max_count = max_count
  # matrix_size <- 6
  # initial_matrix <- matrix(1, nrow = matrix_size, ncol = matrix_size)
  # grid <- initial_matrix
  # row = 4
  # col = 5
  # # island2_rows <- c(7, 7, 7)
  # # island2_cols <- c(7, 8, 9)

  # # initial_matrix[island2_rows, island2_cols] <- 0

  # calculate the distances to the edges
  left_dist <- min(row - 1, col - 1)
  right_dist <- min(NCOL - col, NROW - row)

  # top left corner coordinates of diagonal
  topleft_row <- row - left_dist
  topleft_col <- col - left_dist

  # bottom right corner coordinates of diagonal
  bottomright_row <- row + right_dist
  bottomright_col <- col + right_dist

  # calculate the distances to the edges for top right to bottom left diagonal
  topright_dist <- min(row - 1, NCOL - col)
  bottomleft_dist <- min(NROW - row, col - 1)

  # top right corner coordinates of diagonal
  topright_row <- row - topright_dist
  topright_col <- col + topright_dist

  # bottom left corner coordinates of diagonal
  bottomleft_row <- row + bottomleft_dist
  bottomleft_col <- col - bottomleft_dist

  # coordinates of primary diagonal
  prim_diag_coords <- cbind(
    seq(topleft_row, bottomright_row, by = 1),
    seq(topleft_col, bottomright_col, by = 1)
  )

  # index of starting point in primary diagonal
  prim_index <- which(prim_diag_coords[, 1] == row)

  # coordinates of secondary diagonal
  sec_diag_coords <- cbind(
    seq(bottomleft_row, topright_row, by = -1), # Going from bottomleft row to topright row
    seq(bottomleft_col, topright_col, by = 1) # Going from bottomleft COL to topright COL
    # rev(seq(topright_row, bottomleft_row, by = 1)),  # LATEST WORKING Top Right row to bottom left row then reverse
    # rev(seq(topright_col, bottomleft_col, by = -1))  # LATEST WORKING Top Right COL to bottom left COL then reverse
    # seq(topright_row, bottomleft_row, by = 1),  #  NO REVERSE
    # seq(topright_col, bottomleft_col, by = -1)  #  NO REVERSE
  )

  # index of starting point in secondary diagonal
  sec_index <- which(sec_diag_coords[, 1] == row)

  # primary diagonal (top left to bottom right)
  primary <- grid[prim_diag_coords]

  # secondary diagonal (top right to bottom left)
  secondary <- grid[sec_diag_coords]

  # mmm <- print_diags(
  #   x = grid,
  #   origin_row = row,
  #   origin_col = col,
  #   TL_row= topleft_row,
  #   TL_col =topleft_col,
  #   BR_row = bottomright_row,
  #   BR_col= bottomright_col,
  #   BL_row = bottomleft_row,
  #   BL_col = bottomleft_col,
  #   TR_row = topright_row,
  #   TR_col = topright_col
  #   )
  # plot(terra::rast(mmm))

  return(
    list(
      primary   = list(index = prim_index, diag = primary),
      secondary = list(index = sec_index, diag = secondary)
    )
  )

}

#' Get the distance from a cell to a specific search value in 8 directions
#' Internal function used to get distances in all directions from a given cell (row, column)
#' @param grid matrix, numerics or integers
#' @param row numeric, row numbers
#' @param col, numeric, column number
#' @param search_value numeric, value to search for. Default is 0
#' @param max_count numeric, maximum number of cells to count outwards
#' @noRd
#' @keywords internal
#' @return a numeric vector of length 8, containing the distance from the given cell to search_value in each of the 8 directions
get_cell_distances <- function(grid, row, col, search_value, max_count) {

  # grid dimensions
  NROW = nrow(grid)
  NCOL = ncol(grid)

  # get diagonal vectors for given row, col
  diag_vectors <- get_diag_vectors(
    grid         = grid,
    row          = row,
    col          = col,
    NROW         = NROW,
    NCOL         = NCOL
  )

  # # create vector of cell distances
  cell_distances <- as.integer(
    c(
      get_vector_distances(
        vect         = as.integer(grid[row, ]),
        start_index  = col,
        search_value = search_value,
        max_count    = max_count
      ),
      get_vector_distances(
        vect         = as.integer(grid[,col]),
        start_index  = row,
        search_value = search_value,
        max_count    = max_count
      ),
      get_vector_distances(
        vect         = as.integer(diag_vectors[['primary']][['diag']]),
        start_index =  as.integer(diag_vectors[['primary']][['index']]),
        search_value = search_value,
        max_count    = max_count
      ),
      get_vector_distances(
        vect         = as.integer(diag_vectors[['secondary']][['diag']]),
        start_index =  as.integer(diag_vectors[['secondary']][['index']]),
        search_value = search_value,
        max_count    = max_count
      )
    )
  )

  return(cell_distances)

}

#' Run get_direction_distances, either in parallel or not
#'
#' @param grid matrix
#' @param index_df dataframe
#' @param max_count numeric
#' @param func character indicating the function to apply to each set of distances. Either "mean", "min", "max", or "sum"
#' @param in_parallel logical, whether function should in parallel or not
#' @noRd
#' @keywords internal
#' @return list of summarized distance values for each row in index_df
do_cell_distances <- function(grid, index_df, max_count, func, in_parallel) {

  # match function to R function
  func <- match.fun(func)

  # if in_parallel is TRUE, run function in parallel
  if(in_parallel) {

    # initiate parallel clusters
    clust <- parallel::makeCluster(parallel::detectCores()-1)

    # export variables to clusters and load packages
    parallel::clusterExport(clust, c('get_cell_distances', 'grid', 'index_df', 'max_count', 'get_diag_vectors', 'func',
                                     'get_vector_distances'),
                            envir=environment())
    # system.time({

    # side counts/lengths directions
    dists <- parallel::parLapply(
      cl = clust,
      1:nrow(index_df), function(i) {

        func(
          get_cell_distances(
            grid = grid,
            row  =   index_df$row[i],
            col  =   index_df$col[i],
            search_value = 0,
            max_count = max_count
          ),
          na.rm = TRUE
        )

      })
    # })

    # stop parallel processing clusters
    parallel::stopCluster(clust)

  } else {

    # system.time({
    # get average distance for each cell in matrix
    dists <- lapply(1:nrow(index_df), function(i) {
      # i = 900
      # i = 335

      func(
        get_cell_distances(
          grid = grid,
          row  =   index_df$row[i],
          col  =   index_df$col[i],
          search_value = 0,
          max_count = max_count
        ),
        na.rm = TRUE
      )

    })
    # })
  }

  return(dists)

}

#' Get the distances in 8 directions from a cell to a specified search value
#' For internal use and calculations for get_fetch_directions
#' @param grid matrix
#' @param row numeric, row number to get directions for
#' @param col numeric, column number to get directions for
#' @param search_value numeric, value to search for in each direction
#' @param max_count numeric, maximum number of cells to count outwards
#' @noRd
#' @keywords internal
#' @return list of length 8 containing the distances to search_value
get_direction_distances <- function(grid, row, col, search_value, max_count) {

  # grid = rmat
  # row  =   idx$row[i]
  # col  =   idx$col[i]
  # search_value = 0
  # max_count = max_count

  NROW = nrow(grid)
  NCOL = ncol(grid)

  # get diagonal vectors for given row, col
  diag_vectors <- get_diag_vectors(
    grid         = grid,
    row          = row,
    col          = col,
    NROW         = NROW,
    NCOL         = NCOL
  )

  row_distances <- get_vector_distances(
    vect         = as.integer(grid[row, ]),
    start_index  = col,
    search_value = search_value,
    max_count    = max_count
  )

  col_distances <- get_vector_distances(
    vect         = as.integer(grid[,col]),
    start_index  = row,
    search_value = search_value,
    max_count    = max_count
  )

  diag_top_distances <- get_vector_distances(
    vect         = as.integer(diag_vectors[['primary']][['diag']]),
    start_index  =  as.integer(diag_vectors[['primary']][['index']]),
    search_value = search_value,
    max_count    = max_count
  )

  # col
  # diag_vectors$secondary[col]
  # length(diag_vectors$secondary)
  # min(row, col)
  diag_bottom_distances <- get_vector_distances(
    vect         = as.integer(diag_vectors[['secondary']][['diag']]),
    start_index =  as.integer(diag_vectors[['secondary']][['index']]),
    # vect         = as.integer(diag_vectors$secondary),
    # start_index  = col,
    # start_index  = ifelse(
    #   above_secondary_diag(row, col, NCOL),
    #   col,
    #   row
    #   # row + 1
    #   # (length(as.integer(diag_vectors$secondary)) - col) + 1
    # ),
    # start_index  = length(as.integer(diag_vectors$secondary)) - (row - 1),
    search_value = search_value,
    max_count    = max_count
  )

  # distances in each of the 8 directions
  west <- row_distances[1]
  east <- row_distances[2]

  north <- col_distances[1]
  south <- col_distances[2]

  northwest <- diag_top_distances[1]
  southeast <- diag_top_distances[2]

  northeast <- diag_bottom_distances[1]
  southwest <- diag_bottom_distances[2]

  return(
    list(
      east      = east,
      southeast = southeast,
      south     = south,
      southwest = southwest,
      west      = west,
      northwest = northwest,
      north     = north,
      northeast = northeast
    )
  )

}

#' Run get_direction_distances, either in parallel or not
#'
#' @param grid matrix
#' @param index_df dataframe
#' @param max_count numeric
#' @param in_parallel logical, whether function should in parallel or not
#' @noRd
#' @keywords internal
#' @return list of direction distances for each row index_df
do_dir_distances <- function(grid, index_df, max_count, in_parallel) {

  if(in_parallel) {

    # initiate parallel clusters
    clust <- parallel::makeCluster(  parallel::detectCores()-1)

    # export variables to clusters and load packages
    parallel::clusterExport(clust, c('get_direction_distances', 'grid', 'index_df', 'max_count', 'get_diag_vectors',
                                     'get_vector_distances'),
                            envir=environment())
    # system.time({

    # side counts/lengths directions
    dists <- parallel::parLapply(
      cl = clust,
      1:nrow(index_df), function(i) {

        get_direction_distances(
          grid = grid,
          row  =   index_df$row[i],
          col  =   index_df$col[i],
          search_value = 0,
          max_count = max_count
        )

      })
    # })

    # stop parallel processing clusters
    parallel::stopCluster(clust)

  } else {

    # system.time({
    # get average distance for each cell in matrix
    dists <- lapply(1:nrow(index_df), function(i) {
      # i = 900
      # i = 335

      get_direction_distances(
        grid = grid,
        row  =   index_df$row[i],
        col  =   index_df$col[i],
        search_value = 0,
        max_count = max_count
      )

    })
    # })
  }

  return(dists)

}

#' Find the distance from a given index to the nearest search_value to the left and right of the given index
#' Internal function used for finding distances to the left and right of a given cell in a matrix. This function works by taking a vector,
#' which ends up being the row/column/diagonal values from a given cell, and returns the distance left and right of the cell before the given search_value is found
#'
#' @param vect numeric vector
#' @param start_index numeric, index to start search from
#' @param search_value numeric, value to the left and/or right of the given index to find distances too
#' @param max_count numeric, maximum distance. if no search_value is found or the value is greater than max_count, then max_count is given as the distance
#' @noRd
#' @keywords internal
#' @return numeric vector of length 2,the distance to the left and right of the index, respectively
get_vector_distances <- function(vect, start_index, search_value, max_count) {
  # vect <- row_vector
  # start_index <- col
  # search_value <- 0

  # where are the first indices with the search_value
  cell_stops <- which(vect == search_value)

  # length(which(vect == search_value)) == 0
  # vect == search_value

  # if there are none of the search values in the vector, return a list with the left and right counts as the max_count
  if (!any(vect == search_value)) {
    # if(length(which(vect == search_value)) == 0) {

    return(
      c(max_count, max_count)
    )

  }

  # indices of search values (cells that equalled the search_value) that are to the LEFT of the start_index
  # left_search <- cell_stops[cell_stops <= (start_index + 1)]
  left_search <- cell_stops[cell_stops < start_index]

  # if there are zero indices to the left of the start index, set left to max_count
  if(length(left_search) == 0) {

    # set left to max_count allowed
    left = max_count

  } else {

    # distance to first search value is the current index minus the biggest index in the values found by running which(==search_values)
    left = start_index - left_search[length(left_search)]
    # left = start_index - max(left_search)

    # if left value is greater than the max_count, set left to the max count, otherwise leave it alone
    left <- ifelse(left > max_count, max_count, left)

  }

  # # if left value is greater than the max_count, set left to the max count, otherwise leave it alone
  # left <- ifelse(left > max_count, max_count, left)

  # indices of search values that are to the right of the start_index
  right_search <- cell_stops[cell_stops > start_index]

  # if there are zero indices to the RIGHT of the start index, set RIGHT to max_count
  if(length(right_search) == 0) {

    # set right to max_count allowed
    right = max_count

  } else {
    # distance to first search value is the smallest index in search_values minus current index gives the distance to the nearest search_value in the RIGHT DIRECTION
    right = right_search[1] - start_index
    # right = min(right_search)- start_index

    # if right value is greater than the max_count, set right to the max count, otherwise leave it alone
    right <- ifelse(right > max_count, max_count, right)

  }

  # # if right value is greater than the max_count, set right to the max count, otherwise leave it alone
  # right <- ifelse(right > max_count, max_count, right)
  # list(
  #   left  = left,
  #   right = right
  # )

  return(c(left, right))
  # list(
  #   left  = left,
  #   right = right
  # )
  # )
}

#' Check which half of a matrix a row, column is in along the primary diagonal
#' If you cut a matrix along the secondary diagonal (top left corner to bottom right corner),
#' this condition will check if a row, column is in the top right half or bottom left half of the matrix.
#' @param row numeric, row number
#' @param col numeric, column number
#' @noRd
#' @keywords internal
#' @return logical, TRUE if above the primary diagonal (top right half), otherwise FALSE (bottom left half)
above_priamary_diag <- function(row, col) {
  # If you cut a matrix along the secondary diagonal (top left corner to bottom right corner),
  # this condition will check if a row, column is in the top right half or bottom left half of the matrix.

  # if row is less than or equal to the column, we are in the top right half of matrix
  if (row <= col) {

    return(TRUE)

    # otherwise we are in the bottom left half
  } else {

    return(FALSE)

  }
}

#' Check which half of a matrix a row, column is in along the secondary diagonal
#' If you cut a matrix along the secondary diagonal (bottom left corner to top right corner),
#' this condition will check if a row, column is in the top left half or bottom right half of the matrix.
#' @param row numeric, row number
#' @param col numeric, column number
#' @param NCOL numeric, number of columns in matrix
#' @noRd
#' @keywords internal
#' @return logical, TRUE if above the secondary diagonal (top left half), otherwise FALSE (bottom right half)
above_secondary_diag <- function(row, col, NCOL) {
  # if you cut a matrix along the secondary diagonal (bottom left corner to top right corner),
  # this condition will check if a row, column is in the top left half or bottom right half of the matrix

  # if row + col is less than or equal to the number of columns in the matrix, we are in the top left half
  if(row + col <= NCOL) {

    return(TRUE)

    # otherwise we are in the bottom right half
  } else {

    return(FALSE)

  }
}

#' Create a new matrix with diagonal lines from a cell highlighted/changed
#' Provide the top left, top right, bottom left, and bottom right row and columns representing the diagonals from a given cell and create a new matrix with replaced values.
#'  For internal use.
#'
#' @param x matrix
#' @param origin_row numeric, original row
#' @param origin_col numeric, original column
#' @param TL_row numeric
#' @param TL_col numeric
#' @param BR_row numeric
#' @param BR_col numeric
#' @param BL_row numeric
#' @param BL_col numeric
#' @param TR_row numeric
#' @param TR_col numeric
#' @noRd
#' @keywords internal
#' @return numeric matrix
print_diags <- function(x, origin_row, origin_col,
                        TL_row, TL_col, BR_row, BR_col,
                        BL_row, BL_col, TR_row, TR_col) {
  # x <- grid
  # origin_row = row
  # origin_col = col
  # TL_row= topleft_row
  # TL_col =topleft_col
  # BR_row = bottomright_row
  # BR_col= bottomright_col
  # BL_row = bottomleft_row
  # BL_col = bottomleft_col
  # TR_row = topright_row
  # TR_col = topright_col

  # set primary axis
  x[
    cbind(
      seq(TL_row, BR_row, by = 1),
      seq(TL_col, BR_col, by = 1)
    )
  ]  <- 2

  # set secondary axis
  x[
    cbind(
      # rev(seq(topright_row, bottomleft_row, by = 1)),
      # seq(bottomleft_col, topright_col, by = 1)
      rev(seq(TR_row, BL_row, by = 1)),  # Increasing row indices
      rev(seq(TR_col, BL_col, by = -1))  # Decreasing column indices
      # seq(topright_row, bottomleft_row, by = 1),
      # seq(topright_col, bottomleft_col, by = -1)
      # seq(topright_row, bottomleft_row, by = 1),
      # seq(topright_col, bottomleft_col, by = -1)
    )
  ] <- 2

  x[origin_row, origin_col] <- 9999

  return(x)

}


#' Retrieve length of any diagonal direction for a specific matrix cell
#' @description Given a matrix and the row and column number for the cell of interest, return the number of consecutive non-zero values
#' @param m matrix of interest
#' @param row_index numeric row number of matrix cell to get indices of diagonals for
#' @param col_index numeric column number of matrix cell to get indices of diagonals for
#' @param max_count numeric, maximum number of cells within max_dist
#' @param verbose logical. Whether messages should output to console. Default is FALSE, no messages output.
#' @noRd
#' @keywords internal
#' @return integer of diagonal cell counts or a matrix with diagonal cell counts
diag_dist <- function(
    m,
    row_index  = NULL,
    col_index  = NULL,
    max_count,
    verbose    = FALSE
) {

  # check if col 1
  if(col_index == 1) {

    left_col_start = col_index

  } else {

    left_col_start = col_index - 1

  }

  # sequence of left columns
  left_col_seq <- seq(1, left_col_start, 1)

  # RIGHT COLUMNS
  # check if col 1
  if(col_index == ncol(m)) {

    right_col_start = ncol(m)

  } else {

    right_col_start = col_index + 1

  }

  # sequence of right columns
  right_col_seq <- seq(right_col_start, ncol(m), 1)

  # DOWN ROWS
  # check if row 1
  if(row_index == nrow(m)) {

    down_row_start = row_index

  } else {

    down_row_start = row_index + 1

  }

  # sequence of below rows
  down_row_seq <- seq(down_row_start, nrow(m), 1)

  # UP ROWS
  # check if row 1
  if(row_index == 1) {

    up_row_start = 1

  } else {

    up_row_start = row_index - 1

  }

  # sequence of above rows
  up_row_seq <- seq(up_row_start,  1,  -1)


  # vector diagonal indices
  # DOWN RIGHT
  dr_col = right_col_seq[1:min(length(right_col_seq), length(down_row_seq))]
  dr_row = down_row_seq[1:min(length(right_col_seq), length(down_row_seq))]

  # DOWN LEFT
  dl_col = rev(left_col_seq)[1:min(length(left_col_seq), length(down_row_seq))]
  dl_row = down_row_seq[1:min(length(left_col_seq), length(down_row_seq))]

  # UP RIGHT
  ur_col = right_col_seq[1:min(length(right_col_seq), length(up_row_seq))]
  ur_row = up_row_seq[1:min(length(right_col_seq), length(up_row_seq))]

  # UP LEFT
  ul_col = rev(left_col_seq)[1:min(length(left_col_seq), length(up_row_seq))]
  ul_row = up_row_seq[1:min(length(left_col_seq), length(up_row_seq))]

  # Count vector diagonals
  # DOWN RIGHT
  dr_count <- check_length(
    vect      = m[cbind(dr_row, dr_col)],
    max_count = max_count
  )

  # DOWN LEFT
  dl_count <- check_length(
    vect      = m[cbind(dl_row, dl_col)],
    max_count = max_count
  )

  # UP RIGHT
  ur_count <- check_length(
    vect      = m[cbind(ur_row, ur_col)],
    max_count = max_count
  )

  # UP LEFT
  ul_count <- check_length(
    vect      = m[cbind(ul_row, ul_col)],
    max_count = max_count
  )

  # diagonal count vector
  diag_count <- c(dr_count, dl_count, ur_count, ul_count)

  return(diag_count)

}
# #' Retrieve length of any diagonal direction for a specific matrix cell
# #' @description Given a matrix and the row and column number for the cell of interest, return the number of consecutive non-zero values
# #' @param m matrix of interest
# #' @param row_index numeric row number of matrix cell to get indices of diagonals for
# #' @param col_index numeric column number of matrix cell to get indices of diagonals for
# #' @param vdirect character indicating vertical direction for diagonal. Either "up" or "down". Default is "down".
# #' @param hdirect character indicating horizontal direction for diagonal. Either "left" or "right". Default is "right".
# #' @param max_count numeric, maximum number of cells within max_dist
# #' @param id_col logical. Whether a column stating the direction of the diagonal should be included. Default is FALSE, no direction column is added.
# #' @param as_matrix logical. Whether a matrix with row/column/cell_count/direction is returned or a numeric vector. Default is FALSE, returns a numeric vector of cell counts
# #' @param verbose logical. Whether messages should output to console. Default is FALSE, no messages output.
# #' @return integer of diagonal cell counts or a matrix with diagonal cell counts
# diag_dist <- function(
    #     m,
#     row_index  = NULL,
#     col_index  = NULL,
#     vdirect    = "down",
#     hdirect    = "right",
#     max_count,
#     as_matrix  = FALSE,
#     id_col     = FALSE,
#     verbose    = FALSE
# ) {
#
#   # m <- rmat
#   # hdirect = "right"
#   # vdirect = "down"
#
#   # concatenate direction string
#   direction <- paste0(vdirect, "_", hdirect)
#
#   # system.time(
#   # check if left direction
#   if(grepl("left", hdirect) == TRUE) {
#
#     if(verbose == TRUE) {
#       message(paste0("Sequencing ", hdirect, " diagonals", " - (", direction, ")"))
#     }
#
#     # check if col 1
#     if(col_index == 1) {
#
#       col_start = col_index
#
#     } else {
#
#       col_start = col_index - 1
#
#     }
#
#     # sequence of above columns
#     col_seq <- seq(1, col_start, 1)
#
#
#   } else {    # check if right direction
#
#     # if(grepl("right", hdirect) == TRUE) {
#     if(verbose == TRUE) {
#       message(paste0("Sequencing ", hdirect, " diagonals", " - (", direction, ")"))
#     }
#
#     # check if col 1
#     if(col_index == ncol(m)) {
#
#       col_start = col_index
#
#     } else {
#
#       col_start = col_index + 1
#
#     }
#
#     # sequence of right columns
#     col_seq <- seq(col_start, ncol(m), 1)
#
#   }
# # )
#   # check if above direction
#   if(grepl("up", vdirect) == TRUE) {
#
#     if(verbose == TRUE) {
#       message(paste0("Sequencing ", vdirect, " diagonals", " - (", direction, ")"))
#     }
#
#     # check if row 1
#     if(row_index == 1) {
#
#       row_start = 1
#
#     } else {
#
#       row_start = row_index - 1
#
#     }
#
#     # sequence of above rows
#     row_seq <- seq(row_start,  1,  -1)
#
#   } else {  # check if bottom direction
#
#     # if(grepl("bottom", vdirect) == TRUE) {
#     if(verbose == TRUE) {
#       message(paste0("Sequencing ", vdirect, " diagonals", " - (", direction, ")"))
#     }
#
#     # check if row 1
#     if(row_index == nrow(m)) {
#
#       row_start = row_index
#
#     } else {
#
#       row_start = row_index + 1
#
#     }
#
#     # sequence of below rows
#     row_seq <- seq(row_start, nrow(m), 1)
#
#   }
#
#   # minimum sequence
#   min_seq <- min(length(col_seq), length(row_seq))
#
#   # top left/north west diagonal
#   if(vdirect == "down" & hdirect == "right") {
#
#     # vector diagonal indices
#     diag_row = row_seq[1:min_seq]
#     diag_col = col_seq[1:min_seq]
#
#   } else if(vdirect == "down" & hdirect == "left"){
#
#     # vector diagonal indices
#     diag_row = row_seq[1:min_seq]
#     diag_col = rev(col_seq)[1:min_seq]
#
#     # # diagonal matrix
#     # diag_mat <-
#     #   matrix( cbind(row_index, col_index, check_length(m[cbind(diag_row, diag_col)]), direction),ncol = 4,
#     #     dimnames = list(c("val"), c("row_idx", "col_idx", "cell_count", "direction")) )
#
#     # diagonal dataframe
#     # diag_df <-data.frame(cell_count = check_length(m[cbind(diag_row, diag_col)]))
#
#   } else if(vdirect == "up" & hdirect == "right"){
#
#     # vector diagonal indices
#     diag_row = row_seq[1:min_seq]
#     diag_col = col_seq[1:min_seq]
#
#   } else if(vdirect == "up" & hdirect == "left"){
#
#     # vector diagonal indices
#     diag_row = row_seq[1:min_seq]
#     diag_col = rev(col_seq)[1:min_seq]
#
#   }
#
#   # if a matrix is requested, with or without ID columns (direction)
#   if(as_matrix == TRUE) {
#
#     # add an ID column stating the direction of the diagonal
#     if(id_col == TRUE) {
#
#       # diagonal matrix
#       diag_mat <-
#         matrix(
#           cbind(
#             row_index,
#             col_index,
#             check_length(
#               vect      = m[cbind(diag_row, diag_col)],
#               max_count = max_count
#             ),
#             direction
#           ),
#           ncol     = 4,
#           dimnames = list(c("r"), c("row", "col", "cell_count", "direction"))
#         )
#     } else {
#
#       # diagonal matrix w/o row & column data
#       diag_mat <-
#         matrix(
#           cbind(
#             row_index,
#             col_index,
#             check_length(
#               vect      = m[cbind(diag_row, diag_col)],
#               max_count = max_count
#             ),
#           ),
#           ncol = 3,
#           dimnames = list(c("r"), c("row", "col", "cell_count"))
#         )
#     }
#
#     return(diag_mat)
#
#   } else {   # Default returns an integer with the total diagonal cell count
#
#     # get count of diagonal vector
#     diag_count <- check_length(
#       vect      = m[cbind(diag_row, diag_col)],
#       max_count = max_count
#     )
#
#     return(diag_count)
#
#   }
#
# }

#' Summarize lengths of sequential raster cells diagonals
#' @description Internal function for calculating distance along diagonals for specified indices of a matrix.
#' @param m matrix of interest
#' @param indices_df dataframe of matrix indices from which to calculate lengths from
#' @param func character R function to summarize surrounding vertical/horizontal pixels. One of "mean", "min", "max", or "sum". Default is "mean"
#' @param max_count numeric, maximum cell count. Equal to maximum distance/cell resolution
#' @param ncores numeric, number of cores to use for parallel processing. Default is 2.
#' @param verbose  logical, whether to print messages or not. Default is FALSE.
#' @return list of numeric vectors of diagonal distances
#' @noRd
#' @keywords internal
diag_summary <- function(
    m,
    indices_df,
    func    = "mean",
    max_count,
    ncores  = 2,
    verbose = FALSE
) {
  # get function
  func <- match.fun(func)

  # verbose <- FALSE
  # system.time(
  diags_lst <- lapply(1:nrow(indices_df), function(n) {

    if(verbose == TRUE) {
      message(paste0(n, "/", nrow(indices_df)))
    }

    # diag <- diag_dist(
    #   m          = m,
    #   row_index  = indices_df$row[n],
    #   col_index  = indices_df$col[n],
    #   max_count  = max_count,
    #   verbose    = F
    # ) %>%
    #   func(na.rm = T)

    diag <- func(
              diag_dist(
                m          = m,
                row_index  = indices_df$row[n],
                col_index  = indices_df$col[n],
                max_count  = max_count,
                verbose    = F
                ),
              na.rm = T
              )

  })
  # )
  return(diags_lst)

}
# }) %>%
# func(na.rm = T)
# data.frame(row   = idx$row[n], col   = idx$col[n],cell_count = .,
# direction  = c("up_left", "down_left", "up_right", "down_right"))
# })
# m          = m
# row_index  = indices_df$row[n]
# col_index  = indices_df$col[n]
# vdirect    = diag_dirs[y, 1]
# hdirect    = diag_dirs[y, 2]
# max_count  = max_count
# as_matrix  = FALSE
# id_col     = F
# verbose    = T

# diagonal directions
# diag_dirs <- expand.grid(
#   vdir             = c("up", "down"),
#   hdir             = c("left", "right"),
#   stringsAsFactors = F
# )

#   system.time(
#   # indices list
#   diags_lst <- lapply(1:nrow(indices_df), function(n) {
#
#     if(verbose == TRUE) {
#       message(paste0(n, "/", nrow(indices_df)))
#     }
#
#     diags <- sapply(1:nrow(diag_dirs), function(y) {
#
#       diag_dist(
#         m          = m,
#         row_index  = indices_df$row[n],
#         col_index  = indices_df$col[n],
#         vdirect    = diag_dirs[y, 1],
#         hdirect    = diag_dirs[y, 2],
#         max_count  = max_count,
#         as_matrix  = FALSE,
#         id_col     = F,
#         verbose    = F
#       )
#     }) %>%
#       func(na.rm = T)
#     # data.frame(row   = idx$row[n], col   = idx$col[n],cell_count = .,
#     # direction  = c("up_left", "down_left", "up_right", "down_right"))
#   })
# )
#   return(diags_lst)
# }

# # initiate parallel clusters
# clust <- parallel::makeCluster(ncores)
#
# # export variables to clusters and load packages
# parallel::clusterExport(clust, c('max_count', 'func', 'rmat', 'idx', "diag_dirs",
#                                  'diag_dist', "check_length", "m", "indices_df", "verbose"))
#
# # export packages to clusters and load
# parallel::clusterEvalQ(clust, {
#   library(dplyr)
# })
#
#
# diags <- parallel::parLapply(cl = clust, 1:nrow(indices_df), function(n) {
#     if(verbose == TRUE) {
#       message(paste0(n, "/", nrow(indices_df)))
#     }
#     # parallel::parSapply(cl = clust, 1:nrow(diag_dirs), function(y) {
#       diags <- sapply(1:nrow(diag_dirs), function(y) {
#
#         diag_dist(
#           m          = m,
#           row_index  = indices_df$row[n],
#           col_index  = indices_df$col[n],
#           vdirect    = diag_dirs[y, 1],
#           hdirect    = diag_dirs[y, 2],
#           max_count  = max_count,
#           as_matrix  = FALSE,
#           id_col     = F,
#           verbose    = F
#         )
#       }) %>%
#         func(na.rm = T)
#     })

# stop parallel processing clusters
# parallel::stopCluster(clust)

#' Summarize lengths of sequential raster cells horizontally and vertically
#' @description Internal function for calculating distances horizontally and vertically for specified indices of a matrix.
#' @param m matrix of interest
#' @param indices_df dataframe of matrix indices from which to calculate lengths from
#' @param cell_res numeric, resolution of raster cells. Default is 2000
#' @param max_dist numeric, maximum distance to calculate fetch out to. in meters. Default is 20000m
#' @param func character R function to summarize surrounding vertical/horizontal pixels. One of "mean", "min", "max", or "sum". Default is "mean"
#' @param ncores numeric, number of cores to use for parallel processing. Default is 2.
#' @param verbose logical, whether to print messages or not. Default is FALSE.
#' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster
#' @noRd
#' @keywords internal
#' @return returns a list of summarized side lengths for every water cell, the same length as the number of rows in the input indices_df dataframe
side_summary <- function(
    m,
    indices_df,
    cell_res   = 2000,
    max_dist   = 100000,
    func       = "mean",
    ncores     = 2,
    verbose    = FALSE
) {

  # if smaller number of points, no parallel processing needed
  if(nrow(indices_df) <= 4000) {

    # system.time(
    # side counts/lengths directions
    side_lengths <- lapply(1:nrow(indices_df), function(z) {
      side_count(
        m         = m,
        row_index = indices_df$row[z],
        col_index = indices_df$col[z],
        cell_res  = cell_res,
        max_dist  = max_dist,
        func      = func,
        verbose   = FALSE
      )
    })
    # )
    return(side_lengths)

  } else {

    # medium number of points, limit number of cores will speed up processing
    if(nrow(indices_df) > 4000 & nrow(indices_df) <= 20000) {

      # if more than 4 cores requested and system has more than 4 cores, set ncores to 4
      if(ncores > 4 & parallel::detectCores() > 4 ) {

        ncores <- 4

        # otherwise make ncores 2
      } else {

        ncores <- 2
      }
    }

    # initiate parallel clusters
    clust <- parallel::makeCluster(ncores)

    # export variables to clusters and load packages
    parallel::clusterExport(clust, c('cell_res', 'max_dist', "side_count", "check_cols",
                                     "check_rows", "check_length", "m", "indices_df", "verbose"),
                            envir=environment())

    # side counts/lengths directions
    side_lengths <- parallel::parLapply(cl = clust, 1:nrow(indices_df), function(z) {

      side_count(
        m         = m,
        row_index = indices_df$row[z],
        col_index = indices_df$col[z],
        cell_res  = cell_res,
        max_dist  = max_dist,
        func      = func,
        verbose   = FALSE
      )

    })

    # stop parallel processing clusters
    parallel::stopCluster(clust)

    return(side_lengths)

  }

}


# #' Retrieve indices of diagonal cells from a specific matrix cell in a direction
# #' @description Given a matrix and the row and column number for the cell of interest, return the index values for the diagonal cells in any direction.
# #' @param m matrix of interest
# #' @param row_index numeric row number of matrix cell to get indices of diagonals for
# #' @param col_index numeric column number of matrix cell to get indices of diagonals for
# #' @param direction character indicating which direction get diagonal indices for. Either "top_left", "top_right", "bottom_left", or "bottom_right". Default is "top_left".
# #' @param id_col logical. Whether a column stating the direction of the diagonal should be included. Default is FALSE, no direction column is added.
# #' @param verbose logical. Whether messages should output to console. Default is FALSE, no messages output.
# #' @return dataframe with columns for the row/column indexes of the diagonol cells. Each row of the dataframe indicates a diagonal cell
# index_diag <- function(
    #     m,
#     row_index = NULL,
#     col_index = NULL,
#     direction = "top_left",
#     id_col    = FALSE,
#     verbose   = FALSE
# ) {
#
#   # m <- rmat
#   # row_index <- row_nest$data[[1]]$row
#   # col_index <- row_nest$data[[1]]$col
#
#
#   # top left/north west diagonal
#   if(direction == "top_left") {
#
#     if(verbose == TRUE) {
#       message(paste0("Retrieving top left diaganol indices..."))
#     }
#
#     # check if row 1
#     if(row_index == 1) {
#
#       row_start = 1
#
#     } else {
#
#       row_start = row_index - 1
#
#     }
#
#     # check if col 1
#     if(col_index == 1) {
#
#       col_start = 1
#
#     } else {
#
#       col_start = col_index - 1
#
#     }
#
#     # sequence of northeast columns
#     col_seq <- seq(1, col_start, 1)
#     # col_seq <- seq(1, col_index - 1, 1)
#
#     # sequence of northeast rows
#     row_seq <- seq(row_start,  1,  -1)
#     # row_seq <- seq(row_index - 1,  1,  -1)
#
#     # minimum sequence
#     min_seq <- min(length(col_seq), length(row_seq))
#
#     # indices of diaganols
#     diag_df <- data.frame(
#       diag_row = row_seq[1:min_seq],
#       diag_col = rev(col_seq)[1:min_seq]
#     )
#
#   }
#   # top right/north east diagonal
#   if(direction == "top_right") {
#
#     if(verbose == TRUE) {
#       message(paste0("Retrieving top right diaganol indices..."))
#     }
#
#     # check if row 1
#     if(row_index == 1) {
#
#       row_start = 1
#
#     } else {
#
#       row_start = row_index - 1
#
#     }
#
#     # check if col 1
#     if(col_index == ncol(m)) {
#
#       col_start = col_index
#
#     } else {
#
#       col_start = col_index + 1
#
#     }
#
#     # sequence of northeast columns
#     col_seq <- seq(col_start, ncol(m), 1)
#     # col_seq <- seq(col_index + 1, ncol(m), 1)
#
#     # sequence of northeast rows
#     row_seq <- seq(row_start,  1,  -1)
#
#     # minimum sequence
#     min_seq <- min(length(col_seq), length(row_seq))
#
#     # indices of diaganols
#     diag_df <- data.frame(
#       diag_row = row_seq[1:min_seq],
#       diag_col = col_seq[1:min_seq]
#     )
#
#   }
#
#   # bottom left/south west diagonal
#   if(direction == "bottom_left") {
#
#     if(verbose == TRUE) {
#       message(paste0("Retrieving bottom left diaganol indices..."))
#     }
#
#     # check if row 1
#     if(row_index == nrow(m)) {
#
#       row_start = row_index
#
#     } else {
#
#       row_start = row_index + 1
#
#     }
#
#     # check if col 1
#     if(col_index == 1) {
#
#       col_start = col_index
#
#     } else {
#
#       col_start = col_index - 1
#
#     }
#
#     # sequence of northeast columns
#     col_seq <- seq(1, col_start, 1)
#     # col_seq <- seq(1, col_index - 1, 1)
#
#     # sequence of northeast rows
#     row_seq <- seq(row_start, nrow(m), 1)
#     # row_seq <- seq(row_index + 1, nrow(m), 1)
#
#     # minimum sequence
#     min_seq <- min(length(col_seq), length(row_seq))
#
#     # indices of diagonals
#     diag_row = row_seq[1:min_seq]
#     diag_col = rev(col_seq)[1:min_seq]
#
#     # mvect  <- m[cbind(diag_row, diag_col)]
#
#     check_length(m[cbind(diag_row, diag_col)])
#
#     # idx_mat <-
#     matrix(
#       data = c(row_seq[1:min_seq],  rev(col_seq)[1:min_seq]),
#       nrow = min_seq,
#       ncol = 2
#     )
#
#     # diag_df <-
#     diag_df <- data.frame(
#       diag_row = row_seq[1:min_seq],
#       diag_col = rev(col_seq)[1:min_seq]
#     )
#
#   }
#
#   # bottom right/south east diagonal
#   if(direction == "bottom_right") {
#
#     if(verbose == TRUE) {
#       message(paste0("Retrieving bottom right diaganol indices..."))
#     }
#
#     # check if row 1
#     if(row_index == nrow(m)) {
#
#       row_start = row_index
#
#     } else {
#
#       row_start = row_index + 1
#
#     }
#
#     # check if col 1
#     if(col_index == ncol(m)) {
#
#       col_start = col_index
#
#     } else {
#
#       col_start = col_index + 1
#
#     }
#
#     # sequence of northeast columns
#     col_seq <- seq(col_start, ncol(m), 1)
#     # col_seq <- seq(col_index + 1, ncol(m), 1)
#
#     # sequence of northeast rows
#     row_seq <- seq(row_start, nrow(m), 1)
#     # row_seq <- seq(row_index + 1, nrow(m), 1)
#
#     # minimum sequence
#     min_seq <- min(length(col_seq), length(row_seq))
#
#     # indices of diaganols
#     diag_df <- data.frame(
#       diag_row = row_seq[1:min_seq],
#       diag_col = col_seq[1:min_seq]
#     )
#
#   }
#
#   # if an ID column stating the direction should be added
#   if(id_col == TRUE) {
#
#     diag_df <-
#       diag_df %>%
#       dplyr::mutate(direction = direction)
#
#   }
#
#   return(diag_df)
# }

# #' Count diagonal cells of dataframe by direction
# #' @param m matrix of interest
# #' @param diag_indices dataframe with the indices of matrix diagnols from a reference cell and a "direction" column indicating in which direction does the diagonal extend out
# #' @param cell_res numeric cell resolution. Default is 2000.
# #' @param max_dist numeric max distance to calculate fetch distance. Default is 100000 (100km).
# #' @param verbose logical, whether to print loop messages or not. Default is FALSE.
# #' @return a dataframe with the diagonal cell count from a specified matrix cell, for each specified direction
# #' @importFrom dplyr group_by group_split mutate bind_rows
# count_diag_dir <- function(
    #     m,
#     diag_indices,
#     cell_res   = 2000,
#     max_dist   = 100000,
#     verbose    = FALSE
# ) {
#
#   # split diaganols into groups
#   diag_lst <-
#     diag_indices %>%
#     dplyr::group_by(direction) %>%
#     dplyr::group_split()
#
#   # lapply over each diaganol direction
#   diag_lst <- lapply(1:length(diag_lst), function(x) {
#
#     # count on the diagonals
#     count_diag(
#       m            = m,
#       diag_indices = diag_lst[[x]],
#       as_df        = TRUE,
#       cell_res     = cell_res,
#       max_dist     = max_dist,
#       verbose      = verbose
#     ) %>%
#       dplyr::mutate(
#         direction = diag_lst[[x]]$direction[1]
#       )
#   }) %>%
#     dplyr::bind_rows()
#
#   return(diag_lst)
#
# }

# #' Count the diagnol cells from a matrix index to end of matrix or a stopping value.
# #' @param m matrix to count diagnols
# #' @param diag_indices dataframe with the indices of matrix diagnols from a reference cell
# #' @param stop_value numeric value to stop counting cells if encountered. Default is 0
# #' @param cell_res numeric cell resolution. Default is 2000.
# #' @param max_dist numeric max distance to calculate fetch distance. Default is 100000 (100km).
# #' @param as_df logical, whether return a dataframe or a single count vector. Default is TRUE, returns a dataframe
# #' @param verbose logical, whether to print loop messages or not. Default is TRUE
# #' @return a dataframe or numeric vector indicating the number of diagonal cells before specific value comes up
# count_diag <- function(
    #     m,
#     diag_indices,
#     # row_vect,
#     # col_vect,
#     stop_value = 0,
#     cell_res   = 2000,
#     max_dist   = 100000,
#     as_df      = TRUE,
#     verbose    = TRUE
# ) {
#
#   # max count distance
#   max_count <- floor(max_dist/(cell_res))
#   # max_count <- floor(max_dist/(cell_res*sqrt(2)))
#
#   # start counter
#   counter <- 1
#
#   # loop through diagonal cells until max_count is reached or a 0 is found
#   while (counter <= max_count) {
#
#     # value of nth cell
#     cell_val <- m[row_vect[counter], col_vect[counter]]
#
#     # if cell value is 0, stop loop
#     if (cell_val == stop_value | is.na(cell_val)) {
#
#       if(verbose == TRUE) {
#         message(paste0("--- Stopping diagonal cell count --- "))
#         message(paste0("--- ", stop_value, " or end of matrix reached ---"))
#       }
#
#       break
#
#       # Last iteration message
#     } else if(counter >= max_count) {
#
#       message(paste0("--- Max cell count reached: ", max_count, " ---"))
#
#       # print message until last iteration
#     } else if(counter < max_count) {
#
#       message(paste0("Diagonal: ", counter, " --- (row ", row_vect[counter], ", col ", col_vect[counter], ")"))
#
#     }
#
#     counter <-  counter + 1
#
#   }
#
#   # remove last loop count
#   if (counter >= 1) {
#
#     counter <- counter - 1
#
#   }
#
#   # Whether to return a dataframe or a numberic vector
#   if(as_df == FALSE) {
#
#     return(counter)
#
#   } else {
#
#     # simple dataframe with count value
#     counter_df <- data.frame(
#       cell_count = counter
#     )
#
#     return(counter_df)
#
#   }
#
#   # # loop through diagnols and add to cell counter until a 0 is reached
#   # for (k in 1:nrow(diag_indices)) {
#   #
#   #   if(verbose == TRUE) {
#   #
#   #     message(paste0("Iteration: ", k, " - (row ", diag_indices$diag_row[k], ", col ", diag_indices$diag_col[k], ")"))
#   #
#   #   }
#   #
#   #   # if max number of pixels is reached
#   #   if(counter > max_count) {
#   #
#   #     if(verbose == TRUE) {
#   #       message(paste0("Counter reached max count:"))
#   #       message(paste0("Max count: ", max_count))
#   #       message(paste0("Counter: ", counter))
#   #     }
#   #
#   #     break
#   #
#   #   }
#   #
#   #   # value of kth southwest cell
#   #   cell_val <- m[diag_indices$diag_row[k], diag_indices$diag_col[k]]
#   #
#   #   # if cell value is 0, stop loop
#   #   if (cell_val == stop_value | is.na(cell_val)) {
#   #
#   #     if(verbose == TRUE) {
#   #       message(paste0(stop_value, " detected, stopping counter"))
#   #     }
#   #
#   #     break
#   #   }
#   #
#   #   counter <- counter + 1
#   #
#   # }
#
# }
#' Check boundaries of matrix column
#' @param m matrix of interest
#' @param row_index numeric row index
#' @param col_index numeric col index
#' @noRd
#' @keywords internal
#' @return list of cell indexes
check_cols <- function(m,
                       row_index,
                       col_index
) {

  # check if cells are on the edge of raster
  if(row_index == nrow(m)) {

    # if bottom row
    up_cells   <- 1:(row_index-1)
    down_cells <- (row_index):nrow(m)

  } else if(row_index == 1) {

    # if top row
    up_cells   <- 1:(row_index)
    down_cells <- (row_index+1):nrow(m)

  } else {

    # if any other row
    up_cells   <- 1:(row_index-1)
    down_cells <- (row_index+1):nrow(m)

  }

  updown <- list(up_cells, down_cells)

  return(updown)

}

#' Check boundaries of matrix row
#' @param m matrix of interest
#' @param row_index numeric row index
#' @param col_index numeric col index
#' @noRd
#' @keywords internal
#' @return list of cell indexes
check_rows <- function(m,
                       row_index,
                       col_index
) {

  # check if cells are on the edge of raster
  if(col_index == ncol(m)) {

    # if right most col
    left_cells  <- 1:(col_index-1)
    right_cells <- (col_index):ncol(m)

  } else if(col_index == 1) {

    # if left most col
    left_cells  <- 1:(col_index)
    right_cells <- (col_index+1):ncol(m)

  } else {

    # all other cols
    left_cells  <- 1:(col_index-1)
    right_cells <- (col_index+1):ncol(m)

  }

  leftright <- list(left_cells, right_cells)

  return(leftright)

}
#' Check length of consecutive non-zero cells in vector
#' @param vect numeric vector of 1s and 0s (can include NAs)
#' @param max_count numeric, indicating max number of pixels to count within maximum distance
#' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
#' @noRd
#' @keywords internal
#' @return Count of consecutive non-zero cells from a matrix cell, extending out in a given direction
check_length <- function(
    vect,
    max_count,
    verbose = FALSE
) {

  # vect      = vcol[updown[[2]]]
  # vect      = rev(vcol[updown[[1]]])

  # replace NA values w/ 2
  vect[is.na(vect)] <- 2

  # check for NAs
  na_check <- which(vect == 2)[1]

  # index of first 0 minus 1 gives number of cells before a 0 occurs
  check    <- which(vect == 0)[1]

  # cut off vector at max distance if vector > max pixel count
  if(length(vect) >= max_count) {

    # message(paste0("pass1"))

    # if(verbose == TRUE) {
    #   message(paste0("Trimming vector: ", length(vect), " ---> ", max_count))
    # }

    vect <- vect[1:max_count]

  }

  # if 0 is first value
  if(vect[1] == 0) {

    # message(paste0("pass2"))

    count <- 0

    # if(verbose == TRUE) {
    #   message(paste0("Cell count: ", count))
    # }

    return(count)

  }

  # if 0 is first value
  if(vect[1] == 2) {

    # message(paste0("pass3"))

    count <- as.integer(max_count)

    # count <- NA
    # if(verbose == TRUE) {
    #   message(paste0("Cell count: ", count))
    # }

    return(count)

  }

  # if no NAs found
  if(is.na(na_check)) {

    # message(paste0("pass4"))

    # check that 0s in vector
    if(!is.na(check)) {

      # message(paste0("pass5"))

      # make count index of first zero minus 1
      count <- check - 1

      # convert to integer
      count <- as.integer(count)

      # if(verbose == TRUE) {
      #   message(paste0("Cell count: ", count))
      # }

      return(count)

    }

  }

  # if no 0s found
  if(is.na(check)) {
    # message(paste0("pass6"))

    count <- as.integer(max_count)

    # count <- NA
    # if(verbose == TRUE) {
    #   message(paste0("Cell count: ", count))
    # }

    return(count)

  }

  # if both 0 AND NAs are found
  if(!is.na(check) & !is.na(na_check)) {

    # message(paste0("pass7"))

    # if 0 comes before NA
    if(check < na_check) {

      # message(paste0("pass8"))

      # make count index of first zero minus 1
      count <- check - 1

      # convert to integer
      count <- as.integer(count)

      # if(verbose == TRUE) {
      #   message(paste0("Cell count: ", count))
      # }

      return(count)

      # if NA before 0
    } else {

      # message(paste0("pass9"))

      count <- as.integer(max_count)

      # count <- NA
      # if(verbose == TRUE) {
      #   message(paste0("Cell count: ", count))
      # }

      return(count)

    }

  }


}

# #' Check length of consecutive non-zero cells in vector
# #' @param vect numeric vector of 1s and 0s (can include NAs)
# #' @param max_count numeric, indicating max number of pixels to count within maximum distance
# #' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
# #' @return Count of consecutive non-zero cells from a matrix cell, extending out in a given direction
# check_length <- function(
    #     vect,
#     max_count,
#     verbose = FALSE
# ) {
#
#   # replace NA values w/ 2
#   vect[is.na(vect)] <- 2
#
#   # if NA is first value
#   if(vect[1] == 2) {
#
#     # if number of NAs after first value of NA is longer than the max distance count, make count the max count
#     if(length(which(vect==2)) > max_count) {
#       # message(paste0("if 1"))
#       count <- max_count
#
#       # if number of NAs after first value of NA is shorter than the max distance count, make count NA
#     } else {
#       # message(paste0("if 2"))
#       count <- NA
#
#     }
#
#     # if first cell is a 0, give count of 0
#   } else if(vect[1] == 0) {
#     # message(paste0("else if 1"))
#     count <- 0
#
#     if(verbose == TRUE) {
#       # message(paste0("Cell count: ", count))
#     }
#
#   } else {
#
#     # check for NAs
#     na_check <- which(vect == 2)[1]
#
#     # index of first 0 minus 1 gives number of cells before a 0 occurs
#     check    <- which(vect == 0)[1]
#
#     # if no 0s or NAs are in vector make count length of vector (number of 1s in vector)
#     if(is.na(check) == TRUE & is.na(na_check) == TRUE) {
#       # message(paste0("else, if"))
#       count <- length(vect)
#
#     } else if(is.na(check) == TRUE & is.na(na_check) == FALSE) {
#       # message(paste0("else, else if 1"))
#       # na_check - 1
#       count <- min(length(which(vect==2)), max_count)
#
#     } else if(is.na(check) == FALSE & is.na(na_check) == TRUE) {
#       # message(paste0("else, else if 2"))
#       # make count index of first zero minus 1
#       count <- check - 1
#
#     } else if(is.na(check) == FALSE & is.na(na_check) == FALSE) {
#       # message(paste0("else, else if 3"))
#       # list 0 and NA positions
#       check_lst   <- c(check, na_check)
#
#       # which is closer to point
#       first_check <- which.min(check_lst)
#
#       if(first_check == 1) {
#         # message(paste0("else, else if 3 - if "))
#         count <- check - 1
#
#       } else {
#         # message(paste0("else, else if 3 - else "))
#         count <- min(length(which(vect==2)), max_count)
#
#       }
#
#     }
#
#
#   }
#
#   # convert to integer
#   count <- as.integer(count)
#
#   if(verbose == TRUE) {
#     message(paste0("Cell count: ", count))
#   }
#
#   return(count)
#
# }
#' Check run length encoding output for cell run of zeros
#' @param rle_vect rle object
#' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
#' @noRd
#' @keywords internal
#' @return Count of consecutive non-zero cells from a matrix cell, extending out in a given direction
check_zeros <- function(
    rle_vect,
    verbose = FALSE
) {

  # check if any zeros occur after calculating RLE on a vector
  if(0 %in% unique(rle_vect$values)) {

    if(verbose == TRUE) {
      message(paste0("Calculating consecutive non-zero cell count..."))
    }

    # count of up/down cells before a 0
    count   <- rle_vect$lengths[(which(rle_vect$values == 0)[1]) - 1]

  } else {

    if(verbose == TRUE) {
      message(paste0("Calculating consecutive non-zero cell count..."))
      message(paste0("No zeros found, all cells counted"))
    }

    # count of up/down cells when no zeros
    count   <- rle_vect$lengths

  }

  return(count)

}

#' Count number of consequetive cells around a matrix cell to the up/down/left/right of a matrix cell
#' @param m matrix of interest
#' @param row_index numeric row number of matrix cell to calculate cell count going up/down/left/right of the specified cell
#' @param col_index numeric column number of matrix cell to calculate cell count going up/down/left/right of the specified cell
#' @param cell_res numeric cell resolution. Default is 2000.
#' @param max_dist numeric maximum distance to calculate fetch to. Default is 100000m (100km).
#' @param func character. R function, either, mean, min, max, or sum. Default is mean calculates mean number of cells around cell
#' @param verbose logical, whether to print messages to console. Default is FALSE, no messages print
#' @noRd
#' @keywords internal
#' @return numeric indicating the average cell count in all directions or a dataframe with consequetive cell counts (non-zero cells) up/down/left/right of the specified cell in a matrix
side_count <- function(
    m,
    row_index = NULL,
    col_index = NULL,
    cell_res  = 2000,
    max_dist  = 100000,
    func      = "mean",
    verbose   = FALSE
) {

  # get function
  func <- match.fun(func)

  # max count distance
  max_count <- floor(max_dist/cell_res)

  # columns (up/down)
  # vector of current matrix col
  vcol <- m[, col_index]

  # return index list c(up, down) indices, checks if cells are on the edge of raster
  updown <- check_cols(
    m         = m,
    row_index = row_index,
    col_index = col_index
  )

  # count of up/down cells before a 0, if no zeros occur, returns length of consecutive 1s
  up_count   <- check_length(
    vect      = rev(vcol[updown[[1]]]),
    max_count = max_count
  )

  # up_count
  down_count <- check_length(
    vect      = vcol[updown[[2]]],
    max_count = max_count
  )

  # rows (left/right)
  # vector of current matrix row
  vrow <- m[row_index, ]

  # return index list c(left, right) indices, checks if cells are on the edge of raster
  leftright <- check_rows(
    m         = m,
    row_index = row_index,
    col_index = col_index
  )

  # count of up/down cells before a 0, if no zeros occur, returns length of consecutive 1s
  left_count  <- check_length(
    vect      = rev(vrow[leftright[[1]]]),
    max_count = max_count
  )

  right_count <- check_length(
    vect      = vrow[leftright[[2]]],
    max_count = max_count
  )

  # return average cell count around cell
  # if(as_df == FALSE) {

  # calculate average in all directions
  fcell_count <-  func(c(up_count, down_count, left_count, right_count), na.rm = T)

  return(fcell_count)

  # } else {
  #
  #   # count of cells on up/down/left/right side of cell before a 0 occurs
  #   side_df <- data.frame(
  #     row         = row_index,
  #     col         = col_index,
  #     direction   = c("up", "down", "left", "right"),
  #     cell_count  = c(up_count, down_count, left_count, right_count)
  #   )
  #
  #   # omit NAs
  #   side_df <- na.omit(side_df)
  #
  #   # select cell count and direction columns
  #   side_df <- side_df[c("cell_count", "direction")]
  #
  #   return(side_df)
  #
  # }

}

#' Calculate distances by cell counts
#' @param indices_df dataframe of indices
#' @param cell_res numeric cell resolution. Default is 2000.
#' @param max_dist numeric max distance to calculate fetch distance. Default is 100000 (100km).
#' @param func character R function to apply to distances from each cell. Either mean, min, max, or sum.
#' @param na_remove logical. argument passed to na.rm argument in R func. Default is TRUE, remove NAs
#' @importFrom dplyr group_by summarise
#' @noRd
#' @keywords internal
#' @return dataframe of summarized cell distances
cell_dist <- function(
    indices_df,
    cell_res  = 2000,
    max_dist  = 100000,
    func      = "mean",
    na_remove = TRUE
) {

  # suppress warnings, particularlly dplyr::summarise warning
  suppressWarnings({

    # match function to R function
    func <- match.fun(func)

    # # calculate distance from cell counts
    # cell_distance <-
    #   indices_df %>%
    #   dplyr::group_by(row, col) %>%
    #   dplyr::summarise(
    #     diag_cell_count = func(diag_cell_count, na.rm = na_remove),
    #     side_cell_count = func(side_cell_count, na.rm = na_remove)
    #   )

    # calculate distance from cell counts
    cell_distance <-
      dplyr::summarise(
        dplyr::group_by(indices_df, row, col),
        diag_cell_count = func(diag_cell_count, na.rm = na_remove),
        side_cell_count = func(side_cell_count, na.rm = na_remove)
        )
    # loop through each cell and calculate given function value
    # cell_distance <- lapply(1:nrow(cell_distance), function(i) {
    #   c(cell_distance$row[i], cell_distance$col[i],
    #     func(c(cell_distance$diag_cell_count[i], cell_distance$side_cell_count[i]), na.rm = na_remove)
    #   )
    # }) %>%
    #   do.call(rbind, .) %>%
    #   as.data.frame()
    cell_distance <- as.data.frame(
                          do.call(
                            rbind,
                            lapply(1:nrow(cell_distance), function(i) {
                              c(cell_distance$row[i], cell_distance$col[i],
                                func(c(cell_distance$diag_cell_count[i], cell_distance$side_cell_count[i]), na.rm = na_remove))
                              })
                            )
                        )

    # add cell resolution column
    cell_distance$cell_res <- cell_res

    # add distance column
    cell_distance$distance <- ifelse(cell_distance$V3*cell_distance$cell_res > max_dist, max_dist, cell_distance$V3*cell_distance$cell_res)

    # assign column names to cell_distance dataframe
    names(cell_distance) <- c("row", "col", "cell_count", "cell_res", "distance")

    # calculate distance from cell counts
    # cell_distance <-
    #   indices_df %>%
    #   dplyr::group_by(row, col) %>%
    #   dplyr::summarise(cell_count = func(cell_count, na.rm = na_remove)) %>%
    #   dplyr::ungroup() %>%
    #   dplyr::mutate(
    #     cell_res       = cell_res,
    #     # distance       = cell_count*cell_res
    #     distance       = dplyr::case_when(
    #       cell_count*cell_res > max_dist ~ max_dist,
    #       TRUE                           ~ cell_count*cell_res
    #     )
    #   )

  })

  return(cell_distance)

}

#' Create grid template
#' @description returns a template SpatRaster. All cell values are 1 so the template raster is not actually "empty
#' @param r SpatRaster, SpatVector or SF object with a Projected Coordinate System CRS.
#' @param res numeric resolution for template raster. Default is 2000, generates a template raster with a grid cell resolution of 2000x2000.
#' @importFrom terra rast crs ext vect
#' @importFrom methods is
#' @noRd
#' @keywords internal
#' @return SpatRaster of template grid
make_grid = function(
    r,
    res = 2000
) {

  # if SF object given, convert to SpatVector
  if(methods::is(r, "sf")) {

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

#' Prepare filepath vector/raster data for fetch calculations
#' @description Internal function for checking and aligning raster data to specifications for calculating fetch using get_fetch
#' @param r character, file path to a shp, gpkg, geojson, or tif file in a projected CRS
#' @param water_value numeric, value indicating the value of the water cells within the provided tif file. If a shp, gpkg, or geojson is given, this argument is ignored. Default is NA, such that NA cells are treated as water cells and all other cells are land cells.
#' @param res numeric resolution of the aligned output raster. Default is 2000.
#' @param verbose logical, whether messages should be printed. Default is FALSE, no messages print.
#' @importFrom terra vect rast
#' @importFrom methods is
#' @noRd
#' @keywords internal
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
#' @noRd
#' @keywords internal
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
#' @noRd
#' @keywords internal
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
#' @noRd
#' @keywords internal
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
      terra::resample(
        x      = r,
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

#' Internal r error handling message
#' @noRd
#' @keywords internal
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

#' Internal file path error handling message
#' @noRd
#' @keywords internal
#' @return character error message for invalid file paths
filepath_error_msg = function() {

  paste0(
    "Invalid 'r' file path argument, file path must point to a .shp, .gpkg, .geojson, or .tif file\n"
  )
}
