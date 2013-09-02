#' Douglas Pecker algorithm for line thinning
#' 
#' Implementation of the Douglas-Peucker algorithm for line thinning
#' http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
#' @param points matrix of x and y points
#' @param tol tolerance
simplify_rec <- function(points, tol = 0.01) {
  n <- nrow(points)
  if (n <= 2) return(unique(c(1, n)))
  dist <- with(points, point_line_dist(x, y, x[1], y[1], x[n], y[n]))
  
  if (max(dist, na.rm = T) > tol) {
    furthest <- which.max(dist)
    c(
      simplify_rec(points[1:(furthest), ], tol),
# don't include furthest into list of essential points twice, but base the computation on it
      simplify_rec(points[(furthest):n, ], tol)[-1] + furthest-1 
    )
  } else {
    return(c(1,n))
  }  
}

#' Distance between point and line
#' 
#' Compute distance between point given as (px, py) and line spanned by points (lx1, ly1) and (lx2, ly2).
#' From http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
#' @param px x coordinate of  point outside the 
#' @param py y coordinate of  point
#' @param lx_1, x  coordinate of 1st point spanning a line
#' @param ly_1, y  coordinate of 1st point spanning a line
#' @param lx_2, x  coordinate of 2nd point spanning a line
#' @param ly_2, y  coordinate of 2nd point spanning a line
point_line_dist <- function(px, py, lx_1, ly_1, lx_2, ly_2) {
  abs((lx_2 - lx_1) * (ly_1 - py) - (lx_1 - px) * (ly_2 - ly_1)) /
    sqrt((lx_2 - lx_1) ^ 2 + (ly_2 - ly_1) ^ 2)
}

# Precompute all tolerances so that we can post-select quickly
compute_tol <- function(points, offset = 0) {
  n <- nrow(points)
  if (n <= 2) {
    c()
  } else if (n == 3) {
    with(points,
      point_line_dist(long[2], lat[2], long[1], lat[1], long[3], lat[3]))
  } else {
    dist <- with(points, 
      point_line_dist(long[2:(n-1)], lat[2:(n-1)], long[1], lat[1], long[n], lat[n])
    )
  
    furthest <- which.max(dist)
    if (length(furthest) == 0) browser()
    c(
      compute_tol(points[1:(furthest + 1), ], offset),
      dist[furthest],
      compute_tol(points[(furthest + 1):n, ], furthest + offset)
    )
  }
}