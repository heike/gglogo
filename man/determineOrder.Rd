\name{determineOrder}
\alias{determineOrder}
\title{Determine order in which to pass through a set of points}
\usage{
  determineOrder(x, y)
}
\arguments{
  \item{x}{x coordinates}

  \item{y}{y coordinates}
}
\value{
  data frame with elements x, y, and order
}
\description{
  Greedy algorithm to connect points, with the idea that a
  point is connected by a line with the two points closest
  to each, that haven't yet been connected into the shape.
  Results depend on the starting point.
}

