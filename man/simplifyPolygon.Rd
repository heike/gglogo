\name{simplifyPolygon}
\alias{simplifyPolygon}
\title{Douglas-Peucker algorithm adjusted fo polyons}
\usage{
  simplifyPolygon(points, tol = 1)
}
\arguments{
  \item{points}{matrix of x and y points}

  \item{tol}{tolerance}
}
\description{
  Implementation of the Douglas-Peucker algorithm for line
  thinning
  http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
}

