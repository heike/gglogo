\name{letterToPolygon}
\alias{letterToPolygon}
\title{Convert an image file to a polygon}
\usage{
  letterToPolygon(ch, fontfamily = "Helvetica",
    fontsize = 576, tol = 1, dim = c(480, 480),
    threshold = 0.5, var = "red")
}
\arguments{
  \item{ch}{letter}

  \item{fontfamily}{}

  \item{fontsize}{}

  \item{tol}{tolerance}

  \item{dim}{vector of length two specifying width and
  height (in pixels) of the temporary jpg file created for
  the letter. Defaults to 480 x 480 pixels.}

  \item{threshold}{}

  \item{variable}{one of "red", "green", "blue".}
}
\description{
  Convert an image file to a polygon
}
\examples{
letter <- letterToPolygon("ö", fontfamily="Garamond")
print(qplot(x, y, geom="polygon", data = letter, fill=I("black"), order=order, alpha=I(0.8))+coord_equal())
}

