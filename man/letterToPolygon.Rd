\name{letterToPolygon}
\alias{letterToPolygon}
\title{Convert an image file to a polygon}
\usage{
  letterToPolygon(ch, fontfamily = "Helvetica",
    fontsize = 576, tol = 1)
}
\arguments{
  \item{ch}{letter}

  \item{fontfamily}{}

  \item{fontsize}{}

  \item{tol}{tolerance}
}
\description{
  Convert an image file to a polygon
}
\examples{
letter <- letterToPolygon("ö", fontfamily="Garamond")
print(qplot(x, y, geom="polygon", data = letter, fill=I("black"), order=order, alpha=I(0.8))+coord_equal())
}

