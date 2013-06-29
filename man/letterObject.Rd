\name{letterObject}
\alias{letterObject}
\title{Convert a text element into an R object}
\usage{
  letterObject(ch, fontfamily = "Helvetica",
    fontsize = 576)
}
\arguments{
  \item{ch}{text to be converted, usually just a single
  letter}

  \item{fontfamily}{}

  \item{fontsize}{by default 576. If the resulting string
  exceeds the boundary of the matrix returned, reduced font
  size}
}
\value{
  three dimensional matrix of dimension 480 x 480 x 3 of
  the pixel values, black background and white letter
}
\description{
  Convert a text element into an R object
}
\examples{
plot(letterObject("g", fontfamily="Garamond", fontsize=400))
plot(letterObject("q", fontsize=400))
plot(letterObject("B"))
}

