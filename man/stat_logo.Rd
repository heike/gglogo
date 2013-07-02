\name{stat_logo}
\alias{stat_logo}
\title{calculation of all pieces necessary to plot a logo sequence plot}
\usage{
  stat_logo(mapping = NULL, data = NULL, geom = "logo",
    position = "identity", width = 0.9, drop = "FALSE",
    na.rm = FALSE, ...)
}
\arguments{
  \item{scale}{if "area" (default), all vases have the same
  area (before trimming the tails). If "count", areas are
  scaled proportionally to the number of observations. If
  "width", all vases have the same maximum width.}

  \item{na.rm}{If \code{FALSE} (the default), removes
  missing values with a warning. If \code{TRUE} silently
  removes missing values.}
}
\value{
  A data frame with additional columns:
  \item{density}{density estimate} \item{fivenum}{five
  number summary for boxplots including a list of outliers
  if any} \item{scaled}{density estimate, scaled to maximum
  of 1} \item{count}{density * number of points - probably
  useless} \item{vasewidth}{density scaled for the vase
  plot, according to area, counts or to a constant maximum
  width} \item{n}{number of points} \item{width}{width of
  vase bounding box}
}
\description{
  calculation of all pieces necessary to plot a logo
  sequence plot
}
\examples{
# See geom_logo for examples
# Generate data
data(sequences)
dm2 <- splitSequence(sequences, "peptide")
dm3 <- calcInformation(dm2, pos="position", elems="element", k=21)
ggplot(dm3, aes(x=position, y=elinfo, group=interaction(position, element))) + geom_logo()
}
\seealso{
  \code{\link{geom_vase}} for examples, and
  \code{\link{stat_density}} for examples with data along
  the x axis.
}

