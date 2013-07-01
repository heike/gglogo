\name{stat_logo}
\alias{stat_logo}
\title{calculation of all pieces necessary to plot a logo sequence plot}
\usage{
  stat_logo(mapping = NULL, data = NULL, geom = "logo",
    position = "identity", width = 0.9, drop = "FALSE",
    scale = "area", na.rm = FALSE, ...)
}
\arguments{
  \item{df}{dataframe}

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
\section{Aesthetics}{

}
\examples{
dmlabel <- stat_logo(dm3)
dmlabel$x <- as.numeric(dmlabel$position)-0.4
dmlabel$xmax <- dmlabel$x + 0.8

## very long example - should be in the code not an example
cols <- biovizBase::getBioColor(type="AA_ALPHABET")
#'
#' dm4 <- calcInformation(dm2, pos="position", trt=NULL, elems="element", k=21)
dm4$class <- 1
dmlabel <- stat_logo(dm4)
dmlabel$x <- as.numeric(dmlabel$position)-0.4
dmlabel$xmax <- dmlabel$x+0.8
base <- ggplot(aes(x, y), data=dmlabel) +
  geom_rect(aes(xmin=x, xmax=xmax, ymin=y, ymax=ymax, colour=element, fill=element)) +
  scale_fill_manual(values=cols) +
  scale_colour_manual(values=cols) +
  theme(legend.position="bottom") +
  geom_hline(yintercept=0, colour="grey20") +
  geom_hline(yintercept=-log(1/21, base=2), colour="grey20") +
  scale_x_continuous("position", breaks=1:39, labels=1:39)
ylab("bits")

# use Biovisbase for colors
data(alphabet)
dmletter <- merge(subset(dmlabel, elinfo > 0.25), alphabet, by.x="element", by.y="group")

base + geom_polygon(aes(x=x.x+x.y-0.1, y=y.x+elinfo*y.y, group=interaction(position, element), order=order),
                    alpha=0.9, fill="black", data=dmletter,
                    guides="none") +
  scale_shape_identity() +
  scale_size(range=6*c(0.5, max(dmlabel$freq)), guide="none")
combination of boxplot and 1d kernel density estimate along y axis, for logo plot.
}
\seealso{
  \code{\link{geom_vase}} for examples, and
  \code{\link{stat_density}} for examples with data along
  the x axis.
}

