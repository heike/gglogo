\name{logo}
\alias{logo}
\title{Logo plot}
\usage{
  logo(dm)
}
\arguments{
  \item{dframe}{dataset}

  \item{sequences}{}
}
\description{
  Logo plot
}
\examples{
data(sequences)
dm2 <- splitSequence(sequences, "peptide")
dm3 <- calcInformation(dm2, pos="position", trt="class", elems="element", k=21)
logo(dm=stat_logo(dm3))
dm4 <- calcInformation(dm2, pos="position", trt=NULL, elems="element", k=21)
dm4$class <- 1
logo(dm=stat_logo(dm4)) + facet_wrap(~position, ncol=36)
}

