\name{calcInformation}
\alias{calcInformation}
\title{Compute shannon information based on position and treatment}
\usage{
  calcInformation(dframe, trt = NULL, pos, elems, k = 4)
}
\arguments{
  \item{dframe}{data frame of peptide (or any other)
  sequences and some treatment factors}

  \item{trt}{(vector of) character string(s) of treatment
  information}

  \item{pos}{character string of position}

  \item{elems}{character string of elements}

  \item{k}{alphabet size: 4 for DNA/RNA sequences, 21 for
  standard amino acids}
}
\value{
  extended data frame with additional information of
  shannon info in bits and each elements contribution to
  the total information
}
\description{
  Compute shannon information based on position and
  treatment
}
\examples{
data(sequences)
dm2 <- splitSequence(sequences, "peptide")
dm3 <- calcInformation(dm2, pos="position", trt="class", elems="element", k=21)
# precursor to a logo plot:
library(ggplot2)
library(biovizBase)
qplot(position,  data=dm3, facets=class~., geom="bar", weight=elinfo, fill=element) + scale_fill_manual(values=getBioColor(type="AA_ALPHABET"))
qplot(position,  data=calcInformation(dm2, pos="position", trt=NULL, elems="element", k=21),
geom="bar", weight=elinfo, fill=element) + scale_fill_manual(values=getBioColor(type="AA_ALPHABET"))
}

