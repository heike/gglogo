\name{geom_logo}
\alias{geom_logo}
\title{Sequence logo plots.}
\usage{
  geom_logo(mapping = NULL, data = NULL, stat = "logo",
    position = "stack", ...)
}
\description{
  Sequence logo plots.
}
\section{Aesthetics}{
  \Sexpr[results=rd,stage=build]{ggplot2:::rd_aesthetics("geom",
  "logo")}
}
\examples{
\donttest{
data(sequences)
dm2 <- splitSequence(sequences, "peptide")
dm3 <- calcInformation(dm2, pos="position", elems="element", k=21)
ggplot(dm3, aes(x=position, y=elinfo, group=element)) + geom_logo()
}
}

