\name{geom_logo}
\alias{geom_logo}
\title{Sequence logo plots.}
\usage{
  geom_logo(mapping = NULL, data = NULL, stat = "logo",
    position = "identity", width = 0.9, alpha = 0.9, ...)
}
\description{
  Sequence logo plots.
}
\examples{
\donttest{
data(sequences)
dm2 <- splitSequence(sequences, "peptide")
dm3 <- calcInformation(dm2, pos="position", elems="element", k=21)
library(biovizBase)
cols <- getBioColor(type="AA_ALPHABET")
ggplot(dm3, aes(x=position, y=elinfo, group=element, label=element, fill=element)) + geom_logo() + scale_fill_manual(values=cols)
dm4 <- calcInformation(dm2, pos="position", elems="element", trt="class", k=21)
ggplot(dm4, aes(x=class, y=elinfo, group=element, label=element, fill=element), alpha=0.8) + geom_logo() + scale_fill_manual(values=scales::alpha(cols, 0.8)) + facet_wrap(~position, ncol=18)
}
}

