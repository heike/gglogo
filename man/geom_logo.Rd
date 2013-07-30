\name{geom_logo}
\alias{geom_logo}
\title{Sequence logo plots.}
\usage{
  geom_logo(mapping = NULL, data = NULL, stat = "logo",
    position = "identity", width = 0.9, alpha = 0.25, ...)
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
library(RColorBrewer)
cols <- brewer.pal(10,"Paired")[c(1,2,7,8)]
data(aacids)
dm3b <- merge(dm3, aacids, by.x="element", by.y="AA", all.x=T)
ggplot(dm3b, aes(x=position, y=bits, group=element, label=element, fill=interaction(Polarity, Water))) + geom_logo() + scale_fill_manual(values=cols)
dm4 <- calcInformation(dm2, pos="position", elems="element", trt="class", k=21)
cols2 <- scales::alpha(cols, 0.8)
dm5 <- merge(dm4, aacids, by.x="element", by.y="AA", all.x=T)
ggplot(dm4, aes(x=class, y=bits, group=element, label=element, fill=element), alpha=0.8) + geom_logo() + scale_fill_manual(values=cols) + facet_wrap(~position, ncol=18)
ggplot(dm4, aes(x=position, y=bits, group=element, label=element, fill=element), alpha=0.8) + geom_logo() + scale_fill_manual(values=scales::alpha(cols, 0.8)) + facet_wrap(~class, ncol=1) + theme_bw()
ggplot(dm5, aes(x=class, y=bits, group=element, label=element, fill=interaction(Polarity, Water)), alpha=0.8) + geom_logo() + scale_fill_brewer("Amino-acids properties", palette="Paired") + facet_wrap(~position, ncol=18) + theme(legend.position="bottom") + xlab("")
ggplot(dm5, aes(x=class, y=bits, group=element, label=element, fill=interaction(Water, Polarity)), alpha=0.8) + geom_logo() + scale_fill_manual("Amino acids properties", values=cols) + facet_wrap(~position, ncol=18) + theme_bw() + theme(legend.position="bottom") + xlab("") + ylab("Shannon information in bits")
}
}

