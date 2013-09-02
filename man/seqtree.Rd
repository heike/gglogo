\name{seqtree}
\alias{seqtree}
\title{Cluster sequences according to selected positions}
\usage{
  seqtree(sequences, pos = NULL, k)
}
\arguments{
  \item{sequences,}{}

  \item{pos}{vector of positions for the cluster. If NULL,
  the whole sequence is used}

  \item{k}{number of clusters}
}
\value{
  a data frame containing the clustering in vector cl
}
\description{
  Cluster sequences according to selected positions
}
\examples{
dseq <- seqtree(subset(sequences, class=="gram +")$peptide, pos=c(4,7,24,27), k=4)
dm2 <- splitSequence(dseq, "subseq")
dm3 <- calcInformation(dm2, pos="position", trt="cl", elems="element", k=21, weight="Freq")
library(RColorBrewer)
cols <- brewer.pal(10,"Paired")[c(1,2,7,8)]
data(aacids)
dm3b <- merge(dm3, aacids, by.x="element", by.y="AA", all.x=T)
dm3c <- as.data.frame(xtabs(Freq ~cl, data=dm2)/36)
names(dm3c)[2] <- "width"
dm3b <- merge(dm3b, dm3c, by="cl")
dm3b$width <- dm3b$width/max(dm3b$width)
ggplot(dm3b, aes(x=cl, y=bits, group=element, label=element, fill=Polarity)) + theme_bw() + geom_logo(aes(width=width)) + scale_fill_manual(values=cols) + facet_wrap(~position)
}

