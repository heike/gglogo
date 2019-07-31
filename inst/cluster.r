#' Cluster sequences according to selected positions 
#' 
#' @param sequences, vector of sequences
#' @param pos vector of positions for the cluster. If NULL, the whole sequence is used
#' @param k number of clusters
#' @return a data frame containing the clustering in vector cl
#' @examples
#' dseq <- seqtree(subset(sequences, class=="gram +")$peptide, pos=c(4,7,24,27), k=4)
#' dm2 <- splitSequence(dseq, subseq)
#' dm3 <- calcInformation(dm2, pos = position, trt = cl, elems = element, k = 21, weight ="Freq")
#' library(RColorBrewer)
#' cols <- brewer.pal(10,"Paired")[c(1,2,7,8)]
#' data(aacids)
#' dm3b <- merge(dm3, aacids, by.x="element", by.y="AA", all.x=T)
#' dm3c <- as.data.frame(xtabs(Freq ~cl, data=dm2)/36)
#' names(dm3c)[2] <- "width"
#' dm3b <- merge(dm3b, dm3c, by="cl")
#' dm3b$width <- dm3b$width/max(dm3b$width)
#' ggplot(dm3b, aes(x=cl, y=bits, group=element, label=element, fill=Polarity)) + theme_bw() + geom_logo(aes(width=width)) + scale_fill_manual(values=cols) + facet_wrap(~position) 
seqtree <- function(sequences, pos=NULL, k) {
  Freq <- NA
  
  extract <- function(sequences, pos){
    res <- plyr::ldply(pos, function(x) substr(sequences, x, x))
    unlist(plyr::llply(res, paste, collapse=""), use.names=FALSE)
  }
  seqdist <- function(s1, s2) {
    x1 <- strsplit(as.character(s1), split="")[[1]]
    x2 <- strsplit(as.character(s2), split="")[[1]]
    sum(x1!=x2)/length(x1)
  }
  if (!is.null(pos)) subseq <- extract(sequences, pos)
  else subseq <- sequences
  
  dseq <- as.data.frame(table(subseq))
  dseq <- subset(dseq, Freq > 0)
  n <- length(dseq$subseq)
  d <- matrix(rep(0, n^2), ncol=n, nrow=n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d[j, i] <- seqdist(dseq$subseq[i], dseq$subseq[j])
      d[i, j] <- d[j, i]
    }
  }
  d <- as.dist(d)
  cl <- hclust(d)
  browser()
  dseq$cl <- cutree(cl, k=k)
  dseq
}
