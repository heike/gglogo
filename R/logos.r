#' Reshape data set according to elements in sequences
#' 
#' prepare data set for plotting in a logo
#' @param dframe data frame of peptide (or any other) sequences and some treatment factors
#' @param sequences character string or index for the character vector of (peptide) sequence
#' @export
#' @examples
#' dm2 <- splitSequence(pm, "peptide")
splitSequence <- function(dframe, sequences) {
  seqs <- as.character(dframe[,sequences])
#  require(plyr)
  seqVars <-  data.frame(dframe, ldply(seqs, function(x) unlist(strsplit(x, split=""))))
#  require(reshape2)
  dm <- melt(seqVars, id.vars=names(dframe))
  names(dm) <- c(names(dframe), "position", "element")
  dm
}

#' Compute shannon information based on position and treatment
#' 
#' @param dframe data frame of peptide (or any other) sequences and some treatment factors
#' @param trt (vector of) character string(s) of treatment information
#' @param pos character string of position
#' @param elems character string of elements
#' @param k alphabet size: 4 for DNA/RNA sequences, 21 for standard amino acids
#' @return extended data frame with additional information of shannon info in bits and each elements contribution to the total information
#' @export
#' @examples
#' data(sequences)
#' dm2 <- splitSequence(sequences, "peptide")
#' dm2 <- splitSequence(pm, "peptide")
#' dm3 <- calcInformation(dm2, pos="position", trt="class", elems="element", k=21)
#' # precursor to a logo plot:
#' qplot(position,  data=dm3, facets=class~., geom="bar", weight=elinfo, fill=element) + scale_fill_manual(values=getBioColor(type="AA_ALPHABET"))
calcInformation <- function(dframe, trt, pos, elems, k=4) {

  freqs <- ddply(dframe, c(trt, pos, elems), nrow)
  names(freqs)[ncol(freqs)] <- "freq"
  freqByPos <- ddply(freqs, c(trt, pos), transform, total=sum(freq))
  freqByPos <- ddply(freqByPos, c(trt, pos), transform, info=-sum(freq/total*log(freq/total, base=2)))
  
  freqByPos$info <- -log(1/k, base=2) - with(freqByPos, info)
  freqByPos$elinfo <- with(freqByPos, freq/total*info)
  freqByPos
}

#' Sequence data
#' 
#' @name sequences
#' @title peptide sequence data
#' @description available through biovis redesign contest 2013, see http://www.biovis.net/year/2013/info/redesign-contest
#' published in Wong, B. Nat Methods 7, 889 (2011)
#' @docType data
#' @usage data(sequences)
NULL