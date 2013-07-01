#' Reshape data set according to elements in sequences
#' 
#' prepare data set for plotting in a logo
#' @param dframe data frame of peptide (or any other) sequences and some treatment factors
#' @param sequences character string or index for the character vector of (peptide) sequence
#' @export
#' @examples
#' data(sequences)
#' dm2 <- splitSequence(sequences, "peptide")
splitSequence <- function(dframe, sequences) {
  seqs <- as.character(dframe[,sequences])
  require(plyr)
  seqVars <-  data.frame(dframe, ldply(seqs, function(x) unlist(strsplit(x, split=""))))
  require(reshape2)
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
#' dm3 <- calcInformation(dm2, pos="position", trt="class", elems="element", k=21)
#' # precursor to a logo plot:
#' library(ggplot2)
#' library(biovizBase)
#' qplot(position,  data=dm3, facets=class~., geom="bar", weight=elinfo, fill=element) + scale_fill_manual(elements=getBioColor(type="AA_ALPHABET"))
#' qplot(position,  data=calcInformation(dm2, pos="position", trt=NULL, elems="element", k=21), 
#' geom="bar", weight=elinfo, fill=element) + scale_fill_manual(elements=getBioColor(type="AA_ALPHABET"))
calcInformation <- function(dframe, trt=NULL, pos, elems, k=4) {

  freqs <- ddply(dframe, c(trt, pos, elems), nrow)
  names(freqs)[ncol(freqs)] <- "freq"
  freqByPos <- ddply(freqs, c(trt, pos), transform, total=sum(freq))
  freqByPos <- ddply(freqByPos, c(trt, pos), transform, info=-sum(freq/total*log(freq/total, base=2)))
  
  freqByPos$info <- -log(1/k, base=2) - with(freqByPos, info)
  freqByPos$elinfo <- with(freqByPos, freq/total*info)
  freqByPos
}

#' Logo plot
#' 
#' @param dframe dataset
#' @param sequences
#' @examples
#' data(sequences)
#' dm2 <- splitSequence(sequences, "peptide")
#' dm3 <- calcInformation(dm2, pos="position", trt="class", elems="element", k=21)
#' logo(dm=stat_logo(dm3))
#' dm4 <- calcInformation(dm2, pos="position", trt=NULL, elems="element", k=21)
#' dm4$class <- 1
#' logo(dm=stat_logo(dm4)) + facet_wrap(~position, ncol=36)
logo <- function(dm) {  
  dmlabel <- stat_logo(dm)  
  require(ggplot2)
  library(biovizBase)
  cols <- biovizBase::getBioColor(type="AA_ALPHABET")
  
  
  base <- ggplot(aes(x, y), data=dmlabel) +
    geom_rect(aes(xmin=x, xmax=xmax, ymin=y, ymax=ymax, colour=element, fill=element)) + 
    facet_wrap(~position, ncol=12) + 
    scale_fill_manual(values=cols) + 
    scale_colour_manual(values=cols) + 
    theme(legend.position="bottom") + 
    geom_hline(yintercept=0, colour="grey20") + 
    geom_hline(yintercept=-log(1/21, base=2), colour="grey20") +
    scale_x_continuous("", labels=c("negative","positive"), breaks=c(1,2)) +
    ylab("bits")
  
  # use Biovisbase for colors
  data(alphabet)
  dmletter <- merge(subset(dmlabel, elinfo > 0.25), alphabet, by.x="element", by.y="group")
  
  base + geom_polygon(aes(x=x.x+x.y-0.1, y=y.x+elinfo*y.y, group=interaction(element,class), order=order), 
                      alpha=0.9, fill="black", data=dmletter, 
                      guides="none") + 
    scale_shape_identity() + 
    scale_size(range=6*c(0.5, max(dmlabel$freq)), guide="none")   
}

#' calculation of all pieces necessary to plot a logo sequence plot
#' 
#' @param df dataframe
#' @examples
#' dmlabel <- stat_logo(dm3)
#' dmlabel$x <- as.numeric(dmlabel$position)-0.4
#' dmlabel$xmax <- dmlabel$x + 0.8
#' 
#' ## very long example - should be in the code not an example
#' cols <- biovizBase::getBioColor(type="AA_ALPHABET")
#' #' 
#' #' dm4 <- calcInformation(dm2, pos="position", trt=NULL, elems="element", k=21)
#' dm4$class <- 1
#' dmlabel <- stat_logo(dm4)
#' dmlabel$x <- as.numeric(dmlabel$position)-0.4
#' dmlabel$xmax <- dmlabel$x+0.8
#' base <- ggplot(aes(x, y), data=dmlabel) +
#'   geom_rect(aes(xmin=x, xmax=xmax, ymin=y, ymax=ymax, colour=element, fill=element)) + 
#'   scale_fill_manual(values=cols) + 
#'   scale_colour_manual(values=cols) + 
#'   theme(legend.position="bottom") + 
#'   geom_hline(yintercept=0, colour="grey20") + 
#'   geom_hline(yintercept=-log(1/21, base=2), colour="grey20") +
#'   scale_x_continuous("position", breaks=1:39, labels=1:39)
#' ylab("bits")
#' 
#' # use Biovisbase for colors
#' data(alphabet)
#' dmletter <- merge(subset(dmlabel, elinfo > 0.25), alphabet, by.x="element", by.y="group")
#' 
#' base + geom_polygon(aes(x=x.x+x.y-0.1, y=y.x+elinfo*y.y, group=interaction(position, element), order=order), 
#'                     alpha=0.9, fill="black", data=dmletter, 
#'                     guides="none") + 
#'   scale_shape_identity() + 
#'   scale_size(range=6*c(0.5, max(dmlabel$freq)), guide="none")
stat_logo <- function(dm) {
  # assuming a discrete x variable
  # a numeric y variable
  # a grouping variable 
  # any treatment should either go into facetting (handled outside) or as part of discrete x
  
  dmlabel <- dm[with(dm, order(class, position, freq)),]
  
  dmlabel <- ddply(dmlabel, .(position, class), transform, 
                   ymax = cumsum(elinfo))
  dmlabel$y <- with(dmlabel, ymax-elinfo)
  dmlabel <- ddply(dmlabel, .(position, class), transform, 
                   y1 = max(y))
  dmlabel$y <- with(dmlabel, y-y1)
  dmlabel$ymax <- with(dmlabel, ymax-y1)
  
  dmlabel$x <- 0.6
  dmlabel$xmax <- 1.4
  dmlabel$x[dmlabel$class=="positive"] <- dmlabel$x[dmlabel$class=="positive"]+1
  dmlabel$xmax[dmlabel$class=="positive"] <- dmlabel$xmax[dmlabel$class=="positive"]+1
  dmlabel
}


