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
  levels(dm$position) <- gsub("V"," ", levels(dm$position))
  dm
}

#' Compute shannon information based on position and treatment
#' 
#' @param dframe data frame of peptide (or any other) sequences and some treatment factors
#' @param trt (vector of) character string(s) of treatment information
#' @param pos character string of position
#' @param elems character string of elements
#' @param k alphabet size: 4 for DNA/RNA sequences, 21 for standard amino acids
#' @param weight number of times each sequence is observed, defaults to 1 in case no weight is given
#' @return extended data frame with additional information of shannon info in bits and each elements contribution to the total information
#' @export
#' @examples
#' data(sequences)
#' dm2 <- splitSequence(sequences, "peptide")
#' dm3 <- calcInformation(dm2, pos="position", trt="class", elems="element", k=21)
#' # precursor to a logo plot:
#' library(ggplot2)
#' library(biovizBase)
#' 
calcInformation <- function(dframe, trt=NULL, pos, elems, k=4, weight = NULL, method="shannon") {
  if (is.null(weight)) dframe$wt <- 1
  else dframe$wt <- dframe[,weight]
  
  freqs <- ddply(dframe, c(trt, pos, elems), function(x) sum(x$wt))
## define implicit bindings for variables - not necessary though, since all of the variables do exist
  freq <- NA
  total <- NA
  
  names(freqs)[ncol(freqs)] <- "freq"
  freqByPos <- ddply(freqs, c(trt, pos), transform, total=sum(freq))
  if (method == "shannon") {
    freqByPos <- ddply(freqByPos, c(trt, pos), transform, info=-sum((freq/total*log(freq/total, base=2))[freq>0]))
    freqByPos$info <- -log(1/k, base=2) - with(freqByPos, info)
    freqByPos$bits <- with(freqByPos, freq/total*info)
  }  
  if (method == "relative")
    freqByPos$info <- with(freqByPos, freq/total)
  
  freqByPos
}

#' Logo plot
#' 
#' Simple logo plot of sequences. For more complicated sequence logos, such as with treatment comparisons or subsets see geom_logo.
#' @param sequences vector of text sequences, for which consensus logo is to be shown
#' @return ggplot2 object for simple sequence
#' @export 
#' @examples
#' data(sequences)
#' library(RColorBrewer)
#' cols <- rep(brewer.pal(12, name="Paired"),22)
#' logo(sequences$peptide) + aes(fill=element) + scale_fill_manual(values=cols)
logo <- function(sequences) {  
  position <- NA
  bits <- NA
  element <- NA
  
  dframe <- data.frame(seq=sequences)
  dm2 <- splitSequence(dframe, "seq")
  dm3 <- calcInformation(dm2, pos="position", elems="element", k=length(unique(dm3$element)))
  ggplot(dm3, aes(x=position, y=bits, group=element, label=element)) + geom_logo() 
}


#' Sequence logo plots.
#'
#' @export
#' @examples
#' \donttest{
#' data(sequences)
#' dm2 <- splitSequence(sequences, "peptide")
#' dm3 <- calcInformation(dm2, pos="position", elems="element", k=21)
#' library(biovizBase)
#' cols <- getBioColor(type="AA_ALPHABET")
#' library(RColorBrewer)
#' cols <- brewer.pal(10,"Paired")[c(1,2,7,8)]
#' data(aacids)
#' dm3b <- merge(dm3, aacids, by.x="element", by.y="AA", all.x=T)
#' ggplot(dm3b, aes(x=position, y=bits, group=element, label=element, fill=interaction(Polarity, Water))) + geom_logo() + scale_fill_manual(values=cols)
#' dm4 <- calcInformation(dm2, pos="position", elems="element", trt="class", k=21)
#' cols2 <- scales::alpha(cols, 0.8)
#' dm5 <- merge(dm4, aacids, by.x="element", by.y="AA", all.x=T)
#' ggplot(dm4, aes(x=class, y=bits, group=element, label=element, fill=element), alpha=0.8) + geom_logo() + scale_fill_manual(values=cols) + facet_wrap(~position, ncol=18)
#' ggplot(dm4, aes(x=position, y=bits, group=element, label=element, fill=element), alpha=0.8) + geom_logo() + scale_fill_manual(values=scales::alpha(cols, 0.8)) + facet_wrap(~class, ncol=1) + theme_bw()
#' ggplot(dm5, aes(x=class, y=bits, group=element, label=element, fill=interaction(Polarity, Water)), alpha=0.8) + geom_logo() + scale_fill_brewer("Amino-acids properties", palette="Paired") + facet_wrap(~position, ncol=18) + theme(legend.position="bottom") + xlab("")
#' ggplot(dm5, aes(x=class, y=bits, group=element, label=element, fill=interaction(Water, Polarity)), alpha=0.8) + geom_logo() + scale_fill_manual("Amino acids properties", values=cols) + facet_wrap(~position, ncol=18) + theme_bw() + theme(legend.position="bottom") + xlab("") + ylab("Shannon information in bits")
#' }

geom_logo <- function (mapping = NULL, data = NULL, stat = "logo", position = "identity", width = 0.9, alpha=0.25,
                       ...) {
  GeomLogo$new(mapping = mapping, data = data, stat = stat, 
               position = position, width= width, ...)
}

GeomLogo <- proto(ggplot2:::Geom, {
  objname <- "logo"
  
  reparameterise <- function(., df, params) {
    df
  }
  
  draw <- function(., data = data, scales, coordinates, ...) { 
#    print("draw")
    
#     common <- unique(data.frame(
#       colour = data$colour, 
#       size = data$size, 
#       linetype = data$linetype,
#       fill = "white", #alpha(data$fill, data$alpha),  
#       stringsAsFactors = FALSE
#     ))
    data(alphabet)
    letter <- subset(alphabet, group %in% unique(data$label))
    if (nrow(letter) < 1) {
    #  warning(paste("unrecognized letter in alphabet:", unique(data$label), collapse=","))
      letter <- alphabet[1,]
    }
    data$ROWID <- 1:nrow(data)
    letterpoly <- adply(data, .margins=1, function(x) {
      letter$x <- gglogo:::scaleTo(letter$x, fromRange=c(0,1), toRange=c(x$xmin, x$xmax))
      letter$y <- gglogo:::scaleTo(letter$y, toRange=c(x$ymin, x$ymax))
      letter$group <- interaction(x$ROWID, letter$group)
      letter
    })
#  browser()

    letterpoly$fill <- alpha("black", 0.7) # alpha(letterpoly$fill, letterpoly$alpha) #"white" 
    letterpoly$colour <- NA #"white"
    letterpoly$size <- 0.5    
# for drawing paths
    letterOutline <- letterpoly
    letterOutline$colour <- alpha("black", letterOutline$alpha)
    letterOutline$group <- paste(letterOutline$group, letterOutline$pathGroup, sep=".")
#    browser()
    ggname(.$my_name(), 
           gTree(children=gList(
             GeomRect$draw(data, scales, coordinates, ...),
             GeomPolygon$draw(letterpoly, scales, coordinates, ...) #,
       #      GeomPath$draw(letterOutline, scales, coordinates, ...)
           ))
    )    
  }
  
  guide_geom <- function(.) "polygon"
  
  draw_legend <- function(., data, ...)  {
    data <- aesdefaults(data, .$default_aes(), list(...))
    
    with(data, grobTree(
      rectGrob(gp = gpar(col = colour, fill = alpha(fill, alpha), lty = linetype)),
      linesGrob(gp = gpar(col = colour, lwd = size * .pt, lineend="butt", lty = linetype))
    ))
  }
  
  default_stat <- function(.) StatLogo
  default_pos <- function(.) PositionIdentity
  default_aes <- function(.) aes(weight=1, colour="grey80", fill="white", size=0.1, alpha = NA, shape = 16, linetype = "solid")
  required_aes <- c("x", "y", "group", "label")
  
})

#' calculation of all pieces necessary to plot a logo sequence plot
#' 
#'
#' @param scale if "area" (default), all vases have the same area (before trimming
#'   the tails). If "count", areas are scaled proportionally to the number of
#'   observations. If "width", all vases have the same maximum width.
#' @param na.rm If \code{FALSE} (the default), removes missing values with
#'    a warning. If \code{TRUE} silently removes missing values.
#'
#' @return A data frame with additional columns:
#'   \item{density}{density estimate}
#'   \item{fivenum}{five number summary for boxplots including a list of outliers if any}
#'   \item{scaled}{density estimate, scaled to maximum of 1}
#'   \item{count}{density * number of points - probably useless}
#'   \item{vasewidth}{density scaled for the vase plot, according to area, counts
#'                      or to a constant maximum width}
#'   \item{n}{number of points}
#'   \item{width}{width of vase bounding box}
#' @seealso \code{\link{geom_vase}} for examples, and \code{\link{stat_density}}
#'   for examples with data along the x axis.
#' @export
#' @examples
#' # See geom_logo for examples
#' # Generate data
#' data(sequences)
#' dm2 <- splitSequence(sequences, "peptide")
#' dm3 <- calcInformation(dm2, pos="position", elems="element", k=21)
#' ggplot(dm3, aes(x=position, y=bits, group=interaction(position, element))) + geom_logo()
stat_logo <- function (mapping = NULL, data = NULL, geom = "logo", position = "identity",
                       width = 0.9, drop="FALSE", na.rm = FALSE, ...) {
  StatLogo$new(mapping = mapping, data = data, geom = geom, position = position, width=width, drop=drop, 
               na.rm = na.rm, ...)
}

#' @export
StatLogo <- proto(ggplot2:::Stat, {
  objname <- "logo"
  
  calculate_groups <- function(., data, na.rm = FALSE, width = width, ...) {
#    print("calculate groups")
#     browser()
    data <- remove_missing(data, na.rm, "y", name = "stat_logo", finite = TRUE)
    data <- data[with(data, order(x, y)),]   
    data <- ddply(data, .(x), transform, 
                  ymax = cumsum(y))
    data$ymin <- with(data, ymax-y)
    data <- ddply(data, .(x), transform, 
                  ybase = max(ymin))
    data$ymin <- with(data, ymin-ybase)
    data$ymax <- with(data, ymax-ybase)
 #   browser()
#     xi <- unique(data[,c("x", "width")])
#     xi$w1 <- cumsum(xi$width)-xi$width 
#     xi$w2 <- cumsum(xi$width)
#     xi$w1 <- xi$w1*(max(data$x)-1)/max(xi$w1)
#     xi$w2 <- xi$w2*(max(data$x)-1)/max(xi$w2)
#     data$xmin <- with(data, xi$w1[x]+1)   
#     data$xmax <- with(data, xi$w2[x]+1)   
     data$xmin <- with(data, x - width/2)   
     data$xmax <- with(data, x + width/2)   
    
    .super$calculate_groups(., data, na.rm = na.rm, width = width, ...)
  }
  
  calculate <- function(., data,  scales, binwidth=NULL, origin=NULL, breaks=NULL, width=0.9,
                        na.rm = FALSE, ...) {
#    print("calculate for each group")
    #       browser()
    
    data
  }
  
  default_geom <- function(.) GeomLogo
  required_aes <- c("x", "y", "group", "label")
  
})

#' Cluster sequences according to selected positions 
#' 
#' @param sequences, 
#' @param pos vector of positions for the cluster. If NULL, the whole sequence is used
#' @param k number of clusters
#' @return a data frame containing the clustering in vector cl
#' @examples
#' dseq <- seqtree(subset(sequences, class=="gram +")$peptide, pos=c(4,7,24,27), k=4)
#' dm2 <- splitSequence(dseq, "subseq")
#' dm3 <- calcInformation(dm2, pos="position", trt="cl", elems="element", k=21, weight="Freq")
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
  extract <- function(sequences, pos){
    res <- ldply(pos, function(x) substr(sequences, x, x))
    unlist(llply(res, paste, collapse=""), use.names=FALSE)
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

