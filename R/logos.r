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
#' @export
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
#' ggplot(dm3, aes(x=position, y=elinfo, group=element, label=element, fill=element)) + geom_logo() + scale_fill_manual(values=cols)
#' dm4 <- calcInformation(dm2, pos="position", elems="element", trt="class", k=21)
#' ggplot(dm4, aes(x=class, y=elinfo, group=element, label=element, fill=element), alpha=0.8) + geom_logo() + scale_fill_manual(values=scales::alpha(cols, 0.8)) + facet_wrap(~position, ncol=18)
#' ggplot(dm4, aes(x=position, y=elinfo, group=element, label=element, fill=element), alpha=0.8) + geom_logo() + scale_fill_manual(values=scales::alpha(cols, 0.8)) + facet_wrap(~class, ncol=1) + theme_bw()
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
      warning(paste("unrecognized letter in alphabet:", unique(data$label), collapse=","))
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
  default_aes <- function(.) aes(weight=1, colour="grey20", fill="white", size=0.1, alpha = NA, shape = 16, linetype = "solid")
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
#' ggplot(dm3, aes(x=position, y=elinfo, group=interaction(position, element))) + geom_logo()
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
    # browser()
    data <- remove_missing(data, na.rm, "y", name = "stat_logo", finite = TRUE)
    data <- data[with(data, order(x, y)),]   
    data <- ddply(data, .(x), transform, 
                  ymax = cumsum(y))
    data$ymin <- with(data, ymax-y)
    data <- ddply(data, .(x), transform, 
                  ybase = max(ymin))
    data$ymin <- with(data, ymin-ybase)
    data$ymax <- with(data, ymax-ybase)   
    data$xmin <- with(data, x-width/2)   
    data$xmax <- with(data, x+width/2)   
    
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




