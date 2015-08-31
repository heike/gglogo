#' Reshape data set according to elements in sequences
#' 
#' prepare data set for plotting in a logo
#' @param dframe data frame of peptide (or any other) sequences and some treatment factors
#' @param sequences character string or index for the character vector of (peptide) sequence
#' @export
#' @import plyr
#' @import reshape2
#' @examples
#' data(sequences)
#' dm2 <- splitSequence(sequences, "peptide")
splitSequence <- function(dframe, sequences) {
  seqs <- as.character(dframe[,sequences])
  seqVars <-  data.frame(dframe, ldply(seqs, function(x) unlist(strsplit(x, split=""))))
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
#' @param method either "shannon" or "frequency" for Shannon information or relative frequency of element by position.
#' @return extended data frame with additional information of shannon info in bits and each elements contribution to the total information
#' @export
#' @examples
#' data(sequences)
#' dm2 <- splitSequence(sequences, "peptide")
#' dm3 <- calcInformation(dm2, pos="position", trt="class", elems="element", k=21)
#' # precursor to a logo plot:
#' library(ggplot2)
#' # library(biovizBase)
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
  if (method == "frequency")
    freqByPos$info <- with(freqByPos, freq/total)
  
  freqByPos
}

#' Logo plot
#' 
#' Simple logo plot of sequences. For more complicated sequence logos, such as with treatment comparisons or subsets see geom_logo.
#' @param sequences vector of text sequences, for which consensus logo is to be shown
#' @return ggplot2 object for simple sequence
#' @import ggplot2
#' @export 
#' @examples
#' data(sequences)
#' library(ggplot2)
#' library(RColorBrewer)
#' cols <- rep(brewer.pal(12, name="Paired"),22)
#' logo(sequences$peptide) + aes(fill=element) + scale_fill_manual(values=cols)
logo <- function(sequences) {  
  position <- NA
  bits <- NA
  element <- NA
  
  dframe <- data.frame(seq=sequences)
  dm2 <- splitSequence(dframe, "seq")
  dm3 <- calcInformation(dm2, pos="position", elems="element", k=length(unique(dm2$element)))
  ggplot(dm3, aes(x=position, y=bits, group=element, label=element)) + geom_logo() 
}


#' Sequence logo plots.
#'
#' @param mapping The aesthetic mapping, usually constructed with aes or aes_string. Only needs to be set at the layer level if you are overriding the plot defaults.
#' @param data A layer specific dataset - only needed if you want to override the plot defaults, 
#' @param stat The statistical transformation to use on the data for this layer, 
#' @param position The position adjustment to use for overlappling points on this layer, 
#' @param width maximum width of the letters, defaults to 0.9, 
#' @param alpha amount of alpha blending used for putting letters on top of rectangle, defaults to 0.25,
#' @param ... other arguments passed on to layer. This can include aesthetics whose values you want to set, not map. See layer for more details.
#' @export
#' @importFrom scales alpha
#' @examples
#' \donttest{
#' data(sequences)
#' dm2 <- splitSequence(sequences, "peptide")
#' dm3 <- calcInformation(dm2, pos="position", elems="element", k=21)
#' #library(biovizBase)
#' #cols <- getBioColor(type="AA_ALPHABET")
#' library(RColorBrewer)
#' cols <- brewer.pal(10,"Paired")[c(1,2,7,8)]
#' data(aacids)
#' dm3b <- merge(dm3, aacids, by.x="element", by.y="AA", all.x=T)
#' ggplot(dm3b, aes(x=position, y=bits, group=element, 
#'      label=element, fill=interaction(Polarity, Water))) +
#'      geom_logo() + scale_fill_manual(values=cols)
#' dm4 <- calcInformation(dm2, pos="position", elems="element", trt="class", k=21)
#' cols2 <- scales::alpha(cols, 0.8)
#' dm5 <- merge(dm4, aacids, by.x="element", by.y="AA", all.x=T)
#' ggplot(dm4, aes(x=class, y=bits, group=element, 
#'      label=element, fill=element), alpha=0.8) + 
#'      geom_logo() + scale_fill_manual(values=cols) + facet_wrap(~position, ncol=18)
#' ggplot(dm4, aes(x=position, y=bits, group=element, label=element, fill=element), alpha=0.8) + 
#'      geom_logo() + scale_fill_manual(values=scales::alpha(cols, 0.8)) + 
#'      facet_wrap(~class, ncol=1) + theme_bw()
#' ggplot(dm5, aes(x=class, y=bits, group=element, 
#'      label=element, fill=interaction(Polarity, Water)), alpha=0.8) + 
#'      geom_logo() + scale_fill_brewer("Amino-acids properties", palette="Paired") + 
#'      facet_wrap(~position, ncol=18) + theme(legend.position="bottom") + xlab("")
#' ggplot(dm5, aes(x=class, y=bits, group=element, 
#'      label=element, fill=interaction(Water, Polarity)), alpha=0.8) + 
#'      geom_logo() + scale_fill_manual("Amino acids properties", values=cols) + 
#'      facet_wrap(~position, ncol=18) + theme_bw() + theme(legend.position="bottom") + 
#'      xlab("") + ylab("Shannon information in bits")
#' }

geom_logo <- function (mapping = NULL, data = NULL, stat = "logo", position = "identity", width = 0.9, alpha=0.25,
                       ...) {
  GeomLogo$new(mapping = mapping, data = data, stat = stat, 
               position = position, width= width, ...)
}

#' @import proto
GeomLogo <- proto(ggplot2:::Geom, {
  objname <- "logo"
  
  reparameterise <- function(., df, params) {
    df
  }
  
  draw <- function(., data = data, scales, coordinates, ...) {     
    data(alphabet, envir = environment())
    
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

#' Calculation of all pieces necessary to plot a logo sequence plot
#' 
#'
#' @param mapping The aesthetic mapping, usually constructed with aes or aes_string. Only needs to be set at the layer level if you are overriding the plot defaults.
#' @param data A layer specific dataset - only needed if you want to override the plot defaults, 
#' @param geom The geometric object to use display the data,
#' @param position The position adjustment to use for overlappling points on this layer, 
#' @param width maximum width of the letters, defaults to 0.9, 
#' @param  ... other arguments passed on to layer. This can include aesthetics whose values you want to set, not map. See layer for more details.
#'
#' @return A proto object 
#' @export
#' @examples
#' # See geom_logo for examples
#' # Generate data
#' data(sequences)
#' library(ggplot2)
#' dm2 <- splitSequence(sequences, "peptide")
#' dm3 <- calcInformation(dm2, pos="position", elems="element", k=21)
#' ggplot(dm3, aes(x=position, y=bits, label=element, group=interaction(position, element))) + 
#'      geom_logo()
stat_logo <- function (mapping = NULL, data = NULL, geom = "logo", position = "identity",
                       width = 0.9,  ...) {
  StatLogo$new(mapping = mapping, data = data, geom = geom, position = position, width=width, ...)
}

#' Helper object
#' 
#' shouldn't be used by the user, but has to be exported.
#' @return  proto object
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

