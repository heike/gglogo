#' Convert sequence data to a format suitable for logo plots
#' 
#' @param data data frame with the sequences
#' @param sequences variable containing the sequences
#' @param treatment co-variate(s) used in collecting sequence data
#' @param method either "shannon" or "frequency" for Shannon information or relative frequency of element by position.
#' @param weight numeric variable of weights
#' @return data frame with position, element and information value
#' @export
#' @examples 
#' \donttest{
#' library(ggplot2)
#' data(sequences)
#' 
#' ggplot(data = ggfortify(sequences, peptide, treatment = class)) +
#'   geom_logo(aes(x = class, y = bits, fill = Water, label = element)) + 
#'   facet_wrap(~position)
#'   
#' ggplot(data = ggfortify(sequences, peptide, treatment = class)) +
#'   geom_logo(aes(x = class, y = bits, fill = Polarity, label = element)) + 
#'   facet_wrap(~position, ncol = 18) + 
#'   theme(legend.position = "bottom")
#'}
ggfortify <- function(data, sequences, treatment = NULL, weight = NULL, method = "shannon") {
  aacids <- NULL
  position <- NULL
  element <- NULL
  
  data(aacids, envir = environment())
  
  seqs <- enquo(sequences)
  treatment <- enquo(treatment)
  weight <- enquo(weight)
  
  dm2 <- splitSequence(data, !!seqs)
  
  k <- 4
  if (length(unique(dm2$element)) > 5) k <- 21
  
  dm3 <- calcInformation(dm2, pos = position, elems = element, trt = !!treatment, 
                         weight = !!weight, k = k, method = method)
  
  if (k == 21) # add peptide informatio only for peptides
    dm3 <- merge(dm3, aacids[,-1], by.x="element", by.y="AA", all.x = TRUE)
  
  return(dm3)
}
 
#' Reshape data set according to elements in sequences
#' 
#' prepare data set for plotting in a logo
#' @param dframe data frame of peptide (or any other) sequences and some treatment factors
#' @param sequences column containing the character vector of (peptide) sequence
#' @export
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_at
#' @importFrom dplyr mutate_all
#' @importFrom tidyr unnest
#' @examples
#' data(sequences)
#' dm2 <- splitSequence(sequences, peptide)
splitSequence <- function(dframe, sequences) {
  position <- NULL
  element <- NULL
  
  seqs <- enquo(sequences)
  
  dframe %>%
    mutate_at(vars(!!seqs), as.character) %>%
    mutate(position = list(1:nchar(!!seqs)[1]),
           element = strsplit(!!seqs, split = ""))  %>%
    unnest() %>%
    mutate_at(vars(position, element), as.factor)
}

#' Compute shannon information based on position and treatment
#' 
#' @param dframe data frame of peptide (or any other) sequences and some treatment factors
#' @param trt (vector of) character string(s) of treatment information
#' @param pos variable containing position
#' @param elems variable containing elements
#' @param k alphabet size: 4 for DNA/RNA sequences, 21 for standard amino acids
#' @param weight variable containing number of times each sequence is observed, defaults to 1 in case no weight is given
#' @param method either "shannon" or "frequency" for Shannon information or relative frequency of element by position.
#' @return extended data frame with additional information of shannon info in bits and each elements contribution to the total information
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom rlang quo_is_null
#' @export
#' @examples
#' data(sequences)
#' dm2 <- splitSequence(sequences, peptide)
#' dm3 <- calcInformation(dm2, pos = position, elems = element, trt = class, k = 21)
#' # precursor to a logo plot:
#' library(ggplot2)
#' # library(biovizBase)
#' 
calcInformation <- function(dframe, pos, elems, trt = NULL, weight = NULL, k = 4, method = "shannon") {
  `_default` <- NULL
  freq <- NULL
  total <- NULL
  info <- NULL
  bits <- NULL
  
  weightqo <- enquo(weight)
  if (quo_is_null(enquo(weight))) weightqo <- quo(`_default`)
  
  trtqo <- enquo(trt)
  if (quo_is_null(enquo(trt))) trtqo <- quo(`_default`)
  
  posqo <- enquo(pos)
  elemqo <- enquo(elems)
  
  freqs <- dframe %>%
    mutate(`_default` = 1) %>%
    group_by(!!trtqo, !!posqo, !!elemqo) %>%
    summarise(freq = sum(!!weightqo))
  
  freqByPos <- freqs %>%
    group_by(!!trtqo, !!posqo) %>%
    mutate(total = sum(freq))
  
  if (method == "shannon") {
    freqByPos <- freqByPos %>%
      group_by(!!trtqo, !!posqo) %>%
      mutate(info = -sum((freq/total * log(freq/total, base=2))[freq>0])) %>%
      mutate(info = -log(1/k, base = 2) - info, bits = freq/total * info) %>%
      mutate(info = bits) # why bits and not info?
  } else if (method == "frequency") {
    freqByPos$info <- with(freqByPos, freq/total)
  }
  
  final <- freqByPos %>%
    ungroup()
  
  if (quo_is_null(enquo(trt))) final <- final %>% select(-`_default`)
  
  return(final)
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
  dm2 <- splitSequence(dframe, seq)
  dm3 <- calcInformation(dm2, pos = position, elems = element, k = length(unique(dm2$element)))
  ggplot(dm3, aes(x=position, y=bits, group=element, label=element)) + geom_logo() 
}

#' @rdname stat_logo
#' @return  proto object
#' @export 
#' 
StatLogo <- ggproto("StatLogo", Stat,
  setup_data = function(data, params) {
  #  cat("setup_data in stat logo\n")
    
    data <- remove_missing(data, na.rm=TRUE, "y", name = "stat_logo", finite = TRUE)
    data <- data[with(data, order(PANEL, x, y)),]   
    data$xmin <- with(data, x - params$width/2)   
    data$xmax <- with(data, x + params$width/2)
    data$ymin <- 0
    data$ymax <- data$y

    data
  },                    
  compute_group = function(data, scales, params, na.rm = FALSE, width = 0.9) {
     data
  },
  required_aes = c("x", "y")
)

#' Calculation of all pieces necessary to plot a logo sequence plot
#' 
#'
#' @param mapping The aesthetic mapping, usually constructed with aes or aes_string. Only needs to be set at the layer level if you are overriding the plot defaults.
#' @param data A layer specific dataset - only needed if you want to override the plot defaults, 
#' @param geom The geometric object to use display the data,
#' @param position The position adjustment to use for overlappling points on this layer, 
#' @param show.legend Whether to show the legend or not
#' @param inherit.aes Whether to inherit the aes or not
#' @param width maximum width of the letters, defaults to 0.9, 
#' @param na.rm Whether to remove NAs or not
#' @param  ... other arguments passed on to layer. This can include aesthetics whose values you want to set, not map. See layer for more details.
#'
#' @return A proto object 
#' @export
#' @importFrom grid gpar
#' @examples
#' # See geom_logo for examples
#' # Generate data
#' data(sequences)
#' library(ggplot2)
#' 
#' ggplot(data = ggfortify(sequences, peptide)) + 
#'   geom_logo(aes(x=position, y=bits, label=element, 
#'                 group=interaction(position, element)), 
#'             alpha=0.5)
stat_logo <- function(mapping = NULL, data = NULL, geom = "logo",
                      position = "logo", show.legend = NA, inherit.aes = TRUE, width = 0.9, na.rm = TRUE, ...) {
    layer(
        stat = StatLogo, data = data, mapping = mapping, geom = geom, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(width = width, height = 0.1, na.rm = na.rm, ...)
    )
}

#' @rdname geom_logo
#' @export
#' @importFrom grid grobTree
#' @importFrom plyr adply
GeomLogo <- ggproto("GeomLogo", Geom,
  required_aes = c("x", "y", "group", "label"),
  default_aes = aes(weight = 1, colour = "grey80", fill = "white", size = 0.1, alpha = 0.25, width = 0.9, shape = 16, linetype = "solid"),
  draw_key = draw_key_rect,
  
  draw_panel = function(data, panel_scales, coord, alphabet = NULL) {
    if (is.null(alphabet)) data(alphabet, envir = environment())
    if (is.character(alphabet)) {
      if (length(grep(alphabet, extrafont::fonts())) != 1) stop(sprintf("%s is not a font. Use extrafont::fonts() for an overview of all available fonts.", alphabet))
      
      cat(paste0("creating polygons for font ", alphabet, "\n"))
      alphabet <- createPolygons(c(LETTERS, letters, 0:9), font = alphabet)
    }
#    browser()
    #save(data, file = "testdata.RData")
    #save(panel_scales, file = "panelscales.RData")
    #save(coord, file = "coord.RData")
    
    
    data$ROWID <- 1:nrow(data)
    letterpoly <- adply(data, .margins=1, function(xx) {
      letter <- subset(alphabet, group %in% unique(xx$label))
       if (nrow(letter) < 1) {
      #     warning(paste("unrecognized letter in alphabet:", paste(unique(data$label), collapse=",")))
         letter <- alphabet[1,]
       }
      letter$x <- gglogo:::scaleTo(letter$x, fromRange=c(0,1), toRange=c(xx$xmin, xx$xmax))
      letter$y <- gglogo:::scaleTo(letter$y, toRange=c(xx$ymin, xx$ymax))
      letter$group <- interaction(xx$ROWID, letter$group)
      letter
    })

    # get alpha in rectangle, not the letters
    letterpoly$alpha <- 0.7 # hard re-write
    letterpoly$fill <- alpha("black", 0.7) # alpha(letterpoly$fill, letterpoly$alpha) #"white" 
    letterpoly$colour <- NA #"white"
    letterpoly$size <- 0.5    

    grobTree(
      GeomRect$draw_panel(data, panel_scales, coord),
      GeomPolygon$draw_panel(letterpoly, panel_scales, coord)
    )
  },
  
  draw_legend = function(data, ...)  {
    with(data, grobTree(
      rectGrob(gp = gpar(col = colour, fill = alpha(fill, alpha), lty = linetype)),
      linesGrob(gp = gpar(col = colour, lwd = size * .pt, lineend="butt", lty = linetype))
    ))
  }
)

#' Sequence logo plots.
#'
#' @param mapping The aesthetic mapping, usually constructed with aes or aes_string. Only needs to be set at the layer level if you are overriding the plot defaults.
#' @param data A layer specific dataset - only needed if you want to override the plot defaults, 
#' @param stat The statistical transformation to use on the data for this layer, 
#' @param position The position adjustment to use for overlappling points on this layer, 
#' @param show.legend Whether to show the legend or not
#' @param inherit.aes Whether to inherit the aes or not
#' @param width maximum width of the letters, defaults to 0.9, 
#' @param alpha amount of alpha blending used for putting letters on top of rectangle, defaults to 0.25,
#' @param alphabet Specifies which alphabet is used in rendering the logo. alphabet can be a dataframe (output from createPolygons), a character specifying a font or NULL. If NULL, the default alphabet set is used (based on Helvetica). Use output from `createPolygons` to generate alphabet polygons for a different font.
#' @param na.rm Whether to remove NAs or not
#' @param ... other arguments passed on to layer. This can include aesthetics whose values you want to set, not map. See layer for more details.
#' @export
#' @examples
#' \donttest{
#' library(ggplot2)
#' data(sequences)
#' ggplot(data = ggfortify(sequences, peptide)) +      
#'   geom_logo(aes(x=position, y=bits, group=element, 
#'      label=element, fill=interaction(Polarity, Water)),
#'      alpha = 0.6)  +
#'   scale_fill_brewer(palette="Paired") +
#'   theme(legend.position = "bottom")
#'   
#' ggplot(data = ggfortify(sequences, peptide, treatment = class)) + 
#'   geom_logo(aes(x=class, y=bits, group=element, 
#'      label=element, fill=element)) + 
#'   facet_wrap(~position, ncol=18) +
#'   theme(legend.position = "bottom")
#'   
#' ggplot(data = ggfortify(sequences, peptide, treatment = class)) + 
#'   geom_logo(aes(x=position, y=bits, group=element, label=element, fill=element)) + 
#'   facet_wrap(~class, ncol=1) + theme_bw()
#'   
#' ggplot(data = ggfortify(sequences, peptide, treatment = class)) + 
#'   geom_logo(aes(x=class, y=bits, group=element, 
#'                 label=element, fill=interaction(Polarity, Water))) + 
#'   scale_fill_brewer("Amino-acids properties", palette="Paired") + 
#'   facet_wrap(~position, ncol=18) + 
#'   theme(legend.position="bottom") + 
#'   xlab("") + ylab("Shannon information in bits")
#'   
#' }
geom_logo <- function (mapping = NULL, data = NULL, stat = "logo", position = "logo", 
                       show.legend = NA, inherit.aes = TRUE, width = 0.9, alpha = 0.6,
                       na.rm = TRUE, alphabet=NULL, ...) {
    layer(
        geom = GeomLogo, mapping = mapping,  data = data, stat = stat, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(width = width, na.rm = na.rm, alpha = alpha,  alphabet=alphabet, ...)
    )
}
