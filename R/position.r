#' Logo positioning for overlapping objects on top of one another.
#'
#' \code{position_classic} is stacking objects in an ordered fashion from largest to smallest element,
#' \code{position_logo} reverses the classic order and additionally shifts stacks downward to align the largest objects along their vertical minimum,
#' \code{position_fill} additionally standardises each stack to have unit
#' height.
#'
#' @family position adjustments
#' @seealso See \code{\link{geom_logo}} for
#'   more examples.
#' @export
#' @examples
#' \donttest{
#' library(ggplot2)
#' data(sequences)
#' 
#' # to make the most of comparisons, largest letters ar aligned along their minimum to
#' # work out the main sequence.
#' ggplot(data = ggfortify(sequences, "peptide", treatment = "class")) +
#'   geom_logo(aes(x = class, y = bits, fill = Water, label = element), position="logo") + 
#'   facet_wrap(~position)
#'
#' # in the classic logo plots letters are stacked in an ordered fasahion on top of each other
#' ggplot(data = ggfortify(sequences, "peptide", treatment = "class")) +
#'   geom_logo(aes(x = class, y = bits, fill = Water, label = element), position="classic") + 
#'   facet_wrap(~position)
#' }
position_logo <- function() {
  PositionLogo
}

#' @rdname position_logo
#' @format NULL
#' @usage NULL
#' @export
PositionLogo <- ggproto(
  "PositionLogo", Position,
  # requires one of c("ymax", "y"),
  
  setup_data = function(self, data, params) {
#    browser()
    
#    cat("setup_data in position logo\n")
    data = remove_missing(data, FALSE,
                          c("x", "y", "ymin", "ymax", "xmin", "xmax"), 
                          name = "position_logo")
    
    if (is.null(data$ymax) && is.null(data$y)) {
      message("Missing y and ymax in position = 'stack'. ",
              "Maybe you want position = 'identity'?")
      return(data)
    }
    
    if (!is.null(data$ymin) && !all(data$ymin == 0))
      warning("Stacking not well defined when ymin != 0", call. = FALSE)
    
    data
  },
  
  compute_panel = function(data, params, scales) {
#    browser()
    data <- data[with(data, order(PANEL, x, y)),]   
    data <- ddply(data, .(PANEL, x), transform, 
                  ymax = cumsum(y))
    data$ymin <- with(data, ymax-y)
    data <- ddply(data, .(PANEL, x), transform, 
                  ybase = max(ymin))
    data$ymin <- with(data, ymin-ybase)
    data$ymax <- with(data, ymax-ybase)
    
    #ggplot2:::collide(data, NULL, "position_logo", pos_logo)
    data
  }
)

#' @rdname position_logo
#' @export
position_classic <- function() {
  PositionClassic
}

#' @rdname position_logo
#' @format NULL
#' @usage NULL
#' @export
PositionClassic <- ggproto(
  "PositionClassic", Position,
  # requires one of c("ymax", "y"),
  
  setup_data = function(self, data, params) {
#    browser()
    
#    cat("setup_data in position logo\n")
    data = remove_missing(data, FALSE,
                          c("x", "y", "ymin", "ymax", "xmin", "xmax"), 
                          name = "position_classic")
    
    if (is.null(data$ymax) && is.null(data$y)) {
      message("Missing y and ymax in position = 'classic'. ",
              "Maybe you want position = 'identity'?")
      return(data)
    }
    
    if (!is.null(data$ymin) && !all(data$ymin == 0))
      warning("Stacking not well defined when ymin != 0", call. = FALSE)
    
    data
  },
  
  compute_panel = function(data, params, scales) {
#    browser()
    data <- data[with(data, order(PANEL, x, y)),]   
    data <- ddply(data, .(PANEL, x), transform, 
                  ymin = sum(y) - cumsum(y))
    data$ymax <- with(data, ymin+y)

    #ggplot2:::collide(data, NULL, "position_logo", pos_logo)
    data
  }
)
