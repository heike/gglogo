#' Create a data set of polygons for a set of letters from a specified fontfamily 
#' 
#' @param letters set of characters for which polygons are to be created
#' @param font character describing the name of a font - use fonts() from package extrafont to check on available fonts
#' @param fontsize size of letter to be created - larger means higher resolution, but also bigger result sets
#' @param dim dimensions of the box in which the created letter is supposed to fit
#' @param scale scale the values along the y axis to use result in geom_logo
#' @return data set compatible to work with geom_logo
#' @importFrom magrittr %>%
#' @importFrom purrr map_df
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#' @export
#' @examples
#' \donttest{
#' # check that all letters and digits are nicely shaped:
#' new_alphabet <- createPolygons(c("f", "g", "W", "*"), font="Garamond")
#' 
#' library(ggplot2)
#' qplot(x,y, geom="polygon", data=new_alphabet, facets=~group)
#' }
createPolygons <- function(letters, font, fontsize = 400, dim = c(720, 720), scale = FALSE) {
    region <- group <- x <- y <- pathGroup <- NULL
    
  scale_01 <- function(x) {
    (x - min(x, na.rm=TRUE))/diff(range(x, na.rm=TRUE))
  }

  new_alphabet <- letters %>% map_df(.f = letterToPolygon, fontfamily = font, 
                fontsize = fontsize, dim=dim, .id="region")

  new_alphabet$group <- letters[as.numeric(new_alphabet$region)]
  new_alphabet$pathGroup <- with(new_alphabet, paste(region, group, sep="."))

  if (scale) 
    new_alphabet <- new_alphabet %>% group_by(group) %>% mutate(x = scale_01(x), y=scale_01(y))
  else  
    new_alphabet <- new_alphabet %>% mutate(x = scale_01(x), y=scale_01(y))
  
  new_alphabet %>%
    select(x, y, order, group, pathGroup)
}
