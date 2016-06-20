#' Create a data set of polygons for a set of letters from a specified fontfamily 
#' 
#' @param letters set of characters for which polygons are to be created
#' @param font character describing the name of a font - use fonts() from package extrafont to check on available fonts
#' @param fontsize size of letter to be created - larger means higher resolution, but also bigger result sets
#' @param dim dimensions of the box in which the created letter is supposed to fit
#' @return data set compatible to work with geom_logo
#' @importFrom magrittr %>%
#' @importFrom purrr map_df
#' @export
#' @examples
#' new_alphabet <- createPolygons(c(letters, LETTERS, 0:9), font="Garamond")
#' # check that all letters and digits are nicely shaped:
#' qplot(x,y, geom="polygon", data=new_alphabet, facets=~group)
createPolygons <- function(letters, font, fontsize = 400, dim = c(720, 720)) {
  scale_01 <- function(x) {
    (x - min(x, na.rm=TRUE))/diff(range(x, na.rm=TRUE))
  }
  
  new_alphabet <- letters %>% map_df(.f = letterToPolygon, fontfamily = font, 
                fontsize = fontsize, dim=dim, .id="region")

  new_alphabet$region <- letters[as.numeric(new_alphabet$region)]
  new_alphabet$pathGroup <- with(new_alphabet, paste(region, group, sep="."))
  new_alphabet <- dplyr::rename(new_alphabet, group=region)
  
  new_alphabet <- new_alphabet %>% transform(x = scale_01(x), y=scale_01(y))
  new_alphabet[, c("x", "y", "order", "group", "pathGroup")]
}