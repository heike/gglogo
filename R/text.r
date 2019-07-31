#' Convert a text element into an R object
#' 
#' @param ch text to be converted, usually just a single letter
#' @param fontfamily R has a few default fonts that are always available, such as e.g. Helvetica, Arial, Courier New, and Garamond. Other fonts might be available depending on the platform used. 
#' @param fontsize by default 576. If the resulting string exceeds the boundary of the matrix returned, reduced font size
#' @param dim vector of length two specifying width and height (in pixels) of the temporary jpg file created for the letter. Defaults to 480 x 480 pixels.
#' @return three dimensional matrix of dimension 480 x 480 x 3 of the pixel values, black background and white letter 
#' @importFrom jpeg readJPEG
#' @importFrom grDevices dev.off jpeg
#' @importFrom grid grid.newpage
#' @importFrom grid grid.rect
#' @importFrom grid grid.text
#' @export
#' @examples
#' \donttest{
#' plot(letterObject("g", fontfamily="Garamond", fontsize=400))
#' plot(letterObject("q", fontsize=400))
#' plot(letterObject("B"))
#' }
letterObject <- function(ch, fontfamily="Helvetica", fontsize=576, dim=c(480, 480)) {
  fname <- tempfile(pattern = "file", fileext=".jpg")
  jpeg(filename=fname, width=dim[1], height=dim[2])
  grid.newpage()
  grid.rect(x = 0, y=0, width=3, height=3,
            gp=gpar(fill="black"), draw = TRUE, vp = NULL)
  
  grid.text(ch, 0.5,0.5, gp=gpar(fontsize=fontsize, fontfamily=fontfamily, col="white"))
  dev.off()
  readJPEG(fname)
}

scaleTo <- function(x, fromRange=range(x), toRange=c(0,1)) {
  x <- as.numeric(as.character(x))
  (x-fromRange[1])/diff(fromRange)*diff(toRange) + toRange[1]
}

#' Convert image matrix to a data frame
#' 
#' S3 method to create a data frame from a 3d array. 
#' @param model array as e.g. returned from readJPEG
#' @param data not used by this method
#' @param ... not used by this method
#' @return data frame 
#' \itemize{
#' \item x x coordinate in pixels
#' \item y y coordinate in pixels (usually negative)
#' \item red number vector in (0,1) describing the amount of red of the pixel in an RGB model
#' \item green number vector in (0,1) describing the amount of green of the pixel in an RGB model
#' \item blue number vector in (0,1) describing the amount of blue of the pixel in an RGB model
#' }
#' @export
fortify<-function(model, data, ...){
  UseMethod("fortify")
}

#' @rdname fortify
#' @export
fortify.default <- function(model, data, ...) {
  dims <- dim(model)
  imdf <- adply(model, .margins=1, function(x) x)
  imdf$x <- rep(1:dims[2], length=nrow(imdf)) 
  names(imdf) <- c("y", "red", "green", "blue", "x")
  imdf$y <- -as.numeric(as.character(imdf$y))
  imdf[,c("x", "y", "red", "green", "blue")]
}

#' Determine boundary between foreground and background in an image
#' 
#' @param imdf dataframe describing a pixellated image in x and y. Has to have columns x, y, and var
#' @param var dimension along which foreground and background of a shape in the image are well separated. Usually one of 'red', 'green', or 'blue', but could be extended to any other numerical variable.
#' @param threshold value specifying the cutoff along variable var. Values of var higher than the threshold are considered to belong to the foreground.
#' @return subset of data frame imdf consisting of just boundary points:
#' \itemize{
#' \item x x coordinate in pixels
#' \item y y coordinate in pixels (usually negative)
#' \item red number vector in (0,1) describing the amount of red of the pixel in an RGB model
#' \item green number vector in (0,1) describing the amount of green of the pixel in an RGB model
#' \item blue number vector in (0,1) describing the amount of blue of the pixel in an RGB model
#' }
#' @importFrom stats na.omit
#' @importFrom plyr .
#' @importFrom plyr ddply
#' @export
getOutline <- function(imdf, var="red", threshold=0.5) {
  stopifnot(c("x", "y") %in% names(imdf))
  ## define implicit bindings for variables - not necessary though, 
  ## since all of the variables do exist at this point
  x <- NA
  y <- NA

  
  edgesY <- ddply(imdf, .(y), function(dframe) {
    idx <- which(dframe[, var] > threshold)
    dx <- diff(sort(dframe$x[idx])) 
    nintervals <- sum(dx>1)+1
    jdx <- which(dx > 1)
    start <- idx[c(1, jdx+1)]
    end <- idx[c(jdx, length(idx))]
    dframe[c(start, end),]
  })
  
  edgesX <- ddply(imdf, .(x), function(dframe) {
    idx <- which(dframe[, var] > threshold)
    dx <- diff(sort(-dframe$y[idx])) 
    nintervals <- sum(dx>1)+1
    jdx <- which(dx > 1)
    start <- idx[c(1, jdx+1)]
    end <- idx[c(jdx, length(idx))]
    dframe[c(start, end),]
  })
  
  outline <- na.omit(unique(rbind(edgesX, edgesY)))
  outline[,c("x", "y", "red", "green", "blue")]
}

#' Determine order in which to pass through a set of points
#' 
#' Greedy algorithm to connect points, with the idea that a point is connected by a line with the two points closest to each, that haven't yet been connected into the shape.
#' Results depend on the starting point.
#' @param x x coordinates
#' @param y y coordinates
#' @return vector of indices for ordered traversion along border
#' @export
determineOrder <- function (x, y) {
  determineNext <- function(now, left) {
    x1 <- x[now]
    y1 <- y[now]
    dists <- (x1-x[left])^2 + (y1-y[left])^2
    left[which.min(dists)]
  }
  
  order <- 1
  leftover <- c(1:length(x))[-order]
  now <- order
  while (length(leftover) > 0) {
    now <- determineNext(now, leftover)
    order <- c(order, now)
    leftover <- leftover[-which(leftover==now)]
  }

  order(order)
}

#' Identify different parts of a polygon
#' 
#' @param data is a data frame with coordinates x, y, and order.
#' @param tol numerical tolerance for minimal distance between groups. If this value is not specified, a tolerace is derived from the marginal frequencey break down of observed (squared) distances between consecutive points.
#' @return data frame group variable is added to the input data
#' @export
identifyParts <- function(data, tol = NULL) {
  ## define implicit bindings for variables - not necessary though, 
  ## since all of the variables do exist at this point
  group <- NA
  order <- NA
  data <- data[order(data$order),]
  data$d <- 0
  data$d[-1] <- diff(data$x)^2 + diff(data$y)^2
  
  # find different parts:
  #  idx <- which(data$d > quantile(data$d, probs=0.9))
  if (is.null(tol)) {
    dt <- as.numeric(names(table(data$d)))
    tol <- dt[min(which(diff(dt)>2))]
  }
  idx <- which(data$d > tol)
  data$group <- rep(1:(length(idx)+1),  diff(c(1, idx, nrow(data)+1)))
  dg <- table(data$group)
  idx <- which(dg < 2)
  if (length(idx) > 0) {
    data <- subset(data, !(group %in% idx))
  }
  data
}

determineDirection <- function(x,y) {
  #  positive direction indicates clockwise, negative counter-clockwise order
  sum(diff(x)*(y[-1]+y[-length(y)]))
}

setDirection <- function(polygon, setdir=1) {
  getdir <- determineDirection(polygon$x, polygon$y)
  if (sign(getdir) != setdir) {
    polygon$order <- rev(polygon$order)
  }
  polygon
}

insertIsland <- function(island, main) {
  # main is ordered 1:m, island is ordered 1:n
  # the island can be inserted at any point in the polygon, but we need to make sure the ordering is fixed. 
  # let's take the last spot:
  res <- rbind(main, island, main[nrow(main),])
  res$order <- 1:nrow(res)
  res
}

completePolygon <- function(polygon) {
  polygon <- polygon[order(polygon$order), ]
  polygon <- rbind(polygon, polygon[1,])
  polygon$order <- 1:nrow(polygon)  
  polygon
}


#' Set the orientation of a polygon
#' 
#' @param imdf dataframe describing a pixellated image in x and y. Has to have columns x, y, and group
#' @importFrom plyr llply
#' @return reordered version of data frame imdf consistent with an assumption of group 1 being the main outline and any other groups being cutouts
#' @export
mainPlusIslands <- function(imdf) {
  group <- NA
  
  # assume that first part is the main with additional islands
  main <- completePolygon(unique(subset(imdf, group==1)))
  main <- setDirection(main, 1)
  main <- main[order(main$order),]
  
  lpath2 <- main
  # now make islands of the small groups
  if (max(imdf$group) > 1) {
    islands <- llply(2:max(imdf$group), function(i) {
      l2 <- completePolygon(unique(subset(imdf, group==i)))
      l2 <- setDirection(l2, -1)
      l2[order(l2$order),]    
    })
    
    for (i in 1:length(islands)) {
      lpath2 <- insertIsland(islands[[i]], main=lpath2)
    }
  }
  lpath2
}

#' Douglas-Peucker algorithm adjusted fo polyons
#' 
#' Implementation of the Douglas-Peucker algorithm for line thinning
#' http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
#' @param points matrix of x and y points
#' @param tol tolerance
#' @export
simplifyPolygon <- function(points, tol=1) {
  # thin out polygon in two steps of Douglas-Pecker:
#  source("thin.r")
  #  browser()
  #  n <- nrow(points)
  #  n1 <- n%/%2
  #  res1 <- simplify_rec(points[1:n1,c("x","y")], tol=1)
  #  res2 <- simplify_rec(points[(n1+1):n,c("x","y")], tol=1)
  points[simplify_rec(points, tol=tol),]
}

#' Convert an image file to a polygon
#' 
#' @param ch letter 
#' @param fontfamily R has a few default fonts that are always available, such as e.g. Helvetica, Arial, Courier New, and Garamond. Other fonts might be available depending on the platform used. 
#' @param fontsize by default 576. If the resulting string exceeds the boundary of the matrix returned, reduced font size
#' @param tol tolerance
#' @param dim vector of length two specifying width and height (in pixels) of the temporary jpg file created for the letter. Defaults to 480 x 480 pixels.
#' @param threshold numerical cutoff between 0 and 1
#' @param var one of "red", "green", "blue".
#' @import ggplot2
#' @export
#' @examples
#' library(ggplot2)
#' letter <- letterToPolygon("R", fontfamily="Helvetica")
#' qplot(x, y, geom="polygon", data = letter, fill=I("black"), alpha=I(0.8))+
#'      coord_equal()
letterToPolygon <- function(ch, fontfamily="Helvetica", fontsize=576, tol=1, dim=c(480, 480), threshold=0.5, var="red") {  
  im <- letterObject(ch, fontfamily=fontfamily, fontsize=fontsize, dim=dim)
  imdf <- fortify(im)
  outline <- getOutline(imdf, threshold=threshold, var=var)
    
  outline$order <- determineOrder(outline$x, outline$y)
  letterpath <- identifyParts(outline, tol=5) # puts group into letterpath
  # thin polygons by part
  ## define implicit bindings for variables - not necessary though, 
  ## since all of the variables do exist at this point
  group <- NA
  letterpath2 <- ddply(letterpath, .(group),  simplifyPolygon, tol=tol)
  lpath2 <- mainPlusIslands(letterpath2)

  lpath2
}
