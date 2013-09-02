#' Convert a text element into an R object
#' 
#' @param ch text to be converted, usually just a single letter
#' @param fontfamily R has a few default fonts that are always available, such as e.g. Helvetica, Arial, Courier New, and Garamond. Other fonts might be available depending on the platform used. 
#' @param fontsize by default 576. If the resulting string exceeds the boundary of the matrix returned, reduced font size
#' @param dim vector of length two specifying width and height (in pixels) of the temporary jpg file created for the letter. Defaults to 480 x 480 pixels.
#' @return three dimensional matrix of dimension 480 x 480 x 3 of the pixel values, black background and white letter 
#' @export
#' @examples
#' plot(letterObject("g", fontfamily="Garamond", fontsize=400))
#' plot(letterObject("q", fontsize=400))
#' plot(letterObject("B"))
letterObject <- function(ch, fontfamily="Helvetica", fontsize=576, dim=c(480, 480)) {
  require(ReadImages)
  require(grid)
  fname <- tempfile(pattern = "file", fileext=".jpg")
  jpeg(filename=fname, width=dim[1], height=dim[2])
  grid.newpage()
  grid.rect(x = 0, y=0, width=3, height=3,
            gp=gpar(fill="black"), draw = TRUE, vp = NULL)
  
  grid.text(ch, 0.5,0.5, gp=gpar(fontsize=fontsize, fontfamily=fontfamily, col="white"))
  dev.off()
  read.jpeg(fname)
}

scaleTo <- function(x, fromRange=range(x), toRange=c(0,1)) {
  x <- as.numeric(as.character(x))
  (x-fromRange[1])/diff(fromRange)*diff(toRange) + toRange[1]
}

imageToDFrame <- function(letter) {
  dims <- dim(letter)
  imdf <- adply(letter, .margins=1, function(x) x)
  imdf$x <- rep(1:dims[2], length=nrow(imdf)) 
  names(imdf) <- c("y", "red", "green", "blue", "x")
  imdf$y <- -as.numeric(as.character(imdf$y))
  imdf[,c("x", "y", "red", "green", "blue")]
}

getOutline <- function(imdf, var, threshold) {
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
  
  data.frame(x=x[order], y=y[order], order=1:length(order))
}

identifyParts <- function(letterpath, tol = NULL) {
  stopifnot(c("group") %in% names(letterpath))
  ## define implicit bindings for variables - not necessary though, 
  ## since all of the variables do exist at this point
  group <- NA
  
  letterpath$d <- 0
  letterpath$d[-1] <- diff(letterpath$x)^2 + diff(letterpath$y)^2
  
  # find different parts:
  #  idx <- which(letterpath$d > quantile(letterpath$d, probs=0.9))
  if (is.null(tol)) {
    dt <- as.numeric(names(table(letterpath$d)))
    tol <- dt[min(which(diff(dt)>2))]
  }
  idx <- which(letterpath$d > tol)
  letterpath$group <- rep(1:(length(idx)+1),  diff(c(1, idx, nrow(letterpath)+1)))
  dg <- table(letterpath$group)
  idx <- which(dg < 2)
  if (length(idx) > 0) {
    letterpath <- subset(letterpath, !(group %in% idx))
  }
  letterpath
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

mainPlusIslands <- function(letterpath) {
  # assume that first part is the main with additional islands
  main <- completePolygon(unique(subset(letterpath, group==1)))
  main <- setDirection(main, 1)
  main <- main[order(main$order),]
  
  lpath2 <- main
  # now make islands of the small groups
  if (max(letterpath$group) > 1) {
    islands <- llply(2:max(letterpath$group), function(i) {
      l2 <- completePolygon(unique(subset(letterpath, group==i)))
      l2 <- setDirection(l2, -1)
      l2[order(l2$order),]    
    })
    
    for (i in 1:length(islands)) {
      lpath2 <- insertIsland(islands[[i]], main=lpath2)
    }
  }
  lpath2
}


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
#' @param fontfamily
#' @param fontsize
#' @param tol tolerance
#' @param dim vector of length two specifying width and height (in pixels) of the temporary jpg file created for the letter. Defaults to 480 x 480 pixels.
#' @param threshold
#' @param variable one of "red", "green", "blue".
#' @export
#' @examples
#' letter <- letterToPolygon("ö", fontfamily="Garamond")
#' print(qplot(x, y, geom="polygon", data = letter, fill=I("black"), order=order, alpha=I(0.8))+coord_equal())
letterToPolygon <- function(ch, fontfamily="Helvetica", fontsize=576, tol=1, dim=c(480, 480), threshold=0.5, var="red") {  
  im <- letterObject(ch, fontfamily=fontfamily, fontsize=fontsize, dim=dim)
  imdf <- imageToDFrame(im)
  outline <- getOutline(imdf, threshold=threshold, var=var)

  
    
  letterpath <- determineOrder(outline$x, outline$y)
  letterpath <- identifyParts(letterpath, tol=5) # puts group into letterpath
  # thin polygons by part
  ## define implicit bindings for variables - not necessary though, 
  ## since all of the variables do exist at this point
  group <- NA
  letterpath2 <- ddply(letterpath, .(group),  simplifyPolygon, tol=tol)
  lpath2 <- mainPlusIslands(letterpath2)
#  alphabet <-rbind(alphabet,  data.frame(lpath2, group=ch, region=ch))
#  print(qplot(x, y, geom="polygon", data = lpath2, group=1, fill=I("black"), order=order, alpha=I(0.8))+coord_equal())
  #scan()
  lpath2
}