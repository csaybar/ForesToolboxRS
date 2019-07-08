#' Change detection using the PVts-\eqn{\beta} approach
#'
#' This algorithm will allow to detect disturbances in the forests using
#' all the available Landsat set. In fact, it can also be run with sensors
#' such as MODIS.
#'
#' @section References:
#' Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
#' deforestation detection through using photosynthetic vegetation time
#' series – (PVts-β). Ecological Indicators, 94, 367 379.
#' @section Note:
#' In order to optimise the detections, it is advisable to make a smoothing
#' through the \link[ForesToolboxRS]{smootH} function before detecting changes. The smoothing will
#' allow to eliminate outliers that were not eliminated during the masking
#' of atmospheric artifacts.
#'
#' @param x Matrix, RasterStack, RasterBrick
#' @param start The start of the monitoring time
#' @param end The end of the monitoring time
#' @param threslhold The default threshold is 5 for photosynthetic vegetation,
#' while for indices such as NDVI and EVI the threshold is 3.
#' Please see Tarazona et al. (2018) for more details.
#' @param img The last monitoring image (in case "x" is a matrix)
#' @export
#' @examples
#' library(ForesToolboxRS)
#' library(raster)
#'
#'## Not run:
#'imgs <- as.matrix(stack(list)) # imgs is a matrix
#'dim(imgs)[1]= 1000 is the number of pixels (time series)
#'dim(imgs)[2]=29 is the number of images between 1990 and 2018
#'Change monitoring period 2015-2018. 2015(position 26),2018(position 29)
#'
#'cd <- pvts(x=imgs, start=26, end=29, threshold=5,img=lastimg)
#'plot(cd)
#'
pvts<- function(x,start,end,threshold=5,img) {

  if (is(x, 'matrix')) {

    breakR<-img
    mean<-apply(x[,1:start], 1, mean)
    std<-apply(x[,1:start], 1, sd)
    cd<- ifelse(x[,end] < (mean-threshold*std), 1, 0)
    values(breakR) <- cd; breakR[img<80 | img<0.8]<- 0

  } else if (is(x,'RasterStack') | is(x,'RasterBrick')) {

    breakR <- x[[end]]; x <- as.matrix(x)
    mean<-apply(x[,1:start], 1, mean)
    std<-apply(x[,1:start], 1, sd)
    cd<- ifelse(x[,end] < (mean-threshold*std), 1, 0)
    values(breakR) <- cd; breakR[img<80 | img<0.8]<- 0

  }
  return(breakR)
}
