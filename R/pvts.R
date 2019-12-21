#' Change detection using the PVts-\eqn{\beta} approach
#'
#' @author Yonatan Tarazona
#'
#' This algorithm will allow to detect disturbances in the forests using
#' all the available Landsat set. In fact, it can also be run with sensors
#' such as MODIS.
#'
#' @section References:
#' Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
#' deforestation detection through using photosynthetic vegetation time
#' series (PVts-\eqn{\beta}). Ecological Indicators, 94, 367 379.
#' @section Note:
#' In order to optimise the detections, it is advisable to make a smoothing
#' through the \link[ForesToolboxRS]{smootH} function before detecting changes. The smoothing will
#' allow to eliminate outliers that were not eliminated during the masking
#' of atmospheric artifacts.
#' @importFrom stats sd time ts
#' @importFrom graphics points abline polygon text grid legend plot
#' @importFrom grDevices adjustcolor
#' @importFrom raster values<- as.matrix
#' @param x Vector, Matrix, RasterStack, RasterBrick
#' @param startm The start of the monitoring time
#' @param endm The end of the monitoring time
#' @param threshold The default threshold is 5 for photosynthetic vegetation,
#' while for indices such as NDVI and EVI the threshold is 3.
#' Please see Tarazona et al. (2018) for more details.
#' @param img The image of the position immediately before the monitoring start,
#' i.e. the "start-1" position (in case "x" is a matrix).
#' @param tr The vector of the analysis time range must contain the start time of the
#' time series, the end time and the frequency of the series. For example:
#' tr <- c(1990, 2017, 1) (i.e., the time series starts in 1990, ends in 2017 and
#' has an annual frequency of 1). See \link[stats]{ts} for more details.
#' @param time If it is TRUE the plotting will be with time coordinates
#' (only if the "tr" parameter is within the function).
#' @param vf If the monitoring is with Photosynthetic Vegetation series,
#' then switch to TRUE.
#' @export
#' @examples
#' library(ForesToolboxRS)
#' library(forecast)
#' library(raster)
#'
#' # Example 1.
#' vec <- c(0.86,0.93,0.97,0.91,0.95,0.96,0.91,0.88,0.92,0.89,0.90,0.89,0.91,0.92,0.89,
#'          0.90,0.92,0.84,0.46,0.20,0.27,0.22,0.52,0.63,0.61,0.67,0.64,0.86) # photosynthetic vegetation
#' # time series between 1990 and 2017
#' # We will detect changes in 2008 (position 19)
#' cd <- pvts(x=vec, startm=19, endm=19, threshold= 5)
#'
#' # Example 2.
#' imgs <- system.file("PVts",package = "ForesToolboxRS") %>%
#'   list.files(full.names = TRUE) %>%
#'   stack()
#'
#' # Change monitoring period 2006-2016. Where 2006 is position 7 and 2016 is position 16
#' cd <- pvts(x=imgs, startm=7, endm=16, threshold=3) # EVI indices
#' plot(cd)
#'
pvts<- function(x, startm, endm, threshold=5, img, tr, time=FALSE, vf=FALSE) {

  if (is(x, "vector")) {
    x[x==0 | x== -1] <- NA
    x[summary(x)[7] >= (length(x)-1)] <- 100
    x <- na.interp(x)
    mean<-mean(x[1:(startm-1)]) # mean
    std<-sd(x[1:(startm-1)]) # standard deviation
    li <- mean-threshold*std # lower limit

    if (x[endm] < li) {
      main <- "Breakpoint detected"

      plot.pvts <- function (x , endm) {

        p <- endm
        n <- length(x)/2
        m <- length(x)/20
        # Plot in time coordinates
        if (time) {
          x <- ts(x, start = c(tr[1]), end = c(tr[2]), frequency = tr[3])
          p <- time(x)[endm]
          n <- time(x)[length(x)/2]
          m <- time(x)[length(x)/20]
        }

        plot(x, ylab="Variable", ylim=c(0,1.1), type="l", lwd=1.5, main=main); points(x, pch=20, lwd=1.5, cex=1.5)
        abline(h=li, col="red", lty=2, lwd=2)
        abline(v=c(p - 1/2, p + 1/2), col="blue", lty=3, lwd=2)
        x1<-c(p - 1/2, p - 1/2, p + 1/2, p + 1/2)
        y1<-c(-0.1,1.5,1.5,-0.1); polygon(x1,y1,col=adjustcolor("slateblue1",alpha.f=0.4),border=NA)
        text(n, li-0.05, "Lower limit", col="red", cex = 1.3)
        grid(); legend(m, 0.6, c("Variable","Lower limit","Breakpoint"),inset=.02, cex = 0.9,lty=c(1,2,1),
                       lwd=c(2,2,5), col=c("black","red", "slateblue1"),bty="n")
      }
      breakR <- plot.pvts(x, endm)

    } else stop("Breakpoint not detected")

  } else if (is(x, 'matrix')) {
    breakR<-img

    for (i in 1:dim(x)[1]) {
      x[i,][x[i,] == 0 | x[i,] == -1]<- NA
      x[i,][summary(x[i,])[7] >= (dim(x)[2]-1)] <- 100
      x[i,] <- na.interp(x[i,])
    }

    mean<-apply(x[,1:(startm-1)], 1, mean)
    std<-apply(x[,1:(startm-1)], 1, sd)
    cd<- ifelse(x[,endm] < (mean-threshold*std), 1, 0)
    values(breakR) <- cd

    # Photosynthetic vegetation?
    if (vf) {
      breakR[img<80 | img<0.8 | img<8000]<- 0
    }

  } else if (is(x,'RasterStack') | is(x,'RasterBrick')) {
    img <- x[[startm-1]]
    breakR <- img

    x <- as.matrix(x)
    for (i in 1:dim(x)[1]) {
      x[i,][x[i,] == 0 | x[i,] == -1]<- NA
      x[i,][summary(x[i,])[7] >= (dim(x)[2]-1)] <- 100
      x[i,] <- na.interp(x[i,])
    }

    mean<-apply(x[,1:startm-1], 1, mean)
    std<-apply(x[,1:startm-1], 1, sd)
    cd<- ifelse(x[,endm] < (mean-threshold*std), 1, 0)
    values(breakR) <- cd

    # Photosynthetic vegetation?
    if (vf) {
      breakR[img<80 | img<0.8 | img<8000]<- 0
    }

  } else {

    stop(class(x), ' class is not supported')
  }

  return(breakR)
}
