#' Smoothing time series
#'
#' In order to eliminate outliers in the time series, Hamunyela et al. (2013)
#' smoothing is used (see References).
#'
#' @section References:
#' Hamunyela, E., Verbesselt, J., Roerink, G., & Herold, M. (2013).
#' Trends in spring phenology of western European deciduous forests.
#' Remote Sensing,5(12), 6159-6179.
#'
#' Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
#' deforestation detection through using photosynthetic vegetation time
#' series – (PVts-β). Ecological Indicators, 94, 367 379.
#' @param x Numeric, matrix
#' @param Interp Four interpolation methods are presented, "na.interp",
#' "na.StructTS", "na.approx" and "na.spline". By default is the method
#' "na.interp".
#' @importFrom forecast na.interp
#' @importFrom zoo rollapply
#' @export
#' @examples
#' library(ForesToolboxRS)
#' library(forecast)
#' library(zoo)
#'
#' x <- c(80,78,75,76,79,-100,82,76,81,77,76)
#' smth <- smootH(x)
#' plot(x, type="o", ylab="Reflectance %", xlab="Time")
#' lines(smth, col="blue", type="o")
#'
smootH <- function(x, Interp="na.interp"){

  if (is.vector(x)) {
    x[x==0 | x== -1] <- NA
    x[summary(x)[7] >= (length(x)-1)] <- 100

    # Type of interpolation
    if (Interp=="na.interp") {
      x <- na.interp(x)

    } else if (Interp=="na.StructTS") {
      x <- na.StructTS(x)

    } else if (Interp=="na.approx") {
      x <- na.approx(x)

    } else if (Interp=="na.spline") {
      x <- na.spline(x)

    } else stop("Unsupported interpolation method")

    # We apply Hamunyela Smoothing
    for (j in 2:(length(x)-1)) {
      x[j] <- ifelse(((x[j]-x[j-1]) < -0.01*x[j-1]) & ((x[j]-x[j+1]) < -0.01*x[j+1]),
                      (x[j-1]+x[j+1])/2, x[j])
    }
    np <- x
  } else if (is(x, 'matrix')) {

    for (i in 1:dim(x)[1]) {
      x[i,][x[i,] == 0 | x[i,] == -1]<- NA
      x[i,][summary(x[i,])[7] >= (dim(x)[2]-1)] <- 100

    # Type of interpolation
      if (Interp=="na.interp") {
        x[i,] <- na.interp(x[i,])

      } else if (Interp=="na.StructTS") {
        x[i,] <- na.StructTS(x[i,])

      } else if (Interp=="na.approx") {
        x[i,] <- na.approx(x[i,])

      } else if (Interp=="na.spline") {
        x[i,] <- na.spline(x[i,])

      } else stop("Unsupported interpolation method")
    }
    # We apply Hamunyela Smoothing
    for(i in 1:dim(x)[1]){
      for(j in 2:(dim(x)[2]-1)){
        x[i,][j]<-ifelse(((x[i,][j]-x[i,][j-1])< -0.01*x[i,][j-1]) & ((x[i,][j]-x[i,][j+1])< -0.01*x[i,][j+1]),
                          (x[i,][j-1]+x[i,][j+1])/2,x[i,][j])
      }
    }
    np <- x
  } else if (is(x,'RasterStack') | is(x,'RasterBrick')){
    np <- as.matrix(x)
    for (i in 1:dim(np)[1]) {
      np[i,][np[i,] == 0 | np[i,] == -1]<- NA
      np[i,][summary(np[i,])[7] >= (dim(np)[2]-1)] <- 100

      # Type of interpolation
      if (Interpolation=="na.interp") {
        np[i,] <- na.interp(np[i,])

      } else if (Interpolation=="na.StructTS") {
        np[i,] <- na.StructTS(np[i,])

      } else if (Interpolation=="na.approx") {
        np[i,] <- na.approx(np[i,])

      } else if (Interpolation=="na.spline") {
        np[i,] <- na.spline(np[i,])

      } else stop("Unsupported interpolation method")
    }

    # We apply Hamunyela Smoothing
    for(i in 1:dim(np)[1]){
      for(j in 2:(dim(np)[2]-1)){
        np[i,][j]<-ifelse(((np[i,][j]-np[i,][j-1])< -0.01*np[i,][j-1]) & ((np[i,][j]-np[i,][j+1])< -0.01*np[i,][j+1]),
                          (np[i,][j-1]+np[i,][j+1])/2,np[i,][j])
      }
    }
  }
  else {
    stop(class(x), ' class is not supported')
  }
  return(np)
}
