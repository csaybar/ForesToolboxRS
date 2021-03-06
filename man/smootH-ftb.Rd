% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smootH.R
\name{smootH-ftb}
\alias{smootH-ftb}
\alias{smootH}
\alias{smootH.numeric}
\alias{smootH.matrix}
\alias{smootH.RasterStack}
\alias{smootH.RasterBrick}
\alias{smootH.stars}
\title{Smoothing time series}
\usage{
smootH(x, interp = "na.interp", CLUSTER = NULL, ...)

\method{smootH}{numeric}(x, interp = "na.interp", ...)

\method{smootH}{matrix}(x, interp = "na.interp", CLUSTER = NULL, ...)

\method{smootH}{RasterStack}(x, interp = "na.interp", CLUSTER = NULL,
  ...)

\method{smootH}{RasterBrick}(x, interp = "na.interp", CLUSTER = NULL,
  ...)

\method{smootH}{stars}(x, interp = "na.interp", CLUSTER = NULL, ...)
}
\arguments{
\item{x}{object of class numeric, matrix, RasterStack, RasterBrick or stars}

\item{interp}{Four interpolation methods are presented, "na.interp",
"na.StructTS", "na.approx" and "na.spline". By default is the method
"na.interp".}

\item{CLUSTER}{cluster to use for parallel apply; see \link[parallel]{makeCluster}}

\item{...}{ignored}
}
\description{
In order to eliminate outliers in the time series, Hamunyela et al. (2013)
smoothing is used (see References).
}
\section{References}{

Hamunyela, E., Verbesselt, J., Roerink, G., & Herold, M. (2013).
Trends in spring phenology of western European deciduous forests.
Remote Sensing,5(12), 6159-6179.

Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
deforestation detection through using photosynthetic vegetation time
series (PVts). Ecological Indicators, 94, 367 379.
}

\examples{
library(ForesToolboxRS)
library(stars)

# Example 1
x <- c(80,78,75,76,79,-100,82,76,81,77,76)
smth <- smootH(x)
plot(x, type="o", ylab="Reflectance \%", xlab="Time")
lines(smth, col="blue", type="o")

# Example 2
rasterio <- list(nXOff = 250, nYOff = 250)
pv_serie <- system.file("PVts", package="ForesToolboxRS") \%>\%
  list.files("\\\\.tif$",full.names = TRUE) \%>\%
  read_stars(RasterIO = rasterio)
pv_serie_stars <- smootH(pv_serie)
smoth_mean <- st_apply(pv_serie_stars \%>\% merge,1:2,mean)
simple_mean <- st_apply(pv_serie \%>\% merge,1:2,function(x) mean(x,na.rm=TRUE))
plot(smoth_mean,main = "SmootH mean")
plot(simple_mean,main = "Simple mean")
}
