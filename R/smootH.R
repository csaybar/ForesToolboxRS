#' Smoothing time series
#'
#' In order to eliminate outliers in the time series, Hamunyela et al. (2013)
#' smoothing is used (see References).
#'
#' @name smootH-ftb
#' @section References:
#' Hamunyela, E., Verbesselt, J., Roerink, G., & Herold, M. (2013).
#' Trends in spring phenology of western European deciduous forests.
#' Remote Sensing,5(12), 6159-6179.
#'
#' Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
#' deforestation detection through using photosynthetic vegetation time
#' series – (PVts). Ecological Indicators, 94, 367 379.
#' @param x object of class numeric, matrix, RasterStack, RasterBrick or stars
#' @param CLUSTER cluster to use for parallel apply; see \link[parallel]{makeCluster}
#' @param interp Four interpolation methods are presented, "na.interp",
#' "na.StructTS", "na.approx" and "na.spline". By default is the method
#' "na.interp".
#' @param ... ignored
#' @importFrom stars st_set_dimensions st_apply st_as_stars
#' @importFrom raster stack brick
#' @importFrom parallel parApply
#' @importFrom methods as
#' @importFrom forecast na.interp
#' @importFrom zoo na.spline na.approx
#' @examples
#' library(ForesToolboxRS)
#' library(stars)
#'
#' # Example 1
#' x <- c(80,78,75,76,79,-100,82,76,81,77,76)
#' smth <- smootH(x)
#' plot(x, type="o", ylab="Reflectance %", xlab="Time")
#' lines(smth, col="blue", type="o")
#'
#' # Example 2
#' rasterio <- list(nXOff = 250, nYOff = 250)
#' pv_serie <- system.file("PVts", package="ForesToolboxRS") %>%
#'   list.files("\\.tif$",full.names = TRUE) %>%
#'   read_stars(RasterIO = rasterio)
#' pv_serie_stars <- smootH(pv_serie)
#' smoth_mean <- st_apply(pv_serie_stars %>% merge,1:2,mean)
#' simple_mean <- st_apply(pv_serie %>% merge,1:2,function(x) mean(x,na.rm=TRUE))
#' plot(smoth_mean,main = "SmootH mean")
#' plot(simple_mean,main = "Simple mean")
#' @export
smootH <- function(x, interp="na.interp", CLUSTER=NULL, ...) {
  UseMethod("smootH")
}

#' @name smootH-ftb
#' @export
smootH.numeric  <- function(x, interp="na.interp", ...) {
  na_interpolation(x, interp = interp) %>%
    rcpp_hamunyela()
}

#' @name smootH-ftb
#' @export
smootH.matrix  <- function(x, interp="na.interp", CLUSTER = NULL, ...) {
  # We apply Hamunyela Smoothing
  if (is.null(CLUSTER)) {
    x_na <- apply(x, 1, na_interpolation) %>% t
    apply(x_na, 1, rcpp_hamunyela) %>% t
  } else {
    x_na <- parApply(CLUSTER, x, MARGIN = 1, na_interpolation) %>% t
    parApply(CLUSTER, x_na, MARGIN = 1, rcpp_hamunyela) %>% t
  }
}

#' @name smootH-ftb
#' @export
smootH.RasterStack  <- function(x, interp="na.interp", CLUSTER = NULL, ...) {
  x <- ftb_whatkinditis(x)
  if (length(dim(x)) == 2) {
    x = x %>% merge %>% st_set_dimensions(names = c("x", "y", "band"))
  }
  x %>%
    st_apply(MARGIN = c('x','y'),
             FUN = na_interpolation,
             CLUSTER =  CLUSTER,
             interp = interp) %>%
    st_apply(MARGIN = c('x','y'),
             FUN = rcpp_hamunyela,
             CLUSTER =  CLUSTER) %>%
    st_set_dimensions(names = c('fun',"x","y")) %>%
    split('fun') -> x
  stack(mapply(function(z) as(x[z],'Raster'),seq_len(length(x))))
}


#' @name smootH-ftb
#' @export
smootH.RasterBrick  <- function(x, interp="na.interp", CLUSTER = NULL, ...) {
  x <- ftb_whatkinditis(x)
  if (length(dim(x)) == 2) {
    x = x %>% merge %>% st_set_dimensions(names = c("x", "y", "band"))
  }
  x %>%
    st_apply(MARGIN = c('x','y'),
             FUN = na_interpolation,
             CLUSTER =  CLUSTER,
             interp = interp) %>%
    st_apply(MARGIN = c('x','y'),
             FUN = rcpp_hamunyela,
             CLUSTER =  CLUSTER) %>%
    st_set_dimensions(names = c('fun',"x","y")) %>%
    split('fun') -> x
  brick(mapply(function(z) as(x[z],'Raster'),seq_len(length(x))))
}


#' @name smootH-ftb
#' @export
smootH.stars  <- function(x, interp="na.interp", CLUSTER = NULL, ...) {
  x <- ftb_whatkinditis(x)

  if (length(dim(x)) == 2) {
    x = x %>% merge %>% st_set_dimensions(names = c("x", "y", "band"))
  }

  x %>%
    st_apply(MARGIN = c('x','y'),
             FUN = na_interpolation,
             CLUSTER =  CLUSTER,
             interp = interp) %>%
    st_apply(MARGIN = c('x','y'),
             FUN = rcpp_hamunyela,
             CLUSTER =  CLUSTER) %>%
    st_set_dimensions(names = c('fun',"x","y")) %>%
    split('fun')
}


#' Choose your preference interpolation  method - numeric | integer
#' @noRd
na_interpolation <- function(x, interp = 'na.interp') {
  x[x %in% c(0, -1)] <- NA
  len_x <- length(x) - 1
  x_na_count <- sum(is.na(x)) # count the NA values in numeric vector
  x[x_na_count >= len_x] <- 100

  # Type of interpolation
  if (interp=="na.interp") {
    x <- na.interp(x)

    #} else if (interp=="na.StructTS") {
    #    x <- na.StructTS(x)

  } else if (interp=="na.approx") {
    x <- na.approx(x, na.rm = FALSE)

  } else if (interp=="na.spline") {
    x <- na.spline(x)

  } else stop("Unsupported interpolation method")
  return(x)
}
