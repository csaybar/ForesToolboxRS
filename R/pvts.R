#' Change detection using the PVts-\eqn{\beta} approach
#'
#' This algorithm will allow to detect disturbances in the forests using
#' all the available Landsat set. In fact, it can also be run with sensors
#' such as MODIS.
#'
#' @section References:
#' Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
#' deforestation detection through using photosynthetic vegetation time
#' series – (PVts-b). Ecological Indicators, 94, 367 379.
#' @section Note:
#' In order to optimise the detections, it is advisable to make a smoothing
#' through the \link[ForesToolboxRS]{smootH} function before detecting changes. The smoothing will
#' allow to eliminate outliers that were not eliminated during the masking
#' of atmospheric artifacts.
#'
#' @param x Vector, Matrix, RasterStack, RasterBrick
#' @param startm The start of the monitoring time
#' @param endm The end of the monitoring time (The year in which changes will be detected)
#' @param threshold The default threshold is 5 for photosynthetic vegetation,
#' while for indices such as NDVI and EVI the threshold is 3.
#' Please see Tarazona et al. (2018) for more details.
#' @param img The image of the monitoring start position, i.e. "start" position (in case "x" is a matrix)
#' @param tr The vector of the analysis time range must contain the start time of the
#' time series, the end time and the frequency of the series. For example:
#' tr <- c(1990, 2017, 1) (i.e., the time series starts in 1990, ends in 2017 and
#' has an annual frequency of 1). See \link[stats]{ts} for more details.
#' @param time If it is TRUE the plotting will be with time coordinates.
#' @importFrom stats sd time ts
#' @importFrom graphics points abline polygon text grid legend plot
#' @importFrom grDevices adjustcolor
#' @importFrom raster values<-
#' @examples
#' library(ForesToolboxRS)
#' library(raster)
#' img_list <- system.file("PVts", package="ForesToolboxRS") %>%
#'   list.files("\\.tif$",full.names = TRUE)
#' # Example 1.
#' # photosynthetic vegetation time series between 1990 and 2017.
#' # We will detect changes in 2008 (position 19)
#' vec <- c(0.86, 0.93, 0.97, 0.91, 0.95, 0.96,
#'          0.91, 0.88, 0.92, 0.89, 0.90, 0.89,
#'          0.91, 0.92, 0.89, 0.90, 0.92, 0.84,
#'          0.46, 0.20, 0.27, 0.22, 0.52, 0.63,
#'          0.61, 0.67, 0.64, 0.86)
#' cd <- pvts(x = vec, start = 18, end = 19, threshold = 5)
#' \dontrun{
#' #Example 2.
#' #imgs <- as.matrix(stack(img_list)) # imgs is a matrix
#' #dim(imgs)[1]= 1000 # is the number of pixels (time series)
#' #dim(imgs)[2]=29 # is the number of images between 1990 and 2018
#' ## Change monitoring period 2015-2018. 2015(position 26),2018(position 29)
#' #cd <- pvts(x=imgs, start=26, end=29, threshold=5,img=lastimg)
#' #plot(cd)
#' }
#' @export
pvts <- function(x, startm, endm, threshold = 5, img, tr, time = FALSE) {
  if (is(x, "vector")) {
    mean <- mean(x[1:startm]) # mean
    std <- sd(x[1:startm]) # standard deviation
    li <- mean - threshold * std # lower limit

    if (x[endm] < li) {
      main <- "Breakpoint detected"
      breakR <- plot.pvts(x, endm,tr,time,main,li)
    } else {
      stop("Breakpoint not detected")
    }
  } else if (is(x, "matrix")) {
    breakR <- img
    mean <- apply(x[, 1:startm], 1, mean)
    std <- apply(x[, 1:startm], 1, sd)
    cd <- ifelse(x[, endm] < (mean - threshold * std), 1, 0)
    values(breakR) <- cd
    breakR[img < 80 | img < 0.8] <- 0
  } else if (is(x, "RasterStack") | is(x, "RasterBrick")) {
    breakR <- x[[endm]]
    x <- as.matrix(x)
    mean <- apply(x[, 1:startm], 1, mean)
    std <- apply(x[, 1:startm], 1, sd)
    cd <- ifelse(x[, endm] < (mean - threshold * std), 1, 0)
    values(breakR) <- cd
    breakR[img < 80 | img < 0.8] <- 0
  } else {
    stop(class(x), " class is not supported")
  }

  return(breakR)
}


#' plot pvts
#' @noRd
plot.pvts <- function(x, endm,tr,time,main,li) {
  p <- endm
  n <- length(x) / 2
  m <- length(x) / 20
  # Plot in time coordinates
  if (time) {
    x <- ts(x, start = c(tr[1]), end = c(tr[2]), frequency = tr[3])
    p <- time(x)[endm]
    n <- time(x)[length(x) / 2]
    m <- time(x)[length(x) / 20]
  }

  plot(x, ylab = "Variable", ylim = c(0, 1.1), type = "l", lwd = 1.5, main = main)
  points(x, pch = 20, lwd = 1.5, cex = 1.5)
  abline(h = li, col = "red", lty = 2, lwd = 2)
  abline(v = c(p - 1 / 2, p + 1 / 2), col = "blue", lty = 3, lwd = 2)
  x1 <- c(p - 1 / 2, p - 1 / 2, p + 1 / 2, p + 1 / 2)
  y1 <- c(-0.1, 1.5, 1.5, -0.1)
  polygon(x1, y1, col = adjustcolor("slateblue1", alpha.f = 0.4), border = NA)
  text(n, li - 0.05, "Lower limit", col = "red", cex = 1.3)
  grid()
  legend(m, 0.6, c("Variable", "Lower limit", "Breakpoint"),
         inset = .02, cex = 0.9, lty = c(1, 2, 1),
         lwd = c(2, 2, 5), col = c("black", "red", "slateblue1"), bty = "n"
  )
}
