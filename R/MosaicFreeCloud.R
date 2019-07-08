#' Get the best Available Pixel on an Image Collection
#'
#' A three-step algorithm for creating a mosaic from satellite imagery (see notes).
#'
#' @section Note:
#' The three steps are: \cr
#' \itemize{
#'    \item Chooses the image less affected by atmospheric noise.
#'    \item Calculate the time gap between the image selected in step one to the remaining images.
#'    \item Take the clean pixels  considering the shortest time gap and ending with the image most temporally distant.
#' }
#' @section References:
#' Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
#' deforestation detection through using photosynthetic vegetation time
#' series – (PVts-β). Ecological Indicators, 94, 367 379.
#'
#' @param x filename(character), Raster*  or star object (see \link[stars]{read_stars}).
#' @param time Date or character vector with the dates of the files.
#' @param CLUSTER cluster to use for parallel apply; see \link[parallel]{makeCluster}.
#' @param fit_negative logical; if \code{TRUE}, negative values are replace by NA.
#' @param RasterLayer logical; if \code{TRUE}, return a RasterLayer object.
#' @importFrom methods as
#' @importFrom stats na.omit
#' @importFrom raster merge
#' @importFrom stars st_apply
#' @export
#' @examples
#'library(raster)
#'library(stars)
#'## Creat a mosaic using three landsat images (just one band)
#'
#'img <- "external/data/landsat8_stack.tif"
#'time <- as.Date(c("2016-07-30","2016-08-15","2016-09-16"))
#'mosaic <- MosaicFreeCloud(img, time)
#'plot(mosaic)
#'
MosaicFreeCloud <- function(img, time, CLUSTER = NULL, fit_negative = TRUE,
                                 RasterLayer = TRUE) {

  img <- pvts_whatkinditis(img)

  if (length(dim(img)) == 2) {
    img = img %>% merge %>% st_set_dimensions(names = c("x", "y", "band"))
  }

  # Selecting the best image (1)
  select_img <- img %>%
    st_apply(MARGIN = 3,
             FUN = function(img) sum(img, na.rm = T),
             CLUSTER = CLUSTER) %>%
    unlist %>% order(decreasing = T) %>%
    `[`(1)

  # Estimating the time distance with respect to the best image (2)
  img_ord <- order(abs(time - time[select_img]), decreasing = F)
  img <- img %>% split('band') %>% '['(img_ord)

  # Change negative value by NA (3)
  if (fit_negative) {
    img[img < 0] <- NA
  }

  # Aggregating considering (1) and (2)
  one_img <- img %>% merge %>% st_apply(MARGIN = c("x", "y"),
                                      FUN = function(img) na.omit(img)[1],
                                      CLUSTER = CLUSTER)
  if (RasterLayer) {
    return(as(one_img, "Raster"))
  } else {
    return(one_img)
  }
}

pvts_whatkinditis <-function(img) {
  if (is.character(img)) {
    img = read_stars(img)
  } else if(is(img,'stars')) {
    img
  } else if (is(img,'RasterLayer') | is(img,'RasterStack') | is(img,'RasterBrick')) {
    img = st_as_stars(img)
  } else {
    stop(class(img), ' class is not supported')
  }
  return(img)
}
