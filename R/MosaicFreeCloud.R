#' Get the best Available Pixel on an Image Collection
#'
#' A three-step algorithm for creating a mosaic from satellite imagery (see notes).
#'
#' @name mosaic-free-cloud
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
#' series – (PVts-b). Ecological Indicators, 94, 367 379.
#'
#' @param img Character vector, RasterStack, RasterBrick or 3D/4D star object.
#' @param time Date or character vector with the dates of the files.
#' @param CLUSTER cluster to use for parallel apply; see \link[parallel]{makeCluster}.
#' @param fit_negative logical; if \code{TRUE}, negative values are replace by NA.
#' @param RasterLayer logical; if \code{TRUE}, return a RasterLayer object.
#' @param ... Passed on to \link[stars]{read_stars} parameters.
#' @importFrom methods as
#' @importFrom stats na.omit
#' @importFrom raster merge
#' @importFrom stars st_apply
#' @export
#' @examples
#' library(ForesToolboxRS)
#' library(raster)
#' library(stars)
#' time <- as.Date(c("2016-07-30","2016-08-15","2016-09-16"))
#' rasterio = list(nXOff = 50, nYOff = 50, bands = c(4,3,2))
#'
#' ## Create a mosaic using three landsat images (just one band)
#' img <- system.file("simple_mosaic", package="ForesToolboxRS") %>%
#'   list.files("\\.tif$",full.names = TRUE)
#' free_img <- MosaicFreeCloud(img,time)
#' plot(free_img)
#'
#' ## Create a mosaic using three landsat images (all the bands)
#' img <- system.file("mosaic", package="ForesToolboxRS") %>%
#'   list.files("\\.tif$",full.names = TRUE)
#' free_img <- MosaicFreeCloud(img,time,RasterIO=rasterio)
#' plotRGB(free_img,3,2,1)
#'
#'# Recommend way to use MosaicFreeCloud
#' \dontrun{
#' img <- system.file("mosaic", package="ForesToolboxRS") %>%
#'   list.files("\\.tif$",full.names = TRUE) %>%
#'   read_stars(img,RasterIO=rasterio) %>%
#'   merge
#' CLUSTER = parallel::makeCluster(spec = 4)
#' mosaic <- MosaicFreeCloud(img, time,CLUSTER = CLUSTER)
#' plotRGB(free_img,3,2,1)
#' }
MosaicFreeCloud <- function(img, time, CLUSTER = NULL, fit_negative = TRUE,
                            RasterLayer = TRUE, ...) {
  UseMethod("MosaicFreeCloud")
}

MosaicFreeCloud.default <- function(img, time, CLUSTER = NULL, fit_negative = TRUE,
                                    RasterLayer = TRUE) {
  stop(class(img), " is not supported")
}

#' @name mosaic-free-cloud
#' @export
MosaicFreeCloud.character <- function(img, time, CLUSTER = NULL, fit_negative = TRUE,
                                      RasterLayer = TRUE, ...) {
  x <- read_stars(img, ...)
  img_names <- attr(x, "dimensions")$band$values
  if (is_nD(img) == "3D") {
    x %>%
      merge() %>%
      st_set_dimensions(names = c("x", "y", "time")) %>%
      mosaic_3D(time, CLUSTER, fit_negative) -> x
    if (RasterLayer) {
      return(as(x, "Raster"))
    } else {
      return(x)
    }
  } else if (is_nD(img) == "4D") {
    x %>%
      merge() %>%
      st_set_dimensions(names = c("x", "y", "band", "time")) %>%
      mosaic_4D(time, CLUSTER, fit_negative) -> x
    if (RasterLayer) {
      x <- brick(mapply(function(z) as(x[z], "Raster"), seq_len(length(x))))
      if (!is.null(img_names)) {
        names(x) <- img_names
      }
    }
    return(x)
  }
}

#' @name mosaic-free-cloud
#' @export
MosaicFreeCloud.RasterStack <- function(img, time, CLUSTER = NULL, fit_negative = TRUE,
                                        RasterLayer = TRUE, ...) {
  img %>%
    st_as_stars() %>%
    st_set_dimensions(names = c("x", "y", "time")) %>%
    mosaic_3D(time, CLUSTER, fit_negative) -> img
  if (RasterLayer) {
    return(as(img, "Raster"))
  } else {
    return(img)
  }
}

#' @name mosaic-free-cloud
#' @export
MosaicFreeCloud.RasterBrick <- function(img, time, CLUSTER = NULL, fit_negative = TRUE,
                                        RasterLayer = TRUE, ...) {
  img %>%
    st_as_stars() %>%
    st_set_dimensions(names = c("x", "y", "time")) %>%
    mosaic_3D(time, CLUSTER, fit_negative) -> img
  if (RasterLayer) {
    return(as(img, "Raster"))
  } else {
    return(img)
  }
}

#' @name mosaic-free-cloud
#' @export
MosaicFreeCloud.stars <- function(img, time, CLUSTER = NULL, fit_negative = TRUE,
                                  RasterLayer = TRUE, ...) {
  x_m <- try(merge(img), silent = TRUE)
  if (class(x_m) == "try-error") {
    img
  } else {
    img <- img %>% merge()
  }

  if (is_nD(img) == "3D") {
    img %>%
      st_set_dimensions(names = c("x", "y", "time")) %>%
      mosaic_3D(time, CLUSTER, fit_negative) -> img

    if (RasterLayer) {
      return(as(img, "Raster"))
    } else {
      return(img)
    }
  } else if (is_nD(img) == "4D") {
    img %>%
      st_set_dimensions(names = c("x", "y", "band", "time")) %>%
      mosaic_4D(time, CLUSTER, fit_negative) -> img
    if (RasterLayer) {
      img <- brick(mapply(function(z) as(img[z], "Raster"), seq_len(length(img))))
    }
    return(img)
  }
}


#' mosaic in 4D objects
#' @noRd
mosaic_4D <- function(img, time, CLUSTER = NULL, fit_negative = TRUE) {
  img <- img %>% split("band")
  img_names <- names(img)
  # Selecting the best image (1)
  select_img <- (img[1] * 0 + 1) %>%
    st_apply(
      MARGIN = 3,
      FUN = function(img) sum(img, na.rm = T),
      CLUSTER = CLUSTER
    ) %>%
    unlist() %>%
    order(decreasing = T) %>%
    `[`(1)
  # Estimating the time distance with respect to the best image (2)
  img_ord <- order(abs(time - time[select_img]), decreasing = F)

  my_image <- list()
  for (b in seq_along(img)) {
    img_b <- img[b] %>%
      split("time") %>%
      "["(img_ord)
    # Change negative value by NA (3)
    if (fit_negative) {
      img_b[img_b < 0] <- NA
    }

    # Aggregating considering (1) and (2)
    one_img <- img_b %>%
      merge() %>%
      st_apply(
        MARGIN = c("x", "y"),
        FUN = function(img) na.omit(img)[1],
        CLUSTER = CLUSTER
      )
    my_image[[b]] <- one_img
  }
  stars_merge <- function(x) Reduce("c", x)
  my_image_f <- stars_merge(my_image)
  names(my_image_f) <- img_names
  return(my_image_f)
}

#' mosaic in 3D objects
#' @noRd
mosaic_3D <- function(img, time, CLUSTER = NULL, fit_negative = TRUE,
                      RasterLayer = TRUE) {
  # Selecting the best image (1)
  select_img <- img %>%
    st_apply(
      MARGIN = 3,
      FUN = function(img) sum(img, na.rm = T),
      CLUSTER = CLUSTER
    ) %>%
    unlist() %>%
    order(decreasing = T) %>%
    `[`(1)

  # Estimating the time distance with respect to the best image (2)
  img_ord <- order(abs(time - time[select_img]), decreasing = F)
  img <- img %>%
    split("time") %>%
    "["(img_ord)

  # Change negative value by NA (3)
  if (fit_negative) {
    img[img < 0] <- NA
  }

  # Aggregating considering (1) and (2)
  img %>%
    merge() %>%
    st_apply(
      MARGIN = c("x", "y"),
      FUN = function(img) na.omit(img)[1],
      CLUSTER = CLUSTER
    )
}
