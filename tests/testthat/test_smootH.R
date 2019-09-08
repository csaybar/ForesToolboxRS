library(testthat)
library(ForesToolboxRS)
library(stars)
library(raster)
library(forecast)
context("ForesToolboxRS::smootH")

#tarazona code
smootH_ytarazona <- function(x, interp="na.interp"){

  if (is.vector(x)) {
    x[x==0 | x== -1] <- NA
    x[summary(x)[7] >= (length(x)-1)] <- 100

    # Type of interpolation
    if (interp=="na.interp") {
      x <- na.interp(x)

    } else if (interp=="na.StructTS") {
      x <- na.StructTS(x)

    } else if (interp=="na.approx") {
      x <- na.approx(x)

    } else if (interp=="na.spline") {
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
      if (interp=="na.interp") {
        x[i,] <- na.interp(x[i,])

      } else if (interp=="na.StructTS") {
        x[i,] <- na.StructTS(x[i,])

      } else if (interp=="na.approx") {
        x[i,] <- na.approx(x[i,])

      } else if (interp=="na.spline") {
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
    np <- raster::as.matrix(x)
    for (i in 1:dim(np)[1]) {
      np[i,][np[i,] == 0 | np[i,] == -1]<- NA
      np[i,][summary(np[i,])[7] >= (dim(np)[2]-1)] <- 100

      # Type of interpolation
      if (interp=="na.interp") {
        np[i,] <- na.interp(np[i,])

      } else if (interp=="na.StructTS") {
        np[i,] <- na.StructTS(np[i,])

      } else if (interp=="na.approx") {
        np[i,] <- na.approx(np[i,])

      } else if (interp=="na.spline") {
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

test_that("smootH numeric-na.interp", {
  x <- c(1,2,4,5,6,7,8,9-100,8,9,2,5,10)
  suavizado <- smootH(x)
  expect_equal(mean(suavizado),6.423077, tolerance = 0.00001)
})

test_that("smootH numeric-na.approx", {
  x <- c(1,2,4,5,6,7,8,9-100,8,9,2,5,10)
  suavizado <- smootH(x,"na.approx")
  expect_equal(mean(suavizado),6.423077, tolerance = 0.00001)
})


test_that("smootH numeric-na.spline", {
  x <- c(1,2,4,5,6,7,8,9-100,8,9,2,5,10)
  suavizado <- smootH(x, "na.spline")
  expect_equal(mean(suavizado),6.423077, tolerance = 0.00001)
})

test_that("smootH matrix", {
  data("ForesToolboxRS_dataset")
  pv_serie <- ForesToolboxRS_dataset$PVts
  pv_serie_matrix <- raster::as.matrix(pv_serie)
  ytarazona_results <- smootH_ytarazona(pv_serie_matrix)
  new_approach <- smootH(pv_serie_matrix)
  expect_equal(ytarazona_results,new_approach, check.attributes = FALSE, tolerance = 0.00001)
})


test_that("smootH RasterStack", {
  data("ForesToolboxRS_dataset")
  pv_serie <- ForesToolboxRS_dataset$PVts
  pv_serie_stk <- raster::stack(pv_serie)
  ytarazona_results <- smootH_ytarazona(pv_serie_stk)
  new_approach <- smootH(pv_serie_stk)
  expect_equal(ytarazona_results,raster::as.matrix(new_approach), check.attributes = FALSE, tolerance = 0.00001)
})

test_that("smootH RasterBrick", {
  data("ForesToolboxRS_dataset")
  pv_serie <- ForesToolboxRS_dataset$PVts
  pv_serie_brk <- pv_serie
  ytarazona_results <- smootH_ytarazona(pv_serie_brk)
  new_approach <- smootH(pv_serie_brk)
  expect_equal(ytarazona_results,raster::as.matrix(new_approach), check.attributes = FALSE, tolerance = 0.00001)
})

test_that("smootH stars", {
  data("ForesToolboxRS_dataset")
  pv_serie <- ForesToolboxRS_dataset$PVts
  pv_serie_brk <- smootH(pv_serie)
  pv_serie_stars <- smootH(st_as_stars(pv_serie)) %>% merge
  expect_equal(mean(raster::getValues(raster::mean(pv_serie_brk))),mean(pv_serie_stars$X), tolerance = 0.00001)
})
