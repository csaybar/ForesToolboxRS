library(testthat)
library(ForesToolboxRS)
library(stars)
library(raster)
library(forecast)
context("ForesToolboxRS::MosaicFreeCloud")
rasterio = list(nXOff = 50, nYOff = 50)


test_that("MosaicFreeCloud-character", {
  img <- system.file("simple_mosaic", package="ForesToolboxRS") %>%
    list.files("\\.tif$",full.names = TRUE)
  free_img <- MosaicFreeCloud(img, time, RasterIO=rasterio)
  expect_equal(mean(getValues(free_img),na.rm=TRUE),178.6082,tolerance=0.0001)
})

test_that("MosaicFreeCloud-RasterStack", {
  img <- system.file("simple_mosaic", package="ForesToolboxRS") %>%
    list.files("\\.tif$",full.names = TRUE) %>%
    stack()
  free_img <- MosaicFreeCloud(img, time)
  expect_equal(mean(getValues(free_img),na.rm=TRUE), 192.016, tolerance=0.0001)
})

test_that("MosaicFreeCloud-RasterBrick", {
  img <- system.file("simple_mosaic", package="ForesToolboxRS") %>%
    list.files("\\.tif$",full.names = TRUE) %>%
    stack() %>%
    brick()
  free_img <- MosaicFreeCloud(img, time, RasterIO=rasterio)
  expect_equal(mean(getValues(free_img),na.rm=TRUE),192.016,tolerance=0.0001)
})

test_that("MosaicFreeCloud-stars", {
  img <- system.file("simple_mosaic", package="ForesToolboxRS") %>%
    list.files("\\.tif$",full.names = TRUE) %>%
    read_stars(RasterIO=rasterio)
  free_img <- MosaicFreeCloud(img, time)
  expect_equal(mean(getValues(free_img),na.rm=TRUE), 178.6082,tolerance=0.0001)
})
