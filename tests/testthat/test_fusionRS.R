library(testthat)
library(ForesToolboxRS)
library(stars)
library(raster)
library(forecast)
context("ForesToolboxRS::fusionRS")
test_that("fusionRS", {
  # img <- system.file("simple_mosaic", package="ForesToolboxRS") %>%
  #   list.files("\\.tif$",full.names = TRUE) %>%
  #   stack
  # img <- smootH(img)
  # fusion <- fusionRS(x=img)
  # expect_equal(max(raster::getValues(raster::mean(fusion[[1]]))),1.494531,tolerance=0.0001)
})
