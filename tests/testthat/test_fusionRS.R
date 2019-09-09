library(testthat)
library(ForesToolboxRS)
library(stars)
library(raster)
library(forecast)
context("ForesToolboxRS::fusionRS")
test_that("fusionRS", {
  img <- system.file("mosaic", package="ForesToolboxRS") %>%
    list.files("\\.tif$",full.names = TRUE) %>%
    lapply(brick) %>%
    lapply('[[',1) %>%
    brick
  img <- smootH(img)
  fusion <- fusionRS(x=img)
  expect_equal(max(raster::getValues(raster::mean(fusion[[1]]))),11.85757,tolerance=0.0001)
})
