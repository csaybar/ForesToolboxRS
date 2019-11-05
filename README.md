<img src="https://raw.githubusercontent.com/ytarazona/ForesToolboxRS/master/man/figures/logo.png" align="right" width = 15%/>

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/csaybar/forestoolboxrs?branch=dev&svg=true)](https://ci.appveyor.com/project/csaybar/forestoolboxrs)
[![Travis build
status](https://travis-ci.org/csaybar/ForesToolboxRS.svg?branch=master)](https://travis-ci.org/csaybar/ForesToolboxRS)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Codecov test
coverage](https://codecov.io/gh/csaybar/ForesToolboxRS/branch/master/graph/badge.svg)](https://codecov.io/gh/csaybar/ForesToolboxRS?branch=dev)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#experimental)

## ForesToolboxRS

ForesToolboxRS is a R package that was created to provide a variety of tools and algorithms for processing and analyzing Remote Sensing for Earth Observations, especially for monitoring forest disturbance. All implemented algorithms are based on scientific publications of the highest level. The PVts-β approach, a non-seasonal detection method, is implemented in this package and is capable of reading matrix and raster data. ForesToolboxRs is an initiative whose functions are inspired by the work of [Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical deforestation detection through using photosynthetic vegetation time series – (PVts-β). Ecological Indicators, 94, 367–379.](https://doi.org/10.1016/j.ecolind.2018.07.012)

# Getting Started

## Install

    library(devtools)
    install_github("ytarazona/ForesToolboxRS")

## Some functions

    library(ForesToolboxRS)
    ?MosaicFreeCloud
    ?smootH
    ?pvts
    ?fusionRS
