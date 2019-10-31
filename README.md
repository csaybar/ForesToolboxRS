<img src="https://raw.githubusercontent.com/csaybar/ForesToolboxRS/master/man/figures/logo.png" align="right" width = 15%/>

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

The PVts-β approach is implemented in R with the ForesToolboxRS package.
This algorithm will allow to detect disturbances in the forests using
all the available Landsat set. In fact, it can also be run with sensors
such as MODIS (Moderate Resolution Imaging Spectroradiometer). In
essence, the PVts-β approach is an algorithm that does not model the
seasonal component in time series, instead it only requires calculating
the mean and standard deviation of the time series itself to detect
disturbances. Therefore, the PVts-β approach is a method: i) simple and
intelligent that does not model seasonality, ii) has only one
calibration parameter and iii) can be easily implemented by any standard
user.

For any work you will submit, credits from this algorithm must be given
to:

[Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving
tropical deforestation detection through using photosynthetic vegetation
time series – (PVts-β). Ecological Indicators, 94,
367–379.](https://doi.org/10.1016/j.ecolind.2018.07.012)

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
