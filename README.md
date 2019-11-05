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

ForesToolboxRS es un paquete R que fue creado con la finalidad de proporcionar una serie de herramientas y algoritmos diversos para el procesamiento y análisis de la Teledetección para la Observaciones de la Tierra, especialmente para el monitoreo de la perturbación de bosques. Todos los algoritmos implementados están basados en publicaciones científicas del más alto nivel. El enfoque PVts-β, un método de detección no estacional, está implementado en este paquete y es capaz de leer datos matriciales y raster. ForesToolboxRs es una iniciativa cuyas funciones están inspirados por el trabajo de [Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving
tropical deforestation detection through using photosynthetic vegetation time series – (PVts-β). Ecological Indicators, 94, 367–379.](https://doi.org/10.1016/j.ecolind.2018.07.012)

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
