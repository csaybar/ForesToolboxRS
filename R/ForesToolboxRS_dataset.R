#' ftb_datasets: ForestToolboxRS calculated from Landsat 8 Imagery
#'
#'These images are part of the datasets used in Tarazona et al., 2018 (see notes)
#'covering a deforestation hotspot region ubicated in Madre de Dios, a province
#'in south-eastern Peru.
#'
#' @format ForestToolboxRS_datasets is an object of class list of length 2, each element contain the following information: \cr
#' \itemize{
#'   \item mosaic: A star object with 3 variables that represents the photosynthetic
#'   vegetation values for the dry season of the year 2016 (the values were estimated from Landsat-8).
#'   \item ForestToolboxRS: A star object with 28 variables that represents the year mosaic of the
#'   photosynthetic vegetation for the period 1990-2018(calculated from Landsat-8).
#' }
#' @note
#' Tarazona, Y., Mantas, V. M., & Pereira, A. J. S. C. (2018). Improving tropical
#' deforestation detection through using photosynthetic vegetation time series–(ForestToolboxRS).
#' Ecological indicators, 94, 367-379.
#'
"ForesToolboxRS_dataset"
