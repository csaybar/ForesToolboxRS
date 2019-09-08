#' Is a character, Raster or stars object?
#' @importFrom methods is
#' @importFrom stars read_stars
#' @param x object to evaluate
#' @return a stars object
#' @noRd
ftb_whatkinditis <-function(x) {
  if (is.character(x)) {
    x = read_stars(x)
  } else if(is(x,'stars')) {
    x
  } else if (is(x,'RasterLayer') | is(x,'RasterStack') | is(x,'RasterBrick')) {
    x = st_as_stars(x)
  } else {
    stop(class(x), ' class is not supported')
  }
  return(x)
}
