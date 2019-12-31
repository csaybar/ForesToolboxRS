#' Fusion of images with different observation geometry.
#'
#' This algorithm allows to fusion images coming from different spectral sensors
#' (e.g., optical-optical, optical and SAR or SAR-SAR). Among many of the qualities
#' of this function, it is possible to obtain the contribution (%) of each variable
#' in the fused image.
#'
#' @section References:
#' Tarazona et al...
#'
#' @section Note:
#' Before executing the function, it is recommended that images coming from different
#' sensors or from the same sensor have a co-registration.
#' @importFrom stats prcomp na.omit
#' @importFrom raster getValues as.data.frame brick extent plotRGB
#' @importFrom factoextra get_pca_var
#' @param x Optical image. It could be RasterStack or RasterBrick
#' @param y Radar image. It could be RasterStack or RasterBrick
#' @param na If TRUE the NA values of the images will be omitted from the analysis.
#' @examples
#' library(ForesToolboxRS)
#' library(raster)
#' library(factoextra)
#'
#' # Optical images
#' b1 <- raster(ncol = 100, nrow=100, val = sample(1:2e+15, 10000))
#' b2 <- raster(ncol = 100, nrow=100, val = sample(1:2e+15, 10000))
#' optical <- stack(b1,b2)
#' # Radar images
#' vv <- raster(ncol = 100, nrow=100, val = sample(1:2e+15, 10000))
#' vh <- raster(ncol = 100, nrow=100, val = sample(1:2e+15, 10000))
#' radar <- stack(vv,vh)
#' # Fusion
#' fusion <- fusionRS(x=optical, y=radar)
#' #plotRGB(fusion[[1]], 1,2,3, axes=F, stretch="lin",main ="Fused images")
#'
#' @export
fusionRS <- function(x, y, na = FALSE) {

  if (is(x,'stars')) {
    x <- brick(mapply(function(z) as(x[z],'Raster'),seq_len(length(x))))
  }

  if (is(x, "RasterStack") | is(x, "RasterBrick") & is(y, "RasterStack") | is(y, "RasterBrick")) {

    if (extent(x)==extent(y)){
      img <- stack(x, y)
    } else {
      stop(" The extent of the images are different.")
    }
  } else {
    stop(c(class(x),class(y)), " This classes are not supported yet.")
  }

  df <- as.data.frame(img)

  # Omit NA
  if (na) {
    df <- na.omit(df)
  }

  acp <- prcomp(df, center = TRUE, scale = TRUE) # standardized variables
  var <- summary(acp)$importance[1, ]^2 # Variance
  pov <- summary(acp)$importance[2, ] # Proportion of variance
  varAcum <- summary(acp)$importance[3, ] # Cumulative variance
  corr <- get_pca_var(acp)$cor # Correlation
  contri <- get_pca_var(acp)$contrib # Contribution in %

  for (i in 1:dim(df)[2]) {
    colnames(corr)[i] <- paste("PC", i, sep = "")
    rownames(corr)[i] <- paste("Band", i, sep = "")
    colnames(contri)[i] <- paste("PC", i, sep = "")
    rownames(contri)[i] <- paste("Band", i, sep = "")
  }

  # Storing in a raster the principal components
  acpY <- img
  vr <- getValues(acpY)
  p <- which(!is.na(vr)) # Positions of NA

  # We assign a raster format to each PC
  for (k in 1:dim(df)[2]) {
    acpY[[k]][p] <- acp$x[, k]
  }

  results <- list(FusionRS = acpY, var = var, pov = pov, varAcum = varAcum, corr = corr, contri = contri)
  names(results) <- c(
    "Fused images", "Variance", "Proportion of variance", "Cumulative variance",
    "Correlation", "Contribution in %"
  )
  return(results)
}
