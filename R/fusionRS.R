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
#' @importFrom raster getValues as.data.frame brick
#' @importFrom factoextra get_pca_var
#' @param x RasterStack, RasterBrick
#' @param na If TRUE the NA values of the images will be omitted from the analysis.
#' @examples
#' library(ForesToolboxRS)
#' library(raster)
#' data("ForesToolboxRS_dataset")
#' img <- brick(lapply(ForesToolboxRS_dataset$mosaic, '[[',1))
#' img <- smootH(img)
#' fusion <- fusionRS(x=img)
#' plotRGB(fusion[[1]], axes=FALSE, stretch="lin",main ="Fused images")
#' @export
fusionRS <- function(x, na = TRUE) {
  if (is(x,'stars')) {
    x <- brick(mapply(function(z) as(x[z],'Raster'),seq_len(length(x))))
  }

  if (is(x, "RasterStack") | is(x, "RasterBrick")) {
    df <- as.data.frame(x)

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
    acpY <- x
    vr <- getValues(acpY)
    p <- which(!is.na(vr)) # Positions of NA

    # We assign a raster format to each PC
    for (k in 1:dim(df)[2]) {
      acpY[[k]][p] <- acp$x[, k]
    }
  } else {
    stop(class(x), " class is not supported")
  }

  results <- list(FusionRS = acpY, var = var, pov = pov, varAcum = varAcum, corr = corr, contri = contri)
  names(results) <- c(
    "Fused images", "Variance", "Proportion of variance", "Cumulative variance",
    "Correlation", "Contribution in %"
  )
  return(results)
}
