#' Spectral Mixture Analysis (SMA)
#'
#' The SMA assumes that the energy received, within the field of vision of the
#' remote sensor, can be considered as the sum of the energies received from each
#' dominant endmember. This function addresses a Linear Mixing Model.
#'
#' @section References:
#' Adams, J. B., Smith, M. O., & Gillespie, A. R. (1993). Imaging spectroscopy:
#' Interpretation based on spectral mixture analysis. In C. M. Pieters & P.
#' Englert (Eds.), Remote geochemical analysis: Elements and mineralogical
#' composition. NY: Cambridge Univ. Press 145-166 pp.
#'
#' Shimabukuro, Y.E. and Smith, J., (1991). The least squares mixing models to
#' generate fraction images derived from remote sensing multispectral data.
#' IEEE Transactions on Geoscience and Remote Sensing, 29, pp. 16-21.
#' @section Note:
#' A regression analysis is used to obtain the fractions. In least squares
#' inversion algorithms, the common objective is to estimate abundances that
#' minimize the squared error between the actual spectrum and the estimated spectrum.
#' The values of the fractions will be between 0 and 1.
#' @param img Optical images. It could be RasterStack or RasterBrick.
#' @param endm Endmembers must be a matrix or data.frame and with more than one edmember.
#' Rows represent the endmembers and columns represent the spectral bands.  The number
#' of bands must be greater than the number of endmembers.
#' @importFrom dplyr bind_cols
#' @export
#' @examples
#' \dontrun{
#' library(ForesToolboxRS)
#' library(raster)
#' library(dplyr)
#'
#' #Load an example dataset
#'data(fdata)
#'
#'#Endmembers
#'soil<-c(0.8980,0.8508,0.8704,1.3636,1.6579,1.1420)
#'forest<-c(0.8207,0.7545,0.6548,1.6463,0.9725,0.6673)
#'water<-c(0.9276,0.9570,1.0089,0.6743,0.5220,0.5143)
#'endm <- matrix(c(soil,forest,water), 3, 6, byrow = TRUE, dimnames =
#'list(c("soil", "forest","water"), c("B1", "B2", "B3","B4","B5","B6")))
#'
#'#Unmix the image
#'fractions <- sma(img_ndfi, endm)
#'}
sma <- function(img, endm){

  if (is(img, "RasterStack") | is(img, "RasterBrick")) {
    df <- as.matrix(img)
  } else {
    stop(class(img), " This class is not supported yet.")
  }

  if(dim(df)[2] > dim(endm)[1]){
    if(dim(df)[2] == dim(endm)[2]){

      # Transposed from the endmembers matrix
      M<-t(endm)

      # We calculate fractions through least squares
      lmm <- lapply(1:dim(df)[1], function(i) f<-(tcrossprod((solve(crossprod(M))),M))%*%df[i,])
      n <- bind_cols(lmm)
      val <- t(cbind(n))

      # We estimate the RMSE
      fra_esti <- lapply(1:dim(val)[1], function(i) f <- M%*%val[i,])
      m_df <- lapply(as.list(1:dim(df)[1]), function(x) df[x[1],])
      rmse <- mapply(function(x,y) sqrt(mean((x - y)^2)), m_df, fra_esti)

      # We store the fractions on a raster
      fractions <- cbind(val,rmse)

      sma_fractions <- img[[1:(dim(endm)[1]+1)]]
      for (j in 1:(dim(endm)[1]+1)) {
        sma_fractions[[j]][]<- fractions[,j]
      }

      names(sma_fractions) <- c(as.vector(row.names(endm)), "RMSE")

    } else {
      stop(" The number of values extracted in band should be equal.")
    }
  } else {
    stop(" The number of bands must be greater than the number of endmembers.")
  }

  return(sma_fractions)
}
