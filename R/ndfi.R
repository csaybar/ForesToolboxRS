#' Normalized Difference Fraction Index (NDFI)
#'
#' The NDFI it is sensitive to the state of the canopy cover, and
#' has been successfully applied to monitor forest degradation and deforestation
#' in Peru and Brazil. This index comes from the endmembers Green Vegetation (GV),
#' non-photosynthetic vegetation (NPV), Soil (S) and the reminder is the shade
#' component.
#'
#' @section References:
#' Souza Jr., C.M., Roberts, D.A., Cochrane, M.A., 2005. Combining spectral and
#' spatialinformation to map canopy damage from selective logging and forest
#' fires. Remote Sens. Environ. 98 (2-3), 329-343.
#' @section Note:
#' The fractions will be obtained from the Spectral Mixture Analysis physical model
#' that was implemented in this function.
#' @param img Imagen Landsat TM, ETM+ o OLI. It could be RasterStack or RasterBrick.
#' @param procesLevel Processing level. It is possible to obtain the NDFI from images
#' in surface reflectance (SR) from TM, ETM+ and OLI, or Top of Atmosphere (TOA) values
#' only for Lansat 8 OLI. The default is SR. In addition, for any processing level,
#' the image values must be rescaled between
#' 0 and 10000.
#' @importFrom dplyr bind_cols
#' @examples
#' \dontrun{
#' library(ForesToolboxRS)
#' library(raster)
#' #'Load an example dataset
#' data(fdata)
#'
#' #Endmembers
#' gv<-c(119.0,  475.0,  169.0,  6250.0, 2399.0, 675.0)
#' npv<-c(1514.0, 1597.0, 1421.0, 3053.0, 7707.0, 1975.0)
#' soil<-c(1799.0, 2479.0, 3158.0, 5437.0, 7707.0, 6646.0)
#' cloud<-c(4031.0, 8714.0, 7900.0, 8989.0, 7002.0, 6607.0)
#' endm <- matrix(c(gv,npv,soil,cloud), 4, 6, byrow = TRUE, dimnames =
#' list(c("gv", "npv","soil","cloud"), c("B1", "B2", "B3","B4","B5","B6")))
#'
#' #Unmix the image
#' ndfi <- ndfiSMA(img_ndfi, procesLevel="SR")
#' }
#' @export
ndfiSMA <- function(img, procesLevel="SR"){

  if (is(img, "RasterStack") | is(img, "RasterBrick")) {
    df <- as.matrix(img)
  } else {
    stop(class(img), " This class is not supported yet.")
  }

  if (procesLevel=="SR") {
    # Endmembers
    endm <- matrix(c(119.0,  475.0,  169.0,  6250.0, 2399.0, 675.0,
                     1514.0, 1597.0, 1421.0, 3053.0, 7707.0, 1975.0,
                     1799.0, 2479.0, 3158.0, 5437.0, 7707.0, 6646.0,
                     4031.0, 8714.0, 7900.0, 8989.0, 7002.0, 6607.0),
                   ncol=6, nrow=4, byrow=TRUE) #GV, NPV, Soil, Cloud
    rownames(endm) <- c("gv", "npv", "Soil", "Cloud")

    sma <- function(img, endm){

      if(isTRUE(dim(df)[2] > dim(endm)[1])){
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

          names(sma_fractions) <- c("gv", "npv", "soil", "cloud", "RMSE")

        } else {
          stop(" The number of values extracted in band should be equal.")
        }
      } else {
        stop(" The number of bands must be greater than the number of endmembers.")
      }

      return(sma_fractions)
    }

    fractions <- sma(img, endm)

    gv <- fractions[[1]]; gv[gv<0]<-0 # Green Vegetation
    npv <- fractions[[2]]; npv[npv<0]<-0 # Non Photosynthetic Vegetation
    soil <- fractions[[3]]; soil[soil<0]<-0 # Soil
    cloud <- fractions[[4]]; cloud[cloud<0]<-0 # Cloud

    # gv + npv + soil + cloud
    summed <- gv+npv+soil+cloud
    gvs <- (gv/summed)*100

    # npv + soil + cloud
    npvSoil <- npv+soil+cloud
    # NDFI
    ndfi_index <- (gvs-npvSoil)/(gvs+npvSoil)*100 + 100

  } else if (procesLevel =="TOA") {
    # Endmembers
    endm <- matrix(c(924.6,  933.1,  455.8,  6418.3, 2445.8, 701.6,
                     2319.6, 2055.1, 1707.8, 3221.3, 7753.8, 2001.6,
                     2604.6, 2937.1, 3444.8, 5605.3, 7753.8, 6672.6,
                     4836.6, 9172.1, 8186.8, 9157.3, 7048.8, 6633.6),
                   ncol=6, nrow=4, byrow=TRUE) #GV, NPV, Soil, Cloud
    rownames(endm) <- c("gv", "npv", "Soil", "Cloud")

    sma <- function(img, endm){

      if(isTRUE(dim(df)[2] > dim(endm)[1])){
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

          names(sma_fractions) <- c("gv", "npv", "soil", "cloud", "RMSE")

        } else {
          stop(" The number of values extracted in band should be equal.")
        }
      } else {
        stop(" The number of bands must be greater than the number of endmembers.")
      }

      return(sma_fractions)
    }

    fractions <- sma(img, endm)

    gv <- fractions[[1]]; gv[gv<0]<-0 # Green Vegetation
    npv <- fractions[[2]]; npv[npv<0]<-0 # Non Photosynthetic Vegetation
    soil <- fractions[[3]]; soil[soil<0]<-0 # Soil
    cloud <- fractions[[4]]; cloud[cloud<0]<-0 # Cloud

    # gv + npv + soil + cloud
    summed <- gv+npv+soil+cloud
    gvs <- (gv/summed)*100

    # npv + soil + cloud
    npvSoil <- npv+soil+cloud
    # NDFI
    ndfi_index <- (gvs-npvSoil)/(gvs+npvSoil)*100 + 100

  } else stop(" Processing level not supported.")

  return(ndfi_index)
}

