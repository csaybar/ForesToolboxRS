#' Change detection using the PVts-\eqn{\beta} approach
#'
#' This algorithm will allow to detect disturbances in the forests using
#' all the available Landsat set. In fact, it can also be run with sensors
#' such as MODIS.
#'
#' @section References:
#' Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
#' deforestation detection through using photosynthetic vegetation time
#' series – (PVts-β). Ecological Indicators, 94, 367 379.
#' @section Note:
#' In order to optimise the detections, it is advisable to make a smoothing
#' through the \link[ForesToolboxRS]{smootH} function before detecting changes. The smoothing will
#' allow to eliminate outliers that were not eliminated during the masking
#' of atmospheric artifacts.
#'
#' @param x Matrix, RasterStack, RasterBrick
#' @param start The start of the monitoring time
#' @param end The end of the monitoring time
#' @param threslhold The default threshold is 5 for photosynthetic vegetation,
#' while for indices such as NDVI and EVI the threshold is 3.
#' Please see Tarazona et al. (2018) for more details.
#' @param img The last monitoring image (in case "x" is a matrix)
#' @export
#' @examples
#' library(ForesToolboxRS)
#' library(raster)
#'
#'## Not run:
#'imgs <- as.matrix(stack(list)) # imgs is a matrix
#'dim(imgs)[1]= 1000 is the number of pixels (time series)
#'dim(imgs)[2]=29 is the number of images between 1990 and 2018
#'Change monitoring period 2015-2018. 2015(position 26),2018(position 29)
#'
#'cd <- pvts(x=imgs, start=26, end=29, threshold=5,img=lastimg)
#'plot(cd)
#'
fusionRS<- function(img, ) {

  if (is(img,'RasterStack') | is(img,'RasterBrick')) {

    df <- as.data.frame(img)
    acp<-prcomp(na.omit(df), center = TRUE, scale = TRUE)
    varAcum <- summary(acp)$importance[3,]




  } else if (is(x,'RasterStack') | is(x,'RasterBrick')) {

    breakR <- x[[end]]; x <- as.matrix(x)
    mean<-apply(x[,1:start], 1, mean)
    std<-apply(x[,1:start], 1, sd)
    cd<- ifelse(x[,end] < (mean-threshold*std), 1, 0)
    values(breakR) <- cd; breakR[img<80 | img<0.8]<- 0

  }
  return(breakR)
}


# Convertimos los datos a data.frame
df <- as.data.frame(img)
# ANALISIS EN COMPONENTES PRINCIAPLES -------------------------------------
cor(na.omit(df))
acp<-prcomp(na.omit(df), scale = TRUE, center = TRUE)
summary(acp)
# Ploteo del ACP
fviz_pca_var(acp, col.var = "contrib", bgradient.cols = c("#00AFBB", "#E7B800",
                                                          "#FC4E07"), repel = TRUE
)
# Porcentaje de varianza
eigs <- acp$sdev^2
pv2018 <- c(eigs[1]/sum(eigs),eigs[2]/sum(eigs),eigs[3]/sum(eigs),
            eigs[4]/sum(eigs),eigs[5]/sum(eigs),eigs[6]/sum(eigs),
            eigs[7]/sum(eigs),eigs[8]/sum(eigs))
pv2018_acum<-c(0.6328,0.8350,0.9100,0.97351,0.99485,0.99864,0.99967,1)*100

#Contribuci?n de cada variable a los PC
var <- get_pca_var(acp)
var$contrib
# correlacion
var$cor
# Almacenando en un raster stack
acp2018<-img
vr <- getValues(acp2018)
i <- which(!is.na(vr)) # posiciones de los no NA
# Un bucle para almacenar todas las bandas
for (j in 1:8) {
  acp2018[[j]][i]<- acp$x[,j]
}
# Desvia Standar de un corte de las imagens
shp <- readOGR("H:/PaperSAR-Optico/Datos/shp", "CorteSD")
corte <- crop(acp2018, extent(shp))
sd2018 <- apply(as.matrix(corte), 2, sd)
sd2018_c <- apply(na.omit(as.matrix(acp2018)), 2, sd) # imagen completa

par(mfrow = c(1,2))
plotRGB(img, 3,2,1, axes=T, stretch="lin")
plotRGB(acp2018, 8,8,8, axes=F, stretch="lin",main ="ACP Sentinel-2")
writeRaster(acp2018, "ACP_2018.tif", drivername="GTiff", datatype = "FLT4S")
