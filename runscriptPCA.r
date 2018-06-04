##############  PCA reduction of geomorphological classes in Suriname

### load libraries (install packages if necessary)
library("sp")
library("stats")
library("raster")
library("rgdal")
library("R.matlab")
library("FactoMineR") # may not load running R 3.3 or higher
library("data.table")
library("ade4")
library("graphics")
library("rpart")
library("MASS") # may not load running R 3.3 or higher
library(plotly)
library(factoextra)
library(NbClust)
library(class)
################# ATTENTION!
################# user input required
################# please change paths to relevant raster and shapefiles
################# please alter k-means and nstart to desired values

#### texture analysis inputs
acp320Input <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\input_25112017\\srtmutm_correctedFOTO30-320b1PCA.tif"
acp140Input <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\input_25112017\\srtmutm_correctedFOTO30-140b1PCA.tif"
acp80Input <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\input_25112017\\srtmutm_correctedFOTO30-80b1PCA.tif"

#### all variable inputs
dissectionPath <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\output_15012018\\dissection.tif"
smoothPath <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\output_15012018\\smoothness.tif"
roughPath <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\output_15012018\\roughness.tif"
convexPath <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\output_15012018\\convexity.tif"
altitudePath <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\output_15012018\\altitude.tif"
wetnessPath <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\output_15012018\\wetness.tif"
noisePath <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\output_15012018\\noise.tif"
handPath <-  "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\output_15012018\\hand.tif"

#### k-n variable outputs
outputAcpTextCrit <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\Kmeans_Output\\ACPTexCrit_new.tif"
outputSquk16 <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\Kmeans_Output\\K7_2TexCritnew.tif"


#### extent input. Do not run if you do not have an extent input.
#### an extent will be created further down the code


################# end of user inputs



#Import des données raster (avec les 3 Bandes) des ACP FOTO ayant même étendue (masque : Zone_etude fait dans Arcgis) #
# load foto exported images
ACP30_320 <- brick(acp320Input)

# get extents from smallest image
ACP320extent <- extent(ACP30_320)


ACP30_80 <- setExtent(crop(brick(acp80Input),
                 ACP320extent), ACP30_320)

ACP30_140 <- setExtent(crop(brick(acp140Input),
                            ACP320extent), ACP30_320)


if(!exists("clipExtent")){
  clipExtent <- as(ACP320extent, 'SpatialPolygons')}

#extraction en Matrice suivant l'étendue de la zone d'étude, avec l'outils extract #
extract(ACP30_80, clipExtent) -> Matrice30_80 
extract(ACP30_140, clipExtent) -> Matrice30_140
extract(ACP30_320, clipExtent)-> Matrice30_320 


#Conversion des large list en Matrices
DF30_80<- Map(as.data.frame, Matrice30_80)
rbindlist(DF30_80)-> M30_80

DF30_140 <- Map(as.data.frame, Matrice30_140)
rbindlist(DF30_140)->M30_140

DF30_320 <- Map(as.data.frame, Matrice30_320)
rbindlist(DF30_320)-> M30_320

#Compilation des matrices ACPFoto pour créer la première matrice#
cbind(M30_80, M30_140, M30_320) -> MatrixACP


#import,transformation et création des critères physiques dans la matrice du SRTM90
IndiceDissectionSRTM<- setExtent(crop(brick(dissectionPath),
                                        ACP320extent), ACP30_320)
extract(IndiceDissectionSRTM,clipExtent)-> MatrixDissection
DFDissection<- Map(as.data.frame, MatrixDissection)
rbindlist(DFDissection)->Mdissection

Pente_inf15SRTM <- setExtent(crop(brick(smoothPath),
                                    ACP320extent), ACP30_320)
extract(Pente_inf15SRTM,clipExtent)-> Matricepenteinf15SRTM
DFpenteinf15SRTM<- Map(as.data.frame, Matricepenteinf15SRTM)
rbindlist(DFpenteinf15SRTM)->Mpenteinf15SRTM

Pente_sup30SRTM <- setExtent(crop(brick(roughPath),
                                    ACP320extent), ACP30_320)
extract(Pente_sup30SRTM,clipExtent)-> Matricepentesup30SRTM
DFpentesup30SRTM<- Map(as.data.frame, Matricepentesup30SRTM)
rbindlist(DFpentesup30SRTM)->Mpentesup30SRTM

ConvexitySRTM <- setExtent(crop(brick(convexPath),
                                  ACP320extent), ACP30_320)
extract(ConvexitySRTM,clipExtent)-> MatrixconvexitySRTM
DFconvexitySRTM<- Map(as.data.frame, MatrixconvexitySRTM)
rbindlist(DFconvexitySRTM)->MconvexitySRTM

AltitudeSRTM <- setExtent(crop(brick(altitudePath),
                                 ACP320extent), ACP30_320)
extract(AltitudeSRTM,clipExtent)-> Matricealtitude
DFaltitude<- Map(as.data.frame, Matricealtitude)
rbindlist(DFaltitude)->Maltitude

noiseSRTM <- setExtent(crop(brick(noisePath),
                               ACP320extent), ACP30_320)
extract(noiseSRTM,clipExtent)-> Matrixnoise
DFnoise<- Map(as.data.frame, Matrixnoise)
rbindlist(DFnoise)->Mnoise


wetness <- setExtent(crop(brick(wetnessPath),
                           ACP320extent), ACP30_320)
wetness[is.na(wetness)] <- 0 
extract(wetness,clipExtent) -> MatrixWetness
DFwetness <- Map(as.data.frame, MatrixWetness)
rbindlist(DFwetness) ->MWetness

HANDSRTM <- setExtent(crop(brick(handPath),
                             ACP320extent), ACP30_320)
HANDSRTM[is.na(HANDSRTM)] <- 0 
extract(HANDSRTM,clipExtent) -> MatriceHANDSRTM
DFHANDSRTM <- Map(as.data.frame, MatriceHANDSRTM)
rbindlist(DFHANDSRTM) -> MHANDSRTM


#Compilation des matrices de texture et les critères #  , Mnoise
as.matrix(cbind (MatrixACP, Mpenteinf15SRTM, Mpentesup30SRTM,  MconvexitySRTM, 
                 MHANDSRTM, Maltitude, Mdissection,Mnoise, MWetness )) -> MatrixTexCrit


for (col in 1: ncol(MatrixTexCrit)){
  max <-  quantile(MatrixTexCrit[,col], .99)
  min <-  quantile(MatrixTexCrit[,col], .01)
  
  MatrixTexCrit[,col][MatrixTexCrit[,col] < min] <- min
  MatrixTexCrit[,col][MatrixTexCrit[,col] > max] <- max
}

summary(MatrixTexCrit)


normFun <- function(x){ ((x - min(x)) / (max(x) - min(x)) ) }

MatrixTexCrit <- as.data.frame(lapply(data.frame(MatrixTexCrit), normFun))


dudi.pca(MatrixTexCrit, scannf = FALSE)->ACPTexCrit

barplot(ACPTexCrit$eig/sum(ACPTexCrit$eig))

ncomp=1
sumComp=0
critVal = 85

repeat{
  if(sumComp >= critVal){break}
      sumComp = sumComp + (ACPTexCrit$eig[ncomp]/sum(ACPTexCrit$eig))*100
  ncomp=ncomp+1
  }


#ACP sur la matrice 

dudi.pca(MatrixTexCrit, scannf = FALSE, nf = ncomp)->ACPTexCrit

#créer des raster "squelettes" basés sur la zone d'étude et avec les fenêtres alignées comme la couche 10_27 (3 à faire pour chaque axe)

barplot(ACPTexCrit$eig/sum(ACPTexCrit$eig))
#ACPTexCrit$eig/sum(ACPTexCrit$eig)
#cercle de correlation
s.corcircle(ACPTexCrit$co, xax = 1, yax = 2)

#Classification 
gc()

#PCATextDataframe <-  as.data.frame(lapply(data.frame(ACPTexCrit$li), normFun))

#PCATextDataframe <- PCATextDataframe[sample(nrow(PCATextDataframe),1500),]

#nb <- NbClust(PCATextDataframe, distance = "euclidean", min.nc = 2,
#              max.nc = 8, method = "kmeans")


KnTexCrit<- kmeans(as.data.frame(lapply(data.frame(ACPTexCrit$li), normFun)),
                    5, iter.max = 9999, nstart = round(runif(1, 2, 10)))


cbind( as.data.frame(lapply(data.frame(ACPTexCrit$li), normFun)), "labels" =
matrix(KnTexCrit1$cluster, ncol = 1) ) -> matrix.PCA.kmeans


plot_ly(matrix.PCA.kmeans[sample(nrow(matrix.PCA.kmeans), 1000),], x= ~Axis1, y= ~Axis2, z= ~Axis3,
        color = ~labels, colors = c('#b9936c', '#3e4444', '#dac292', '#82b74b',  '#e6e2d3'))


outputSquk16 <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\Kmeans_Output\\K5_HW_TexCritnew2501_8.tif"
rasterize (clipExtent, ACP30_80)->Squk16TexCrit

values(Squk16TexCrit)[!is.na(values(Squk16TexCrit))]<- KnTexCrit$cluster
writeRaster(Squk16TexCrit, filename = outputSquk16, format = "GTiff", overwrite=TRUE)

#test1<- Squk16TexCrit
#test1[test1 == 5,] <- matrix.LHR.Kmeans$labels
#outputtest <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\Kmeans_Output\\test1.tif"
#writeRaster(test1, filename = outputtest, format = "GTiff", overwrite=TRUE)
#####################################################################
########################  treat data ###################################

#add lables to texture matrix

PCAMatrixLabels     <- matrix.PCA.kmeans

for (i in 1:5){
  
  x <- PCAMatrixLabels[PCAMatrixLabels$labels == i, ]
  
  assign( paste0('PCATextDataframe', i) ,
  x[sample(nrow(x),500),])
 
}

PCATextDataframe <- rbind(PCATextDataframe1,PCATextDataframe2,PCATextDataframe3,PCATextDataframe4,PCATextDataframe5)
PCATextDataframe <- PCATextDataframe[sample(nrow(PCATextDataframe), nrow(PCATextDataframe)),]

#nb <- NbClust(PCATextDataframe[,1:8], distance = "euclidean", min.nc = 3,
#              max.nc = 8, method = "kmeans")

#plot_ly(PCATextDataframe, x= ~Axis1, y= ~Axis2, z= ~Axis3,
#        color = ~labels, colors = c('#b9936c', '#3e4444', '#dac292', '#82b74b',  '#e6e2d3'))


matrix.Means <- matrix(,ncol =ncol(PCAMatrixLabels) -1, nrow = 5)
matrix.SD <- matrix(,ncol=ncol(PCAMatrixLabels) -1, nrow = 5)

for(i in 1:nrow(matrix.Means)){
  for(j in 1:ncol(matrix.Means)){
    matrix.Means[i,j] <- mean(PCAMatrixLabels[PCAMatrixLabels[,9] == i, j] )
    matrix.SD[i,j] <- sd(PCAMatrixLabels[PCAMatrixLabels[,9] == i, j] )
  }
}


  
  for(j in 1:5){
  assign(paste0("trainingSD.",j),  
                     PCAMatrixLabels[(PCAMatrixLabels[,1] >= (matrix.Means[j,1] - matrix.SD[j,1])  &
                     PCAMatrixLabels[,1] <= (matrix.Means[j,1] + matrix.SD[j,1])) &
                                         
                    (PCAMatrixLabels[,2] >= (matrix.Means[j,2] - matrix.SD[j,2])  &
                     PCAMatrixLabels[,2] <= (matrix.Means[j,2] + matrix.SD[j,2])) &
                                          
                    (PCAMatrixLabels[,3] >= (matrix.Means[j,3] - matrix.SD[j,3])  &
                     PCAMatrixLabels[,3] <= (matrix.Means[j,3] + matrix.SD[j,3]))
                    ,c(1,2,3,9)])
         
  assign(paste0("trainingSD.",j),
         cbind(eval(as.name(paste0("trainingSD.",j) )), "labelsEst" =
               matrix(j, nrow = nrow(eval(as.name(paste0("trainingSD.",j) ))), ncol = 1 )))
         
  }

trainingSD <- rbind(trainingSD.1, trainingSD.2, trainingSD.3, trainingSD.4, trainingSD.5)


minQuant <- 0.3 
maxQuant <- 0.7


for(j in 1:5){
 assign(paste0("trainingQuant.", j),
  PCAMatrixLabels[  (PCAMatrixLabels[,1] >= (quantile(PCAMatrixLabels[PCAMatrixLabels[,9] == j,1], minQuant))  &
                     PCAMatrixLabels[,1] <= (quantile(PCAMatrixLabels[PCAMatrixLabels[,9] == j,1], maxQuant))) &
                                    
                    (PCAMatrixLabels[,2] >= (quantile(PCAMatrixLabels[PCAMatrixLabels[,9] == j,2], minQuant))  &
                     PCAMatrixLabels[,2] <= (quantile(PCAMatrixLabels[PCAMatrixLabels[,9] == j,2], maxQuant))) &
                                    
                    (PCAMatrixLabels[,3] >= (quantile(PCAMatrixLabels[PCAMatrixLabels[,9] == j,3], minQuant))  &
                     PCAMatrixLabels[,3] <= (quantile(PCAMatrixLabels[PCAMatrixLabels[,9] == j,3], maxQuant))) 
                     ,c(1,2,3,9)])
  
  assign(paste0("trainingQuant.",j),
         cbind(eval(as.name(paste0("trainingQuant.",j) )), "labelsEst" =
                 matrix(j, nrow = nrow(eval(as.name(paste0("trainingQuant.",j) ))), ncol = 1 )))
  
  }


trainingQuant <- rbind(trainingQuant.1, trainingQuant.2, trainingQuant.3, trainingQuant.4, trainingQuant.5)






######################################################################
########################## knn ###################################

  knnlabelsSD <- as.integer(
    
    knn( trainingSD[,1:3], PCAMatrixLabels[,1:3], trainingSD[,5], k = 30)
    
    )

knnlabelsQuant <- as.integer(
  
  knn( trainingQuant[,1:3], PCAMatrixLabels[,1:3], trainingQuant[,5], k = 30)
  
)


rasterize (clipExtent, ACP30_80)-> knnRaster
values(knnRaster)[!is.na(values(knnRaster))]<- knnlabelsQuant
outputKNN <- "V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\Kmeans_Output\\knn_Quant.tif"
writeRaster(knnRaster, filename = outputKNN, format = "GTiff", overwrite=TRUE)

  
######################################################################
######################## neuralnet ###################################


trainsSample <- 2000
trainingNN <- trainingSD[1:trainsSample,]
threshold <- 0.01

nHidden = c(10,10)

NNformula <- as.formula(paste0("labelsEst ~ ", paste0("Axis",1:3, collapse = "+")))

Nnet <- neuralnet(formula = NNformula, data = trainingNN, hidden = c(nHidden),
                  algorithm = 'rprop+',
                  learningrate.limit = list(min = c(1e-10), max = c(0.01)),
                  learningrate.factor = list(minus = c(0.5), plus = c(1.2)),
                  err.fct = "sse",
                  act.fct = "tanh",
                  threshold = threshold,
                  lifesign = "full",
                  lifesign.step = 500,
                  stepmax = 3e05)

######################################################################
###################### random forest #################################



######################################################################
##################### tensorflow #####################################

