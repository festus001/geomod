#' @param
#' @return
#' @export
#' Step 1: Load the libraries
library(rgdal);
library(MASS);
library(ggplot2);
library(maptools);
library(spatstat);
library(rgdal);
library(MASS);
library(Hmisc);
library(survival);
library(splines2);
library(gstat);
library(nlme);
library(nnet);
library(lattice)

#' Step2: Set up working directory
setwd("D:/Academic/MSC1/Data/NgenoF.")
#getwd()


#' Step3: Import the data and explore data structures
data=read.table("SPT32.txt", header=T);#memory.limit(size=4500)
str(data)
Spline_SPT<- function (Borehole,row, path)
  data_path <- paste(path, "/", Borehole,row, ".txt", sep="")
x<-c(data$X)
y<-c(data$Y)
Depth<-c(data$Depth)
SPT<-c(data$SPT)
plot(Depth,SPT,col = "red", pch=19, main = paste0("",format(row)))
spline(x,y=NULL, n = 14*length(1), method = "hyman",xmin=1,xmax=14,ties=mean)
#lines(spline(Depth,SPT,n=14,method="hyman",xmin=1,xmax=14,ties=mean))
#out<-spline(Depth, SPT, n = 14, method = "hyman",xmin = 1, xmax = 14, ties = mean)
lines(smooth.spline(Depth, SPT), col = 2)


#' R code for interpolating geotechnical data
#' Step 4: Converted the data into spatial data-frame
#coordinates(data)=~X+Y;
coordinates(data)=c("X","Y")
proj4string(data)=CRS("+proj=utm +ellps=WGS84")
str(data)

#' Step 5: Import the environmental correlates and convert them into spatial dataframe
predictors=readGDAL("lithology.asc")
predictors$PC2=readGDAL("dem.asc")$band1
predictors$PC3=readGDAL("morphology.asc")$band1
predictors$PC1=predictors$band1
predictors$band1=NULL
object.size(predictors)
proj4string(predictors)=CRS("+proj=utm +ellps=WGS84")
newpredictors <- spTransform(predictors, CRS("+proj=utm +ellps=WGS84"))
str(predictors)
summary(predictors)
predictors

#' Step 6 - Overlay the environmental correlates
#overlay the predictors
predictors.ov=over(data,predictors)
data$PC1=predictors.ov$PC1
data$PC2=predictors.ov$PC2
data$PC3=predictors.ov$PC3
str(predictors)
predictors.ov
summary(predictors.ov)

#' Step 7 - Soil Property Prediction Model (Regression)
por.lme=lm((SPT)~(PC1+PC2+PC3), data)
cor(fitted(por.lme),data$SPT) # The  variation accounted by large-scale trends
summary(por.lme) # To obtain summary results of the fitted model
hist(residuals(por.lme)) # To check for normal distribution of the residuals
plot(fitted(por.lme)~data$SPT) # To obtain a plot of the comparison of model fit
predictors$PredSPT=predict(por.lme, predictors)
plot(predictors$PredSPT)

#' Step 8 - Variogram analysis of the autocorrelation
por.rev=variogram(residuals(por.lme)~1, data,alpha=c(0,90,180,270))
plot(por.rev)
por.rvgm=fitted(por.rev, vgm(nugget=100,model="Gau",range=3000, sill=200))
plot(por.rev, por.rvgm, plot.nu=F)
por.rk=krige(residuals(por.lme)~1, data, predictors, por.rvgm)
str(por.rk)

spplot(por.rk)#spplot(por.rk["var1.var"],scales=list(draw=T))

#' Step 9 - To export the output to other GIS software
writeGDAL(por.rk["var1.pred"], "residnew.mpr","ILWIS")#writeGDAL(predictors["PredPI2"])



