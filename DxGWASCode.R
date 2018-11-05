setwd("/Users/jmshook/Desktop")
source("https://bioconductor.org/biocLite.R")
biocLite("Biobase", "multtest")
library(Biobase)
library(multtest)
install.packages('Rcpp', dependencies = TRUE)
library(Rcpp)
install.packages('gplots')
install.packages('LDheatmap')
install.packages('genetics')
install.packages('compiler')
install.packages('scatterplot3d')
library(gplots)
library(LDheatmap)
library(genetics)
library(compiler)
library(scatterplot3d)
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
source("http://www.zzlab.net/GAPIT/emma.txt")
GD <- read.table("all_GD.txt", header = T) #genotypic data file
myGM <- read.table("all_GM.txt",header = T)#genetic map file
myGM$chr <- as.numeric(substr(myGM$chr,3,4))#Convert chromosome number
################
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("bp488001.csv",header = T)#phenotypic data
y <- Y[names(Y)=='bp488001']
line = Y[names(Y) %in% c('acid', 'bp488001')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/bp488001")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("bpmvall.csv",header = T)#phenotypic data
y <- Y[names(Y)=='bpmvall']
line = Y[names(Y) %in% c('acid', 'bpmvall')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/bpmvall")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("bsr491584.csv",header = T)#phenotypic data
y <- Y[names(Y)=='bsr491584']
line = Y[names(Y) %in% c('acid', 'bsr491584')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/bsr491584")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("bsrcode491584.csv",header = T)#phenotypic data
y <- Y[names(Y)=='bsrcode491584']
line = Y[names(Y) %in% c('acid', 'bsrcode491584')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/bsrcode491584")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("bsrcode492477.csv",header = T)#phenotypic data
y <- Y[names(Y)=='bsrcode492477']
line = Y[names(Y) %in% c('acid', 'bsrcode492477')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/bsrcode492477")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
##########
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("bsrall.csv",header = T)#phenotypic data
y <- Y[names(Y)=='bsrall']
line = Y[names(Y) %in% c('acid', 'bsrall')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/bsrall")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("bsrcodeall.csv",header = T)#phenotypic data
y <- Y[names(Y)=='bsrcodeall']
line = Y[names(Y) %in% c('acid', 'bsrcodeall')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/bsrcodeall")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("crds.csv",header = T)#phenotypic data
y <- Y[names(Y)=='crds']
line = Y[names(Y) %in% c('acid', 'crds')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/crds")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("crdr.csv",header = T)#phenotypic data
y <- Y[names(Y)=='crdr']
line = Y[names(Y) %in% c('acid', 'crdr')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/crdr")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("Fe026.csv",header = T)#phenotypic data
y <- Y[names(Y)=='Fe026']
line = Y[names(Y) %in% c('acid', 'Fe026')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/Fe026")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("Fe32perc.csv",header = T)#phenotypic data
y <- Y[names(Y)=='Fe32perc']
line = Y[names(Y) %in% c('acid', 'Fe32perc')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/Fe32perc")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("pmv.csv",header = T)#phenotypic data
y <- Y[names(Y)=='pmv']
line = Y[names(Y) %in% c('acid', 'pmv')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/pmv")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("Fe2.csv",header = T)#phenotypic data
y <- Y[names(Y)=='Fe2']
line = Y[names(Y) %in% c('acid', 'Fe2')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/Fe2")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("Fe2.10001.csv",header = T)#phenotypic data
y <- Y[names(Y)=='Fe2.10001']
line = Y[names(Y) %in% c('acid', 'Fe2.10001')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/Fe2.10001")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("Fe2.10002.csv",header = T)#phenotypic data
y <- Y[names(Y)=='Fe2.10002']
line = Y[names(Y) %in% c('acid', 'Fe2.10002')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/Fe2.10002")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("Fe2.10003.csv",header = T)#phenotypic data
y <- Y[names(Y)=='Fe2.10003']
line = Y[names(Y) %in% c('acid', 'Fe2.10003')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/Fe2.10003")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
#############
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("Fe2.10004.csv",header = T)#phenotypic data
y <- Y[names(Y)=='Fe2.10004']
line = Y[names(Y) %in% c('acid', 'Fe2.10004')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/Fe2.10004")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1']
line = Y[names(Y) %in% c('acid', 'PRR1')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.10001.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.10001']
line = Y[names(Y) %in% c('acid', 'PRR1.10001')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.10001")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.10002.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.10002']
line = Y[names(Y) %in% c('acid', 'PRR1.10002')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.10002")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.10003.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.10003']
line = Y[names(Y) %in% c('acid', 'PRR1.10003')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.10003")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.10004.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.10004']
line = Y[names(Y) %in% c('acid', 'PRR1.10004')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.10004")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.11002.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.11002']
line = Y[names(Y) %in% c('acid', 'PRR1.11002')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.11002")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.11003.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.11003']
line = Y[names(Y) %in% c('acid', 'PRR1.11003')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.11003")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.488001.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.488001']
line = Y[names(Y) %in% c('acid', 'PRR1.488001')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.488001")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.461592.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.461592']
line = Y[names(Y) %in% c('acid', 'PRR1.461592')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.461592")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.492577.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.492577']
line = Y[names(Y) %in% c('acid', 'PRR1.492577')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.492577")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR1.492990.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR1.492990']
line = Y[names(Y) %in% c('acid', 'PRR1.492990')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR1.492990")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR2.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR2']
line = Y[names(Y) %in% c('acid', 'PRR2')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR2")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR3.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR3']
line = Y[names(Y) %in% c('acid', 'PRR3')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR3")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR3.491592.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR3.491592']
line = Y[names(Y) %in% c('acid', 'PRR3.491592')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR3.491592")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR3.492577.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR3.492577']
line = Y[names(Y) %in% c('acid', 'PRR3.492577')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR3.492577")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR3.492758.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR3.492758']
line = Y[names(Y) %in% c('acid', 'PRR3.492758')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR3.492758")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("Prr3.492990.csv",header = T)#phenotypic data
y <- Y[names(Y)=='Prr3.492990']
line = Y[names(Y) %in% c('acid', 'Prr3.492990')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/Prr3.492990")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR4.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR4']
line = Y[names(Y) %in% c('acid', 'PRR4')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR4")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR4.492758.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR4.492758']
line = Y[names(Y) %in% c('acid', 'PRR4.492758')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR4.492758")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR4.492990.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR4.492990']
line = Y[names(Y) %in% c('acid', 'PRR4.492990')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR4.492990")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR5.492990.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR5.492990']
line = Y[names(Y) %in% c('acid', 'PRR5.492990')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR5.492990")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR5.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR5']
line = Y[names(Y) %in% c('acid', 'PRR5')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR5")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR7.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR7']
line = Y[names(Y) %in% c('acid', 'PRR7')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR7")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR7.491404.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR7.491404']
line = Y[names(Y) %in% c('acid', 'PRR7.491404')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR7.491404")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR7.491592.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR7.491592']
line = Y[names(Y) %in% c('acid', 'PRR7.491592')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR7.491592")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR7.492448.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR7.492448']
line = Y[names(Y) %in% c('acid', 'PRR7.492448')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR7.492448")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR7.492990.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR7.492990']
line = Y[names(Y) %in% c('acid', 'PRR7.492990')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR7.492990")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR10.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR10']
line = Y[names(Y) %in% c('acid', 'PRR10')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR10")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR12.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR12']
line = Y[names(Y) %in% c('acid', 'PRR12')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR12")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR17.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR17']
line = Y[names(Y) %in% c('acid', 'PRR17')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR17")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR17.491404.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR17.491404']
line = Y[names(Y) %in% c('acid', 'PRR17.491404')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR17.491404")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR17.492448.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR17.492448']
line = Y[names(Y) %in% c('acid', 'PRR17.492448')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR17.492448")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR17.492990.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR17.492990']
line = Y[names(Y) %in% c('acid', 'PRR17.492990')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR17.492990")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR20.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR20']
line = Y[names(Y) %in% c('acid', 'PRR20')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR20")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR25.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR25']
line = Y[names(Y) %in% c('acid', 'PRR25')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR25")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR25.491404.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR25.491404']
line = Y[names(Y) %in% c('acid', 'PRR25.491404')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR25.491404")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR25.491592.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR25.491592']
line = Y[names(Y) %in% c('acid', 'PRR25.491592')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR25.491592")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR25.492488.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR25.492488']
line = Y[names(Y) %in% c('acid', 'PRR25.492488')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR25.492488")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR25.492990.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR25.492990']
line = Y[names(Y) %in% c('acid', 'PRR25.492990')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR25.492990")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PRR30T.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PRR30T']
line = Y[names(Y) %in% c('acid', 'PRR30T')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PRR30T")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PYU.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PYU']
line = Y[names(Y) %in% c('acid', 'PYU')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PYU")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PYU.11002.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PYU.11002']
line = Y[names(Y) %in% c('acid', 'PYU.11002')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PYU.11002")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("PYU.11003.csv",header = T)#phenotypic data
y <- Y[names(Y)=='PYU.11003']
line = Y[names(Y) %in% c('acid', 'PYU.11003')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/PYU.11003")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("Northern.Stem.Canker.csv",header = T)#phenotypic data
y <- Y[names(Y)=='Northern.Stem.Canker']
line = Y[names(Y) %in% c('acid', 'Northern.Stem.Canker')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/Northern.Stem.Canker")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("NSC.491493.csv",header = T)#phenotypic data
y <- Y[names(Y)=='NSC.491493']
line = Y[names(Y) %in% c('acid', 'NSC.491493')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/NSC.491493")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SMV.G1.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SMV.G1']
line = Y[names(Y) %in% c('acid', 'SMV.G1')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SMV.G1")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SMV.G3.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SMV.G3']
line = Y[names(Y) %in% c('acid', 'SMV.G3')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SMV.G3")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SMV.G6.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SMV.G6']
line = Y[names(Y) %in% c('acid', 'SMV.G6')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SMV.G6")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SMV.G7.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SMV.G7']
line = Y[names(Y) %in% c('acid', 'SMV.G7')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SMV.G7")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SDS.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SDS']
line = Y[names(Y) %in% c('acid', 'SDS')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SDS")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SDS.493025.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SDS.493025']
line = Y[names(Y) %in% c('acid', 'SDS.493025')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SDS.493025")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SDS.493026.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SDS.493026']
line = Y[names(Y) %in% c('acid', 'SDS.493026')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SDS.493026")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN1.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN1']
line = Y[names(Y) %in% c('acid', 'SCN1')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN1")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
#########
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN1FI.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN1FI']
line = Y[names(Y) %in% c('acid', 'SCN1FI')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN1FI")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN1.494409.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN1.494409']
line = Y[names(Y) %in% c('acid', 'SCN1.494409')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN1.494409")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN1FI.494409.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN1FI.494409']
line = Y[names(Y) %in% c('acid', 'SCN1FI.494409')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN1FI.494409")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN2.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN2']
line = Y[names(Y) %in% c('acid', 'SCN2')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN2")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN2FI.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN2FI']
line = Y[names(Y) %in% c('acid', 'SCN2FI')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN2FI")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3']
line = Y[names(Y) %in% c('acid', 'SCN3')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3FI.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3FI']
line = Y[names(Y) %in% c('acid', 'SCN3FI')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3FI")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3code.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3code']
line = Y[names(Y) %in% c('acid', 'SCN3code')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3code")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3.83.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3.83']
line = Y[names(Y) %in% c('acid', 'SCN3.83')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3.83")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3code.83.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3code.83']
line = Y[names(Y) %in% c('acid', 'SCN3code.83')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3code.83")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3.491576.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3code.491576']
line = Y[names(Y) %in% c('acid', 'SCN3code.491576')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3code.491576")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3.491577.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3code.491577']
line = Y[names(Y) %in% c('acid', 'SCN3code.491577')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3code.491577")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3.492579.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3code.492579']
line = Y[names(Y) %in% c('acid', 'SCN3code.492579')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3code.492579")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3.491576.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3.491576']
line = Y[names(Y) %in% c('acid', 'SCN3.491576')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3.491576")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3.491577.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3.491577']
line = Y[names(Y) %in% c('acid', 'SCN3.491577')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3.491577")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN3.492579.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN3.492579']
line = Y[names(Y) %in% c('acid', 'SCN3.492579')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN3.492579")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN4.492579.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN4.492579']
line = Y[names(Y) %in% c('acid', 'SCN4.492579')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN4.492579")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
######
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
# Y <- read.csv("SCN4.492579.csv",header = T)#phenotypic data
# y <- Y[names(Y)=='SCN4code.492579']
# line = Y[names(Y) %in% c('acid', 'SCN4code.492579')]
# line1 <- line[(line$acid %in% GD$PI),]
# myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
# myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
# g= nrow(myY1)
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN4code.492579")
# myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN4.83.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN4code.83']
line = Y[names(Y) %in% c('acid', 'SCN4code.83')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN4code.83")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN4.83.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN4.83']
line = Y[names(Y) %in% c('acid', 'SCN4.83')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN4.83")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN4.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN4code']
line = Y[names(Y) %in% c('acid', 'SCN4code')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN4code")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN4.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN4']
line = Y[names(Y) %in% c('acid', 'SCN4')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN4")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN5code.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN5code']
line = Y[names(Y) %in% c('acid', 'SCN5code')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN5code")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN5.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN5']
line = Y[names(Y) %in% c('acid', 'SCN5')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN5")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN5FI.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN5FI']
line = Y[names(Y) %in% c('acid', 'SCN5FI')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN5FI")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05FI)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN5.83.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN5code.83']
line = Y[names(Y) %in% c('acid', 'SCN5code.83')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN5code.83")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN5.83.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN5.83']
line = Y[names(Y) %in% c('acid', 'SCN5.83')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN5.83")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
# Y <- read.csv("SCN5.491576.csv",header = T)#phenotypic data
# y <- Y[names(Y)=='SCN5code.491576']
# line = Y[names(Y) %in% c('acid', 'SCN5code.491576')]
# line1 <- line[(line$acid %in% GD$PI),]
# myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
# myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
# g= nrow(myY1)
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN5code.491576")
# myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
# Y <- read.csv("SCN5.491576.csv",header = T)#phenotypic data
# y <- Y[names(Y)=='SCN5.491576']
# line = Y[names(Y) %in% c('acid', 'SCN5.491576')]
# line1 <- line[(line$acid %in% GD$PI),]
# myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
# myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
# g= nrow(myY1)
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN5.491576")
# myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
# Y <- read.csv("SCN5.491577.csv",header = T)#phenotypic data
# y <- Y[names(Y)=='SCN5code.491577']
# line = Y[names(Y) %in% c('acid', 'SCN5code.491577')]
# line1 <- line[(line$acid %in% GD$PI),]
# myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
# myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
# g= nrow(myY1)
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN5code.491577")
# myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
# Y <- read.csv("SCN5.491577.csv",header = T)#phenotypic data
# y <- Y[names(Y)=='SCN5.491577']
# line = Y[names(Y) %in% c('acid', 'SCN5.491577')]
# line1 <- line[(line$acid %in% GD$PI),]
# myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
# myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
# g= nrow(myY1)
# setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN5.491577")
# myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN14.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN14code']
line = Y[names(Y) %in% c('acid', 'SCN14code')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN14code")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN14.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN14']
line = Y[names(Y) %in% c('acid', 'SCN14')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN14")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
###############
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN14.491576.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN14code.491576']
line = Y[names(Y) %in% c('acid', 'SCN14code.491576')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN14code.491576")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN14.491576.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN14.491576']
line = Y[names(Y) %in% c('acid', 'SCN14.491576')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN14.491576")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN14.491577.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN14code.491577']
line = Y[names(Y) %in% c('acid', 'SCN14code.491577')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN14code.491577")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("SCN14.491577.csv",header = T)#phenotypic data
y <- Y[names(Y)=='SCN14.491577']
line = Y[names(Y) %in% c('acid', 'SCN14.491577')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/SCN14.491577")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("saltreact.csv",header = T)#phenotypic data
y <- Y[names(Y)=='saltreact']
line = Y[names(Y) %in% c('acid', 'saltreact')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/saltreact")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("CHLOROSIS.csv",header = T)#phenotypic data
y <- Y[names(Y)=='CHLOROSIS']
line = Y[names(Y) %in% c('acid', 'CHLOROSIS')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/CHLOROSIS")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx")
Y <- read.csv("RUST.csv",header = T)#phenotypic data
y <- Y[names(Y)=='RUST']
line = Y[names(Y) %in% c('acid', 'RUST')]
line1 <- line[(line$acid %in% GD$PI),]
myY1 <- line[which(line$acid %in% GD$PI),]#select lines have been genotyped
myGD <- GD[which(GD$PI %in% myY1$acid),]#selection lines have observation available
g= nrow(myY1)
setwd("/Users/jmshook/Desktop/Beagle GWAS/Dx/RUST")
myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = myGM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
