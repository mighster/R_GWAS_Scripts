#Package Control
install.packages('dplyr')
install.packages('Rcpp', dependencies = TRUE)
install.packages('gplots')
install.packages('LDheatmap')
install.packages('genetics')
install.packages('scatterplot3d')
install.packages('synbreed')
install.packages('compiler')
library(dplyr)
library(compiler)
library(synbreed)
library(Rcpp)
library(gplots)
library(LDheatmap)
library(genetics)
library(scatterplot3d)
source("https://bioconductor.org/biocLite.R")
biocLite("Biobase")
biocLite("multtest")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
source("http://www.zzlab.net/GAPIT/emma.txt")
library(multtest)
setwd('C:/Users/jmshook/Desktop/KGF/PhenotypeFile/')
results_dir <- ("C:/Users/jmshook/Desktop/KGF/Results")

Y <- read.csv('C:/Users/jmshook/Desktop/KGF/PhenotypeFile/KGF_AdjustedBLUPsAllDays_thinned_Oct19.csv')
Y <- read.csv('C:/Users/jmshook/Desktop/KGF/PhenotypeFile/composite_traits.csv')

#remove BAD genotypes PI594438 (273) PI495832 (202) PI467333A (382) PI473899 (390) PI506528 (204) PI445845 (367)
Y <- Y %>%
  filter(!((Entry == "273")|(Entry == "202")|(Entry  == "182")|( df$Entry  == "190")|( df$Entry  == "204")|( df$Entry  == "267" )|( df$Entry  == "9" )|( df$Entry  == "98" )))

GM  <- read.table('C:/Users/jmshook/Desktop/KGF/PhenotypeFile/298_nmis0.1_maf.05_GM.txt',head=TRUE)
GD <- read.table('C:/Users/jmshook/Desktop/KGF/PhenotypeFile/298_nmis0.1_maf.05_GD.txt',head=TRUE)
y <- "TRL"

colnames(Y)

df_merged <- inner_join(Y, GD, by = 'PI')
df_merged[1:10,1:10]
str(df_merged)
unique(df_merged$PI)
#DATA<-df_merged[-c(21,23,25),-c(2:127)]
Y[1:10,5:31]


Y$Shape<-as.factor(Y$Shape)

Y

str(Y)

x=120
colnames(Y)[119]
x=6
for (x in 5:20){
  y <- colnames(Y)[x] # obtain the name of the trait. All files have phenotypes in second column (starting from 1, not 0)
  line = Y[names(Y) %in% c('PI', y)] # obtain list of acession identifiers found in the phenotype file
  myY <- line[which(line$PI %in% GD$PI),] # select lines that have been genotyped
  myGD <- as.data.frame(GD[which(line$PI %in% GD$PI),]) #reduce GD file to genotypes specific to this test
  g= nrow(myY) # number of accessions in current GWAS test... to be used as an argument in the "GAPIT" function
  dir.create(paste(results_dir,y)) # create a folder within the "Results" directory to store the output of the GWAS test
  setwd(paste(results_dir,y)) # change the directory to output the GAPIT results into the folder just created
  myGAPIT <- GAPIT(Y = myY, GD = GD, GM = GM, PCA.total =3, group.from =g, group.to = g,SNP.fraction = 1, SNP.MAF = 0.05)
}

setwd("/Users/jmshook/Desktop/KGF/Oct31Results")
files <- list.files(pattern = "\\.GWAS.Results.csv$", recursive = TRUE)
data2=lapply(files, read.table, header=T, sep=",")
for (i in 1:length(data2)){data2[[i]]<-cbind(data2[[i]],files[i])}
data_rbind <- do.call("rbind", data2) 
colnames(data_rbind)[c(1:10)]<-c("SNP", "Chromosome", "Position", "P.value", "maf", "nobs", "Rsquare.of.Model.without.SNP", "Rsquare.of.Model.with.SNP", "`FDR_Adjusted_P-values`","file.name")
signif <- data_rbind[(data_rbind[,9]<0.05),]
write.csv(signif,'SignificantSNPsForKevin_10312018.csv')

write.csv(myGD, 'C:/Users/jmshook/Desktop/KGF/PhenotypeFile/myGD.csv')

head(data2)
