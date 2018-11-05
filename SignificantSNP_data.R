#Install packages
install.packages(c("fields","RColorBrewer","mapplots"))
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
library(dplyr)
install.packages('adegenet')
library(adegenet)
install.packages("lmap")

#Parallel computing
library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl<-makeCluster(4,type="SOCK")
registerDoSNOW(cl)


##########################################################################################
########################################################################
##########################################################################################

#Import new df
df<-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data Stuff/GENO_GWAS_Merged_299.csv", sep=",", header=T, check.names = FALSE)
df<-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data Stuff/GWAS_Merged_299.csv", sep=",", header=T, check.names = FALSE)

df[1:5,1:15]
df1 <- df

str(df)
#These are significant SNPs identified by the GWAS taken from John 
#search for them and copy their nucleotide column
ss715579841 <- df[grepl("ss715579841", df$ID),]
ss715583518 <- df[grepl("ss715583518", df$ID),]
ss715585038 <- df[grepl("ss715585038", df$ID),]
ss715586566 <- df[grepl("ss715586566", df$ID),]
ss715590058 <- df[grepl("ss715590058", df$ID),]
ss715598986 <- df[grepl("ss715598986", df$ID),]
ss715603553 <- df[grepl("ss715603553", df$ID),]
ss715606432 <- df[grepl("ss715606432", df$ID),]
ss715606449 <- df[grepl("ss715606449", df$ID),]
ss715611836 <- df[grepl("ss715611836", df$ID),]
ss715611871 <- df[grepl("ss715611871", df$ID),]
ss715612051 <- df[grepl("ss715612051", df$ID),]
ss715612130 <- df[grepl("ss715612130", df$ID),]
ss715612131 <- df[grepl("ss715612131", df$ID),]
ss715612159 <- df[grepl("ss715612159", df$ID),]
ss715612222 <- df[grepl("ss715612222", df$ID),]
ss715617353 <- df[grepl("ss715617353", df$ID),]
ss715617534 <- df[grepl("ss715617534", df$ID),]
ss715617605 <- df[grepl("ss715617605", df$ID),]
ss715618282 <- df[grepl("ss715618282", df$ID),]
ss715618290 <- df[grepl("ss715618290", df$ID),]
ss715625627 <- df[grepl("ss715625627", df$ID),]
ss715623707 <- df[grepl("ss715623707", df$ID),]
ss715629364 <- df[grepl("ss715629364", df$ID),]
ss715625627 <- df[grepl("ss715625627", df$ID),]
ss715623707 <- df[grepl("ss715623707", df$ID),]
ss715631578 <- df[grepl("ss715631578", df$ID),]
ss715631926 <- df[grepl("ss715631926", df$ID),]
ss715631927 <- df[grepl("ss715631927", df$ID),]
ss715636162 <- df[grepl("ss715636162", df$ID),]
ss715636682 <- df[grepl("ss715636682", df$ID),]
ss715636812 <- df[grepl("ss715636812", df$ID),]

#bind them together
SNPs <- rbind(ss715579841,ss715583518,ss715585038,ss715586566,ss715590058,ss715598986,ss715603553,ss715606432,ss715606449,ss715611836,ss715611871,ss715612051,ss715612130,ss715612131,ss715612159,ss715612222,ss715617353,ss715617534,ss715617605,ss715618282,ss715618290,ss715625627,ss715623707,ss715629364,ss715631578,ss715631926,ss715631927,ss715636162,ss715636682,ss715636812)
SNPs <- SNPs[,c(1,5:length(SNPs))]
SNPs[1:10,1:20]

#Transcribe the dataframe
tSNPs <- t(SNPs)

write.csv(SNPs,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/GWAS/significantSNPs_ATGCdata.csv", row.names = T)
write.csv(tSNPs,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/GWAS/significantSNPs_ATGCdata_transposed_1.csv", row.names = T)

SNPsEntry<-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/GWAS/significantSNPs_ATGCdata_transposed.csv", sep=",", header=T, check.names = FALSE)
BigData <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Data/BlueSteelData_Sept25.csv", sep=",", header=T, check.names = FALSE)

#remove the PI number names
SNPsEntry <- SNPsEntry[1:300,2:32]

#merge the two df's together like VLOOKUP
Merged_df <- left_join(BigData, SNPsEntry, by="Entry")
Merged_df[1:10,1:88]

write.csv(Merged_df,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/GWAS/BigData_significantSNPs_ATGCdata.csv", row.names = T)

SNPsEntry[1:10,1:10]
#entry numbers in order
Entries <- (SNPsEntry[1:292,1])

##########################################################################################
########################################################################
##########################################################################################

#Import new df
GD<-read.table("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data Stuff/GWAS_GD.txt", sep="\t", header=T)
GD[1:5,1:15]

t_GD <- t(GD)
t_GD[1:300,1:5]
str(t_GD)
t_GD[,1]

str(GD)
#These are significant SNPs identified by the GWAS taken from John 
#search for them and copy their nucleotide column

ss715579841 <- GD[grepl("ss715579841", GD$ID),]

t_GD

ss715583518 <- GD$ss715583518
ss715585038 <- GD$ss715585038
ss715586566 <- GD$ss715586566
ss715590058 <- GD$ss715590058
ss715598986 <- GD$ss715598986
ss715603553 <- GD$ss715603553
ss715606432 <- GD$ss715606432
ss715606449 <- GD$ss715606449
ss715611836 <- GD$ss715611836
ss715611871 <- GD$ss715611871
ss715612051 <- GD$ss715612051
ss715612130 <- GD$ss715612130
ss715612131 <- GD$ss715612131
ss715612159 <- GD$ss715612159
ss715612222 <- GD$ss715612222
ss715617353 <- GD$ss715617353
ss715617534 <- GD$ss715617534
ss715617605 <- GD$ss715617605
ss715618282 <- GD$ss715618282
ss715618290 <- GD$ss715618290
ss715625627 <- GD$ss715625627
ss715623707 <- GD$ss715623707
ss715629364 <- GD$ss715629364
ss715625627 <- GD$ss715625627
ss715623707 <- GD$ss715623707
ss715631578 <- GD$ss715631578
ss715631926 <- GD$ss715631926
ss715631927 <- GD$ss715631927
ss715636162 <- GD$ss715636162
ss715636682 <- GD$ss715636682
ss715636812 <- GD$ss715636812

onlySNPs <- cbind(ss715579841,ss715583518,ss715585038,ss715586566,ss715590058,ss715598986,ss715603553,ss715606432,ss715606449,ss715611836,ss715611871,ss715612051,ss715612130,ss715612131,ss715612159,ss715612222,ss715617353,ss715617534,ss715617605,ss715618282,ss715618290,ss715625627,ss715623707,ss715629364,ss715631578,ss715631926,ss715631927,ss715636162,ss715636682,ss715636812)

#Convert GD into matrix form 
geno = as.matrix(onlySNPs[,2:31])

#Compute genetic distance using Nei Distance
gdist = Gdist(geno, method = 1)

neinan <- nei.dist(geno)
