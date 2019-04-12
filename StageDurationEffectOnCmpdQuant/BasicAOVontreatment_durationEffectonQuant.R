#Load libraries
library("MetaboAnalystR")
library("broom")
library("tidyr")
library("dplyr")
library("data.table")
library(effsize)
library(ggplot2)
library(readxl)


###load metabs data
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "~/Desktop/UCSD-SALK-UW/NOAA/metabolomics/data/for metaboanalyst/339117_data4MetaboAnalyst.csv", "rowu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-RemoveMissingPercent(mSet, percent=0.5)
mSet<-ImputeVar(mSet, method="min")
dataset <- data.frame(mSet$dataSet)
ncmpds <- length(grep("orig.", colnames(dataset)))-1

data1 <- data.frame(mSet$dataSet)
data1 <- data1[,c(6,grep("^procr.", colnames(data1)))]
data1$sample <- rownames(data1)
data1 <- data1[,c(ncol(data1),2:ncol(data1)-1)]
data1$pH <- ifelse(grepl("igh", substr(data1$cls,1,4)), "High", "Low")
data1$DO <- ifelse(grepl("igh",substr(data1$cls,6,10)), "High", "Low")
data1 <- data1[,c(1,2,654,655,3:653)]
colnames(data1) <- c("mx_sample","treatment", "pH", "DO", as.list(substring(colnames(data1[,grep("procr.",colnames(data1))]),7)))

STACKED_data <- tidyr::gather(data1, "analyte", "quant", 5:ncol(data1))
#STACKED_data$analyte <- for (i in 1:nrow(STACKED_data))
ifelse(!grepl("[A-Za-z]", STACKED_data$analyte[i]), paste("unk_",STACKED_data$analyte[i]),STACKED_data$analyte[i])



#count compounds after removing those with noise abundance
length(unique(STACKED_data$analyte))
#651

#ggplot(STACKED_data[grep("pipe", STACKED_data$analyte),], aes(treatment, quant)) + geom_violin(aes(fill= treatment)) + 
#geom_boxplot(aes(fill= treatment), width = 0.1) + theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 4)) +
#facet_wrap(~analyte, scales = "free")


###load lipids data

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "~/Desktop/UCSD-SALK-UW/NOAA/metabolomics/ARTool_analysis/GenLipidsAll.csv", "rowu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-RemoveMissingPercent(mSet, percent=0.5)
mSet<-ImputeVar(mSet, method="min")
dataset <- data.frame(mSet$dataSet)
ncmpds <- length(grep("orig.", colnames(dataset)))-1

data2 <- data.frame(mSet$dataSet)
data2 <- data2[,c(6,grep("^procr.", colnames(data2)))]
data2$sample <- rownames(data2)
data2 <- data2[,c(ncol(data2),2:ncol(data2)-1)]
data2$pH <- ifelse(grepl("igh", substr(data2$cls,1,4)), "High", "Low")
data2$DO <- ifelse(grepl("igh",substr(data2$cls,6,10)), "High", "Low")
ncol(data2)
data2 <- data2[,c(1,2,3099, 3100, 3:3098)]
colnames(data2) <- c("mx_sample","treatment", "pH", "DO", as.list(substring(colnames(data2[,grep("procr.",colnames(data2))]),7)))

STACKED_data2 <- tidyr::gather(data2, "analyte", "quant", 5:ncol(data2))

#count compounds after removing those with noise abundance
length(unique(STACKED_data2$analyte))
#3096



################
###STATISTICS###
################

#Check if there are effects from stage duration or days on treatment

duration <- read.csv("~/Desktop/UCSD-SALK-UW/NOAA/metabolomics/data/crab_duration_info.csv", stringsAsFactors = FALSE)

library(readxl)
sampleidinfo <- read_xlsx("~/Desktop/UCSD-SALK-UW/NOAA/metabolomics/data/SampleID_CrabID_info.xlsx")

colnames(STACKED_data)[1] <- "Sample.Number"

STACKED_data <- merge(STACKED_data,sampleidinfo, by = "Sample.Number")
STACKED_data <- merge(STACKED_data,duration, by = "crabID")

colnames(STACKED_data2)[1] <- "mx.sample"
STACKED_data2 <- merge(STACKED_data2,sampleidinfo, by = "mx.sample")
STACKED_data2 <- merge(STACKED_data2,duration, by = "crabID")


#test for M duration effect on quant
fit <- aov(quant ~ M, data = STACKED_data)
summary(fit)
#     Df    Sum Sq   Mean Sq F value Pr(>F)
#M               1 4.164e+08 4.164e+08   0.001  0.972
#Residuals   39058 1.308e+16 3.349e+11  


#test for J duration effect on quant
fit <- aov(quant ~ J1, data = STACKED_data)
summary(fit)
#              Df    Sum Sq   Mean Sq F value Pr(>F)
#J1              1 1.244e+10 1.244e+10   0.037  0.847
#Residuals   39058 1.308e+16 3.349e+11 

#test for MtoJ2 duration effect on quant
fit <- aov(quant ~ MtoJ2, data = STACKED_data)
summary(fit)
#Df    Sum Sq   Mean Sq F value Pr(>F)
#MtoJ2           1 3.407e+10 3.407e+10   0.102   0.75
#Residuals   39058 1.308e+16 3.349e+11   

#test for startTOfreeze duration effect on quant
fit <- aov(quant ~ startTOfreeze, data = STACKED_data)
summary(fit)

#Df    Sum Sq   Mean Sq F value Pr(>F)
#startTOfreeze     1 3.407e+10 3.407e+10   0.102   0.75
#Residuals     39058 1.308e+16 3.349e+11          

#test for treatment effect on quant
fit <- aov(quant ~ treatment.x, data = STACKED_data)
summary(fit)

#Df    Sum Sq   Mean Sq F value Pr(>F)
#treatment.x     3 1.329e+11 4.429e+10   0.132  0.941
#Residuals   39056 1.308e+16 3.349e+11 

### for lipids


#test for M duration effect on quant
fit <- aov(quant ~ M, data = STACKED_data2)
summary(fit)
#Df    Sum Sq   Mean Sq F value Pr(>F)
#M                1 3.349e+10 3.349e+10   1.699  0.192
#Residuals   185758 3.662e+15 1.971e+10     

#test for J duration effect on quant
fit <- aov(quant ~ J1, data = STACKED_data2)
summary(fit)
#Df    Sum Sq   Mean Sq F value Pr(>F)  
#J1               1 1.272e+11 1.272e+11    6.45 0.0111 *
#  Residuals   185758 3.662e+15 1.971e+10     

#test for MtoJ2 duration effect on quant
fit <- aov(quant ~ MtoJ2, data = STACKED_data2)
summary(fit)
#Df    Sum Sq   Mean Sq F value  Pr(>F)   
#MtoJ2            1 1.562e+11 1.562e+11   7.924 0.00488 **
#  Residuals   185758 3.662e+15 1.971e+10   

#test for startTOfreeze duration (days on experiment) effect on quant
fit <- aov(quant ~ startTOfreeze, data = STACKED_data2)
summary(fit)
#Df    Sum Sq   Mean Sq F value  Pr(>F)   
#startTOfreeze      1 1.562e+11 1.562e+11   7.924 0.00488 **
#  Residuals     185758 3.662e+15 1.971e+10                   


#test for treatment effect on quant
fit <- aov(quant ~ cls, data = STACKED_data2)
summary(fit)
#Df    Sum Sq   Mean Sq F value  Pr(>F)    
#cls              3 6.259e+11 2.086e+11   10.59 5.9e-07 ***
#  Residuals   185756 3.661e+15 1.971e+10 


