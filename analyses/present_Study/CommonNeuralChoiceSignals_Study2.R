library(R.matlab)
library(languageR)
library(lme4)
library(MASS)
library(ggplot2) #cookbook for R/ Graphs
library(memisc)
library(reshape)
library(reshape2) #melt and cast -> restructure and aggregate data
library(data.table)
library(psych)
library(doBy)
library(heplots)
library(plyr) #necessary for ddply
library(matrixStats) 
library(foreign) 
library(Hmisc)
library(lmerTest)
library (stringr)
library(gdata)
library(Rmisc)
library(effects)
library(RColorBrewer)
library(sjPlot)
library(buildmer)

# get to the right place:
fileLoc <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(fileLoc) # go to script location first
setwd("../..") 

basepath <- getwd()


didLmerConverge = function(lmerModel){
  relativeMaxGradient=signif(max(abs(with(lmerModel@optinfo$derivs,solve(Hessian,gradient)))),3)
  if (relativeMaxGradient < 0.001) {
    cat(sprintf("\tThe relative maximum gradient of %s is less than our 0.001 criterion.\n\tYou can safely ignore any warnings about a claimed convergence failure.\n\n", relativeMaxGradient))
  }
  else {
    cat(sprintf("The relative maximum gradient of %s exceeds our 0.001 criterion.\nThis looks like a real convergence failure; maybe try simplifying your model?\n\n", relativeMaxGradient))
  }
}

datapath <- '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASS_EEG/Data/Export/'

### load data ###
# choice data
input_file = paste0(datapath, 'allSubDataTable505.xls')
a1 = read.xls(input_file)

##### load in EEG data ################

# 
input_file = input_file = paste0(datapath, 'CPPPisauro-700-200.mat')
CPP = readMat(input_file)

CPP = CPP$CPPPisauro.700.200 #--> da muss man die Variable angeben, als die man in Matlab gespeichert hat
CPP = as.data.frame(CPP)
colnames(CPP)[c( 1, 2, 3, 4, 5, 
                      6, 7, 8, 9, 10, 
                      11, 12, 13, 14, 15, 
                      16, 17, 18, 19, 20, 
                      21, 22, 23, 24, 25, 
                      26, 27, 28, 29, 30,
                      31, 32, 33, 34, 35,
                      36, 37, 38, 39, 40,
                      41, 42, 43, 44, 45,
                      46, 47, 48, 49, 50,
                      51, 52, 53, 54, 55,
                      56, 57, 58, 59, 60,
                      61, 62, 63, 64, 65)]=c('FP1', 'FPz', 'FP2', 'AF3', 
                                             'AFz', 'AF4', 'F7', 'F5', 'F3', 
                                             'F1','Fz', 'F2', 'F4', 'F6', 
                                             'F8', 'FT9', 'FT7', 'FC5', 'FC3', 
                                             'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 
                                             'FT8', 'FT10', 'T7', 'C5', 'C3', 
                                             'C1','C2', 'C4', 'C6', 'T8', 
                                             'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 
                                             'CPz','CP2', 'CP4', 'CP6', 'TP8', 
                                             'TP10', 'P7', 'P5', 'P3', 'P1', 
                                             'Pz', 'P2', 'P4', 'P6', 'P8', 
                                             'PO3', 'POz', 'PO4', 'O1', 'Oz', 
                                             'O2', 'LO1', 'IO1', 'IO2', 'LO2', 
                                             'Cz')

# convert variables and set contrasts
a1$SubNum <- factor(a1$SubNum)
a1$sVD <- scale(a1$VD/10, scale=FALSE, center=TRUE)
a1$cRT <- scale(a1$RT, scale=FALSE, center=TRUE)


a1s <- a1[!a1$SubNum==5000 & !a1$SubNum== 5022, ]
a1s$SubNum <- factor(a1s$SubNum)

a1s$CPPPz <- CPP$Pz

# lower activity for longer RTs
print (summary(CPPRTmod <- lmer(CPPPz~ cRT+(cRT|SubNum), a1s[!a1s$Choice==-1,], 
                                REML=FALSE)))

sv1_max <- svd(getME(CPPRTmod, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

# higher activity for larger VD
print (summary(CPPRTmodsv <- lmer(CPPPz~ sVD+(sVD|SubNum), a1s[!a1s$Choice==-1,], 
                                  REML=FALSE)))

sv1_max <- svd(getME(CPPRTmod, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1) 

tab_model(CPPRTmod, CPPRTmodsv,show.stat = TRUE, col.order= c("est", "ci", "stat", "p"))


