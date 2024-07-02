library(R.matlab)
library(languageR)
library(lme4)
library(MASS)
library(ggplot2) #cookbook for R/ Graphs
library(hexbin)
library(memisc)
library(reshape)
library(reshape2) #melt and cast -> restructure and aggregate data
library(data.table)
library(coin) #for permutation tests
library(psych)
library(doBy)
library(heplots)
library(plyr) #necessary for ddply
library(matrixStats) 
library(foreign) 
library(Hmisc)
library(stringr)
library(gridExtra)
library(grid)
library(gdata)
library(effects)
library(ggExtra)
library(jtools)
library(sjPlot)


setwd("~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Analyses/R")


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


didLmerConverge = function(lmerModel){
  relativeMaxGradient=signif(max(abs(with(lmerModel@optinfo$derivs,solve(Hessian,gradient)))),3)
  if (relativeMaxGradient < 0.001) {
    cat(sprintf("\tThe relative maximum gradient of %s is less than our 0.001 criterion.\n\tYou can safely ignore any warnings about a claimed convergence failure.\n\n", relativeMaxGradient))
  }
  else {
    cat(sprintf("The relative maximum gradient of %s exceeds our 0.001 criterion.\nThis looks like a real convergence failure; maybe try simplifying your model?\n\n", relativeMaxGradient))
  }
}

# for behavioral analyses
input_file = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/Export/allSubDataTable.xls'
a1 = read.xls(input_file)

atmp <- a1
# for CPP analyses
input_file = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/Export/Pzrout.xls'
Pzrdat = read.xls(input_file, blank.lines.skip=FALSE) # without the blank lines thing, the matrices don't match!

# exclude 1 participant with incomplete data
Pzrdat <- Pzrdat[!a1$SubNum==5028,]
a1 <- a1[!a1$SubNum==5028,]

# for plotting relationship between PCs and otehr variables
input_file = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/Export/PCscores_new.xls'
PCdat = read.xls(input_file)



## data and variable prep
str(a1)

a1$SubNum <- as.factor(a1$SubNum)
a1$isrightchosen <- a1$Choice1 -1
a1$RightvsLeft <- (a1$Valuer-a1$Valuel)/10
a1$sVD <- scale(a1$MaxvMinBid/10, scale=FALSE, center=TRUE)
a1$savV <- scale(a1$AvBid/10, scale=FALSE, center=TRUE)
a1$sAnx <- scale(a1$Anxious, scale=FALSE, center=TRUE) 
a1$sLiking <- scale(a1$Liking, scale=FALSE, center=TRUE) 
a1$sConf <- scale(a1$Confident, scale=FALSE, center=TRUE) 
a1$iszero <- rep(0, length(a1$SubNum))
a1$RT <- a1$RTeval1*1000
a1$cRT <- scale(a1$RTeval1, scale=FALSE, center=TRUE)

acf <-  a1[!is.na(a1$Liking),]



# generate Fig 2B top
pVDIST <- ggplot(data=a1, aes(x=MaxvMinBid, y= AvBid))+ geom_point(shape="O", color="#C0C0C0")+theme_bw(12)+
  xlab("value difference") + ylab("average set value") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="bottom")

pdf("~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Docu/Figures/value_dist_w_hist.pdf", width = 5, height = 5)#, units = 'cm', res = 200, compression = 'lzw'
ggMarginal(pVDIST, type = "histogram", fill="#C0C0C0", xparams = list(size=0.1), yparams = list(size=0.1))
dev.off()



### Simple behavior (also to generate fig 2B, center)
# choices
print (summary(ChoiceMod <- glmer(isrightchosen~ RightvsLeft+ savV+(RightvsLeft|SubNum), data = acf[!acf$Choice1==-1,], 
                                  family = binomial))) 

sv1_max <- svd(getME(ChoiceMod, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


eff_df <- effect("RightvsLeft", ChoiceMod,  xlevels=list(RightvsLeft=seq(min(acf$RightvsLeft, na.rm=TRUE), max(acf$RightvsLeft, na.rm=TRUE), 0.01)))
contmain <- as.data.frame(eff_df)
pChoiceAcc <- ggplot(data=contmain, aes(x=RightvsLeft, y=fit)) +theme_bw(12)+ ylim(0,1)+ geom_ribbon(data=contmain, aes(x=RightvsLeft, max = upper, min = lower),alpha=0.3, inherit.aes = FALSE)+ geom_line(color="#000000", size=1)+
  xlab("Right - Left Value") + ylab("Probability Right Chosen") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+geom_vline(xintercept = 0, linetype=2, size=0.2) + geom_hline(yintercept = 0.5, linetype=2, size=0.2)+ theme(legend.position="bottom")

#RT
print (summary(RTmod0 <- lmer(RT~ sVD+ savV+(sVD+ savV|SubNum), acf[!acf$Choice1==-1,], 
                              REML=FALSE)))

sv1_max <- svd(getME(RTmod0, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

# generate Table S1
tab_model(ChoiceMod,RTmod0, transform = NULL,show.stat = TRUE, col.order= c("est", "ci", "stat", "p"))

eff_df <- effect("sVD", RTmod0,  xlevels=list(sVD=seq(min(acf$sVD, na.rm=TRUE), max(acf$sVD, na.rm=TRUE), 0.1)))
contmain <- as.data.frame(eff_df)
contmain$sVD <- contmain$sVD - min(contmain$sVD, na.rm=TRUE)

contmainC <- contmain
contmainC$ValType <- rep("VD", length(contmain$fit))
contmainC <- rename(contmainC, c("sVD"="Value"))


eff_df <- effect("savV", RTmod0,  xlevels=list(savV=seq(min(acf$savV, na.rm=TRUE), max(acf$savV, na.rm=TRUE), 0.1)))
contmain <- as.data.frame(eff_df)
contmain$savV <- contmain$savV - min(contmain$savV, na.rm=TRUE)

contmainuC <- contmain
contmainuC$ValType <- rep("OV", length(contmain$fit))
contmainuC <- rename(contmainuC, c("savV"="Value"))

contmainall <- rbind(contmainC, contmainuC)

contmainall$ValType <- ordered(contmainall$ValType, levels=c("VD", "OV"))


pRT <- ggplot(data=contmainall, aes(x=Value, y=fit, linetype= ValType))+theme_bw(12)+ geom_ribbon(data=contmainall, aes(x=Value, max = upper, min = lower, fill = ValType),alpha=0.2, inherit.aes = FALSE) + geom_line()+
  xlab("Value") + ylab("RT")+scale_fill_manual(values=c("#000000","#C0C0C0C0")) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position=c(0.8, 0.8))


pdf("~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Docu/Figures/ChoiceRT.pdf", width = 4, height = 2)
multiplot(pChoiceAcc,pRT, cols=2)
dev.off()



# generate Fig2B bottom
##### relationships Anxiety, Confidence and Liking and OV, VD, just for plot ################

print (summary(AnxMod0 <- lmer(Anxious ~ sVD+ savV + I(savV*savV)+(sVD+ savV |SubNum), a1[!a1$Choice1==-1,], 
                               REML=FALSE)))

sv1_max <- svd(getME(AnxMod0, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


eff_df <- Effect("sVD", AnxMod0,  xlevels=list(sVD=seq(min(a1$sVD, na.rm=TRUE), max(a1$sVD, na.rm=TRUE), 0.01)))
contmain <- as.data.frame(eff_df)
contmain$sVD <- contmain$sVD - min(contmain$sVD, na.rm=TRUE)

contmainC <- contmain
contmainC$ValType <- rep("Choice", length(contmain$fit))
contmainC <- rename(contmainC, c("sVD"="Value"))


eff_df <- Effect("savV", AnxMod0,  xlevels=list(savV=seq(min(a1$savV, na.rm=TRUE), max(a1$savV, na.rm=TRUE), 0.01)))
contmain <- as.data.frame(eff_df)
contmain$savV <- contmain$savV - min(contmain$savV, na.rm=TRUE)


contmainuC <- contmain
contmainuC$ValType <- rep("Valuation", length(contmain$fit))
contmainuC <- rename(contmainuC, c("savV"="Value"))

contmainall <- rbind(contmainC, contmainuC)

pPCAnx <- ggplot(data=contmainall, aes(x=Value, y=fit, color=ValType)) + geom_line()+scale_colour_manual(name="Option",values=c("#005F96","#C55A11"))+theme_bw(12)+ geom_ribbon(data=contmainall, aes(x=Value, max = upper, min = lower, fill = ValType),alpha=0.3, inherit.aes = FALSE)+scale_y_continuous(limits=c(-1,5.5), breaks=c(1,5))+
  xlab("Value") + ylab("Anxiety")+scale_fill_manual(name="Option",values=c("#005F96","#C55A11")) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")+coord_cartesian(ylim = c(1, 5)) 

### Confidence

print (summary(ConfMod0 <- lmer(Confident ~ sVD+ savV + I(savV*savV)+(sVD+ savV |SubNum), a1[!a1$Choice1==-1,], 
                               REML=FALSE)))

sv1_max <- svd(getME(AnxMod0, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


eff_df <- Effect("sVD", ConfMod0,  xlevels=list(sVD=seq(min(a1$sVD, na.rm=TRUE), max(a1$sVD, na.rm=TRUE), 0.01)))
contmain <- as.data.frame(eff_df)
contmain$sVD <- contmain$sVD - min(contmain$sVD, na.rm=TRUE)

contmainC <- contmain
contmainC$ValType <- rep("Choice", length(contmain$fit))
contmainC <- rename(contmainC, c("sVD"="Value"))


eff_df <- Effect("savV", ConfMod0,  xlevels=list(savV=seq(min(a1$savV, na.rm=TRUE), max(a1$savV, na.rm=TRUE), 0.01)))
contmain <- as.data.frame(eff_df)
contmain$savV <- contmain$savV - min(contmain$savV, na.rm=TRUE)

contmainuC <- contmain
contmainuC$ValType <- rep("Valuation", length(contmain$fit))
contmainuC <- rename(contmainuC, c("savV"="Value"))

contmainall <- rbind(contmainC, contmainuC)
pPCConf <- ggplot(data=contmainall, aes(x=Value, y=fit, color=ValType)) + geom_line()+scale_colour_manual(name="Option",values=c("#005F96","#C55A11"))+theme_bw(12)+ geom_ribbon(data=contmainall, aes(x=Value, max = upper, min = lower, fill = ValType),alpha=0.3, inherit.aes = FALSE)+scale_y_continuous( breaks=c(1,5))+#limits=c(-1,5.5),
  xlab("Value") + ylab("Confidence")+scale_fill_manual(name="Option",values=c("#005F96","#C55A11")) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")+coord_cartesian(ylim = c(1, 5)) 


## Liking
print (summary(LikMod0 <- lmer(Liking ~ sVD+ savV + I(savV*savV)+(sVD+ savV |SubNum), a1[!a1$Choice1==-1,], 
                               REML=FALSE)))

sv1_max <- svd(getME(LikMod0, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


eff_df <- Effect("sVD", LikMod0,  xlevels=list(sVD=seq(min(a1$sVD, na.rm=TRUE), max(a1$sVD, na.rm=TRUE), 0.01)))
contmain <- as.data.frame(eff_df)
contmain$sVD <- contmain$sVD - min(contmain$sVD, na.rm=TRUE)

contmainC <- contmain
contmainC$ValType <- rep("Choice", length(contmain$fit))
contmainC <- rename(contmainC, c("sVD"="Value"))


eff_df <- Effect("savV", LikMod0,  xlevels=list(savV=seq(min(a1$savV, na.rm=TRUE), max(a1$savV, na.rm=TRUE), 0.01)))
contmain <- as.data.frame(eff_df)
contmain$savV <- contmain$savV - min(contmain$savV, na.rm=TRUE)

contmainuC <- contmain
contmainuC$ValType <- rep("Valuation", length(contmain$fit))
contmainuC <- rename(contmainuC, c("savV"="Value"))

contmainall <- rbind(contmainC, contmainuC)


pPCLik <- ggplot(data=contmainall, aes(x=Value, y=fit, color=ValType)) + geom_line()+scale_colour_manual(name="Option",values=c("#005F96","#C55A11"))+theme_bw(12)+ geom_ribbon(data=contmainall, aes(x=Value, max = upper, min = lower, fill = ValType),alpha=0.3, inherit.aes = FALSE)+scale_y_continuous(limits=c(-1,5.5), breaks=c(1,5))+
  xlab("Value") + ylab("Liking")+scale_fill_manual(name="Option",values=c("#005F96","#C55A11")) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")+coord_cartesian(ylim = c(1, 5)) 

pdf("~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Docu/Figures/Anx_Conf_Lik_VD_OV.pdf", width = 5, height = 4)#, units = 'cm', res = 200, compression = 'lzw'
multiplot(pPCAnx, pPCConf, pPCLik, cols = 3)
dev.off()


### CPP analysis

input_file = '~/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/BASB_EEG/Data/Export/CPPPisauro700200.mat' # where is your file?
CPPPisauro = readMat(input_file)

CPPPisauro = CPPPisauro$CPPPisauro700200 #--> specify variable as saved in Matlab
CPPPisauro = as.data.frame(CPPPisauro) # convert to data frame
# specify column names
colnames(CPPPisauro)[c( 1, 2, 3, 4, 5, 
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
                      61)]=c('FP1', 'FPz', 'FP2', 'AF3', 
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
                                             'O2',  
                                             'Cz')

CPPPisauro <- CPPPisauro[!atmp$SubNum==5028,]

a1$CPPPz <- CPPPisauro$Pz


# lower activity for longer RTs
print (summary(CPPRTmod <- lmer(CPPPz~ cRT+(1|SubNum), a1[!a1$Choice1==-1,], 
                                      REML=FALSE)))

sv1_max <- svd(getME(CPPRTmod, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

# higher activity for larger VD
print (summary(CPPRTmodsv <- lmer(CPPPz~ sVD+(1|SubNum), a1[!a1$Choice1==-1,], 
                                  REML=FALSE)))

sv1_max <- svd(getME(CPPRTmod, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1) 

tab_model(CPPRTmod, CPPRTmodsv,show.stat = TRUE, col.order= c("est", "ci", "stat", "p"))
