# CLM Genomic Factor Analysis

# Conduct EFA in the odd autosomes
setwd("../GSEM/Format_GWAS/")

library(dplyr)
library(GenomicSEM)
library(Matrix)
library(stats)
library(tidyr)
library(gdata)
library(data.table)
library(plyr)
library(e1071)
library(readr)
library(lavaan)
library(stringr)
library(splitstackshape)
library(R.utils)

#vector of munged summary statisitcs
traits<-c("AF.sumstats.gz","HF.sumstats.gz","Stroke.sumstats.gz","PAD.sumstats.gz","CAD.sumstats.gz","VHD.sumstats.gz","NAFLD2.sumstats.gz","ALT.sumstats.gz","AST.sumstats.gz","GGT.sumstats.gz","ALP.sumstats.gz","BMI.sumstats.gz","WC.sumstats.gz","BFP.sumstats.gz","T2D.sumstats.gz","HbA1c.sumstats.gz","HTN.sumstats.gz","nonHDL.sumstats.gz","HDL.sumstats.gz","logTG.sumstats.gz","CRP.sumstats.gz")

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(0.12,0.05,0.06,0.01,0.16,0.08,0.01,NA,NA,NA,NA,NA,NA,NA,0.10,NA,0.27,NA,NA,NA,NA)

#vector of population prevalences
population.prev<-c(0.02,0.02,0.01,0.02,0.08,0.08,0.25,NA,NA,NA,NA,NA,NA,NA,0.09,NA,0.27,NA,NA,NA,NA)

#the folder of LD scores
ld<-"eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"eur_w_ld_chr/"

#name the traits
trait.names<-c("AF","HF","Stroke","PAD","CAD","VHD","NAFLD2","ALT","AST","GGT","ALP","BMI","WC","BFP","T2D","HbA1c","HTN","nonHDL","HDL","logTG","CRP")


#run LDSC - odd
LDSCoutput_odd<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names,
                     ldsc.log="CLM_LDSCoutput_ODDchrEFA_21traits",select = "ODD")

#optional command to save the output as a .RData file for later use
save(LDSCoutput_odd,file="CLM_LDSCoutput_ODDchrEFA_21traits.RData")

#load("CLM_LDSCoutput_ODDchrEFA_21traits.RData")
#save as separate object for smoothing
S_LD<-as.matrix((nearPD(LDSCoutput_odd$S, corr = FALSE))$mat)

#look at ratio and  difference of before and after smoothing
S_LD/LDSCoutput_odd$S
S_LD-LDSCoutput_odd$S

##############################################
######code to calculate number of factors needed
##############################################

eig <- eigen(cov2cor(S_LD), only.values = TRUE)
eig <- eig$values; nk <- length(eig); k <- 1:nk
criteria <- mean(eig)
# Compute the proportion of total variance explained by the eigenvalues.
for (r in eig) {
  cat(round(r / sum(eig), digits = 3), "")
}
# Compute the Kaiser rule.
nkaiser <- sum(eig >= rep(criteria, nk))
# Compute the acceleration factor.
aparallel <- rep(criteria, length(eig))
pred.eig <- af <- rep(NA, nk)
for (j in 2:(nk - 1)) {
  if (eig[j - 1] >= aparallel[j - 1]) {
    af[j] <- (eig[j + 1] - 2 * eig[j]) + eig[j - 1]
  }
}
naf <- which(af == max(af, na.rm = TRUE), TRUE) - 1
# Compute the optimal coordinates.
proportion <- eig/sum(eig)
cumulative <- proportion
for (i in 2:nk) cumulative[i] = cumulative[i - 1] + proportion[i]
proportion[proportion < 0] <- 0
cond1 <- TRUE
cond2 <- TRUE
i <- 0
pred.eig <- af <- rep(NA, nk)
while ((cond1 == TRUE) && (cond2 == TRUE) && (i < nk)) {
  i <- i + 1
  ind <- k[c(i + 1, nk)]
  vp.p <- lm(eig[c(i + 1, nk)] ~ ind)
  vp.prec <- pred.eig[i] <- sum(c(1, i) * coef(vp.p))
  cond1 <- (eig[i] >= vp.prec)
  cond2 <- (eig[i] >= aparallel[i])
  noc <- i - 1
}

# Print results.
cat("Number of factors to retain according to the Kaiser rule:", nkaiser, "\n",
    "Number of factors to retain according to the acceleration factor:", naf, "\n",
    "Number of factors to retain according to the optimal coordinates:", noc, "\n")


#function to automatically write model based on EFA results
#the arguments to the function are: 
##Loadings: the matrix of EFA lodaings
#S_LD: the LDSC genetic covariance matrix
#cutoff: the EFA loadings cutoff
#fix_resid: whether to apply constraint on all variables to keep residual variances above 0
#bifactor: whether to specify a bifactor model
#mustload: whether all variables should load on at least one factor, even if they don't meet EFA cutoff
#common: whether to specify a common factor model 

write.model<-function(Loadings,S_LD,cutoff,fix_resid=TRUE,bifactor=FALSE,mustload=FALSE,common=FALSE){
  Model<-""
  if(common == TRUE){
    for(f in 1){
      u<-1
      Model1<-""
      for(i in 1:nrow(S_LD)){
        if(u == 1){
          linestart<-paste("F", f, "=~",  colnames(S_LD)[i], sep = "")
          u<-u+1
          linemid<-""
        }else{
          linemid<-paste(linemid, " + ", colnames(S_LD)[i], sep = "")
        }
      }
    }
    Model<-paste(Model,linestart, linemid, " \n ", sep="")
  }else{
    if(mustload == TRUE){
      Mins<-apply(abs(Loadings),1,max)
      for(i in 1:nrow(Loadings)){
        for(f in 1:ncol(Loadings)){
          if(Mins[i] == Loadings[i,f]){
            Loadings[i,f]<-cutoff+.01
          }
        }
      }
    }
    
    
    for(f in 1:ncol(Loadings)){
      u<-1
      Model1<-""
      for(i in 1:nrow(Loadings)){
        if(abs(Loadings[i,f]) > cutoff){
          if(u == 1){
            linestart<-paste("F", f, "=~NA*",  colnames(S_LD)[i], sep = "")
            u<-u+1
            linemid<-""
          }else{
            linemid<-paste(linemid, " + ", colnames(S_LD)[i], sep = "")
          }
        }
      }
      Model<-paste(Model,linestart, linemid, " \n ", sep="")
      linestart<-""
      linemid<-""
    }
    
    if(bifactor==TRUE){
      Model_bi<-""
      u<-1
      for(i in 1:ncol(S_LD)){
        b<-grepl(colnames(S_LD)[i], Model)
        if(b == TRUE){
          if(u == 1){
            linestart_bi<-paste("Common_F", "=~",  colnames(S_LD)[i], sep = "")
            u<-u+1
            linemid_bi<-""
          }else{
            linemid_bi<-paste(linemid_bi, " + ", colnames(S_LD)[i], sep = "")
          }
        }
      }
      Model_bi<-paste(linestart_bi,linemid_bi," \n ", sep="")
      r<-1
      Factor_bi<-""
      for(i in 1:ncol(Loadings)){
        Factor_bi_start<-paste("Common_F~~0*", "F", i, " \n ", sep = "")
        Factor_bi<-paste(Factor_bi,Factor_bi_start,sep="")
      }
      
      Modelsat<-""
      for (i in 1:(ncol(Loadings)-1)) {
        linestartc <- paste("", "F", i, "~~0*F", i+1, sep = "")
        if (ncol(Loadings)-i >= 2) {
          linemidc <- ""
          for (j in (i+2):ncol(Loadings)) {
            linemidc <- paste("", linemidc, "F", i, "~~0*F", j, " \n ", sep="")
            
          }
        } else {linemidc <- ""}
        Modelsat <- paste(Modelsat, linestartc, " \n ", linemidc, sep = "")
      } 
      
      Model<-paste(Model,Model_bi,Factor_bi,Modelsat)
    }
    
    
  }
  
  if(fix_resid==TRUE){
    Model3<-""
    #create unique combination of letters for residual variance parameter labels
    n<-combn(letters,4)[,sample(1:14000, ncol(S_LD), replace=FALSE)]
    for(i in 1:ncol(S_LD)){
      if(grepl(colnames(S_LD)[i],Model) == TRUE){
        linestart3a <- paste(colnames(S_LD)[i], " ~~ ",  paste(n[,i],collapse=""), "*", colnames(S_LD)[i], sep = "")
        linestart3b <- paste(paste(n[,i],collapse=""), " > .0001", sep = "")
        Model3<-paste(Model3, linestart3a, " \n ", linestart3b, " \n ", sep = "")
      }
    }
    Model<-paste(Model,Model3,sep="")
  }
  
  
  return(Model)
}

##5-factor model - Kaiser
EFA_5<-factanal(factors = 5, covmat = S_LD, n.obs = 2, rotation = "promax")
print(EFA_5, sort=TRUE)
print(EFA_5,digits=3,cutoff=.30,sort=TRUE)
Loadings_5<-(EFA_5$loadings[1:ncol(S_LD),1:5])
Model_5<-write.model(Loadings_5,S_LD,.30,mustload=TRUE)

##4-factor model - Kaiser
EFA_4<-factanal(factors = 4, covmat = S_LD, n.obs = 2, rotation = "promax")
print(EFA_4, sort=TRUE)
print(EFA_4,digits=3,cutoff=.30,sort=TRUE)
Loadings_4<-(EFA_4$loadings[1:ncol(S_LD),1:4])
Model_4<-write.model(Loadings_4,S_LD,.30,mustload=TRUE)

##3-factor model - Kaiser
EFA_3<-factanal(factors = 3, covmat = S_LD, n.obs = 2, rotation = "promax")
print(EFA_3, sort=TRUE)
print(EFA_3,digits=3,cutoff=.30,sort=TRUE)
Loadings_3<-(EFA_3$loadings[1:ncol(S_LD),1:3])
Model_3<-write.model(Loadings_3,S_LD,.30,mustload=TRUE)

######## CFA

#run LDSC - even
LDSCoutput_even<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names,
                      ldsc.log="CLM_LDSCoutput_EVEN_CFA_21traits",select = "EVEN")

#optional command to save the output as a .RData file for later use
save(LDSCoutput_even,file="CLM_LDSCoutput_EVEN_CFA_21traits.RData")

#load("CLM_LDSCoutput_EVEN_CFA_21traits.RData")

# run 5-factor CFA models
# .30 loadings
CFAeven_5 <- usermodel(LDSCoutput_even, estimation = "DWLS", model = Model_5, 
                       CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAeven_5

# run 4-factor CFA models
# .30 loadings
CFAeven_4 <- usermodel(LDSCoutput_even, estimation = "DWLS", model = Model_4, 
                       CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAeven_4

# run 3-factor CFA models
# .30 loadings
CFAeven_3 <- usermodel(LDSCoutput_even, estimation = "DWLS", model = Model_3, 
                       CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAeven_3

####### GENOME-WIDE STUFF
#run LDSC - all autosomes
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names,
                 ldsc.log="CLM_LDSCoutput_ALL_21traits")


#optional command to save the output as a .RData file for later use
save(LDSCoutput,file="CLM_LDSCoutput_ALL_21traits.RData")

#load("CLM_LDSCoutput_ALL_21traits.RData")
# Genome-wide common factor model
CommonFactor<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")

CommonFactor

# genome wide CFA - 6 factor model

# genome wide CFA - 5 factor model
# .30 loadings
CFAall_5 <- usermodel(LDSCoutput, estimation = "DWLS", model = Model_5, 
                      CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAall_5

# genome wide CFA - 4 factor model
# .30 loadings
CFAall_4 <- usermodel(LDSCoutput, estimation = "DWLS", model = Model_4, 
                      CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAall_4
#$modelfit
#chisq  df p_chisq      AIC       CFI      SRMR
#df 45479.42 182       0 45577.42 0.9519347 0.0757914

# genome wide CFA - 3 factor model
# .30 loadings
CFAall_3 <- usermodel(LDSCoutput, estimation = "DWLS", model = Model_3, 
                      CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFAall_3



# bifactor model with 4 factors - are the loadings for F1 near 0?
bifactor_model <- "BF=~NA*AF + HF + Stroke + PAD + CAD + VHD + NAFLD2 + ALT + AST + GGT + ALP + BMI + WC + BFP + T2D + HbA1c + HTN + nonHDL + HDL + logTG + CRP
F1=~NA*AF + HF + Stroke + PAD + CAD + VHD + HbA1c + HTN
F2=~NA*BMI + WC + BFP + T2D + CRP
F3=~NA*NAFLD2 + ALT + AST + GGT + ALP + T2D
F4=~NA*nonHDL + HDL + logTG
BF~~0*F1 
BF~~0*F2
BF~~0*F3
BF~~0*F4
F1 ~~ F2
F1 ~~ F3
F1 ~~ F4
F2 ~~ F3
F2 ~~ F4
F3 ~~ F4
HF ~~ a*HF
a > .0001"


bifactor <- usermodel(LDSCoutput, estimation = "DWLS", model = bifactor_model, 
                      CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
bifactor



# remove any loadings on F1-4 that are below .30
bifactor_model_above30 <- "BF=~NA*AF + HF + Stroke + PAD + CAD + VHD + NAFLD2 + ALT + AST + GGT + ALP + BMI + WC + BFP + T2D + HbA1c + HTN + nonHDL + HDL + logTG + CRP
F1=~NA*AF + HF + Stroke + PAD + CAD + VHD + HTN
F2=~NA*BMI + WC + BFP
F3=~NA*NAFLD2 + ALT + AST + GGT
F4=~a*nonHDL + a*logTG
BF~~0*F1 
BF~~0*F2
BF~~0*F3
BF~~0*F4
F1 ~~ F2
F1 ~~ F3
F1 ~~ F4
F2 ~~ F3
F2 ~~ F4
F3 ~~ F4
HF ~~ b*HF
b > .001"


bifactor_above30 <- usermodel(LDSCoutput, estimation = "DWLS", model = bifactor_model_above30, 
                              CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
bifactor_above30


