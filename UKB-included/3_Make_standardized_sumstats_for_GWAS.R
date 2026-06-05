setwd("../GSEM/RAw_GWAS/")

# Load packages required
library(devtools)
library(GenomicSEM)
library(data.table)

### Specify the arugments for sumstats funtion
# make sure files and trait.names match order and names that were used for the LDSC step

trait.names<-c("AF","HF","Stroke","PAD","CAD","VHD","NAFLD2","ALT","AST","GGT","ALP","BMI","WC","BFP","T2D","HbA1c","HTN","nonHDL","HDL","logTG","CRP")

files=c("AF_EUR_format.txt.gz","HF_EUR_format.txt.gz","Stroke_EUR_format.txt.gz","PAD_EUR_format.txt.gz","CAD_EUR_format.txt.gz","VHD_EUR_format.txt.gz","NAFLD2_EUR_format.txt.gz","ALT_EUR_format.txt.gz","AST_EUR_format.txt.gz","GGT_EUR_format.txt.gz","ALP_EUR_format.txt.gz","BMI_EUR_format.txt.gz","WC_EUR_format.txt.gz","BFP_EUR_format.txt.gz","T2D_EUR_format.txt.gz","HbA1c_EUR_format.txt.gz","HTN_EUR_format.txt.gz","nonHDL_EUR_format.txt.gz","HDL_EUR_format.txt.gz","logTG_EUR_format.txt.gz","CRP_EUR_format.txt.gz")

ref="reference.1000G.maf.0.005.txt" 

se.logit=c(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)

OLS=c(F,F,F,F,F,F,F,T,T,T,T,T,T,T,F,T,F,T,T,T,T)

linprob=c(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)

N=c(NA,NA,NA,NA,NA,NA,NA,752839,394642,752766,753010,694649,394642,394642,NA,146806,NA,926571,1244580,1253277,575531)

prepped_sumstats <- sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=OLS,linprob=linprob,
                                             N=N,info.filter=0.6,maf.filter=0.01,keep.indel=FALSE,parallel=TRUE,cores=NULL)

fwrite(prepped_sumstats,file="../GSEM/Format_GWAS/prepped_sumstats_CLM_21traits.txt", row.names = FALSE, quote = FALSE) # save all sumstats into one file
