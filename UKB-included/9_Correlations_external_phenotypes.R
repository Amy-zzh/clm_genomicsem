setwd("../GSEM/External/Format_GWAS/")

library(dplyr)
library(GenomicSEM)
library(Matrix)
library(stats)
#install.packages("corrplot")
library(corrplot)
#install.packages("qgraph")
library(qgraph)
library(tidyr)
library(gdata)
library(data.table)

traits<-c("GF.sumstats.gz","F1.sumstats.gz","F2.sumstats.gz","F3.sumstats.gz","F4.sumstats.gz","PL.sumstats.gz","FI.sumstats.gz","BMR.sumstats.gz","INC.sumstats.gz","EDU.sumstats.gz","INSOM.sumstats.gz","SD.sumstats.gz","COFF.sumstats.gz","PA.sumstats.gz","WALK.sumstats.gz","LST.sumstats.gz","CIG.sumstats.gz","SI.sumstats.gz","ALC.sumstats.gz","ANX.sumstats.gz","MDD.sumstats.gz","PD.sumstats.gz","AD.sumstats.gz","ALS.sumstats.gz","GMV.sumstats.gz","WMH.sumstats.gz","BA.sumstats.gz","INT.sumstats.gz","FIS.sumstats.gz","RT.sumstats.gz","OCC.sumstats.gz","LONE.sumstats.gz","CKD.sumstats.gz","KF.sumstats.gz","NEP.sumstats.gz","PCOS.sumstats.gz","ENDO.sumstats.gz","AAM.sumstats.gz","AANM.sumstats.gz","CD.sumstats.gz","UC.sumstats.gz","IBD.sumstats.gz","IBS.sumstats.gz","CC.sumstats.gz","PUD.sumstats.gz","ASTH.sumstats.gz","COPD.sumstats.gz","PNE.sumstats.gz","COVID.sumstats.gz","FEV1.sumstats.gz","FVC.sumstats.gz","FF.sumstats.gz","RA.sumstats.gz","OST.sumstats.gz","Gout.sumstats.gz","HGS.sumstats.gz","SSTI.sumstats.gz","UTI.sumstats.gz","CI.sumstats.gz","AI.sumstats.gz","SEP.sumstats.gz","WBCC.sumstats.gz","RBCC.sumstats.gz","PCOS.sumstats.gz","NC.sumstats.gz","LC.sumstats.gz","TP.sumstats.gz","IL6.sumstats.gz","VD.sumstats.gz","EST.sumstats.gz","TEST.sumstats.gz","ALB.sumstats.gz","FERR.sumstats.gz","CALC.sumstats.gz","Urate.sumstats.gz","SOD.sumstats.gz","POT.sumstats.gz","eGFR.sumstats.gz","UACR.sumstats.gz","BUN.sumstats.gz")
trait.names<-c("GF","F1","F2","F3","F4","PL","FI","BMR","INC","EDU","INSOM","SD","COFF","PA","WALK","LST","CIG","SI","ALC","ANX","MDD","PD","AD","ALS","GMV","WMH","BA","INT","FIS","RT","OCC","LONE","CKD","KF","NEP","PCOS","ENDO","AAM","AANM","CD","UC","IBD","IBS","CC","PUD","ASTH","COPD","PNE","COVID","FEV1","FVC","FF","RA","OST","Gout","HGS","SSTI","UTI","CI","AI","SEP","WBCC","RBCC","PCOS","NC","LC","TP","IL6","VD","EST","TEST","ALB","FERR","CALC","Urate","SOD","POT","eGFR","UACR","BUN")

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.30,NA,NA,NA,NA,NA,NA,NA,NA,0.02,0.04,0.04,0.14,0.20,NA,NA,NA,NA,NA,NA,NA,0.18,0.09,0.02,0.01,0.02,0.05,NA,NA,0.30,0.27,0.42,0.11,0.04,0.004,0.21,0.06,0.08,0.01,NA,NA,NA,0.23,0.26,0.05,NA,0.04,0.05,0.01,0.18,0.03,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
#vector of population prevalences
population.prev<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.30,NA,NA,NA,NA,NA,NA,NA,NA,0.02,0.04,0.001,0.05,0.0001,NA,NA,NA,NA,NA,NA,NA,0.18,0.13,0.02,0.01,0.06,0.01,NA,NA,0.001,0.002,0.003,0.10,0.04,0.004,0.07,0.12,0.01,0.01,NA,NA,NA,0.01,0.39,0.04,NA,0.14,0.05,0.01,0.18,0.005,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)

#the folder of LD scores
ld<-"../GSEM/Format_GWAS/eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"../GSEM/Format_GWAS/eur_w_ld_chr/"

#LDSC
LDSCoutput_pairs <- ldsc(traits = traits,
                         sample.prev = sample.prev,
                         population.prev = population.prev,
                         ld = ld, wld = wld,
                         trait.names = trait.names,
                         ldsc.log = "CRHM_LDSCoutput_BF")  


#optional command to save the output as a .RData file for later use
save(LDSCoutput_pairs,file="CRHM_LDSCoutput_BF.RData")

#load("CRHM_LDSCoutput_BF.RData")

## MAKE DATA TABLE OF LDSC output
# convert to correlation matrix
cov_matrix<-as.matrix(LDSCoutput_pairs$S)
rownames(cov_matrix) <- colnames(cov_matrix)

# make cov table
cov_dat <- data.frame(col = rep(colnames(cov_matrix), each = nrow(cov_matrix)), 
                      row = rep(rownames(cov_matrix), ncol(cov_matrix)), 
                      value = as.vector(cov_matrix))
names(cov_dat) <- c("Trait 1", "Trait 2", "Genetic Covariance")
head(cov_dat)

# pull the SEs of the sampling covariance matrix
k<-nrow(LDSCoutput_pairs$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput_pairs$V))
rownames(SE) <- rownames(cov_matrix)
colnames(SE) <- colnames(cov_matrix)

SE_cov <- data.frame(col = rep(colnames(SE), each = nrow(SE)), 
                     row = rep(rownames(SE), ncol(SE)), 
                     value = as.vector(SE))
names(SE_cov) <- c("Trait 1", "Trait 2", "Genetic Covariance SE")
head(SE_cov)

#obtain the p-values and create a matrix and a data table of these in same format as gen cov info
library(magrittr)

p<-2*pnorm(abs(LDSCoutput_pairs$S/SE),lower.tail=FALSE)
rownames(p) <- rownames(cov_matrix)
colnames(p) <- colnames(cov_matrix)
upperTriangle(p) <- lowerTriangle(p, byrow=TRUE) # add pvalues to upper triangle of matrix

p_fdr <- matrix(p.adjust(as.vector(p), method='fdr'),ncol=80)
rownames(p_fdr) <- rownames(cov_matrix)
colnames(p_fdr) <- colnames(cov_matrix)
upperTriangle(p_fdr) <- lowerTriangle(p_fdr, byrow=TRUE) # add pvalues to upper triangle of matrix

p0_dat <- data.frame(col = rep(colnames(p), each = nrow(p)), 
                     row = rep(rownames(p), ncol(p)), 
                     value = as.vector(p))
names(p0_dat) <- c("Trait 1", "Trait 2", "P0 value")
head(p0_dat)

p_dat <- data.frame(col = rep(colnames(p_fdr), each = nrow(p_fdr)), 
                    row = rep(rownames(p_fdr), ncol(p_fdr)), 
                    value = as.vector(p_fdr))
names(p_dat) <- c("Trait 1", "Trait 2", "P value")
head(p_dat)

# Calculate & create correlation matrix using same code as function that prints the log file
ratio <- tcrossprod(1 / sqrt(diag(LDSCoutput_pairs$S)))
S_Stand <- LDSCoutput_pairs$S * ratio

rownames(S_Stand) <- colnames(S_Stand)
corr_dat <- data.frame(col = rep(colnames(S_Stand), each = nrow(S_Stand)), 
                       row = rep(rownames(S_Stand), ncol(S_Stand)), 
                       value = as.vector(S_Stand))
names(corr_dat) <- c("Trait 1", "Trait 2", "Genetic Correlation")
head(corr_dat)

#calculate the ratio of the rescaled and original S matrices
scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)

#rescale the sampling correlation matrix by the appropriate diagonals
V_Stand <- LDSCoutput_pairs$V * tcrossprod(scaleO)

#enter SEs from diagonal of standardized V to get genetic correlation SE's matching the log file
r<-nrow(LDSCoutput_pairs$S)
SE_Stand<-matrix(0, r, r)
SE_Stand[lower.tri(SE_Stand,diag=TRUE)] <-sqrt(diag(V_Stand))

colnames(SE_Stand) <- colnames(S_Stand)
rownames(SE_Stand) <- rownames(S_Stand)

corrSE_dat <- data.frame(col = rep(colnames(SE_Stand), each = nrow(SE_Stand)), 
                         row = rep(rownames(SE_Stand), ncol(SE_Stand)), 
                         value = as.vector(SE_Stand))
names(corrSE_dat) <- c("Trait 1", "Trait 2", "Genetic Correlation SE")
head(corrSE_dat)

# make a table with all details
all_dat <- cov_dat
all_dat$`Genetic Covariance SE` <- SE_cov$`Genetic Covariance SE`
all_dat$`Genetic Correlation` <- corr_dat$`Genetic Correlation`
all_dat$`Genetic Correlation SE` <- corrSE_dat$`Genetic Correlation SE`
all_dat$`P value` <- p_dat$`P value`
all_dat$`P0 value` <- p0_dat$`P0 value`
head(all_dat)

# remove all diagonal rows (i.e. trait with itself)
all_dat <- subset(all_dat, all_dat$`Trait 1` != all_dat$`Trait 2`)
# remove the replicated upper diagonal rows (so that each pairwise correlation is only included once in the table)
all_dat <- subset(all_dat, all_dat$`Genetic Covariance SE`!=0 & all_dat$`Genetic Correlation SE` !=0)
write.csv(all_dat, file = "CRHM_GeneticCorrelations_BF.csv")

