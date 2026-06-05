## CLM PROJECT LDSC

setwd("../GSEM/non-UKB/Format_GWAS/")

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

list.files("../GSEM/Format_GWAS")
#vector of munged summary statisitcs
traits<-c("AF.sumstats.gz","HF.sumstats.gz","Stroke.sumstats.gz","PAD.sumstats.gz","CAD.sumstats.gz","VHD.sumstats.gz","NAFLD.sumstats.gz","ALT.sumstats.gz","AST.sumstats.gz","GGT.sumstats.gz","ALP.sumstats.gz","BMI.sumstats.gz","WC.sumstats.gz","BFP.sumstats.gz","T2D.sumstats.gz","HbA1c.sumstats.gz","HTN.sumstats.gz","nonHDL.sumstats.gz","HDL.sumstats.gz","logTG.sumstats.gz","CRP.sumstats.gz")

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(0.20,0.08,0.13,0.01,0.11,0.08,0.01,NA,NA,NA,NA,NA,NA,NA,0.17,NA,0.31,NA,NA,NA,NA)

#vector of population prevalences
population.prev<-c(0.20,0.08,0.13,0.01,0.11,0.08,0.01,NA,NA,NA,NA,NA,NA,NA,0.17,NA,0.31,NA,NA,NA,NA)

#the folder of LD scores
ld<-"eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"eur_w_ld_chr/"

#name the traits
trait.names<-c("AF","HF","Stroke","PAD","CAD","VHD","NAFLD","ALT","AST","GGT","ALP","BMI","WC","BFP","T2D","HbA1c","HTN","nonHDL","HDL","logTG","CRP")

#run LDSC
LDSCoutput_all<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names,
                     ldsc.log="CLM_LDSCoutput_ALLchrALLphenos")

#optional command to save the output as a .RData file for later use
save(LDSCoutput_all,file="CLM_LDSCoutput_ALLchrALLphenos.RData")

#load("CLM_LDSCoutput_ALLchrALLphenos.RData")

## MAKE DATA TABLE OF LDSC output
# convert to correlation matrix
cov_matrix<-as.matrix(LDSCoutput_all$S)
rownames(cov_matrix) <- colnames(cov_matrix)

# make cov table
cov_dat <- data.frame(col = rep(colnames(cov_matrix), each = nrow(cov_matrix)), 
                      row = rep(rownames(cov_matrix), ncol(cov_matrix)), 
                      value = as.vector(cov_matrix))
names(cov_dat) <- c("Trait 1", "Trait 2", "Genetic Covariance")
head(cov_dat)

# pull the SEs of the sampling covariance matrix
k<-nrow(LDSCoutput_all$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput_all$V))
rownames(SE) <- rownames(cov_matrix)
colnames(SE) <- colnames(cov_matrix)

SE_cov <- data.frame(col = rep(colnames(SE), each = nrow(SE)), 
                     row = rep(rownames(SE), ncol(SE)), 
                     value = as.vector(SE))
names(SE_cov) <- c("Trait 1", "Trait 2", "Genetic Covariance SE")
head(SE_cov)

#obtain the p-values and create a matrix and a data table of these in same format as gen cov info
library(magrittr)

p<-2*pnorm(abs(LDSCoutput_all$S/SE),lower.tail=FALSE)
rownames(p) <- rownames(cov_matrix)
colnames(p) <- colnames(cov_matrix)
upperTriangle(p) <- lowerTriangle(p, byrow=TRUE) # add pvalues to upper triangle of matrix

p_fdr <- matrix(p.adjust(as.vector(p), method='fdr'),ncol=21)
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
ratio <- tcrossprod(1 / sqrt(diag(LDSCoutput_all$S)))
S_Stand <- LDSCoutput_all$S * ratio

rownames(S_Stand) <- colnames(S_Stand)
corr_dat <- data.frame(col = rep(colnames(S_Stand), each = nrow(S_Stand)), 
                       row = rep(rownames(S_Stand), ncol(S_Stand)), 
                       value = as.vector(S_Stand))
names(corr_dat) <- c("Trait 1", "Trait 2", "Genetic Correlation")
head(corr_dat)

#calculate the ratio of the rescaled and original S matrices
scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)

#rescale the sampling correlation matrix by the appropriate diagonals
V_Stand <- LDSCoutput_all$V * tcrossprod(scaleO)

#enter SEs from diagonal of standardized V to get genetic correlation SE's matching the log file
r<-nrow(LDSCoutput_all$S)
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
write.csv(all_dat, file = "CLM_GeneticCorrelations_AllPhenotypes.csv")

# calculate the mean genetic correlation for each traits with all other traits
means <- numeric(length(trait.names))
for(i in 1:length(trait.names))
{
  trait <- trait.names[i]
  trait_subset <- subset(all_dat, all_dat$`Trait 1` == trait | all_dat$`Trait 2` == trait)
  means[i] <- mean(trait_subset$`Genetic Correlation`)
} # loop to calculate mean correlation with all other traits for each trait

mean_dat <- data.frame(trait.names, means) # create table of mean genetic correlations to use to establish overall comparisons
write.csv(mean_dat, file="CLM_MeanGeneticCorrelations.csv")

# Make a subset table of genetic correlations that are highly multicollinear (i.e. rg >.9)
collinear <- subset(all_dat, all_dat$`Genetic Correlation` >= 0.9)
write.csv(collinear, file="CLM_HighCollinearity.csv")

#check ratio difference before and after smoothing
S_LD<-as.matrix((nearPD(LDSCoutput_all$S, corr = FALSE))$mat)

#look at ratio and  difference of before and after smoothing
S_LD/LDSCoutput_all$S
S_LD-LDSCoutput_all$S # looks ok only very minimal smoothing was required

