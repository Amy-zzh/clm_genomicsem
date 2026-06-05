setwd("../GSEM/non-UKB/Format_GWAS/BF_results/")

require(stringr)
require(data.table)
require(RColorBrewer)
require(lattice)
require(dplyr)
require(rlist)
library(R.utils)
#install.packages("qqman")
library(qqman)

# save .tab files from array jobs into separate directories for each latent factor
# Here we provide an example for the General Factor but we repeated for each of the other latent factors too
tempfreq<-list.files(pattern="*.tab")
tempfreq<-str_sort(tempfreq, numeric = TRUE)

Data<-rbindlist(lapply(tempfreq, fread, header=TRUE,fill=TRUE))
rm(tempfreq)
head(Data)

# QC
table(is.na(Data$SNP))
# none
table(duplicated(Data$SNP))
# none
table(Data$error)
# none
table(Data$warning)
#warn <- subset(Data, Data$warning != 0) # there were 290 SNPs with negative lv variances, 1 SNP with negative ov variance & 1 SNP with non positive definite cov matrix of latent vars

#rm(warn)
#rm(sub)
# there doesn't seem to be an obvious reason as to why SNPs are erroneous so we remove SNPs with warnings - 63392 SNPs in total removed
Data <- subset(Data, Data$warning == 0)
table(Data$Z_smooth)
# 0 SNP removed
Data <- subset(Data, Data$Z_smooth <=1.96)
# Number of final SNPs:
# 2070948 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 38932 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 6349 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_CLMGWAS_BF_sumstats.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_BF_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_CLMGWAS_BF_topQsnps.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

Data_short <- Data %>% select(SNP, CHR, BP) # subset to only req columns for clumping to save mem
## Calculate regions around Q SNPs to remove
# add two columns to the Q SNPlist that provide a range of +/- 1,000,000 base pairs around each of the position points to get the range of values to filter on for each Q hit
Q_SNPs$max <- Q_SNPs$BP + 1000000
Q_SNPs$min <- Q_SNPs$BP - 1000000
head(Q_SNPs)
Q_SNPs_max = as.list(Q_SNPs$max)
Q_SNPs_min = as.list(Q_SNPs$min)

nestedDT <- list()

for (i in 1:nrow(Q_SNPs)) {
  CHROM <- Q_SNPs$CHR[i]
  MAX <- Q_SNPs$max[i]
  MIN <- Q_SNPs$min[i]
  temp <- Data_short[CHR == CHROM & BP <= MAX & BP >= MIN]
  nestedDT[[i]] <- temp
}

names(nestedDT) <- Q_SNPs$SNP
Q_all<-rbindlist(nestedDT, use.names=T)
# remove duplicate SNPs
table(duplicated(Q_all$SNP))
Q_all<-subset(Q_all, !(duplicated(Q_all$SNP)))
fwrite(Q_all, file = "bifactor_CLMGWAS_BF_clumpedQsnps.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
# remove clumped Q SNPs from main file and save


# Remove Q SNPs from data and save for post GWAS analyses where don't want to include Q signal in results
Data_noQ<-anti_join(Data, Q_all, by='SNP')
table(Data_noQ$Q_SNP_pval <= 1e-08) # check Q is definitely removed

##Calculate Effective Sample Size for non Q sample
#restrict to MAF of 40% and 10%
CorrelatedFactors<-subset(Data_noQ, Data_noQ$MAF <= .4 & Data_noQ$MAF >= .1)

N_hat<-mean(1/((2*CorrelatedFactors$MAF*(1-CorrelatedFactors$MAF))*CorrelatedFactors$SE^2))
N_hat
# Bifactor = 156599

Data_noQ$N <- 156599
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 1e-8)
# Bifactor: 2903 SNPs

fwrite(Data_noQ, file = "bifactor_CLMGWAS_BF_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_BF_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_CLMGWAS_BF_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_CLMGWAS_BF_noQ_postGWAS.txt")

# MUNGE FORMATTED DATA FOR GSEM
files<-fread("bifactor_CLMGWAS_BF_noQ_postGWAS.txt.gz")

files <- files %>% select(CHR, BP, SNP, A1, A2, BETA, SE, P, MAF, N) # select columns needed for GSEM
names(files) <- c("CHR", "POS", "SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N") # rename selected columns for easy recognition in subsequent packages
head(files)

fwrite(files, file = "bifactor_CLMGWAS_BF_forsumstat.txt", 
       row.names = FALSE, quote = FALSE, col.names = TRUE,
       sep = " ")
gzip("bifactor_CLMGWAS_BF_forsumstat.txt")

files<-c("../GSEM/non-UKB/Format_GWAS/BF_results/bifactor_CLMGWAS_BF_forsumstat.txt.gz")
trait.names<-c("BF")
#sample size: NA if SNP-specific Neff used, Neff if sum of cohorts Neffs used or binary single cohort; total N if one samPAe continuous 
N=c(156599)
#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

##########################################################################
setwd("../GSEM/non-UKB/Format_GWAS/F1_results/")

require(stringr)
require(data.table)
require(RColorBrewer)
require(lattice)
require(dplyr)
require(rlist)
library(R.utils)
#install.packages("qqman")
library(qqman)

# save .tab files from array jobs into separate directories for each latent factor
# Here we provide an example for the General Factor but we repeated for each of the other latent factors too
tempfreq<-list.files(pattern="*.tab")
tempfreq<-str_sort(tempfreq, numeric = TRUE)

Data<-rbindlist(lapply(tempfreq, fread, header=TRUE,fill=TRUE))
rm(tempfreq)
head(Data)

# QC
table(is.na(Data$SNP))
# none
table(duplicated(Data$SNP))
# none
table(Data$error)
# none
table(Data$warning)
#warn <- subset(Data, Data$warning != 0) # there were 290 SNPs with negative lv variances, 1 SNP with negative ov variance & 1 SNP with non positive definite cov matrix of latent vars

#rm(warn)
#rm(sub)
# there doesn't seem to be an obvious reason as to why SNPs are erroneous so we remove SNPs with warnings - 63392 SNPs in total removed
Data <- subset(Data, Data$warning == 0)
nrow(Data)
table(Data$Z_smooth)
# 0 SNP removed
Data <- subset(Data, Data$Z_smooth <=1.96)
nrow(Data)
# Number of final SNPs:
# 2070948 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 4786 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 851 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_CLMGWAS_F1_sumstats.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_F1_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_CLMGWAS_F1_topQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

Data_short <- Data %>% select(SNP, CHR, BP) # subset to only req columns for clumping to save mem
## Calculate regions around Q SNPs to remove
# add two columns to the Q SNPlist that provide a range of +/- 1,000,000 base pairs around each of the position points to get the range of values to filter on for each Q hit
Q_SNPs$max <- Q_SNPs$BP + 1000000
Q_SNPs$min <- Q_SNPs$BP - 1000000
head(Q_SNPs)
Q_SNPs_max = as.list(Q_SNPs$max)
Q_SNPs_min = as.list(Q_SNPs$min)

nestedDT <- list()

for (i in 1:nrow(Q_SNPs)) {
  CHROM <- Q_SNPs$CHR[i]
  MAX <- Q_SNPs$max[i]
  MIN <- Q_SNPs$min[i]
  temp <- Data_short[CHR == CHROM & BP <= MAX & BP >= MIN]
  nestedDT[[i]] <- temp
}

names(nestedDT) <- Q_SNPs$SNP
Q_all<-rbindlist(nestedDT, use.names=T)
# remove duplicate SNPs
table(duplicated(Q_all$SNP))
Q_all<-subset(Q_all, !(duplicated(Q_all$SNP)))
fwrite(Q_all, file = "bifactor_CLMGWAS_F1_clumpedQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)
# remove clumped Q SNPs from main file and save


# Remove Q SNPs from data and save for post GWAS analyses where don't want to include Q signal in results
Data_noQ<-anti_join(Data, Q_all, by='SNP')
table(Data_noQ$Q_SNP_pval <= 1e-08) # check Q is definitely removed
# Bifactor: 1828981 SNPs

##Calculate Effective Sample Size for non Q sample
#restrict to MAF of 40% and 10%
CorrelatedFactors<-subset(Data_noQ, Data_noQ$MAF <= .4 & Data_noQ$MAF >= .1)

N_hat<-mean(1/((2*CorrelatedFactors$MAF*(1-CorrelatedFactors$MAF))*CorrelatedFactors$SE^2))
N_hat
# Bifactor = 171580

Data_noQ$N <- 171580
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 5e-8)
# Bifactor: 2064 SNPs

fwrite(Data_noQ, file = "bifactor_CLMGWAS_F1_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_F1_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_CLMGWAS_F1_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_CLMGWAS_F1_noQ_postGWAS.txt")

# MUNGE FORMATTED DATA FOR GSEM
files<-fread("bifactor_CLMGWAS_F1_noQ_postGWAS.txt.gz")

files <- files %>% select(CHR, BP, SNP, A1, A2, BETA, SE, P, MAF, N) # select columns needed for GSEM
names(files) <- c("CHR", "POS", "SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N") # rename selected columns for easy recognition in subsequent packages
head(files)

fwrite(files, file = "bifactor_CLMGWAS_F1_forsumstat.txt", 
       row.names = FALSE, quote = FALSE, col.names = TRUE,
       sep = " ")
gzip("bifactor_CLMGWAS_F1_forsumstat.txt")

files<-c("bifactor_CLMGWAS_F1_forsumstat.txt.gz")
trait.names<-c("F1")
#sample size: NA if SNP-specific Neff used, Neff if sum of cohorts Neffs used or binary single cohort; total N if one samPAe continuous 
N=c(171580)
#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

##########################################################################
setwd("../GSEM/non-UKB/Format_GWAS/F2_results/")

require(stringr)
require(data.table)
require(RColorBrewer)
require(lattice)
require(dplyr)
require(rlist)
library(R.utils)
#install.packages("qqman")
library(qqman)

# save .tab files from array jobs into separate directories for each latent factor
# Here we provide an example for the General Factor but we repeated for each of the other latent factors too
tempfreq<-list.files(pattern="*.tab")
tempfreq<-str_sort(tempfreq, numeric = TRUE)

Data<-rbindlist(lapply(tempfreq, fread, header=TRUE,fill=TRUE))
rm(tempfreq)
head(Data)

# QC
table(is.na(Data$SNP))
# none
table(duplicated(Data$SNP))
# none
table(Data$error)
# none
table(Data$warning)
#warn <- subset(Data, Data$warning != 0) # there were 290 SNPs with negative lv variances, 1 SNP with negative ov variance & 1 SNP with non positive definite cov matrix of latent vars

#rm(warn)
#rm(sub)
# there doesn't seem to be an obvious reason as to why SNPs are erroneous so we remove SNPs with warnings - 63392 SNPs in total removed
Data <- subset(Data, Data$warning == 0)
nrow(Data)
table(Data$Z_smooth)
# 0 SNP removed
Data <- subset(Data, Data$Z_smooth <=1.96)
nrow(Data)
# Number of final SNPs:
# 2070948 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 42113 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 5186 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_CLMGWAS_F2_sumstats.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_F2_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_CLMGWAS_F2_topQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

Data_short <- Data %>% select(SNP, CHR, BP) # subset to only req columns for clumping to save mem
## Calculate regions around Q SNPs to remove
# add two columns to the Q SNPlist that provide a range of +/- 1,000,000 base pairs around each of the position points to get the range of values to filter on for each Q hit
Q_SNPs$max <- Q_SNPs$BP + 1000000
Q_SNPs$min <- Q_SNPs$BP - 1000000
head(Q_SNPs)
Q_SNPs_max = as.list(Q_SNPs$max)
Q_SNPs_min = as.list(Q_SNPs$min)

nestedDT <- list()

for (i in 1:nrow(Q_SNPs)) {
  CHROM <- Q_SNPs$CHR[i]
  MAX <- Q_SNPs$max[i]
  MIN <- Q_SNPs$min[i]
  temp <- Data_short[CHR == CHROM & BP <= MAX & BP >= MIN]
  nestedDT[[i]] <- temp
}

names(nestedDT) <- Q_SNPs$SNP
Q_all<-rbindlist(nestedDT, use.names=T)
# remove duplicate SNPs
table(duplicated(Q_all$SNP))
Q_all<-subset(Q_all, !(duplicated(Q_all$SNP)))
fwrite(Q_all, file = "bifactor_CLMGWAS_F2_clumpedQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)
# remove clumped Q SNPs from main file and save


# Remove Q SNPs from data and save for post GWAS analyses where don't want to include Q signal in results
Data_noQ<-anti_join(Data, Q_all, by='SNP')
table(Data_noQ$Q_SNP_pval <= 1e-08) # check Q is definitely removed
# Bifactor: 2070948 SNPs

##Calculate Effective Sample Size for non Q sample
#restrict to MAF of 40% and 10%
CorrelatedFactors<-subset(Data_noQ, Data_noQ$MAF <= .4 & Data_noQ$MAF >= .1)

N_hat<-mean(1/((2*CorrelatedFactors$MAF*(1-CorrelatedFactors$MAF))*CorrelatedFactors$SE^2))
N_hat
# Bifactor = 423057

Data_noQ$N <- 423057
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 5e-8)
# Bifactor: 11662  SNPs

fwrite(Data_noQ, file = "bifactor_CLMGWAS_F2_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_F2_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_CLMGWAS_F2_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_CLMGWAS_F2_noQ_postGWAS.txt")

# MUNGE FORMATTED DATA FOR GSEM
files<-fread("bifactor_CLMGWAS_F2_noQ_postGWAS.txt.gz")

files <- files %>% select(CHR, BP, SNP, A1, A2, BETA, SE, P, MAF, N) # select columns needed for GSEM
names(files) <- c("CHR", "POS", "SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N") # rename selected columns for easy recognition in subsequent packages
head(files)

fwrite(files, file = "bifactor_CLMGWAS_F2_forsumstat.txt", 
       row.names = FALSE, quote = FALSE, col.names = TRUE,
       sep = " ")
gzip("bifactor_CLMGWAS_F2_forsumstat.txt")

files<-c("bifactor_CLMGWAS_F2_forsumstat.txt.gz")
trait.names<-c("F2")
#sample size: NA if SNP-specific Neff used, Neff if sum of cohorts Neffs used or binary single cohort; total N if one samPAe continuous 
N=c(423057)
#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

##########################################################################
setwd("../GSEM/non-UKB/Format_GWAS/F3_results/")

require(stringr)
require(data.table)
require(RColorBrewer)
require(lattice)
require(dplyr)
require(rlist)
library(R.utils)
#install.packages("qqman")
library(qqman)

# save .tab files from array jobs into separate directories for each latent factor
# Here we provide an example for the General Factor but we repeated for each of the other latent factors too
tempfreq<-list.files(pattern="*.tab")
tempfreq<-str_sort(tempfreq, numeric = TRUE)

Data<-rbindlist(lapply(tempfreq, fread, header=TRUE,fill=TRUE))
rm(tempfreq)
head(Data)

# QC
table(is.na(Data$SNP))
# none
table(duplicated(Data$SNP))
# none
table(Data$error)
# none
table(Data$warning)
#warn <- subset(Data, Data$warning != 0) # there were 290 SNPs with negative lv variances, 1 SNP with negative ov variance & 1 SNP with non positive definite cov matrix of latent vars

#rm(warn)
#rm(sub)
# there doesn't seem to be an obvious reason as to why SNPs are erroneous so we remove SNPs with warnings - 63392 SNPs in total removed
Data <- subset(Data, Data$warning == 0)
nrow(Data)
table(Data$Z_smooth)
# 0 SNP removed
Data <- subset(Data, Data$Z_smooth <=1.96)
nrow(Data)
# Number of final SNPs:
# 2070948 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 7910 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 2099 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_CLMGWAS_F3_sumstats.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_F3_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_CLMGWAS_F3_topQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

Data_short <- Data %>% select(SNP, CHR, BP) # subset to only req columns for clumping to save mem
## Calculate regions around Q SNPs to remove
# add two columns to the Q SNPlist that provide a range of +/- 1,000,000 base pairs around each of the position points to get the range of values to filter on for each Q hit
Q_SNPs$max <- Q_SNPs$BP + 1000000
Q_SNPs$min <- Q_SNPs$BP - 1000000
head(Q_SNPs)
Q_SNPs_max = as.list(Q_SNPs$max)
Q_SNPs_min = as.list(Q_SNPs$min)

nestedDT <- list()

for (i in 1:nrow(Q_SNPs)) {
  CHROM <- Q_SNPs$CHR[i]
  MAX <- Q_SNPs$max[i]
  MIN <- Q_SNPs$min[i]
  temp <- Data_short[CHR == CHROM & BP <= MAX & BP >= MIN]
  nestedDT[[i]] <- temp
}

names(nestedDT) <- Q_SNPs$SNP
Q_all<-rbindlist(nestedDT, use.names=T)
# remove duplicate SNPs
table(duplicated(Q_all$SNP))
Q_all<-subset(Q_all, !(duplicated(Q_all$SNP)))
fwrite(Q_all, file = "bifactor_CLMGWAS_F3_clumpedQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)
# remove clumped Q SNPs from main file and save

# Remove Q SNPs from data and save for post GWAS analyses where don't want to include Q signal in results
Data_noQ<-anti_join(Data, Q_all, by='SNP')
table(Data_noQ$Q_SNP_pval <= 1e-08) # check Q is definitely removed
# Bifactor: 2696844 SNPs

##Calculate Effective Sample Size for non Q sample
#restrict to MAF of 40% and 10%
CorrelatedFactors<-subset(Data_noQ, Data_noQ$MAF <= .4 & Data_noQ$MAF >= .1)

N_hat<-mean(1/((2*CorrelatedFactors$MAF*(1-CorrelatedFactors$MAF))*CorrelatedFactors$SE^2))
N_hat
# Bifactor = 53594

Data_noQ$N <- 53594
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 5e-8)
# Bifactor: 1718 SNPs

fwrite(Data_noQ, file = "bifactor_CLMGWAS_F3_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_F3_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_CLMGWAS_F3_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_CLMGWAS_F3_noQ_postGWAS.txt")

# MUNGE FORMATTED DATA FOR GSEM
files<-fread("bifactor_CLMGWAS_F3_noQ_postGWAS.txt.gz")

files <- files %>% select(CHR, BP, SNP, A1, A2, BETA, SE, P, MAF, N) # select columns needed for GSEM
names(files) <- c("CHR", "POS", "SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N") # rename selected columns for easy recognition in subsequent packages
head(files)

fwrite(files, file = "bifactor_CLMGWAS_F3_forsumstat.txt", 
       row.names = FALSE, quote = FALSE, col.names = TRUE,
       sep = " ")
gzip("bifactor_CLMGWAS_F3_forsumstat.txt")

files<-c("bifactor_CLMGWAS_F3_forsumstat.txt.gz")
trait.names<-c("F3")
#sample size: NA if SNP-specific Neff used, Neff if sum of cohorts Neffs used or binary single cohort; total N if one samPAe continuous 
N=c(53594)
#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

##########################################################################
setwd("../GSEM/non-UKB/Format_GWAS/F4_results/")

require(stringr)
require(data.table)
require(RColorBrewer)
require(lattice)
require(dplyr)
require(rlist)
library(R.utils)
#install.packages("qqman")
library(qqman)

# save .tab files from array jobs into separate directories for each latent factor
# Here we provide an example for the General Factor but we repeated for each of the other latent factors too
tempfreq<-list.files(pattern="*.tab")
tempfreq<-str_sort(tempfreq, numeric = TRUE)

Data<-rbindlist(lapply(tempfreq, fread, header=TRUE,fill=TRUE))
rm(tempfreq)
head(Data)

# QC
table(is.na(Data$SNP))
# none
table(duplicated(Data$SNP))
# none
table(Data$error)
# none
table(Data$warning)
#warn <- subset(Data, Data$warning != 0) # there were 290 SNPs with negative lv variances, 1 SNP with negative ov variance & 1 SNP with non positive definite cov matrix of latent vars

#rm(warn)
#rm(sub)
# there doesn't seem to be an obvious reason as to why SNPs are erroneous so we remove SNPs with warnings - 63392 SNPs in total removed
Data <- subset(Data, Data$warning == 0)
nrow(Data)
table(Data$Z_smooth)
# 0 SNP removed
Data <- subset(Data, Data$Z_smooth <=1.96)
nrow(Data)
# Number of final SNPs:
# 2070948 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 10182 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 4827 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_CLMGWAS_F4_sumstats.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_F4_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_CLMGWAS_F4_topQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

Data_short <- Data %>% select(SNP, CHR, BP) # subset to only req columns for clumping to save mem
## Calculate regions around Q SNPs to remove
# add two columns to the Q SNPlist that provide a range of +/- 1,000,000 base pairs around each of the position points to get the range of values to filter on for each Q hit
Q_SNPs$max <- Q_SNPs$BP + 1000000
Q_SNPs$min <- Q_SNPs$BP - 1000000
head(Q_SNPs)
Q_SNPs_max = as.list(Q_SNPs$max)
Q_SNPs_min = as.list(Q_SNPs$min)

nestedDT <- list()

for (i in 1:nrow(Q_SNPs)) {
  CHROM <- Q_SNPs$CHR[i]
  MAX <- Q_SNPs$max[i]
  MIN <- Q_SNPs$min[i]
  temp <- Data_short[CHR == CHROM & BP <= MAX & BP >= MIN]
  nestedDT[[i]] <- temp
}

names(nestedDT) <- Q_SNPs$SNP
Q_all<-rbindlist(nestedDT, use.names=T)
# remove duplicate SNPs
table(duplicated(Q_all$SNP))
Q_all<-subset(Q_all, !(duplicated(Q_all$SNP)))
fwrite(Q_all, file = "bifactor_CLMGWAS_F4_clumpedQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)
# remove clumped Q SNPs from main file and save


# Remove Q SNPs from data and save for post GWAS analyses where don't want to include Q signal in results
Data_noQ<-anti_join(Data, Q_all, by='SNP')
table(Data_noQ$Q_SNP_pval <= 1e-08) # check Q is definitely removed
# Bifactor: 1812165 SNPs

##Calculate Effective Sample Size for non Q sample
#restrict to MAF of 40% and 10%
CorrelatedFactors<-subset(Data_noQ, Data_noQ$MAF <= .4 & Data_noQ$MAF >= .1)

N_hat<-mean(1/((2*CorrelatedFactors$MAF*(1-CorrelatedFactors$MAF))*CorrelatedFactors$SE^2))
N_hat
# Bifactor = 1001422

Data_noQ$N <- 1001422
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 5e-8)
# Bifactor: 11579 SNPs

fwrite(Data_noQ, file = "bifactor_CLMGWAS_F4_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_CLMGWAS_F4_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_CLMGWAS_F4_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_CLMGWAS_F4_noQ_postGWAS.txt")

# MUNGE FORMATTED DATA FOR GSEM
files<-fread("bifactor_CLMGWAS_F4_noQ_postGWAS.txt.gz")

files <- files %>% select(CHR, BP, SNP, A1, A2, BETA, SE, P, MAF, N) # select columns needed for GSEM
names(files) <- c("CHR", "POS", "SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N") # rename selected columns for easy recognition in subsequent packages
head(files)

fwrite(files, file = "bifactor_CLMGWAS_F4_forsumstat.txt", 
       row.names = FALSE, quote = FALSE, col.names = TRUE,
       sep = " ")
gzip("bifactor_CLMGWAS_F4_forsumstat.txt")

files<-c("bifactor_CLMGWAS_F4_forsumstat.txt.gz")
trait.names<-c("F4")
#sample size: NA if SNP-specific Neff used, Neff if sum of cohorts Neffs used or binary single cohort; total N if one samPAe continuous 
N=c(1001422)
#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)



###################################################################################
setwd("../GSEM/non-UKB/Format_GWAS/")

traits<-c("BF.sumstats.gz","F1.sumstats.gz","F2.sumstats.gz","F3.sumstats.gz","F4.sumstats.gz")
trait.names<-c("BF","F1","F2","F3","F4")

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(NA,NA,NA,NA,NA)
#vector of population prevalences
population.prev<-c(NA,NA,NA,NA,NA)

#the folder of LD scores
ld<-"../GSEM/Format_GWAS/eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"../GSEM/Format_GWAS/eur_w_ld_chr/"

# 运行LDSC
LDSCoutput_pairs <- ldsc(traits = traits,
                         sample.prev = sample.prev,
                         population.prev = population.prev,
                         ld = ld, wld = wld,
                         trait.names = trait.names,
                         ldsc.log = "CLM_LDSCoutput_BF")  


#optional command to save the output as a .RData file for later use
save(LDSCoutput_pairs,file="CLM_LDSCoutput_BF.RData")
