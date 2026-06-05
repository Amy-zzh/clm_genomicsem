setwd("../GSEM/Format_GWAS/BF_results/")

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
Data <- subset(Data, Data$Z_smooth == 0)
# Number of final SNPs:
# 4319864 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 225729 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 31330 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_frailtyGWAS_BF_sumstats.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_BF_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_frailtyGWAS_BF_topQsnps.txt", 
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
fwrite(Q_all, file = "bifactor_frailtyGWAS_BF_clumpedQsnps.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
# remove clumped Q SNPs from main file and save
# Bifactor: removed 2,556,852 SNPs

# Remove Q SNPs from data and save for post GWAS analyses where don't want to include Q signal in results
Data_noQ<-anti_join(Data, Q_all, by='SNP')
table(Data_noQ$Q_SNP_pval <= 1e-08) # check Q is definitely removed

##Calculate Effective Sample Size for non Q sample
#restrict to MAF of 40% and 10%
CorrelatedFactors<-subset(Data_noQ, Data_noQ$MAF <= .4 & Data_noQ$MAF >= .1)

N_hat<-mean(1/((2*CorrelatedFactors$MAF*(1-CorrelatedFactors$MAF))*CorrelatedFactors$SE^2))
N_hat
# Bifactor = 2762458

Data_noQ$N <- 2762458
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 5e-8)
# Bifactor: 3235 SNPs

fwrite(Data_noQ, file = "bifactor_frailtyGWAS_BF_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_BF_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_frailtyGWAS_BF_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_frailtyGWAS_BF_noQ_postGWAS.txt")


##########################################################################
setwd("../GSEM/Format_GWAS/F1_results/")

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
Data <- subset(Data, Data$Z_smooth == 0)
nrow(Data)
# Number of final SNPs:
# 4319864 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 26678 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 11487 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_frailtyGWAS_F1_sumstats.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_F1_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_frailtyGWAS_F1_topQsnps.txt", 
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
fwrite(Q_all, file = "bifactor_frailtyGWAS_F1_clumpedQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)
# remove clumped Q SNPs from main file and save


# Remove Q SNPs from data and save for post GWAS analyses where don't want to include Q signal in results
Data_noQ<-anti_join(Data, Q_all, by='SNP')
table(Data_noQ$Q_SNP_pval <= 1e-08) # check Q is definitely removed
# Bifactor: 3357471 SNPs

##Calculate Effective Sample Size for non Q sample
#restrict to MAF of 40% and 10%
CorrelatedFactors<-subset(Data_noQ, Data_noQ$MAF <= .4 & Data_noQ$MAF >= .1)

N_hat<-mean(1/((2*CorrelatedFactors$MAF*(1-CorrelatedFactors$MAF))*CorrelatedFactors$SE^2))
N_hat
# Bifactor = 784730

Data_noQ$N <- 784730
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 5e-8)
# Bifactor: 11579 SNPs

fwrite(Data_noQ, file = "bifactor_frailtyGWAS_F1_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_F1_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_frailtyGWAS_F1_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_frailtyGWAS_F1_noQ_postGWAS.txt")

##########################################################################
setwd("../GSEM/Format_GWAS/F2_results/")

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
Data <- subset(Data, Data$Z_smooth == 0)
nrow(Data)
# Number of final SNPs:
# 4319864 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 42113 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 5186 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_frailtyGWAS_F2_sumstats.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_F2_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_frailtyGWAS_F2_topQsnps.txt", 
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
fwrite(Q_all, file = "bifactor_frailtyGWAS_F2_clumpedQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)
# remove clumped Q SNPs from main file and save


# Remove Q SNPs from data and save for post GWAS analyses where don't want to include Q signal in results
Data_noQ<-anti_join(Data, Q_all, by='SNP')
table(Data_noQ$Q_SNP_pval <= 1e-08) # check Q is definitely removed
# Bifactor: 2840479 SNPs

##Calculate Effective Sample Size for non Q sample
#restrict to MAF of 40% and 10%
CorrelatedFactors<-subset(Data_noQ, Data_noQ$MAF <= .4 & Data_noQ$MAF >= .1)

N_hat<-mean(1/((2*CorrelatedFactors$MAF*(1-CorrelatedFactors$MAF))*CorrelatedFactors$SE^2))
N_hat
# Bifactor = 582909

Data_noQ$N <- 582909
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 5e-8)
# Bifactor: 9292 SNPs

fwrite(Data_noQ, file = "bifactor_frailtyGWAS_F2_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_F2_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_frailtyGWAS_F2_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_frailtyGWAS_F2_noQ_postGWAS.txt")

##########################################################################
setwd("../GSEM/Format_GWAS/F3_results/")

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
Data <- subset(Data, Data$Z_smooth == 0)
nrow(Data)
# Number of final SNPs:
# 4319864 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 66464 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 12529 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_frailtyGWAS_F3_sumstats.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_F3_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_frailtyGWAS_F3_topQsnps.txt", 
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
fwrite(Q_all, file = "bifactor_frailtyGWAS_F3_clumpedQsnps.txt", 
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
# Bifactor = 614401

Data_noQ$N <- 614401
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 5e-8)
# Bifactor: 1718 SNPs

fwrite(Data_noQ, file = "bifactor_frailtyGWAS_F3_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_F3_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_frailtyGWAS_F3_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_frailtyGWAS_F3_noQ_postGWAS.txt")

##########################################################################
setwd("../GSEM/Format_GWAS/F4_results/")

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
Data <- subset(Data, Data$Z_smooth == 0)
nrow(Data)
# Number of final SNPs:
# 4319864 SNPs

# Assess Q 
Qsig <- 5e-8/5
table(Data$Q_SNP_pval <= 1e-08)
# Bifactor: 32595 SNPs

# How many of the sig SNPs are Qsig too?
sig <- subset(Data, Data$Pval_Estimate <= 5e-8)
table(sig$Q_SNP_pval <= 1e-08)
# Bifactor = 14340 SNPs

#save file of full sumstats (i.e. still with Q SNPs in)
fwrite(Data, file = "bifactor_frailtyGWAS_F4_sumstats.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_F4_sumstats.txt")

# clump out SNPs around Qsig portions
Q_SNPs <- subset(Data, Data$Q_SNP_pval <= 1e-08, select=c(SNP, CHR, BP)) # make list of Qsig SNPs and positional info
head(Q_SNPs)
fwrite(Q_SNPs, file = "bifactor_frailtyGWAS_F4_topQsnps.txt", 
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
fwrite(Q_all, file = "bifactor_frailtyGWAS_F4_clumpedQsnps.txt", 
       row.names = FALSE, quote = FALSE, col.names = T)
# remove clumped Q SNPs from main file and save


# Remove Q SNPs from data and save for post GWAS analyses where don't want to include Q signal in results
Data_noQ<-anti_join(Data, Q_all, by='SNP')
table(Data_noQ$Q_SNP_pval <= 1e-08) # check Q is definitely removed
# Bifactor: 3485811 SNPs

##Calculate Effective Sample Size for non Q sample
#restrict to MAF of 40% and 10%
CorrelatedFactors<-subset(Data_noQ, Data_noQ$MAF <= .4 & Data_noQ$MAF >= .1)

N_hat<-mean(1/((2*CorrelatedFactors$MAF*(1-CorrelatedFactors$MAF))*CorrelatedFactors$SE^2))
N_hat
# Bifactor = 1418556

Data_noQ$N <- 1418556
head(Data_noQ)
rm(CorrelatedFactors)

# Assess p values
range(Data_noQ$Pval_Estimate,na.rm=T)
table(Data_noQ$Pval_Estimate <= 5e-8)
# Bifactor: 11579 SNPs

fwrite(Data_noQ, file = "bifactor_frailtyGWAS_F4_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)

gzip("bifactor_frailtyGWAS_F4_sumstats_noQ.txt")

# save GWAS sumstats for PRS without Q
head(Data_noQ)
Data_postGWAS <- Data_noQ %>% select(SNP, CHR, BP, MAF, A1, A2, est, SE, Pval_Estimate, N)
names(Data_postGWAS) <- c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","P","N")

fwrite(Data_postGWAS, file = "bifactor_frailtyGWAS_F4_noQ_postGWAS.txt", 
            row.names = FALSE, quote = FALSE, col.names = T)
gzip("bifactor_frailtyGWAS_F4_noQ_postGWAS.txt")