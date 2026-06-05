#### Gene mapping & prioritisation script

### SMR results
library(tidyverse)
library(data.table)
library(stringr)
library(dplyr)
library(stats)
library(UpSetR)
library(Hmisc)
library(tidyr)

### Read in the sig SMR results and make gene lists for each latent factor

## BF has 4015 sig SMR; F1 = 2171 sig; F2 = 4036; F3 = 22845 sig; F4 = 8266 sig; F5 = 13772 sig; F6 = 2950 sig

SMR <- fread("../GSEM/SMR/Results/smr_sig_hits_multi.csv")
head(SMR)

table(SMR$Factor)

# Make a subset for each latent factor, remove NAs in gene column and remove gene duplicates
SMR_BF <- subset(SMR, Factor == "BF", select=c(Gene))
SMR_BF <- subset(SMR_BF, !is.na(SMR_BF$Gene))
SMR_BF <- subset(SMR_BF, !(duplicated(SMR_BF$Gene))) 
nrow(SMR_BF) # 492 genes
SMR_BF$SMR <- 1

SMR_F1 <- subset(SMR, Factor == "F1", select=c(Gene))
SMR_F1 <- subset(SMR_F1, !is.na(SMR_F1$Gene))
SMR_F1 <- subset(SMR_F1, !(duplicated(SMR_F1$Gene)))
nrow(SMR_F1) # 2386 genes
SMR_F1$SMR <- 1

SMR_F2 <- subset(SMR, Factor == "F2", select=c(Gene))
SMR_F2 <- subset(SMR_F2, !is.na(SMR_F2$Gene))
SMR_F2 <- subset(SMR_F2, !(duplicated(SMR_F2$Gene)))
nrow(SMR_F2) # 1856 genes
SMR_F2$SMR <- 1

SMR_F3 <- subset(SMR, Factor == "F3", select=c(Gene))
SMR_F3 <- subset(SMR_F3, !is.na(SMR_F3$Gene))
SMR_F3 <- subset(SMR_F3, !(duplicated(SMR_F3$Gene)))
nrow(SMR_F3) # 850 genes
SMR_F3$SMR <- 1

SMR_F4 <- subset(SMR, Factor == "F4", select=c(Gene))
SMR_F4 <- subset(SMR_F4, !is.na(SMR_F4$Gene))
SMR_F4 <- subset(SMR_F4, !(duplicated(SMR_F4$Gene)))
nrow(SMR_F4) # 2091 genes
SMR_F4$SMR <- 1


rm(SMR)

## Get MAGMA results
MAGMA_BF <- fread("../GSEM/BF_FUMA/BF_FUMA/magma.genes.out")
nrow(MAGMA_BF)
MAGMA_BF <- subset(MAGMA_BF, P <= 0.05/2504, select = c(SYMBOL)) 
MAGMA_BF$MAGMA <- 1
names(MAGMA_BF) <- c("Gene", "MAGMA")
MAGMA_BF <- subset(MAGMA_BF, !is.na(MAGMA_BF$Gene))
MAGMA_BF <- subset(MAGMA_BF, !(duplicated(MAGMA_BF$Gene)))
nrow(MAGMA_BF) # 136 genes

MAGMA_F1 <- fread("../GSEM/F1_FUMA/F1_FUMA/magma.genes.out")
nrow(MAGMA_F1)
MAGMA_F1 <- subset(MAGMA_F1, P <= 0.05/12470, select = c(SYMBOL)) 
MAGMA_F1$MAGMA <- 1
names(MAGMA_F1) <- c("Gene", "MAGMA")
MAGMA_F1 <- subset(MAGMA_F1, !is.na(MAGMA_F1$Gene))
MAGMA_F1 <- subset(MAGMA_F1, !(duplicated(MAGMA_F1$Gene)))
nrow(MAGMA_F1) # 473 genes

MAGMA_F2 <- fread("../GSEM/F2_FUMA/F2_FUMA/magma.genes.out")
nrow(MAGMA_F2)
MAGMA_F2 <- subset(MAGMA_F2, P <= 0.05/9879, select = c(SYMBOL)) 
MAGMA_F2$MAGMA <- 1
names(MAGMA_F2) <- c("Gene", "MAGMA")
MAGMA_F2 <- subset(MAGMA_F2, !is.na(MAGMA_F2$Gene))
MAGMA_F2 <- subset(MAGMA_F2, !(duplicated(MAGMA_F2$Gene)))
nrow(MAGMA_F2) # 499 genes

MAGMA_F3 <- fread("../GSEM/F3_FUMA/F3_FUMA/magma.genes.out")
nrow(MAGMA_F3)
MAGMA_F3 <- subset(MAGMA_F3, P <= 0.05/7948, select = c(SYMBOL)) 
MAGMA_F3$MAGMA <- 1
names(MAGMA_F3) <- c("Gene", "MAGMA")
MAGMA_F3 <- subset(MAGMA_F3, !is.na(MAGMA_F3$Gene))
MAGMA_F3 <- subset(MAGMA_F3, !(duplicated(MAGMA_F3$Gene)))
nrow(MAGMA_F3) # 157 genes

MAGMA_F4 <- fread("../GSEM/F4_FUMA/F4_FUMA/magma.genes.out")
nrow(MAGMA_F4)
MAGMA_F4 <- subset(MAGMA_F4, P <= 0.05/11827, select = c(SYMBOL)) 
MAGMA_F4$MAGMA <- 1
names(MAGMA_F4) <- c("Gene", "MAGMA")
MAGMA_F4 <- subset(MAGMA_F4, !is.na(MAGMA_F4$Gene))
MAGMA_F4 <- subset(MAGMA_F4, !(duplicated(MAGMA_F4$Gene)))
nrow(MAGMA_F4) # 589 genes

# Position, eQTL & CI mapping
FUMA_BF <- fread("../GSEM/BF_FUMA/BF_FUMA/genes.txt")
head(FUMA_BF)
POS_BF <- subset(FUMA_BF, posMapSNPs != 0, select = c(symbol, posMapSNPs))
names(POS_BF) <- c("Gene", "Position")
POS_BF$Position <- 1
POS_BF <- subset(POS_BF, !is.na(POS_BF$Gene))
POS_BF <- subset(POS_BF, !(duplicated(POS_BF$Gene)))
nrow(POS_BF) # 54 genes

QTL_BF <- subset(FUMA_BF, eqtlMapSNPs != 0, select = c(symbol, eqtlMapSNPs))
names(QTL_BF) <- c("Gene", "eQTL")
QTL_BF$eQTL <- 1
QTL_BF <- subset(QTL_BF, !is.na(QTL_BF$Gene))
QTL_BF <- subset(QTL_BF, !(duplicated(QTL_BF$Gene)))
nrow(QTL_BF) # 129 genes

CI_BF <- subset(FUMA_BF, ciMap == "Yes", select = c(symbol, ciMap))
names(CI_BF) <- c("Gene", "CI")
CI_BF$CI <- 1
CI_BF <- subset(CI_BF, !is.na(CI_BF$Gene))
CI_BF <- subset(CI_BF, !(duplicated(CI_BF$Gene)))
nrow(CI_BF) # 380 genes


FUMA_F1 <- fread("../GSEM/F1_FUMA/F1_FUMA/genes.txt")
head(FUMA_F1)
POS_F1 <- subset(FUMA_F1, posMapSNPs != 0, select = c(symbol, posMapSNPs))
names(POS_F1) <- c("Gene", "Position")
POS_F1$Position <- 1
POS_F1 <- subset(POS_F1, !is.na(POS_F1$Gene))
POS_F1 <- subset(POS_F1, !(duplicated(POS_F1$Gene)))
nrow(POS_F1) # 417 genes

QTL_F1 <- subset(FUMA_F1, eqtlMapSNPs != 0, select = c(symbol, eqtlMapSNPs))
names(QTL_F1) <- c("Gene", "eQTL")
QTL_F1$eQTL <- 1
QTL_F1 <- subset(QTL_F1, !is.na(QTL_F1$Gene))
QTL_F1 <- subset(QTL_F1, !(duplicated(QTL_F1$Gene)))
nrow(QTL_F1) # 1434 genes

CI_F1 <- subset(FUMA_F1, ciMap == "Yes", select = c(symbol, ciMap))
names(CI_F1) <- c("Gene", "CI")
CI_F1$CI <- 1
CI_F1 <- subset(CI_F1, !is.na(CI_F1$Gene))
CI_F1 <- subset(CI_F1, !(duplicated(CI_F1$Gene)))
nrow(CI_F1) # 1588 genes



FUMA_F2 <- fread("../GSEM/F2_FUMA/F2_FUMA/genes.txt")
head(FUMA_F2)
POS_F2 <- subset(FUMA_F2, posMapSNPs != 0, select = c(symbol, posMapSNPs))
names(POS_F2) <- c("Gene", "Position")
POS_F2$Position <- 1
POS_F2 <- subset(POS_F2, !is.na(POS_F2$Gene))
POS_F2 <- subset(POS_F2, !(duplicated(POS_F2$Gene)))
nrow(POS_F2) # 439 genes

QTL_F2 <- subset(FUMA_F2, eqtlMapSNPs != 0, select = c(symbol, eqtlMapSNPs))
names(QTL_F2) <- c("Gene", "eQTL")
QTL_F2$eQTL <- 1
QTL_F2 <- subset(QTL_F2, !is.na(QTL_F2$Gene))
QTL_F2 <- subset(QTL_F2, !(duplicated(QTL_F2$Gene)))
nrow(QTL_F2) # 1475 genes

CI_F2 <- subset(FUMA_F2, ciMap == "Yes", select = c(symbol, ciMap))
names(CI_F2) <- c("Gene", "CI")
CI_F2$CI <- 1
CI_F2 <- subset(CI_F2, !is.na(CI_F2$Gene))
CI_F2 <- subset(CI_F2, !(duplicated(CI_F2$Gene)))
nrow(CI_F2) # 1563 genes


FUMA_F3 <- fread("../GSEM/F3_FUMA/F3_FUMA/genes.txt")
head(FUMA_F3)
POS_F3 <- subset(FUMA_F3, posMapSNPs != 0, select = c(symbol, posMapSNPs))
names(POS_F3) <- c("Gene", "Position")
POS_F3$Position <- 1
POS_F3 <- subset(POS_F3, !is.na(POS_F3$Gene))
POS_F3 <- subset(POS_F3, !(duplicated(POS_F3$Gene)))
nrow(POS_F3) # 94 genes

QTL_F3 <- subset(FUMA_F3, eqtlMapSNPs != 0, select = c(symbol, eqtlMapSNPs))
names(QTL_F3) <- c("Gene", "eQTL")
QTL_F3$eQTL <- 1
QTL_F3 <- subset(QTL_F3, !is.na(QTL_F3$Gene))
QTL_F3 <- subset(QTL_F3, !(duplicated(QTL_F3$Gene)))
nrow(QTL_F3) # 249 genes

CI_F3 <- subset(FUMA_F3, ciMap == "Yes", select = c(symbol, ciMap))
names(CI_F3) <- c("Gene", "CI")
CI_F3$CI <- 1
CI_F3 <- subset(CI_F3, !is.na(CI_F3$Gene))
CI_F3 <- subset(CI_F3, !(duplicated(CI_F3$Gene)))
nrow(CI_F3) # 428 genes



FUMA_F4 <- fread("../GSEM/F4_FUMA/F4_FUMA/genes.txt")
head(FUMA_F4)
POS_F4 <- subset(FUMA_F4, posMapSNPs != 0, select = c(symbol, posMapSNPs))
names(POS_F4) <- c("Gene", "Position")
POS_F4$Position <- 1
POS_F4 <- subset(POS_F4, !is.na(POS_F4$Gene))
POS_F4 <- subset(POS_F4, !(duplicated(POS_F4$Gene)))
nrow(POS_F4) # 543 genes

QTL_F4 <- subset(FUMA_F4, eqtlMapSNPs != 0, select = c(symbol, eqtlMapSNPs))
names(QTL_F4) <- c("Gene", "eQTL")
QTL_F4$eQTL <- 1
QTL_F4 <- subset(QTL_F4, !is.na(QTL_F4$Gene))
QTL_F4 <- subset(QTL_F4, !(duplicated(QTL_F4$Gene)))
nrow(QTL_F4) # 1370 genes

CI_F4 <- subset(FUMA_F4, ciMap == "Yes", select = c(symbol, ciMap))
names(CI_F4) <- c("Gene", "CI")
CI_F4$CI <- 1
CI_F4 <- subset(CI_F4, !is.na(CI_F4$Gene))
CI_F4 <- subset(CI_F4, !(duplicated(CI_F4$Gene)))
nrow(CI_F4) # 2043 genes



# Make full joined tables for each factor
setwd("../GSEM/Format_GWAS")

BF_all <- SMR_BF %>%
  full_join(MAGMA_BF, by="Gene") %>%
  full_join(POS_BF, by="Gene") %>%
  full_join(QTL_BF, by="Gene") %>%
  full_join(CI_BF, by="Gene")

BF_all$CI <- as.numeric(BF_all$CI)
BF_all <- BF_all %>% replace_na(list(SMR = 0, MAGMA = 0, Position = 0, eQTL = 0, CI = 0))
nrow(BF_all) # 886 genes

write_csv(BF_all, "CLM_BF_mapped_genes_all_methods.csv")

BF_all$TotalMappings <- BF_all$SMR + BF_all$MAGMA + BF_all$Position + BF_all$eQTL + BF_all$CI
BF_prioritised <- subset(BF_all, TotalMappings >= 3)
nrow(BF_prioritised) # 67 genes

write_csv(BF_prioritised, "CLM_BF_prioritised_genelist.csv")


# Make full joined tables for each factor
F1_all <- SMR_F1 %>%
  full_join(MAGMA_F1, by="Gene") %>%
  full_join(POS_F1, by="Gene") %>%
  full_join(QTL_F1, by="Gene") %>%
  full_join(CI_F1, by="Gene")

F1_all$CI <- as.numeric(F1_all$CI)
F1_all <- F1_all %>% replace_na(list(SMR = 0, MAGMA = 0, Position = 0, eQTL = 0, CI = 0))
nrow(F1_all) # 4217 genes

write_csv(F1_all, "CLM_F1_mapped_genes_all_methods.csv")

F1_all$TotalMappings <- F1_all$SMR + F1_all$MAGMA + F1_all$Position + F1_all$eQTL + F1_all$CI
F1_prioritised <- subset(F1_all, TotalMappings >= 3) 
nrow(F1_prioritised) # 514 genes

write_csv(F1_prioritised, "CLM_F1_prioritised_genelist.csv")

# Make full joined tables for each factor
F2_all <- SMR_F2 %>%
  full_join(MAGMA_F2, by="Gene") %>%
  full_join(POS_F2, by="Gene") %>%
  full_join(QTL_F2, by="Gene") %>%
  full_join(CI_F2, by="Gene")

F2_all$CI <- as.numeric(F2_all$CI)
F2_all <- F2_all %>% replace_na(list(SMR = 0, MAGMA = 0, Position = 0, eQTL = 0, CI = 0))
nrow(F2_all) # 3809 genes

write_csv(F2_all, "CLM_F2_mapped_genes_all_methods.csv")

F2_all$TotalMappings <- F2_all$SMR + F2_all$MAGMA + F2_all$Position + F2_all$eQTL + F2_all$CI
F2_prioritised <- subset(F2_all, TotalMappings >= 3) # 532 genes
nrow(F2_prioritised)

write_csv(F2_prioritised, "CLM_F2_prioritised_genelist.csv")


# Make full joined tables for each factor
F3_all <- SMR_F3 %>%
  full_join(MAGMA_F3, by="Gene") %>%
  full_join(POS_F3, by="Gene") %>%
  full_join(QTL_F3, by="Gene") %>%
  full_join(CI_F3, by="Gene")

F3_all$CI <- as.numeric(F3_all$CI)
F3_all <- F3_all %>% replace_na(list(SMR = 0, MAGMA = 0, Position = 0, eQTL = 0, CI = 0))
nrow(F3_all) # 1279 genes

write_csv(F3_all, "CLM_F3_mapped_genes_all_methods.csv")

F3_all$TotalMappings <- F3_all$SMR + F3_all$MAGMA + F3_all$Position + F3_all$eQTL + F3_all$CI
F3_prioritised <- subset(F3_all, TotalMappings >= 3) 
nrow(F3_prioritised) # 116 genes

write_csv(F3_prioritised, "CLM_F3_prioritised_genelist.csv")


# Make full joined tables for each factor
F4_all <- SMR_F4 %>%
  full_join(MAGMA_F4, by="Gene") %>%
  full_join(POS_F4, by="Gene") %>%
  full_join(QTL_F4, by="Gene") %>%
  full_join(CI_F4, by="Gene")

F4_all$CI <- as.numeric(F4_all$CI)
F4_all <- F4_all %>% replace_na(list(SMR = 0, MAGMA = 0, Position = 0, eQTL = 0, CI = 0))
nrow(F4_all) # 4109 genes

write_csv(F4_all, "CLM_F4_mapped_genes_all_methods.csv")

F4_all$TotalMappings <- F4_all$SMR + F4_all$MAGMA + F4_all$Position + F4_all$eQTL + F4_all$CI
F4_prioritised <- subset(F4_all, TotalMappings >= 3)
nrow(F4_prioritised) # 650 genes

write_csv(F4_prioritised, "CLM_F4_prioritised_genelist.csv")




