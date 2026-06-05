
setwd("../GSEM/non-UKB/Format_GWAS/")

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


load("CLM_LDSCoutput_ALLchrALLphenos.RData")

bifactor_model <- "BF=~NA*AF + HF + Stroke + PAD + CAD + VHD + NAFLD + ALT + AST + GGT + ALP + BMI + WC + BFP + T2D + HbA1c + HTN + nonHDL + HDL + logTG + CRP
F1=~NA*AF + HF + Stroke + PAD + CAD + VHD + HTN
F2=~NA*BMI + WC + BFP
F3=~NA*NAFLD + ALT + AST + GGT
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
F3 ~~ F4"


bifactor <- usermodel(LDSCoutput_all, estimation = "DWLS", model = bifactor_model, 
                              CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
bifactor

