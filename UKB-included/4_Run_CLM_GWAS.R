## Run the GWAS of the bifactor model on the Colorado Research Computing Cluster

#set the working directory
#loading necessary R packages
require(gdata)
require(lavaan)
require(doParallel)
require(stringr)
require(Matrix)
require(R.utils)
require(dplyr)
require(utils)

require(GenomicSEM)
require(data.table)

starttime=Sys.time()
paste("Starttime:")
starttime

#Sys.sleep(60)

paste("Loading data:")
Sys.time()-starttime

#source("userGWAS_sansMPI.R")

#load the LDSCoutput
load("CLM_LDSCoutput_ALL_21traits.RData")

paste("Taking diag")
Sys.time()-starttime

#load and subset the sumstats
sumstats<-fread("prepped_sumstats_CLM_21traits.txt")

#print message:
paste("Data loaded")
Sys.time()-starttime

#get environment variables from SLURM corresponding to array task ID and ntasks:
args <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[1])
#cores <- as.numeric(args[2])

#split sumstats into chunks:
snpj <- 1:nrow(sumstats)
chunks <- split(snpj, ceiling(seq_along(snpj)/4384))
chunk_idx <- job_id
chunk <- chunks[[chunk_idx]]


#write in the model
model <- "BF=~NA*AF + HF + Stroke + PAD + CAD + VHD + NAFLD2 + ALT + AST + GGT + ALP + BMI + WC + BFP + T2D + HbA1c + HTN + nonHDL + HDL + logTG + CRP
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
BF ~ SNP
F1 ~ SNP
F2 ~ SNP
F3 ~ SNP
F4 ~ SNP
AF ~~ b*AF
b > .001
HF ~~ c*HF
c > .001
Stroke ~~ d*Stroke
d > .001
PAD ~~ e*PAD
e > .001
CAD ~~ f*CAD
f > .001
VHD ~~ h*VHD
h > .001
NAFLD2 ~~ i*NAFLD2
i > .001
ALT ~~ j*ALT
j > .001
AST ~~ k*AST
k > .001
GGT ~~ l*GGT
l > .001
ALP ~~ m*ALP
m > .001
BMI ~~ n*BMI
n > .001
WC ~~ o*WC
o > .001
BFP ~~ p*BFP
p > .001
T2D ~~ q*T2D
q > .001
HbA1c ~~ r*HbA1c
r > .001
HTN ~~ s*HTN
s > .001
nonHDL ~~ t*nonHDL
t > .001
HDL ~~ u*HDL
u > .001
logTG ~~ v*logTG
v > .001
CRP ~~ w*CRP
w > .001"

sub=c("BF~SNP","F1~SNP","F2~SNP","F3~SNP","F4~SNP")

paste("Running GWAS")
Sys.time()-starttime

results<-userGWAS(covstruc=LDSCoutput_all, SNPs=sumstats[chunk,], model=model, parallel=T, sub=sub, cores=12, GC="none",smooth_check=TRUE, printwarn=TRUE,toler=1e-50,fix_measurement=TRUE,Q_SNP=TRUE)

paste("Writing GWAS output")
Sys.time()-starttime

fwrite(results[[1]], file = paste0("../BF_results/CLM21_bifactor_BF_", chunk_idx, ".tab"), sep="\t", row.names=F, quote=F)
fwrite(results[[2]], file = paste0("../F1_results/CLM21_bifactor_F1_", chunk_idx, ".tab"), sep="\t", row.names=F, quote=F)
fwrite(results[[3]], file = paste0("../F2_results/CLM21_bifactor_F2_", chunk_idx, ".tab"), sep="\t", row.names=F, quote=F)
fwrite(results[[4]], file = paste0("../F3_results/CLM21_bifactor_F3_", chunk_idx, ".tab"), sep="\t", row.names=F, quote=F)
fwrite(results[[5]], file = paste0("../F4_results/CLM21_bifactor_F4_", chunk_idx, ".tab"), sep="\t", row.names=F, quote=F)
