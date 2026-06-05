# bifactor frailty S-GSEM enrichment analysis to get results from Stratified-GSEM
library(GenomicSEM)
library(data.table)
library(dplyr)
setwd("../GSEM/Format_GWAS/")
load("bifactor_CLM_GSEM_sldsc.RData")
colnames(GSEM_sldsc$S)
for(i in 1:length(GSEM_sldsc$S)){
  colnames(GSEM_sldsc$S[[i]])<-c("AF","HF","Stroke","PAD","CAD","VHD","NAFLD2","ALT","AST","GGT","ALP","BMI","WC","BFP","T2D","HbA1c","HTN","nonHDL","HDL","logTG","CRP")
  print(i)
}

# had to constrain all indicators to have variance above 0 for the stratified model

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
HF ~~ b*HF
b > .001"

params<-c("BF~~BF","F1~~F1", "F2~~F2", "F3~~F3", "F4~~F4")

std.lv=TRUE

#estimate enrichment using the enrich function
enrich<-enrich(s_covstruc=GSEM_sldsc,model=model,params=params,std.lv=std.lv)
save(enrich, file= "bifactor_CLM_SGSEM_enrichment.Rdata")

# Tables
load("bifactor_CLM_SGSEM_enrichment.Rdata")

Cahoy_Names<-read.csv("../GSEM/Format_GWAS/Cahoy_1000Gv3_ldscores/Cahoy_1000Gv3_ldscores/Cahoy.ldcts.csv",header=F)
Chromatin_Names<-read.csv("../GSEM/Format_GWAS/Multi_tissue_chromatin_1000Gv3_ldscores/Multi_tissue_chromatin_1000Gv3_ldscores/Multi_tissue_chromatin.ldcts.csv",header=F)
Gene_Expr_Names<-read.csv("../GSEM/Format_GWAS/Multi_tissue_gene_expr_1000Gv3_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/Multi_tissue_gene_expr.ldcts.csv",header=F)



BF<-data.frame(enrich[[1]])
idx1 <- match(BF$Annotation, Chromatin_Names$V2)
BF$Annotation <- ifelse(
  is.na(idx1),
  BF$Annotation,                 # 未匹配 → 保持原来的
  Chromatin_Names$V1[idx1]        # 匹配 → 替换成 V1
)
idx2 <- match(BF$Annotation, Cahoy_Names$V2)
BF$Annotation <- ifelse(
  is.na(idx2),
  BF$Annotation,                 # 未匹配 → 保持原来的
  Cahoy_Names$V1[idx2]        # 匹配 → 替换成 V1
)
idx3 <- match(BF$Annotation, Gene_Expr_Names$V2)
BF$Annotation <- ifelse(
  is.na(idx3),
  BF$Annotation,                 # 未匹配 → 保持原来的
  Gene_Expr_Names$V1[idx3]        # 匹配 → 替换成 V1
)
View(BF)
BF$Annotation <- sub("L2$", "", BF$Annotation)
BF$FDR <- p.adjust(BF$Enrichment_p_value, method = "BH")
BF$FDR_sig <- ifelse(BF$FDR < 0.05, "yes", "no")
BF$BONF <- p.adjust(BF$Enrichment_p_value, method = "bonferroni")
BF$BONF_sig <- ifelse(BF$BONF < 0.05, "yes", "no")
BF <- BF[ !grepl("control", BF$Annotation, ignore.case = TRUE), ]
BF <- BF[ , !(names(BF) %in% c("lhs", "op", "rhs", "Error", "Warning")) ]
fwrite(BF,"BF_enrichment.csv")



F1<-data.frame(enrich[[2]])
idx1 <- match(F1$Annotation, Chromatin_Names$V2)
F1$Annotation <- ifelse(
  is.na(idx1),
  F1$Annotation,                 # 未匹配 → 保持原来的
  Chromatin_Names$V1[idx1]        # 匹配 → 替换成 V1
)
idx2 <- match(F1$Annotation, Cahoy_Names$V2)
F1$Annotation <- ifelse(
  is.na(idx2),
  F1$Annotation,                 # 未匹配 → 保持原来的
  Cahoy_Names$V1[idx2]        # 匹配 → 替换成 V1
)
idx3 <- match(F1$Annotation, Gene_Expr_Names$V2)
F1$Annotation <- ifelse(
  is.na(idx3),
  F1$Annotation,                 # 未匹配 → 保持原来的
  Gene_Expr_Names$V1[idx3]        # 匹配 → 替换成 V1
)
View(F1)
F1$Annotation <- sub("L2$", "", F1$Annotation)
F1$FDR <- p.adjust(F1$Enrichment_p_value, method = "BH")
F1$FDR_sig <- ifelse(F1$FDR < 0.05, "yes", "no")
F1$BONF <- p.adjust(F1$Enrichment_p_value, method = "bonferroni")
F1$BONF_sig <- ifelse(F1$BONF < 0.05, "yes", "no")
F1 <- F1[ !grepl("control", F1$Annotation, ignore.case = TRUE), ]
F1 <- F1[ , !(names(F1) %in% c("lhs", "op", "rhs", "Error", "Warning")) ]
fwrite(F1,"F1_enrichment.csv")



F2<-data.frame(enrich[[3]])
idx1 <- match(F2$Annotation, Chromatin_Names$V2)
F2$Annotation <- ifelse(
  is.na(idx1),
  F2$Annotation,                 # 未匹配 → 保持原来的
  Chromatin_Names$V1[idx1]        # 匹配 → 替换成 V1
)
idx2 <- match(F2$Annotation, Cahoy_Names$V2)
F2$Annotation <- ifelse(
  is.na(idx2),
  F2$Annotation,                 # 未匹配 → 保持原来的
  Cahoy_Names$V1[idx2]        # 匹配 → 替换成 V1
)
idx3 <- match(F2$Annotation, Gene_Expr_Names$V2)
F2$Annotation <- ifelse(
  is.na(idx3),
  F2$Annotation,                 # 未匹配 → 保持原来的
  Gene_Expr_Names$V1[idx3]        # 匹配 → 替换成 V1
)
View(F2)
F2$Annotation <- sub("L2$", "", F2$Annotation)
F2$FDR <- p.adjust(F2$Enrichment_p_value, method = "BH")
F2$FDR_sig <- ifelse(F2$FDR < 0.05, "yes", "no")
F2$BONF <- p.adjust(F2$Enrichment_p_value, method = "bonferroni")
F2$BONF_sig <- ifelse(F2$BONF < 0.05, "yes", "no")
F2 <- F2[ !grepl("control", F2$Annotation, ignore.case = TRUE), ]
F2 <- F2[ , !(names(F2) %in% c("lhs", "op", "rhs", "Error", "Warning")) ]
fwrite(F2,"F2_enrichment.csv")





F3<-data.frame(enrich[[4]])
idx1 <- match(F3$Annotation, Chromatin_Names$V2)
F3$Annotation <- ifelse(
  is.na(idx1),
  F3$Annotation,                 # 未匹配 → 保持原来的
  Chromatin_Names$V1[idx1]        # 匹配 → 替换成 V1
)
idx2 <- match(F3$Annotation, Cahoy_Names$V2)
F3$Annotation <- ifelse(
  is.na(idx2),
  F3$Annotation,                 # 未匹配 → 保持原来的
  Cahoy_Names$V1[idx2]        # 匹配 → 替换成 V1
)
idx3 <- match(F3$Annotation, Gene_Expr_Names$V2)
F3$Annotation <- ifelse(
  is.na(idx3),
  F3$Annotation,                 # 未匹配 → 保持原来的
  Gene_Expr_Names$V1[idx3]        # 匹配 → 替换成 V1
)
View(F3)
F3$Annotation <- sub("L2$", "", F3$Annotation)
F3$FDR <- p.adjust(F3$Enrichment_p_value, method = "BH")
F3$FDR_sig <- ifelse(F3$FDR < 0.05, "yes", "no")
F3$BONF <- p.adjust(F3$Enrichment_p_value, method = "bonferroni")
F3$BONF_sig <- ifelse(F3$BONF < 0.05, "yes", "no")
F3 <- F3[ !grepl("control", F3$Annotation, ignore.case = TRUE), ]
F3 <- F3[ , !(names(F3) %in% c("lhs", "op", "rhs", "Error", "Warning")) ]
fwrite(F3,"F3_enrichment.csv")




F4<-data.frame(enrich[[5]])
idx1 <- match(F4$Annotation, Chromatin_Names$V2)
F4$Annotation <- ifelse(
  is.na(idx1),
  F4$Annotation,                 # 未匹配 → 保持原来的
  Chromatin_Names$V1[idx1]        # 匹配 → 替换成 V1
)
idx2 <- match(F4$Annotation, Cahoy_Names$V2)
F4$Annotation <- ifelse(
  is.na(idx2),
  F4$Annotation,                 # 未匹配 → 保持原来的
  Cahoy_Names$V1[idx2]        # 匹配 → 替换成 V1
)
idx3 <- match(F4$Annotation, Gene_Expr_Names$V2)
F4$Annotation <- ifelse(
  is.na(idx3),
  F4$Annotation,                 # 未匹配 → 保持原来的
  Gene_Expr_Names$V1[idx3]        # 匹配 → 替换成 V1
)
View(F4)
F4$Annotation <- sub("L2$", "", F4$Annotation)
F4$FDR <- p.adjust(F4$Enrichment_p_value, method = "BH")
F4$FDR_sig <- ifelse(F4$FDR < 0.05, "yes", "no")
F4$BONF <- p.adjust(F4$Enrichment_p_value, method = "bonferroni")
F4$BONF_sig <- ifelse(F4$BONF < 0.05, "yes", "no")
F4 <- F4[ !grepl("control", F4$Annotation, ignore.case = TRUE), ]
F4 <- F4[ , !(names(F4) %in% c("lhs", "op", "rhs", "Error", "Warning")) ]
fwrite(F4,"F4_enrichment.csv")

#remove annotations with high levels of smoothing
#Data<-subset(Data, Data$Z_smooth < 1.96)

#add a column specifying factor name for plotting
BF$Factor <- "GF"
F1$Factor <- "F1"
F2$Factor <- "F2"
F3$Factor <- "F3"
F4$Factor <- "F4"

# make a combined table of all the annotations for Common Factor & residuals that pass the QC

head(BF)
head(F1)
head(F2)
head(F3)
head(F4)

BF_edit <- BF %>%
  dplyr::rename(
    BF_enrichment = Enrichment,
    BF_SE = Enrichment_SE,
    BF_P = Enrichment_p_value
  )


F1_edit <- F1 %>%
  dplyr::rename(
    F1_enrichment = Enrichment,
    F1_SE = Enrichment_SE,
    F1_P = Enrichment_p_value
  )


F2_edit <- F2 %>%
  dplyr::rename(
    F2_enrichment = Enrichment,
    F2_SE = Enrichment_SE,
    F2_P = Enrichment_p_value
  )


F3_edit <- F3 %>%
  dplyr::rename(
    F3_enrichment = Enrichment,
    F3_SE = Enrichment_SE,
    F3_P = Enrichment_p_value
  )

F4_edit <- F4 %>%
  dplyr::rename(
    F4_enrichment = Enrichment,
    F4_SE = Enrichment_SE,
    F4_P = Enrichment_p_value
  )
# cbind on Names & V1
full_table <- cbind(BF_edit, F1_edit, F2_edit, F3_edit, F4_edit)
head(full_table)
write.csv(full_table, file = "CLMSGSEM_Enrichmentresults_FullTable.csv")

#rbind factors 1 to 3 to make the data long
enrich_long <- rbind(BF,F1,F2,F3,F4)
head(enrich_long)
write.csv(enrich_long, file = "CLMSGSEM_Enrichmentresults_LongTable.csv")
enrich_group <- enrich_long

#Cardiovascular
data_Cardiovascular<-subset(enrich_group,Annotation %in% c("Artery_Aorta","Artery_Coronary","Artery_Tibial","Heart_Atrial_Appendage","Heart_Left_Ventricle"))
data_Cardiovascular$Group<- "Cardiovascular"
head(data_Cardiovascular)

data_Digestive<-subset(enrich_group,Annotation %in% c("Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Liver","Pancreas","Small_Intestine_Terminal_Ileum","Stomach"))
data_Digestive$Group<- "Digestive"
head(data_Digestive)

data_Endocrine<-subset(enrich_group,Annotation %in% c("Minor_Salivary_Gland","Adipose_Subcutaneous","Adipose_Visceral_(Omentum)","Adrenal_Gland","Breast_Mammary_Tissue","Ovary","Pituitary","Prostate","Testis","Thyroid"))
data_Endocrine$Group<- "Exo-/Endocrine"
head(data_Endocrine)

data_Immune<-subset(enrich_group,Annotation %in% c("Cells_EBV-transformed_lymphocytes","Spleen","Whole_Blood"))
data_Immune$Group<- "Hemic and Immune"
head(data_Immune)

data_Integumentay<-subset(enrich_group,Annotation %in% c("Cells_Transformed_fibroblasts","Skin_Not_Sun_Exposed_(Suprapubic)","Skin_Sun_Exposed_(Lower_leg)"))
data_Integumentay$Group<- "Integumentay"
head(data_Integumentay)

data_Musculoskeletal<-subset(enrich_group,Annotation %in% c("Muscle_Skeletal"))
data_Musculoskeletal$Group<- "Musculoskeletal"
head(data_Musculoskeletal)

data_Nervous<-subset(enrich_group,Annotation %in% c("Brain_Amygdala","Brain_Anterior_cingulate_cortex_(BA24)","Brain_Caudate_(basal_ganglia)","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_(BA9)","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_(basal_ganglia)","Brain_Putamen_(basal_ganglia)","Brain_Spinal_cord_(cervical_c-1)","Brain_Substantia_nigra","Nerve_Tibial"))
data_Nervous$Group<- "Nervous"
head(data_Nervous)

data_Respiratory<-subset(enrich_group,Annotation %in% c("Lung"))
data_Respiratory$Group<- "Respiratory"
head(data_Respiratory)

data_Urogenital<-subset(enrich_group,Annotation %in% c("Bladder","Uterus","Vagina","Cervix_Ectocervix","Cervix_Endocervix","Fallopian_Tube","Kidney_Cortex"))
data_Urogenital$Group<- "Urogenital"
head(data_Urogenital)

data<-rbind(data_Cardiovascular,data_Digestive,data_Endocrine,data_Immune,data_Integumentay,data_Musculoskeletal,data_Nervous,data_Respiratory,data_Urogenital)
write.csv(data, file = "CLMSGSEM_all_forplot.csv")