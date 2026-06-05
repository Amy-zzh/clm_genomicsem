## Written by Jonny Flint

##This document provides code and methods used for analyzing polygenic risk scores (PRS) and structural equation modeling (SEM) within the context of genome-wide association studies (GWAS). The code includes steps for regression modeling, elastic net regression, quintile analysis, and cognitive change modeling.

# Packages and Data Loading
library(dplyr)
library(ggplot2)
library(tidyverse)
library(psych)
library(Rcpp)
library(MASS)
library(lavaan)
library(nnet)
library(glmnet)
library(boot)
library(broom)
library(survival)
library(pscl)
library(data.table)
library(caret)
###read data
CLMall=fread("../UKB_data/UKB_CLM_PRS.csv")

#Ensure that the input data ( in this case CLMall and other datasets) is cleaned and preprocessed according to the documented protocol and genotyping process.
#Regression Models for Univariate Analysis
#Creating a Multivariable Predictor
#The multivariable predictor is generated using the following steps:
  
CLMallmultinoCov <- lm(formula = CLM ~ z_bf_total + z_f1_total + 
                         z_f2_total + z_f3_total_r + z_f4_total, 
                          data = CLMall)
summary(CLMallmultinoCov)

CLMall$multi <- predict(CLMallmultinoCov, newdata = CLMall)
CLMall$multiscale <- scale(CLMall$multi)
#Regression Models with Covariates
modelmulticov <- lm(formula = CLM ~ age + sex + GPC1 + GPC2 + GPC3 + 
                      GPC4 + GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, 
                    data = CLMall)
summary(modelmulticov)

CLMAllmodelmulti <- lm(formula = CLM ~ multiscale + age + sex + GPC1 + 
                         GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, 
                       data = CLMall)
summary(CLMAllmodelmulti)
# Regression Models with Individual PRS Predictors
# The following loop iteratively fits models for each PRS predictor and extracts the results:
  
predictors <- c("z_bf_total", "z_f1_total", "z_f2_total", 
                  "z_f3_total_r", "z_f4_total", "multiscale")
#################CLM
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CLM ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CLM ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
               family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)
                                                
write.csv(combined_results, file = "../GSEM/Format_GWAS/CLM_regression_results.csv")




predictors <- c("z_bf_total", "multiscale")

#################CVD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CVD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CVD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                 GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
               family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/CVD_final_regression_results.csv")


#################MASLD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("MASLD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(MASLD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                 GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
               family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/MASLD_final_regression_results.csv")

#################AF_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("AF_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(AF_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                 GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
               family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/AF_final_regression_results.csv")


#################CAD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CAD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CAD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/CAD_final_regression_results.csv")


#################PAD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("PAD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(PAD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/PAD_final_regression_results.csv")

#################Stroke_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Stroke_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(Stroke_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/Stroke_final_regression_results.csv")

#################HF_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("HF_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(HF_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/HF_final_regression_results.csv")


#################VHD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("VHD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(VHD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/VHD_final_regression_results.csv")


#################HS_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("HS_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(HS_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/HS_final_regression_results.csv")

#################Dyslipid_final
CLMall$Dyslipid_final[CLMall$Dyslipid_final == 2] <- 1

fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Dyslipid_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(Dyslipid_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/Dyslipid_final_regression_results.csv")

#################T2D_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("T2D_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(T2D_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/T2D_final_regression_results.csv")

#################Obesity_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Obesity_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(Obesity_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/Obesity_final_regression_results.csv")



#################HTN_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("HTN_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(HTN_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/HTN_final_regression_results.csv")

fwrite(CLMall, file = "../UKB_data/UKB_CLM_PRS_addMulti.csv")








################Age<60
###read data
CLMall=fread("../UKB_data/UKB_CLM_PRS.csv")
CLMall <- CLMall[age < 60]

#Ensure that the input data ( in this case CLMall and other datasets) is cleaned and preprocessed according to the documented protocol and genotyping process.
#Regression Models for Univariate Analysis
#Creating a Multivariable Predictor
#The multivariable predictor is generated using the following steps:

CLMallmultinoCov <- lm(formula = CLM ~ z_bf_total + z_f1_total + 
                         z_f2_total + z_f3_total_r + z_f4_total, 
                       data = CLMall)
summary(CLMallmultinoCov)

CLMall$multi <- predict(CLMallmultinoCov, newdata = CLMall)
CLMall$multiscale <- scale(CLMall$multi)
#Regression Models with Covariates
modelmulticov <- lm(formula = CLM ~ age + sex + GPC1 + GPC2 + GPC3 + 
                      GPC4 + GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, 
                    data = CLMall)
summary(modelmulticov)

CLMAllmodelmulti <- lm(formula = CLM ~ multiscale + age + sex + GPC1 + 
                         GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, 
                       data = CLMall)
summary(CLMAllmodelmulti)
# Regression Models with Individual PRS Predictors
# The following loop iteratively fits models for each PRS predictor and extracts the results:

predictors <- c("z_bf_total", "z_f1_total", "z_f2_total", 
                "z_f3_total_r", "z_f4_total", "multiscale")
#################CLM
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CLM ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CLM ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/CLM_regression_results_under60.csv")




predictors <- c("z_bf_total", "multiscale")

#################CVD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CVD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CVD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/CVD_final_regression_results_under60.csv")


#################MASLD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("MASLD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(MASLD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/MASLD_final_regression_results_under60.csv")

#################AF_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("AF_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(AF_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/AF_final_regression_results_under60.csv")


#################CAD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CAD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CAD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/CAD_final_regression_results_under60.csv")


#################PAD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("PAD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(PAD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/PAD_final_regression_results_under60.csv")

#################Stroke_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Stroke_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(Stroke_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/Stroke_final_regression_results_under60.csv")

#################HF_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("HF_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(HF_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/HF_final_regression_results_under60.csv")


#################VHD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("VHD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(VHD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/VHD_final_regression_results_under60.csv")


#################HS_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("HS_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(HS_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/HS_final_regression_results_under60.csv")

#################Dyslipid_final
CLMall$Dyslipid_final[CLMall$Dyslipid_final == 2] <- 1

fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Dyslipid_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(Dyslipid_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/Dyslipid_final_regression_results_under60.csv")

#################T2D_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("T2D_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(T2D_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/T2D_final_regression_results_under60.csv")

#################Obesity_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Obesity_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(Obesity_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/Obesity_final_regression_results_under60.csv")



#################HTN_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("HTN_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(HTN_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/HTN_final_regression_results_under60.csv")

fwrite(CLMall, file = "../UKB_data/UKB_CLM_PRS_addMulti.csv")








################Age≥60
###read data
CLMall=fread("../UKB_data/UKB_CLM_PRS.csv")
CLMall <- CLMall[age >= 60]

#Ensure that the input data ( in this case CLMall and other datasets) is cleaned and preprocessed according to the documented protocol and genotyping process.
#Regression Models for Univariate Analysis
#Creating a Multivariable Predictor
#The multivariable predictor is generated using the following steps:

CLMallmultinoCov <- lm(formula = CLM ~ z_bf_total + z_f1_total + 
                         z_f2_total + z_f3_total_r + z_f4_total, 
                       data = CLMall)
summary(CLMallmultinoCov)

CLMall$multi <- predict(CLMallmultinoCov, newdata = CLMall)
CLMall$multiscale <- scale(CLMall$multi)
#Regression Models with Covariates
modelmulticov <- lm(formula = CLM ~ age + sex + GPC1 + GPC2 + GPC3 + 
                      GPC4 + GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, 
                    data = CLMall)
summary(modelmulticov)

CLMAllmodelmulti <- lm(formula = CLM ~ multiscale + age + sex + GPC1 + 
                         GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, 
                       data = CLMall)
summary(CLMAllmodelmulti)
# Regression Models with Individual PRS Predictors
# The following loop iteratively fits models for each PRS predictor and extracts the results:

predictors <- c("z_bf_total", "z_f1_total", "z_f2_total", 
                "z_f3_total_r", "z_f4_total", "multiscale")
#################CLM
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CLM ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CLM ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/CLM_regression_results_above60.csv")




predictors <- c("z_bf_total", "multiscale")

#################CVD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CVD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CVD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/CVD_final_regression_results_above60.csv")


#################MASLD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("MASLD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(MASLD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/MASLD_final_regression_results_above60.csv")

#################AF_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("AF_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(AF_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/AF_final_regression_results_above60.csv")


#################CAD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CAD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CAD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/CAD_final_regression_results_above60.csv")


#################PAD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("PAD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(PAD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/PAD_final_regression_results_above60.csv")

#################Stroke_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Stroke_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(Stroke_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/Stroke_final_regression_results_above60.csv")

#################HF_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("HF_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(HF_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/HF_final_regression_results_above60.csv")


#################VHD_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("VHD_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(VHD_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/VHD_final_regression_results_above60.csv")


#################HS_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("HS_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(HS_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/HS_final_regression_results_above60.csv")

#################Dyslipid_final
CLMall$Dyslipid_final[CLMall$Dyslipid_final == 2] <- 1

fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Dyslipid_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(Dyslipid_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/Dyslipid_final_regression_results_above60.csv")

#################T2D_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("T2D_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(T2D_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/T2D_final_regression_results_above60.csv")

#################Obesity_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Obesity_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(Obesity_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/Obesity_final_regression_results_above60.csv")



#################HTN_final
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("HTN_final ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(HTN_final ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/HTN_final_regression_results_above60.csv")

fwrite(CLMall, file = "../UKB_data/UKB_CLM_PRS_addMulti.csv")






#######################single disease
###read data
CLMall=fread("../UKB_data/UKB_CLM_PRS_marcov.csv")


predictors <- c("z_af_total", "z_vhd_total", "z_cad_total", 
                "z_pad_total", "z_stroke_total", "z_hf_total", "z_masld_total", "z_t2d_total")
#################CLM
fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("CLM ~", predictor, 
                              "+ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + GPC5 + GPC6 + 
                              GPC7 + GPC8 + GPC9 + GPC10"))
  model <- glm(formula, data = CLMall, 
               family = binomial(link = "logit"))
  r2_full <- pR2(model)["McFadden"]
  r2_base <- pR2(cov_model)["McFadden"]
  r_squared_perc <- (r2_full - r2_base) / r2_base * 100
  
  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_perc = r_squared_perc)
  
  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

CLMcov <- glm(CLM ~ sex + age + GPC1 + GPC2 + GPC3 + GPC4 + 
                GPC5 + GPC6 + GPC7 + GPC8 + GPC9 + GPC10, data = CLMall, 
              family = binomial(link = "logit"))

results_list <- lapply(predictors, fit_model_extract_results, cov_model = CLMcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

combined_results$OR=exp(combined_results$estimate)
combined_results$OR_LCI=exp(combined_results$`CI Lower`)
combined_results$OR_UCI=exp(combined_results$`CI Upper`)

write.csv(combined_results, file = "../GSEM/Format_GWAS/CLM_regression_results.csv")
