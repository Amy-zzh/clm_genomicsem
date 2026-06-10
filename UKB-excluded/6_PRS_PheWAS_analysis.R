library(data.table)
df1 <- fread("../UKB_data/UKB_CLM_forPheWAS.csv")
head(df1)
df2 <- fread("../UKB_data/phenotype_phecode.csv")
df2[] <- lapply(df2, function(x) if(is.logical(x)) as.integer(x) else x)
names(df2) <- paste0("phecode_", gsub("\\.", "", names(df2)))
names(df2)[1] <- "id"

df1$id <- as.integer(df1$id)
df2$id <- as.integer(df2$id)
data <- merge(df1, df2, by.x='id', by.y='id', all.x=TRUE, all.y=TRUE)

phenotypes <- colnames(data)[grepl("phecode_", colnames(data))]
n_phecode <- length(phenotypes)

phenotypes_keep <- phenotypes[sapply(data[, phenotypes, with=FALSE], function(x) {
  n_case <- sum(x == 1, na.rm = TRUE)
  n_control <- sum(x == 0, na.rm = TRUE)
  n_case >= 200 & n_control >= 1
})]
n_phecode_keep <- length(phecodes_keep)


z_95 <- qnorm(0.975,mean=0,sd=1)
idx <- 1
results <- NULL
for(phecode in phecodes_keep){
  
  phe_subset <- data[!is.na(data[[phecode]]), ]
  if(length(unique(phe_subset[[phecode]])) < 2) next
  covariates <- c("clm_prs", "age", paste0("GPC",1:10))
  if("sex" %in% colnames(phe_subset) && length(unique(na.omit(phe_subset$sex))) >= 2){
    covariates <- c("sex", covariates)
  }
  
  prsmodel_syntax <- as.formula(paste0(phecode, " ~ ", paste(covariates, collapse=" + ")))
  
  # glm 回归
  prs_model <- glm(prsmodel_syntax, family=binomial(), data = phe_subset, na.action = na.omit)
  n_observation <- nobs(prs_model)
  beta <- as.numeric(coef(summary(prs_model))[,'Estimate'][2])
  SE <- as.numeric(coef(summary(prs_model))[,'Std. Error'][2])
  beta_ci_lower <- beta - (z_95 * SE)
  beta_ci_upper <- beta + (z_95 * SE)
  Z <- as.numeric(coef(summary(prs_model))[,'z value'][2])
  P <- as.numeric(coef(summary(prs_model))[,'Pr(>|z|)'][2])
  OR <- exp(beta)
  OR_ci_lower <- exp(beta_ci_lower)
  OR_ci_upper <- exp(beta_ci_upper)
  
  result <- cbind(phecode, n_observation, beta, SE, beta_ci_lower, beta_ci_upper, Z, P, OR, OR_ci_lower, OR_ci_upper)
  
  idx <- idx + 1
  cat(paste0("[", idx, "/", n_phecode_keep, "]\n"))
  cat(result)
  cat('\n')
  results <- rbind(results, result)
}


results <- as.data.frame(results)
results <- setNames(results, c("Phecode", "n_obs", "Beta", "SE", "Beta_CI_lower", "Beta_CI_upper", "Z-score", "P-value", 
                               "OR", "OR_CI_lower", "OR_CI_upper"))
head(results)

map=fread("../UKB_data/phecode_definitions1.2.csv")
head(map)

data_all <- merge(results, map, by.x='Phecode', by.y='phecode_phecode', all.x=TRUE, all.y=FALSE)
first_upper <- function(s) {
  paste0(toupper(substr(s,1,1)), substr(s,2,nchar(s)))
}

data_all$category=first_upper(data_all$category)
head(data_all)
data_all <- data_all[data_all$category != "NULL", ]
data_all$`P-value` <- as.numeric(data_all$`P-value`)
data_all <- data_all[order(-data_all$`P-value`), ]
data_remove <- data_all[data_all$`P-value` ==0, ]

write.csv(data_remove, "../GSEM/removed_Pvalue_lt_0.csv", row.names = FALSE)

data_all <- data_all[data_all$`P-value` !=0 | is.na(data_all$`P-value`), ]
write.csv(data_all, file="../GSEM/CLM_multiPRS_PheWAS_results.csv")
