library(Seurat)
library(seismicGWAS)
library(SingleCellExperiment)
library(data.table)

setwd("../GSEM/")


###############################seismicGWAS
rds_file <- "MASLD_snRNA_seq_seurat_v4.rds"
af <- readRDS(rds_file)

af_sce <- Seurat::as.SingleCellExperiment(af)

my_ce11_type_col <- "Cell_type_broad"

af_sscore <-calc_specificity(af_sce,ct_label_col=my_ce11_type_col)

bf_gwas_file <-"../GSEM/BF_FUMA/BF_FUMA/magma.genes.out"

bf_magma_gwas <- data.table::fread(bf_gwas_file,header = TRUE)
bf_magma_gwas <- as.data.frame(bf_magma_gwas)

bf_gwas_results <- get_ct_trait_associations(
  sscore =af_sscore,
  magma = bf_magma_gwas,
  magma_gene_col="SYMBOL"
)
bf_gwas_results

p_assoc <- plot_top_associations(bf_gwas_results,limit=20)
print(p_assoc)

target_ct <- NULL 
top_row_idx<-which.min(bf_gwas_results$FDR)
target_ct<-bf_gwas_results$cell_type[top_row_idx]

magma_for_inf <- bf_magma_gwas
magma_for_inf$GENE<-magma_for_inf$SYMBOL

tryCatch({
  if (!(target_ct %in% colnames(af_sscore))){
  stop(paste("False",target_ct,"do not exist"))
}
  inf_genes <-find_inf_genes(target_ct,af_sscore, magma_for_inf)
  if(nrow(inf_genes)>0){
    message("Plotting...")
    p_genes <-plot_inf_genes(inf_genes,gene_col = 'gene',num_labels=20)
    print(p_genes)
    message(">>>Success!")
  } else{
    message(">>>Warning")
  }
  },error =function(e){
  message("False, details:",e)
  message("Please check",target_ct,"exist or not)
})


