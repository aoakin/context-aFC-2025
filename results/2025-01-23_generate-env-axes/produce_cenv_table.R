library(rtracklayer)
library(dplyr)
library(data.table)

gtf <- readGFF('../../dat/Homo_sapiens.GRCh38.113.chr.gtf')
genes_protein_lncRNA <- gtf %>%
  dplyr::filter(type == "gene") %>% 
  dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
  dplyr::pull("gene_id")

expr <- read.table('../../dat/MedianGeneCount_GTEx_v8.tsv', sep="\t", header=TRUE) #rows: genes, cols: tissues # nolint: line_length_linter.
afcn_files <- list.files('../../dat/afcn')
afcn_files <- file.path('../../dat/afcn/', afcn_files)
afcn_tissues <- gsub("^.*/([^.]+)\\..*$", "\\1", afcn_files)
expr <- expr %>%
  select(Name, Description, any_of(afcn_tissues)) # drop tissues(indivs) without afc
expr <- expr[rowSums(expr[-c(1,2)]) > 10,] # filter out genes with low expression (<= 10 counts across all samples)
expr$Name <-  gsub("\\.\\d+$", "", expr$Name) # remove .# at the end of every gene name
expr <- expr[expr$Name %in% genes_protein_lncRNA,]
# expr <- expr[rowSums(expr[,3:ncol(expr)] != 0) > 0,] # remove columns only containing zeroes
gene_descriptions <- setNames(expr$Description, expr$Name)
tissue_types <- colnames(expr)[-c(1:2)]
t_expr <- transpose(expr)
names(t_expr) <- t_expr[1,] # make gene IDs the column names
t_expr <- t_expr[-(1:2),]       # remove gene id and description
# head(t_expr)[,c(0:5)] # view first 5 row,col in expr. data
# t_expr <- sapply(t_expr, as.numeric)
t_expr <- as.data.frame(sapply(t_expr, as.numeric))
log2p1 <- function(x) {
  return (log2(x+1))
}
t_expr <- data.frame(lapply(t_expr, log2p1))
rownames(t_expr) <- tissue_types
# pc <- prcomp(t_expr, scale = TRUE, center = TRUE) # zero centered, unit variance BEFORE analysis
pc <- prcomp(t_expr, scale = TRUE, center = TRUE) # zero centered, variance not normalized BEFORE analysis | why? different genes need not be weighted equally
n_hits <- 20
total_pc_axes <- ncol(pc$rotation)

pc_dat <- data.frame(pc$x)
rownames(pc_dat) <- tissue_types

cenv_df <- pc_dat
for (axis in colnames(pc_dat)) {
  pc_axis_coords <- pc_dat[,axis]
  cenv_i <- ecdf(pc_axis_coords)
  cenv_df[,axis] <- sapply(X = pc_dat[,axis],FUN = cenv_i)
}
head(cenv_df[1:5])
write.table(tibble::rownames_to_column(cenv_df, "TissueType"), '../../dat/GTEx_cenv_data_median_counts.tsv', sep='\t', row.names = FALSE)


# Generate tables with less PCs to avoid square matrix problem
pc_eigenval <- (pc$sdev)^2
total_pc_var <- sum(pc_eigenval)

for (i in 1:total_pc_axes) {
  cumvar <- sum(pc_eigenval[1:i])/total_pc_var
  if (cumvar >= 0.60 && !exists("cenv_cutoff_60p")) {
    cenv_cutoff_60p <- i
  }
  if (cumvar >= 0.70 && !exists("cenv_cutoff_70p")) {
    cenv_cutoff_70p <- i
  }
  if (cumvar >= 0.90 && !exists("cenv_cutoff_90p")) {
    cenv_cutoff_90p <- i
  }
  if (cumvar >= 0.95 && !exists("cenv_cutoff_95p")) {
    cenv_cutoff_95p <- i
  }
  if (cumvar >= 0.99 && !exists("cenv_cutoff_99p")) {
    cenv_cutoff_99p <- i
  }
}
write.table(tibble::rownames_to_column(cenv_df[1:cenv_cutoff_60p], "TissueType"), '../../dat/GTEx_cenv_data_median_counts_60p.tsv', sep='\t', row.names = FALSE)
write.table(tibble::rownames_to_column(cenv_df[1:cenv_cutoff_70p], "TissueType"), '../../dat/GTEx_cenv_data_median_counts_70p.tsv', sep='\t', row.names = FALSE)
write.table(tibble::rownames_to_column(cenv_df[1:cenv_cutoff_90p], "TissueType"), '../../dat/GTEx_cenv_data_median_counts_90p.tsv', sep='\t', row.names = FALSE)
write.table(tibble::rownames_to_column(cenv_df[1:cenv_cutoff_95p], "TissueType"), '../../dat/GTEx_cenv_data_median_counts_95p.tsv', sep='\t', row.names = FALSE)
write.table(tibble::rownames_to_column(cenv_df[1:cenv_cutoff_99p], "TissueType"), '../../dat/GTEx_cenv_data_median_counts_99p.tsv', sep='\t', row.names = FALSE)