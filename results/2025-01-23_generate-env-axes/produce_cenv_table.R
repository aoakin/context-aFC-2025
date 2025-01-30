library(rtracklayer)
library(dplyr)
library(data.table)

gtf <- readGFF('../../dat/Homo_sapiens.GRCh38.113.chr.gtf')
genes_protein_lncRNA <- gtf %>%
  dplyr::filter(type == "gene") %>% 
  dplyr::filter(gene_biotype %in% c("protein_coding", "lncRNA")) %>%
  dplyr::pull("gene_id")

expr <- read.table('../../dat/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct', skip=2, sep="\t", header=TRUE) #rows: genes, cols: tissues # nolint: line_length_linter.
expr$Name <-  gsub("\\.\\d+$", "", expr$Name) # remove .# at the end of every gene name
expr <- expr[expr$Name %in% genes_protein_lncRNA,]
expr <- expr[rowSums(expr[,3:ncol(expr)] != 0) > 0,] # remove columns only containing zeroes
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
pc <- prcomp(t_expr, scale = TRUE, center = TRUE) # zero centered, unit variance BEFORE analysis
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

write.table(tibble::rownames_to_column(cenv_df, "TissueType"), '../../dat/GTEx_cenv_data.tsv', sep='\t', row.names = FALSE)