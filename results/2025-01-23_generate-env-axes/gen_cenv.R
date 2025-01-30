# working through the analysis

library(ggplot2)
library(ggfortify)
library(data.table)
library(corrplot)
library(rmcorr)
library(rtracklayer)
library(dplyr)
library(tibble)

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



# Experiment 1: remove blood and pancreas; all genes
# expr <- subset(expr, select=-c(`Whole.Blood`, `Pancreas`))
## Result: 73% of expression variance primarily explained by mitochondria

# Experiment 2: all tissues; remove mitochondrial genes
# expr <- expr[!grepl("^MT-", expr$Description), ]
## Result: PC1 blood vs not blood, PC2 pancreas vs not pancreas

# Experiment 3: remove blood and pancreas; remove mitochondral genes
# expr <- subset(expr, select=-c(`Whole.Blood`, `Pancreas`))
# expr <- expr[!grepl("^MT-", expr$Description), ]
## Result: slighly more reasonable!

# Experiment 4: remove blood, pancreas, and pituitary; remove mitochondrial genes
# expr <- subset(expr, select=-c(`Whole.Blood`, `Pancreas`, `Pituitary`))
# expr <- expr[!grepl("^MT-", expr$Description), ]


tissue_types <- colnames(expr)[-c(1:2)]
table(duplicated(expr$Name)) # check if any duplicates; answer: all FALSE
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
head(t_expr)[,c(1:5)] # view first 5 row,col in expr. data
# expr <- t(expr) # cols
# expr_mat <- (t(expr[,-c(1:2)])

pc <- prcomp(t_expr, scale = TRUE, center = TRUE) # zero centered, unit variance BEFORE analysis
# pc <- prcomp(t_expr) # zero centered, unit variance BEFORE analysis

str(t_expr)
str(pc)


n_hits <- 20
total_pc_axes <- ncol(pc$rotation)


hit_df_rownames <- paste0("PC", 1:total_pc_axes)
hit_df_colnames <- paste0("hit_", 1:n_hits)
high_hits_df <- data.frame(matrix(NA, 
                  nrow = length(hit_df_rownames), 
                  ncol = length(hit_df_colnames),
                  dimnames = list(hit_df_rownames, hit_df_colnames)))


for (axis in colnames(pc$rotation)) {
  gene_eigenvec <- pc$rotation[,axis]
  top_genes_ens_id <- names(head(sort(abs(gene_eigenvec), decreasing=TRUE), n_hits))
  top_genes <- unname(gene_descriptions[top_genes_ens_id]) # this mess gets the top n genes
  top_gene_eigenvec <- gene_eigenvec[top_genes_ens_id]
  gene_eigenvec_direction <- ifelse(top_gene_eigenvec > 0, "+", "-")
  top_genes_and_direction <- paste(top_genes, gene_eigenvec_direction, sep = "/")
  high_hits_df[axis, ] <- top_genes_and_direction
}


# abs(pc$rotation[,1])

high_hits_df['PC1',]
high_hits_df['PC2',]

write.table(x = tibble::rownames_to_column(high_hits_df, "PCaxis"), file='./highest_loaded_eigengenes.tsv', sep='\t', row.names = FALSE)

pc


# pc$rotation
# top_hits_cenv1 <- tail(sort(pc$rotation[,1]), 10) # top 10 loaded genes on env xis 1 
# # cenv 1
# gene_descriptions[names(top_hits_cenv1)]

plot(pc) # TODO figure out how to plot pc axies
# head(t_expr)

pc_eigenval = (pc$sdev)^2
total_pc_var <- sum(pc_eigenval)


for (i in 1:length(pc_eigenval)) {
  cumulative_var = sum(pc_eigenval[1:i])/total_pc_var
  cat(sprintf("Cumulative variance at PC%s: %s\n", i, cumulative_var))
}

# Cumulative variance thresholds - not normalized
# > 90.0%: PC6   [~91.38%]
# > 95.0%: PC10  [~95.12%]
# > 99.0%: PC20  [~99.06%]
# > 99.7%: PC28  [~99.74%]

# gene_descriptions['ENSG00000210082']

autoplot(pc, x = 1, y = 2, data = t_expr, label = TRUE, label.label = rownames(t_expr)) + theme_bw()

min(pc$x[,1])

pc_dat <- data.frame(pc$x)
nrow(pc_dat) == length(tissue_types) # check if number of rows is equal to number of tissue types
rownames(pc_dat) <- tissue_types

# pc_dat

subset(pc_dat[order(pc_dat$PC1, decreasing = TRUE),], select=c('PC1')) # view tissues in order on pc axis

# tissue_types

rownames(pc_dat)[which.max(pc_dat$PC1)] 

# pc_dat$PC1

j = 1


cenv_df <- pc_dat

for (axis in colnames(pc_dat)) {
  pc_axis_coords <- pc_dat[,axis]
  cenv_i <- ecdf(pc_axis_coords)
  cenv_df[,axis] <- sapply(X = pc_dat[,axis],FUN = cenv_i)
}

plot(sort(cenv_df$PC51))



cenv <- ecdf(pc_dat$PC7)
plot(cenv)

# pc_dat$PC1

cenv(pc_dat$PC7)

# 
# # how correlated are pancreatic and whole blood expression?
# 
# gtex_corr_mat <- cor(subset(expr, select= -c(`Name`,`Description`)))
# testRes = cor.mtest(subset(expr, select= -c(`Name`,`Description`)), conf.level = 0.95)
# 
# 
# # gtex_rmcorr_mat <- rmcorr_mat(participant = )
# 
# ## specialized the insignificant value according to the significant level
# corrplot(gtex_corr_mat, p.mat = testRes$p, sig.level = 0.10, order = 'hclust', addrect = 2)
# 
# my.colramp <- colorRampPalette(c("black", "white","brown4"))
# 
# corrplot(gtex_corr_mat,
#          tl.col = 'black',
#          tl.cex=0.5, 
#          diag = TRUE, 
#          order = 'hclust', 
#          p.mat = testRes$p, 
#          sig.level = 0.10,
#          addrect = 18,
#          method = "color",
#          col = my.colramp(100))
# 
# testRes$p
# 
# subset(expr, select= -c(`Name`,`Description`))
# 
# summary(lm(Whole.Blood ~ `Brain...Hippocampus`, data=expr))
