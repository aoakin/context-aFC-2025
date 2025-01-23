library(ggplot2)

expr <- read.table('../../dat/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct', skip=2, sep="\t")

prcomp(expr)

