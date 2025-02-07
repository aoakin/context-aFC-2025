library(dplyr)
library(tibble)
library(tidyr)
library(data.table)

cenv <- read.table('../../dat/GTEx_cenv_data_median_counts.tsv', header=TRUE)
cenv <- tibble::column_to_rownames(cenv, "TissueType")
afcn_files <- list.files('../../dat/afcn')
afcn_files <- file.path('../../dat/afcn/', afcn_files)
afcn_tissues <- gsub("^.*/([^.]+)\\..*$", "\\1", afcn_files)

afcn_gv_tbl <- data.frame(matrix(vector(), 0, 3,
					dimnames = list(c(), c("gene_id", "variant_id", "tissue"))),
					stringsAsFactors=FALSE)

i = 1
for (file in afcn_files) {
	tissue_afc <- read.csv(file)
	tissue_afc <- tissue_afc[,c('gene_id', 'variant_id', 'log2_aFC')]
	tissue_afc$tissue <- afcn_tissues[i]
	afcn_gv_tbl <- rbind(afcn_gv_tbl, tissue_afc)
	i = i+1
}

cenv_colnames <- colnames(cenv)
print(head(cenv[1:5]))
print(head(afcn_gv_tbl))


for (i in cenv_colnames) {
	afcn_gv_tbl[,i] <- cenv[afcn_gv_tbl$tissue, i]
	# afcn_gv_tbl$i <- cenv[afcn_gv_tbl$tissue == rownames(cenv),i]
}

afcn_gv_tbl <- afcn_gv_tbl %>% drop_na() # drop eQTLs with n/a

print(as.data.table(afcn_gv_tbl)[,1:6][sample(.N,6)])

unique_gv_combos <- afcn_gv_tbl %>% distinct(gene_id, variant_id)

pval_tbl <- data.frame(matrix(vector(), 0, 2+ncol(cenv),
						dimnames = list(c(), c("gene_id", "variant_id", paste0("pval_",colnames(cenv))))),
						stringsAsFactors=FALSE)
effectsize_tbl <- data.frame(matrix(vector(), 0, 2+ncol(cenv),
						dimnames = list(c(), c("gene_id", "variant_id", paste0("effsize_",colnames(cenv))))),
						stringsAsFactors=FALSE)
corr_tbl <- data.frame(matrix(vector(), 0, 2+ncol(cenv),
						dimnames = list(c(), c("gene_id", "variant_id", paste0("corr_",colnames(cenv))))),
						stringsAsFactors=FALSE)


cat("Generating p-val_table\n")
# pb <- txtProgressBar(min = 0, max = nrow(unique_gv_combos), initial = 0, style=3)  # progress bar for generating correlations

for (i in 1:nrow(unique_gv_combos)) {
	gv_combo = unique_gv_combos[i, ]
	cat("G: ", gv_combo$gene_id, " V: ", gv_combo$variant_id, "|", i, "out of", nrow(unique_gv_combos), "\n")
	gv_afc_cenv <- afcn_gv_tbl %>% 
		filter(gene_id == gv_combo$gene_id, variant_id == gv_combo$variant_id) %>%
		select(-c('gene_id','variant_id'))
	gv_afc_cenv <- column_to_rownames(gv_afc_cenv, "tissue")
	gv_afc <- gv_afc_cenv$log2_aFC
	gv_cenv <- gv_afc_cenv %>% select(-log2_aFC)
	fit_list <- list()
	p.val_list <- list()
	effectsize_list <- list()
	corr_list <- list()
	for (j in names(gv_cenv)) {
		fit_list[[j]] <- lm(gv_afc ~ get(j), gv_cenv)
		p.val_list[j] <- summary(fit_list[[j]])$coefficients[,'Pr(>|t|)'][['get(j)']]
		effectsize_list[j]  <- summary(fit_list[[j]])$coefficients[2,1]
		corr_list[j] <- summary(fit_list[[j]])$r.squared
		# err=summary(fit)$coefficients[2,2] coef + c(-1,1)*err*qt(0.975, 42)
	}
	pval_tbl[nrow(pval_tbl) + 1,] <- append(list(gv_combo$gene_id, gv_combo$variant_id), p.val_list)
	effectsize_tbl[nrow(effectsize_tbl) + 1,] <- append(list(gv_combo$gene_id, gv_combo$variant_id), effectsize_list)
	corr_tbl[nrow(corr_tbl) + 1,] <- append(list(gv_combo$gene_id, gv_combo$variant_id), corr_list)
	# setTxtProgressBar(pb,i) # update progress bar
}
# close(pb)

cat('Saving P-val table\n')
print(as_tibble(pval_tbl))
write.csv(pval_tbl, "./pval_table_afc_cenv_median_counts.tsv", sep='\t')

cat('Saving Effect size table\n')
print(as_tibble(effectsize_tbl))
write.csv(effectsize_tbl, "./es_table_afc_cenv_median_counts.tsv", sep='\t')

cat('Saving Correlation table\n')
print(as_tibble(corr_tbl))
write.csv(corr_tbl, "./corr_table_afc_cenv_median_counts.tsv", sep='\t')



# gv_combo = unique_gv_combos[1,]
# print(gv_combo$gene_id)
# print(gv_combo$variant_id)

# gv_afc_cenv <- afcn_gv_tbl %>% 
# 	filter(gene_id == gv_combo$gene_id, variant_id == gv_combo$variant_id) %>%
# 	select(-c('gene_id','variant_id'))

# gv_afc_cenv <- column_to_rownames(gv_afc_cenv, "tissue")


# gv_afc <- gv_afc_cenv$log2_aFC
# gv_cenv <- gv_afc_cenv %>% select(-log2_aFC)

# print(head(gv_afc_cenv)[,1:6])
# print(head(gv_afc))
# print(head(gv_cenv)[,1:6])
# # print(as.data.table(gv_afc)[sample(.N,6)])
# # print(as.data.table(gv_cenv)[,1:6][sample(.N,6)])
# corr_list = list()
# p.val_list = list()
# for (j in names(gv_cenv)) {
# 	corr_list[[j]] <- lm(gv_afc ~ get(j), gv_cenv)
# 	p.val_list[j] <- summary(corr_list[[j]])$coefficients[,'Pr(>|t|)'][['get(j)']]
# }

# corr_tbl <- data.frame(matrix(vector(), 0, 2+ncol(gv_cenv),
# 					dimnames = list(c(), c("gene_id", "variant_id", paste0("corr",colnames(gv_cenv))))),
# 				 	stringsAsFactors=FALSE)



# corr_tbl[nrow(corr_tbl) + 1,] <- append(list(gv_combo$gene_id, gv_combo$variant_id), p.val_list)

# print(as_tibble(corr_tbl))

# # afc

# # print(afcn_gv_tbl[afcn_gv_tbl$gene_id == 'ENSG00000197580' & afcn_gv_tbl$variant_id == 'chr11_112138402_A_T_b38',])

# # print(afcn_gv_tbl[afcn_gv_tbl['gene_id' == unique_gv_combos[1,'gene_id'] & afcn_gv_tbl['variant_id' == unique_gv_combos[1,'varint_id']]],])