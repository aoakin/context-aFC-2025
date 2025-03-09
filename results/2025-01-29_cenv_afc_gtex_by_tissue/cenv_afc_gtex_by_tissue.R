library(dplyr)
library(tibble)
library(tidyr)
library(data.table)
library(ShapleyValue)
library(tictoc)
# library(fitdistrplus)

# Generate coefficient vector
coef_vec_generator <- function(fit){
	# returns vector of coefficients
	##  will return NA for coefficients where 0 \in 95% CI
	s <- summary(fit)
	coef_vec <- s$coefficients[-1,1] # drops intercept
	total_terms <- length(s$coefficients[-1,1])
	confint_tbl <- data.frame(confint(fit)[-1,])
	colnames(confint_tbl)=c('lower_95ci', 'upper_95ci')
	contains_zero_vec <- confint_tbl %>%
	mutate(contains_zero = between(rep(0,total_terms), lower_95ci, upper_95ci)) %>%
	pull(contains_zero)
	coef_vec[contains_zero_vec] <-  NA
	# return(contains_zero_vec)
	return(coef_vec)
}
# Generate full model p-value
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

InvNormTransform <- function(x) {
	return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))) #-0.5 accounts for missing data? Credit: 
}

# cenv <- data.table::fread('../../dat/GTEx_cenv_data_median_counts.tsv', header=TRUE) # all of data
cenv <- data.table::fread('../../dat/GTEx_cenv_data_median_counts_70p.tsv', header=TRUE) # 90% of data
cenv <- tibble::column_to_rownames(cenv, "TissueType")
cat("Total cenv columns:", ncol(cenv), "\n")
print(as_tibble(cenv[]))
afcn_files <- list.files('../../dat/afcn')
afcn_files <- file.path('../../dat/afcn/', afcn_files)
afcn_tissues <- gsub("^.*/([^.]+)\\..*$", "\\1", afcn_files)
cenv <- cenv[rownames(cenv) %in% afcn_tissues, ] # drop cenv values not in afc
print(table(afcn_tissues %in% rownames(cenv)))

afcn_gv_tbl <- data.frame(matrix(vector(), 0, 3,
					dimnames = list(c(), c("gene_id", "variant_id", "tissue"))),
					stringsAsFactors=FALSE)

i = 1
for (file in afcn_files) {
	tissue_afc <- data.table::fread(file)
	tissue_afc <- tissue_afc[,c('gene_id', 'variant_id', 'log2_aFC')]
	tissue_afc$tissue <- afcn_tissues[i]
	afcn_gv_tbl <- rbind(afcn_gv_tbl, tissue_afc)
	i = i+1
}

cenv_colnames <- colnames(cenv)
print(cenv_colnames)
print(head(cenv[1:5]))
print(head(afcn_gv_tbl))


for (i in cenv_colnames) {
	afcn_gv_tbl[,i] <- cenv[afcn_gv_tbl$tissue, i]
	# afcn_gv_tbl$i <- cenv[afcn_gv_tbl$tissue == rownames(cenv),i]
}

afcn_gv_tbl <- afcn_gv_tbl %>% drop_na() # drop eQTLs with n/a
afcn_gv_tbl$log2_aFC <- InvNormTransform(afcn_gv_tbl$log2_aFC)

print(as.data.table(afcn_gv_tbl)[,1:6][sample(.N,6)])

unique_gv_combos <- afcn_gv_tbl %>% distinct(gene_id, variant_id)
print("Total unique GV combos considered:")
print(nrow(unique_gv_combos))

pval_tbl <- data.frame(matrix(vector(), 0, 2+ncol(cenv),
						dimnames = list(c(), c("gene_id", "variant_id", paste0("pval_",cenv_colnames)))),
						stringsAsFactors=FALSE)
qval_tbl <- data.frame(matrix(vector(), 0, 2+ncol(cenv),
						dimnames = list(c(), c("gene_id", "variant_id", paste0("qval_",cenv_colnames)))),
						stringsAsFactors=FALSE)
effectsize_tbl <- data.frame(matrix(vector(), 0, 2+ncol(cenv),
						dimnames = list(c(), c("gene_id", "variant_id", paste0("effsize_",cenv_colnames)))),
						stringsAsFactors=FALSE)
corr_tbl <- data.frame(matrix(vector(), 0, 2+ncol(cenv),
						dimnames = list(c(), c("gene_id", "variant_id", paste0("corr_",cenv_colnames)))),
						stringsAsFactors=FALSE)


cenv_p.val_lst <- list()
cenv_q.val_lst <- list()
effectsize_lst <- list()
cenv_fit_p <- c()
corr_lst <- list()

preFDR_total_sig_hit <- 0

# afc_cenv_formula <- paste("log2_aFC ~", paste0(cenv_colnames, collapse=" + "))

cat("Generating p-val_table\n")
# pb <- txtProgressBar(min = 0, max = nrow(unique_gv_combos), initial = 0, style=3)  # progress bar for generating correlations
# for (i in 5786:5787) { # TESTING

for (i in 1:nrow(unique_gv_combos)) {
	gv_combo = unique_gv_combos[i, ]
	cat("G: ", gv_combo$gene_id, " V: ", gv_combo$variant_id, "|", i, "out of", nrow(unique_gv_combos))
	gv_afc_cenv <- afcn_gv_tbl %>% 
		filter(gene_id == gv_combo$gene_id, variant_id == gv_combo$variant_id) %>%
		dplyr::select(-c(gene_id,variant_id))
	gv_afc_cenv <- column_to_rownames(gv_afc_cenv, "tissue")
	# fit_list <- list()
	# p.val_list <- list()
	# q.val_list <- list()
	# effectsize_list <- list()
	# corr_list <- list()
	
	cenv_fit <- lm(log2_aFC ~ ., data=gv_afc_cenv)
	cenv_fit_p[i] <- lmp(cenv_fit)
	if (is.nan(cenv_fit_p[i])) {cenv_fit_p[i] <- 999}
	
	# print(lmp(cenv_fit))
	effectsize_lst[[i]] <- coef_vec_generator(cenv_fit)
	if (cenv_fit_p[i] < 0.05) {
		preFDR_total_sig_hit = preFDR_total_sig_hit + 1
		cat("*\tpre-FDR hit #", preFDR_total_sig_hit, "\n")
	} else {
		cat("\n")
	}

	
	# s <- summary(cenv_fit)
	# print(s)
	# cenv_p.val_lst[[i]] <- s$coefficients[,'Pr(>|t|)'][-1] # remove intercept
	# cenv_q.val_lst[[i]] <- p.adjust(cenv_p.val_lst[[i]], method="BH")
	# effectsize_lst[[i]] <- s$coefficients[-1,1]
	# print('generating shapley values')
	# shapley_lst[[i]] <- as.numeric(shapleyvalue(gv_afc, gv_cenv)['Shapley Value',]) #COSTLY

	# for (j in names(gv_cenv)) {
	# 	fit_list[[j]] <- lm(gv_afc ~ get(j), gv_cenv)
	# 	p.val_list[j] <-  summary(fit_list[[j]])$coefficients[,'Pr(>|t|)'][['get(j)']]
	# 	q.val_list[j] <- p.adjust(p.val_list[[j]]) # not using mutoss::BH bc it's 1.5x slower on avg
	# 	effectsize_list[j] <- summary(fit_list[[j]])$coefficients[2,1]
	# 	corr_list[j] <- summary(fit_list[[j]])$r.squared
	# }
	# pval_tbl[nrow(pval_tbl) + 1,] <- append(list(gv_combo$gene_id, gv_combo$variant_id), p.val_list)
	# qval_tbl[nrow(qval_tbl) + 1,] <- append(list(gv_combo$gene_id, gv_combo$variant_id), q.val_list)
	# effectsize_tbl[nrow(effectsize_tbl) + 1,] <- append(list(gv_combo$gene_id, gv_combo$variant_id), effectsize_list)
	# corr_tbl[nrow(corr_tbl) + 1,] <- append(list(gv_combo$gene_id, gv_combo$variant_id), corr_list)
}
# print(as_tibble(gv_afc_cenv))
# p.val_list
# print(p.val_list)
# print(cenv_p.val_lst[1])
# print(class(cenv_p.val_lst[[1]]))
# pval_tbl <- data.table::rbindlist(cenv_p.val_lst)
# pval_tbl <- as.data.frame(do.call(rbind,cenv_p.val_lst))
# qval_tbl <- as.data.frame(do.call(rbind,cenv_q.val_lst))
# effectsize_tbl <- as.data.frame(do.call(rbind, effectsize_lst))
# # corr_tbl <- as.data.frame(do.call(rbind, corr_lst))
# shapley_tbl <- as.data.frame(do.call(rbind, shapley_lst))

es_tbl <- as.data.frame(do.call(rbind, effectsize_lst))
report_tbl <- as.data.frame(do.call(cbind, list(unique_gv_combos, cenv_fit_p, es_tbl)))
colnames(report_tbl) <- c('gene_id', 'variant_id', 'cenv_fit_p', cenv_colnames)
report_tbl <- report_tbl[report_tbl$cenv_fit_p != 999,]
cenv_fit_q <- p.adjust(report_tbl$cenv_fit_p, method="BH")
report_tbl['cenv_fit_q'] <- cenv_fit_q
report_tbl <- report_tbl %>% relocate(cenv_fit_q, .after = cenv_fit_p)
# print(pval_tbl)
# print(qval_tbl)
# print(effectsize_tbl)
# print(shapley_tbl)
# print(corr_tbl)
# quit(save="no")
# close(pb)

# print(as_tibble(report_tbl[5786,]))

cat('Saving report table\n')
print(as_tibble(report_tbl))
cat("Total significant eQTLs (post-FDR):", sum(cenv_fit_q < 0.05), "\n")
write.csv(report_tbl, "./report_table_afc_cenv_median_counts.csv", row.names=FALSE)

# cat('Saving P-val table\n')
# print(as_tibble(pval_tbl))
# write.csv(pval_tbl, "./pval_table_afc_cenv_median_counts.csv", row.names=FALSE)

# cat('Saving Q-val table\n')
# print(as_tibble(qval_tbl))
# write.csv(qval_tbl, "./qval_table_afc_cenv_median_counts.csv", row.names=FALSE)

# cat('Saving Effect size table\n')
# print(as_tibble(effectsize_tbl))
# write.csv(effectsize_tbl, "./es_table_afc_cenv_median_counts.csv", row.names=FALSE)

# cat('Saving Correlation table\n')
# print(as_tibble(corr_tbl))
# write.csv(corr_tbl, "./corr_table_afc_cenv_median_counts.csv", row.names=FALSE)

# cat('Saving Shapley Contribution table\n')
# print(as_tibble(shapley_tbl))
# write.csv(shapley_tbl, "./shapley_table_afc_cenv_median_counts.csv", row.names=FALSE)



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