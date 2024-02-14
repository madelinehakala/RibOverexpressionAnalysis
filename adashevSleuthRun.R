library("sleuth")

### Tj run
tj_metadata <- read.table("/home/madeline/research/adashev_paper_analysis/tj_metadata.txt", header = TRUE)
tj_metadata
tj_metadata$condition <- factor(tj_metadata$condition, ref = "unsorted")
levels(tj_metadata$condition) <- list(unsorted = "unsorted", sorted = "sorted")
levels(tj_metadata$condition) 
adashev_tj_so <- sleuth_prep(tj_metadata, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)
adashev_tj_so <- sleuth_fit(adashev_tj_so, ~condition, 'full')
adashev_tj_so <- sleuth_fit(adashev_tj_so, ~1, 'reduced')
adashev_tj_so <- sleuth_lrt(adashev_tj_so, 'reduced', 'full')
models(adashev_tj_so)
tj_adashev_table <- sleuth_results(adashev_tj_so, 'reduced:full','lrt', show_all = FALSE)
View(tj_adashev_table)
tj_adashev_table$symbol <- ""
tj_adashev_tablet$gene_id <- ""
for (i in 1:nrow(tj_adashev_table)) { 
  id <- toString(tj_adashev_table[i,'target_id'])
  row <- which(fbtr_to_fbgn_2$transcript_ID == id ,arr.ind=TRUE)
  if (toString(fbtr_to_fbgn_2[row, 'gene_symbol']) != "")
    tj_adashev_table[i,'symbol'] <- toString(fbtr_to_fbgn_2[row, 'gene_symbol'])
}
tj_adashev_table
tj_adashev_table_q_0.05 <- dplyr::filter(tj_adashev_table, qval <= 0.05)
View(tj_adashev_table_q_0.05)
tj_adashev_table_p_0.05 <- dplyr::filter(tj_adashev_table, pval <= 0.05)
View(tj_adashev_table_p_0.05)
hist(tj_adashev_table$pval)
plot_bootstrap(adashev_tj_so, "FBtr0070195", units = "tpm", color_by = "condition")
adashev_tj_so_wt <- sleuth_prep(tj_metadata, 
                                full_model = ~condition, 
                                read_bootstrap_tpm = TRUE,
                                extra_bootstrap_summary = TRUE,
                                transformation_function = function(x) log2(x + 0.5)) 
adashev_tj_so_wt <- sleuth_fit(adashev_tj_so_wt)
models(adashev_tj_so_wt)
adashev_tj_so_wt <- sleuth_wt(adashev_tj_so_wt, which_beta = 'conditionz_sorted')
tj_adashev_table_wt <- sleuth_results(adashev_tj_so_wt, test = 'conditionz_sorted','wt', show_all = FALSE)
View(tj_adashev_table_wt)
tj_adashev_all <- left_join(tj_adashev_table, tj_adashev_table_wt, 
                            by = "target_id", suffix = c(".lrt", ".wt"))
View(tj_adashev_all)
tj_adashev_all <- select(tj_adashev_all, -c(test_stat, rss, degrees_free, mean_obs.lrt, var_obs.lrt, tech_var.lrt, sigma_sq.lrt,smooth_sigma_sq.lrt, final_sigma_sq.lrt, pval.wt, qval.wt, se_b, mean_obs.wt, var_obs.wt, tech_var.wt, sigma_sq.wt, smooth_sigma_sq.wt, final_sigma_sq.wt))
tj_adashev_all <- relocate(tj_adashev_all, symbol, .before = target_id)
fwrite(tj_adashev_all, "tj_adashev_all.csv")
tj_adashev_p_0.05_b_1 <- filter(tj_adashev_all, (pval.lrt <= 0.05 & b > 1))
View(tj_adashev_p_0.05_b_1)
fwrite(tj_adashev_p_0.05_b_1, "tj_adashev_p_0.05_b_1.csv")

### Eya run
eya_metadata <- read.table("/home/madeline/research/adashev_paper_analysis/eya_metadata.txt", header = TRUE)
eya_metadata
adashev_eya_so <- sleuth_prep(eya_metadata, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)
adashev_eya_so <- sleuth_fit(adashev_eya_so, ~condition, 'full')
adashev_eya_so <- sleuth_fit(adashev_eya_so, ~1, 'reduced')
adashev_eya_so <- sleuth_lrt(adashev_eya_so, 'reduced', 'full')
models(adashev_eya_so)
eya_adashev_table <- sleuth_results(adashev_eya_so, 'reduced:full','lrt', show_all = FALSE)
View(eya_adashev_table)
eya_adashev_table$symbol <- ""
eya_adashev_tablet$gene_id <- ""
fbtr_to_fbgn_2 <- fread("fbtr_to_fbgn2021_06.tsv")
for (i in 1:nrow(eya_adashev_table)) { 
  id <- toString(eya_adashev_table[i,'target_id'])
  row <- which(fbtr_to_fbgn_2$transcript_ID == id ,arr.ind=TRUE)
  if (toString(fbtr_to_fbgn_2[row, 'gene_symbol']) != "")
    eya_adashev_table[i,'symbol'] <- toString(fbtr_to_fbgn_2[row, 'gene_symbol'])
}
eya_adashev_table
adashev_eya_so_wt <- sleuth_prep(eya_metadata, 
                                 full_model = ~condition, 
                                 read_bootstrap_tpm = TRUE,
                                 extra_bootstrap_summary = TRUE,
                                 transformation_function = function(x) log2(x + 0.5)) 
adashev_eya_so_wt <- sleuth_fit(adashev_eya_so_wt)
models(adashev_eya_so_wt)
adashev_eya_so_wt <- sleuth_wt(adashev_eya_so_wt, which_beta = 'conditionz_sorted')
eya_adashev_table_wt <- sleuth_results(adashev_eya_so_wt, test = 'conditionz_sorted','wt', show_all = FALSE)
View(eya_adashev_table_wt)
eya_adashev_all <- left_join(eya_adashev_table, eya_adashev_table_wt, 
                             by = "target_id", suffix = c(".lrt", ".wt"))
View(eya_adashev_all)
eya_adashev_all <- select(eya_adashev_all, -c(test_stat, rss, degrees_free, mean_obs.lrt, var_obs.lrt, tech_var.lrt, sigma_sq.lrt,smooth_sigma_sq.lrt, final_sigma_sq.lrt, pval.wt, qval.wt, se_b, mean_obs.wt, var_obs.wt, tech_var.wt, sigma_sq.wt, smooth_sigma_sq.wt, final_sigma_sq.wt))
eya_adashev_all <- relocate(eya_adashev_all, symbol, .before = target_id)
fwrite(eya_adashev_all, "eya_adashev_all.csv")
eya_adashev_p_0.05_b_1 <- filter(eya_adashev_all, (pval.lrt <= 0.05 & b > 1))
View(eya_adashev_p_0.05_b_1)
fwrite(eya_adashev_p_0.05_b_1, "eya_adashev_p_0.05_b_1.csv")
                                 
