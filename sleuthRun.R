library(sleuth)
library(data.table)
library(tidyverse)
library(patchwork)

stab_oe = fread("oe_kallisto_results.txt")
annot = fread("Dm_BDGP6.32.cdna_annotation.txt")
annot = rename(annot, target_id=transcript)
so_oe = sleuth_prep(stab_oe, ~condition, target_mapping = annot, transform_fun_counts= function(x) log2(x + 0.5), extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)
so_oe = sleuth_fit(so_oe, ~condition, 'full')
so_oe = sleuth_fit(so_oe, ~1, 'reduced')
so_oe = sleuth_lrt(so_oe, 'reduced', 'full')
models(so_oe)
so_oe = sleuth_wt(so_oe,which_beta = 'conditionRib', which_model='full')
oe_lrt = sleuth_results(so_oe, 'reduced:full', 'lrt', show_all = FALSE)
oe_wald = sleuth_results(so_oe, 'conditionRib', 'wt', show_all = FALSE)

oe_results = inner_join(oe_lrt,oe_wald,by=c("target_id","gene","gene_symbol","description")) %>% rename(transcript=target_id, oe_LRT_pval=pval.x, oe_LRT_qval=qval.x, oe_test_stat=test_stat, oe_Wald_pval=pval.y, oe_Wald_qval=qval.y, oe_b=b, oe_se_b=se_b) %>% select(transcript, gene, gene_symbol, description, oe_LRT_pval, oe_LRT_qval, oe_test_stat, oe_Wald_pval, oe_Wald_qval, oe_b, oe_se_b)
fwrite(oe_results, "Mierisch_lab_oe_Rib_w1118_results.csv")

fdr_results = filter(oe_results, oe_LRT_qval < 0.05, oe_b > 0.5 | oe_b < -0.5)
dim(fdr_results)[1]
fwrite(fdr_results, "Mierisch_lab_FDR0.05_opp_direction_FC>0.5_oe_Rib_w1118_results.csv")


fdr_results = filter(oe_results, oe_LRT_qval < 0.1, oe_b > 0.5 | oe_b < -0.5)
dim(fdr_results)[1]
fwrite(fdr_results, "Mierisch_lab_FDR0.1_opp_direction_FC>0.5_oe_Rib_w1118_results.csv")


