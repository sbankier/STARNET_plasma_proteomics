suppressMessages({
      library(UpSetR)
})

#plot correlations between protein and gene expression for all proteins
corr_count <- read.csv('~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/total_STARNET_pro_exp-match_correlations_intersection_q0.05.tsv', header=TRUE, sep='\t')
corr_plot <- corr_count[, -c(1, 2)]
png("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/total_STARNET_pro_exp-match_correlations_q0.05_upset.png", 
      units="cm", width=23, height=15, res=300)
upset(corr_plot, 
      nsets = 7, order.by = "freq", sets.bar.color = "#56B4E9", 
      mainbar.y.label = "Number of associated proteins (5% FDR)", 
      sets.x.label = "Associations per tissue", 
      text.scale = c(1.75, 1.75, 1.75, 1.5, 1.5, 1.75))
dev.off()

#plot results from multiple linear regression
lr_count <- read.csv("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/STARNET_LM_Olink_GE_tissues_intersection_5FDR.tsv", header=TRUE, sep="\t")
lr_plot <- lr_count[, -c(1, 2)]
png("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/STARNET_LM_Olink_GE_tissues_UpSet_5FDR.png", 
      units="cm", width=23, height=15, res=300)
upset(lr_plot, 
      nsets = 7, order.by = "freq", sets.bar.color = "#56B4E9", 
      mainbar.y.label = "Number of associated proteins (5% FDR)", 
      sets.x.label = "Associations per tissue", 
      text.scale = c(1.75, 1.75, 1.75, 1.5, 1.5, 1.75))
dev.off()

#plot overlapping QTLs
qt_count <- read.csv("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_cis-pQTL_tissue-intersection_count_20180383_20181235_cases_combined_unique.tsv", header=TRUE, sep="\t")
qt_plot <- qt_count[, -c(1, 2)]
png("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/visualisation/QTL_overlap_UpSet_5FDR.png", 
      units="cm", width=23, height=15, res=300)
upset(qt_plot, 
      nsets = 7, order.by = "freq", sets.bar.color = "#56B4E9", 
      mainbar.y.label = "Number of unqiue proteins (5% FDR)", 
      sets.x.label = "e/pQTLs per tissue", 
      text.scale = c(1.75, 1.75, 1.75, 1.5, 1.5, 1.75))
dev.off()

#plot QTLs in LD
qt_count <- read.csv("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_cis-pQTL_tissue-intersection_count_20180383_20181235_cases_combined_unique_LD.tsv", header=TRUE, sep="\t")
qt_plot <- qt_count[, -c(1, 2)]
png("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/visualisation/QTL_LD_UpSet_5FDR.png", 
      units="cm", width=23, height=15, res=300)
upset(qt_plot, 
      nsets = 7, order.by = "freq", sets.bar.color = "#56B4E9", 
      mainbar.y.label = "Number of unqiue proteins (5% FDR)", 
      sets.x.label = "e/pQTLs per tissue", 
      text.scale = c(1.75, 1.75, 1.75, 1.5, 1.5, 1.75))
dev.off()

#plot colocalised QTLs
coloc_count <- read.csv("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/total_STARNET_Olink_cis-SNPs_500Kb_cases_combined_cis-pQTLs_FDR5_eQTL_coloc_PP0.5_intersection.tsv", header=TRUE, sep="\t")
coloc_plot <- coloc_count[, -c(1, 2)]
png("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/QTL_coloc_UpSet_PP0.5.png", 
      units="cm", width=23, height=15, res=300)
upset(coloc_plot, 
      nsets = 7, order.by = "freq", sets.bar.color = "#56B4E9", 
      mainbar.y.label = "Number of unqiue proteins (PP > 0.5)", 
      sets.x.label = "Colocalisations per tissue", 
      text.scale = c(1.75, 1.75, 1.6, 1.5, 1.5, 1.75))
dev.off()

#plot correlations between tissue specific GRNs and proteins
spec_count <- read.csv("~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_tissue_specific_intersection_5FDR.tsv", header=TRUE, sep="\t")
spec_plot <- spec_count[, -c(1, 2)]
png("~/sean/INTRePID/STARNET/gene_expression/eigengenes/figures/Olink_cases_combined_adjusted_eigengene_correlations_tissue_specific_5FDR_UpSet.png", 
      units="cm", width=23, height=15, res=300)
upset(spec_plot, 
      nsets = 7, order.by = "freq", sets.bar.color = "#56B4E9", 
      mainbar.y.label = "Number of proteins (5% FDR)", 
      sets.x.label = "Proteins per tissue", 
      text.scale = c(1.75, 1.75, 1.75, 1.5, 1.5, 1.75))
dev.off()

#plot which tissues are represented in correlations between cross-tissue GRNs and proteins
cross_count <- read.csv("~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_cross_tissue_intersection_networks_5FDR.tsv", header=TRUE, sep="\t")
cross_plot <- cross_count[, -c(1, 2)]
png("~/sean/INTRePID/STARNET/gene_expression/eigengenes/figures/Olink_cases_combined_adjusted_eigengene_correlations_cross_tissue_intersection_networks_5FDR_UpSet.png", 
      units="cm", width=23, height=15, res=300)
upset(cross_plot, 
      nsets = length(cross_plot), order.by = "freq", sets.bar.color = "#56B4E9", 
      mainbar.y.label = "Number of networks (5% FDR)", 
      sets.x.label = "Networks per tissue", 
      text.scale = c(1.75, 1.75, 1.75, 1.5, 1.5, 1.75))
dev.off()

#plot correlations between cross-tissue GRNs and proteins when protein is in the GRN
cross_t_count <- read.csv("~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_cross_tissue_intersection_5FDR_GRN_present.tsv", header=TRUE, sep="\t")
cross_t_plot <- cross_t_count[, -c(1, 2)]
png("~/sean/INTRePID/STARNET/gene_expression/eigengenes/figures/Olink_cases_combined_adjusted_eigengene_correlations_cross_tissue_5FDR_GRN_present_UpSet.png", 
      units="cm", width=23, height=15, res=300)
upset(cross_t_plot, 
      nsets = length(cross_t_plot), order.by = "freq", sets.bar.color = "#56B4E9", 
      mainbar.y.label = "Number of proteins (5% FDR)", 
      sets.x.label = "Proteins per tissue", 
      text.scale = c(1.75, 1.75, 1.75, 1.5, 1.5, 1.75))
dev.off()

#plot network-trait associations
net_count <- read.csv("~/sean/INTRePID/STARNET/gene_expression/eigengenes/protein_cross_t_unique_genes_5FDR_module_intersection.tsv", header=TRUE, sep='\t', check.names=FALSE)
net_plot <- net_count[, -c(1)]
png("~/sean/INTRePID/STARNET/gene_expression/eigengenes/figures/protein_cross_t_unique_genes_5FDR_module_intersection_UpSet.png", 
      units="cm", width=23, height=15, res=300)
upset(net_plot, 
      nsets = length(net_plot), order.by = "freq", sets.bar.color = "#56B4E9", 
      mainbar.y.label = "Number of networks (5% FDR)", 
      sets.x.label = "Number of associations", 
      text.scale = c(1.75, 1.75, 1.75, 1.5, 1.5, 1.75))
dev.off()