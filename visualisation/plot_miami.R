suppressMessages({
	library(dplyr)
	library(data.table)
	library(gwaRs)
	library(ggplot2)
	library(patchwork)
})

#import overlapping e/pQTLs
epq <- fread("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_FDR5_eQTL_overlap.tsv", data.table=FALSE)

#import STARNET plasma cis-pQTLs and protein reference
pq <- fread('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_P1.tsv', data.table=FALSE)
pro_ref <- fread("~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_geneposition.tsv", data.table=FALSE)

#import SNP reference data  
snpref <- fread('~/sean/INTRePID/STARNET/proteins/pQTLs/STARNET_Olink_cis-SNPs_500Kb.tsv', data.table=FALSE)
snpref <- snpref[c('rs_number', 'chromosome', 'position')]
snpref <- snpref %>% distinct(rs_number, .keep_all = TRUE)

#annotate pQTLs
pq_mer <- merge(pq, snpref, by='rs_number', how='inner')
pq_mer <- pq_mer[c('Gene_name', 'rs_number', 'chromosome', 'position', 'p-value')]
colnames(pq_mer) <- c('Gene_name', 'SNP', 'CHR', 'BP', 'P')
pq_mer$Trait <- "Plasma"

#function to replace SNP ID with gene name for lead SNPs
rep_lead_snp <- function(df) {
  df <- df %>%
    group_by(Gene_name) %>%
    mutate(SNP = ifelse(row_number() == which.min(P), Gene_name, SNP)) %>%
    ungroup()

  return(df)
}

#annotate lead SNPs for pQTLs
pq_res <- rep_lead_snp(pq_mer)

#get pQTL 5% FDR threshold
pq_filt <- filter(pq, FDR <= 0.05) %>% arrange(desc(FDR))
pq_thrd <- pq_filt[1,8]

#import STARNET cis-eQTLs across tissues
eq_res <- list()
eq_thrds <- list()
eq_top <- list()
tis_lab <- c("LIV", "SKLM", "AOR", "MAM", "VAF", "SF", "Blood")
for (tis in tis_lab){
  eq <- fread(paste0("~/sean/INTRePID/STARNET/gene_expression/eQTLs/", tis, "_STARNET_cases_matched_QTLs_cis-eQTLs_500Kb.tsv"), sep="\t", header=TRUE, check.names=TRUE, quote="") %>% as.data.frame()
  
  #annotate eQTLs
  eq_mer <- merge(eq, snpref, by='rs_number', how='inner')
  eq_mer <- eq_mer[c('Gene_name', 'rs_number', 'chromosome', 'position', 'p.value')]
  colnames(eq_mer) <- c('Gene_name', 'SNP', 'CHR', 'BP', 'P')
  eq_mer$Trait <- tis

  #replace SNP ID with gene ID for top genes
  eq_ano <- rep_lead_snp(eq_mer)
  eq_res[[tis]] <- eq_ano

  #get top shared e/pQTLs
  epq_t <- filter(epq, tissue == tis) %>% 
  	arrange(desc(eQTL_FDR)) %>%
  	distinct(Gene_name, .keep_all = TRUE)
  epq_top <- head(epq_t, 10)
  eq_top[[tis]] <- epq_top$Gene_name

  #get pQTL 5% FDR threshold
  eq_filt <- filter(eq, FDR <= 0.05) %>% arrange(desc(FDR))
  eq_thrd <- eq_filt[1,5]
  eq_thrds[[tis]] <- eq_thrd
} 

#obtain SNP data for plotting
plot_data <- list()
for (tis in tis_lab){
	pq_lab <- select(pq_res, -Gene_name)
	eq_lab <- select(eq_res[[tis]], -Gene_name)
	t_data <- rbind(pq_lab, eq_lab)
	plot_data[[tis]] <- t_data %>% mutate(CHR = as.integer(CHR)) %>% arrange(CHR)
}

#plot Miami plots for different STARNET tissues
ano_liv <- eq_top[["LIV"]]
m_liv <- mirrored_man_plot(plot_data[["LIV"]], trait1 = "Plasma", trait2 = "LIV",
	genomewideline_trait1 = -log10(pq_thrd), genomewideline_trait2 = -log10(eq_thrds[["LIV"]]),
	title = NULL, trait2_chromCols = c("burlywood2", "burlywood4"), annotateSNP = ano_liv)

ano_sklm <- eq_top[["SKLM"]]
m_sklm <- mirrored_man_plot(plot_data[["SKLM"]], trait1 = "Plasma", trait2 = "SKLM",
	genomewideline_trait1 = -log10(pq_thrd), genomewideline_trait2 = -log10(eq_thrds[["SKLM"]]),
	title = NULL, trait2_chromCols = c("seagreen2", "seagreen4"), annotateSNP = ano_sklm)

ano_aor <- eq_top[["AOR"]]
m_aor <- mirrored_man_plot(plot_data[["AOR"]], trait1 = "Plasma", trait2 = "AOR",
	genomewideline_trait1 = -log10(pq_thrd), genomewideline_trait2 = -log10(eq_thrds[["AOR"]]),
	title = NULL, trait2_chromCols = c("firebrick2", "firebrick4"), annotateSNP = ano_aor)

ano_mam <- eq_top[["MAM"]]
m_mam <- mirrored_man_plot(plot_data[["MAM"]], trait1 = "Plasma", trait2 = "MAM",
	genomewideline_trait1 = -log10(pq_thrd), genomewideline_trait2 = -log10(eq_thrds[["MAM"]]),
	title = NULL, trait2_chromCols = c("darkgoldenrod2", "darkgoldenrod4"), annotateSNP = ano_mam)

ano_vaf <- eq_top[["VAF"]]
m_vaf <- mirrored_man_plot(plot_data[["VAF"]], trait1 = "Plasma", trait2 = "VAF",
	genomewideline_trait1 = -log10(pq_thrd), genomewideline_trait2 = -log10(eq_thrds[["VAF"]]),
	title = NULL, trait2_chromCols = c("darkorchid2", "darkorchid4"), annotateSNP = ano_vaf)

ano_sf <- eq_top[["SF"]]
m_sf <- mirrored_man_plot(plot_data[["SF"]], trait1 = "Plasma", trait2 = "SF",
	genomewideline_trait1 = -log10(pq_thrd), genomewideline_trait2 = -log10(eq_thrds[["SF"]]),
	title = NULL, trait2_chromCols = c("deeppink2", "deeppink4"), annotateSNP = ano_sf)

ano_blood <- eq_top[["Blood"]]
m_blood <- mirrored_man_plot(plot_data[["Blood"]], trait1 = "Plasma", trait2 = "Blood",
	genomewideline_trait1 = -log10(pq_thrd), genomewideline_trait2 = -log10(eq_thrds[["Blood"]]),
	title = NULL, trait2_chromCols = c("deepskyblue2", "deepskyblue4"), annotateSNP = ano_blood)

#obtain patchwork of all Miami plots
miami_patch <- m_liv/(m_sklm|m_blood)/(m_aor|m_mam)/(m_sf|m_vaf)

#export patchwork
ggsave(filename="~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/visualisation/miami_plts/STARNET_cases_cis-QTL_miami_patchwork.png",
	plot=miami_patch, units="cm", width=35, height=45, dpi=300)