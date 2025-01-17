suppressMessages({
	library(data.table)
	library(dplyr)
	library(tidyverse)
	library(coloc, lib='~/temp/Rlib')
})

#set script paramters
tis <- commandArgs(trailingOnly=TRUE)[1]
thrd <- 100000

#import full summary statistics
pq_filename <- paste("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_P1.tsv", sep="")
pq <- fread(pq_filename, check.names=TRUE)
eq_filename <- paste("~/sean/INTRePID/STARNET/gene_expression/eQTLs/", tis, "_STARNET_cases_matched_QTLs_cis-eQTLs_500Kb.tsv", sep="")
eq <- fread(eq_filename, check.names=TRUE)

#import shared STARNET p/eQTLs and filter for tissue
qtl <- fread("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_eQTL_overlap_LD.tsv", check.names=TRUE)
qtl <- filter(qtl, tissue == tis)
qtl <- filter(qtl, Gene_name %in% eq$Gene_name)

#import reference datasets
snp_info <- fread("~/sean/INTRePID/STARNET/proteins/pQTLs/STARNET_Olink_cis-SNPs_500Kb.tsv", check.names=TRUE)
gene_ref <- fread("~/sean/resources/gene_annotation_GRCh37_biomart.tsv", check.names=TRUE)
ol_ref <- fread("~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv", check.names=TRUE)
maf <- fread("~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_MAFs.tsv", check.names=TRUE)

#Obtain SNP info for eQTLs and pQTLs
snp_info <- snp_info %>% distinct(rs_number, .keep_all= TRUE)
pq_ref <- merge(pq, select(snp_info, -c('Gene_name')), by='rs_number')
eq_ref <- merge(eq, select(snp_info, -c('Gene_name')), by='rs_number')

#run coloc as a function
run_coloc <- function(gene) {
	#obtain transcription start site for gene and SNPs within set threshold
	print(gene)
	tss <- filter(gene_ref, Gene.name==gene)[[1,2]]
	pq_id <- filter(ol_ref, Gene_name==gene)[[1,1]]
	eq_id <- filter(ol_ref, Gene_name==gene)[[1,4]]
	up <- tss+thrd
	dn <- tss-thrd
	pq_filt <- filter(pq_ref, OlinkID == pq_id & position < up & position > dn)
	eq_filt <- filter(eq_ref, Gene.stable.ID == eq_id & position < up & position > dn)

	#calculate beta variance as square of standard error
	pq_filt$betavar <- pq_filt$beta_se^2
	eq_filt$betavar <- eq_filt$beta_se^2

	#import common SNP MAFs
	pq_maf <- merge(pq_filt, maf, by="rs_number")
	eq_maf <- merge(eq_filt, maf, by="rs_number")

	#select for coloc input
	pq_sel <- subset(pq_maf, select = c(beta, betavar, rs_number, position, MAF))
	eq_sel <- subset(eq_maf, select = c(beta, betavar, rs_number, position, MAF))

	#remove duplicate SNPs and missing values
	pq_na <- pq_sel[!duplicated(pq_sel$rs_number), ]
	eq_na <- eq_sel[!duplicated(eq_sel$rs_number), ]
	pq_na <- na.omit(pq_na)
	eq_na <- na.omit(eq_na)

	#obtain and filter for common SNPs
	com_snps <- intersect(pq_na$rs_number, eq_na$rs_number)
	pq_stats <- filter(pq_na, rs_number %in% com_snps)
	eq_stats <- filter(eq_na, rs_number %in% com_snps)

	#obtain data lists for input to coloc
	pq_data <- list(beta=pq_stats$beta, varbeta=pq_stats$betavar, N=461, type="quant", snp=pq_stats$rs_number, position=pq_stats$position, MAF=pq_stats$MAF)
	eq_data <- list(beta=eq_stats$beta, varbeta=eq_stats$betavar, N=523, type="quant", snp=eq_stats$rs_number, position=eq_stats$position, MAF=eq_stats$MAF)

	#run coloc for single variant assumption with sensitivity analysis
	coloc_res <- coloc.abf(dataset1=pq_data, dataset2=eq_data)
	return(coloc_res$summary)
}

#Run coloc for all genes in dataset
gene_list <- unique(qtl$Gene_name)
coloc_list <- sapply(gene_list, run_coloc)
coloc_df <- t(data.frame(coloc_list))

#export output
out_filename <- paste("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/", tis, "_STARNET_Olink_cis-SNPs_500Kb_cases_combined_cis-pQTLs_LD_eQTL_coloc.tsv", sep="")
write.table(coloc_df, file=out_filename, quote=FALSE, sep='\t', row.names=TRUE)