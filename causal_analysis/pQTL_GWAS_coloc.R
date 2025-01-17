suppressMessages({
    library(data.table)
    library(dplyr)
    library(coloc, lib='~/temp/Rlib')
})

#set analysis parameters
p <- commandArgs(trailingOnly=TRUE)[1]
thrd <- 100000

#import GWAS reference
gw_ref <- fread("~/sean/INTRePID/GWAS/GWAS_ref.csv", data.table=FALSE)
gw_info <- filter(gw_ref, phenotype==p)
gw_n <- gw_info[1, 2]
gw_fn <- gw_info[1, 3]

#import full GWAS summary statistics
gwas <- fread(gw_fn, check.names=TRUE, data.table=FALSE, 
    select=c("rs_number", "chromosome", "base_pair_location", "beta", "standard_error"))

#format GWAS summary statistics
cn <- c("rs_number", "chr", "POS", "beta", "se")
colnames(gwas) <- cn

#import STARNET cis-pQTLs
pq <- fread("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_P1.tsv", check.names=TRUE, data.table=FALSE)

#import protein and SNP reference
snp_ref <- fread('~/sean/INTRePID/STARNET/proteins/pQTLs/STARNET_Olink_cis-SNPs_500Kb.tsv', check.names=TRUE, data.table=FALSE)
pro_ref <- fread('~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_geneposition.tsv', check.names=TRUE, data.table=FALSE)
maf <- fread("~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_MAFs.tsv", check.names=TRUE, data.table=FALSE)
g_names <- subset(pro_ref, select=c("Gene_name"))

#import and format overlap results
ov_fn <- paste("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/GWAS/STARNET_Olink_20180383_20181235_cases_combined_cis-pQTLs_5FDR_", p, "_overlap.tsv", sep='')
pq_gwas <- fread(ov_fn, check.names=TRUE, data.table=FALSE)
pq_gwas_ref <- merge(pq_gwas, g_names, by="Gene_name")

#Function for coloc
run_coloc <- function(gene, gw, n) {
    #select SNPs within defined threshold of cis-protein
    print(gene)
    tss <- filter(pro_ref, Gene_name==gene)[[1,7]]
    chr <- filter(pro_ref, Gene_name==gene)[[1,6]]
    up <- tss+thrd
    dn <- tss-thrd
    snps <- filter(snp_ref, Gene_name == gene & chromosome == chr & position < up & position > dn) 
    pq_filt <- merge(pq, snps, by="rs_number")
    gw_filt <- merge(gw, snps, by="rs_number")

    #calculate beta variance as square of standard error
    pq_filt$betavar <- pq_filt$beta_se^2
    gw_filt$betavar <- gw_filt$se^2

    #obtain MAF
    pq_maf <- merge(pq_filt, maf, by="rs_number")
    gw_maf <- merge(gw_filt, maf, by="rs_number")

    #select for coloc input
    pq_sel <- subset(pq_maf, select = c(beta, betavar, rs_number, position, MAF))
    gw_sel <- subset(gw_maf, select = c(beta, betavar, rs_number, position, MAF))

    #remove duplicate SNPs and missing values
    pq_na <- pq_sel[!duplicated(pq_sel$rs_number), ]
    gw_na <- gw_sel[!duplicated(gw_sel$rs_number), ]
    pq_na <- na.omit(pq_na)
    gw_na <- na.omit(gw_na)

    #obtain and filter for common SNPs
    com_snps <- intersect(pq_na$rs_number, gw_na$rs_number)
    pq_stats <- filter(pq_na, rs_number %in% com_snps)
    gw_stats <- filter(gw_na, rs_number %in% com_snps)

    #obtain data lists for input to coloc
    pq_data <- list(beta=pq_stats$beta, varbeta=pq_stats$betavar, N=461, type="quant", snp=pq_stats$rs_number, position=pq_stats$position, MAF=pq_stats$MAF)
    gw_data <- list(beta=gw_stats$beta, varbeta=gw_stats$betavar, N=n, type="quant", snp=gw_stats$rs_number, position=gw_stats$position, MAF=gw_stats$MAF)

    #run coloc for single variant assumption with sensitivity analysis
    coloc_res <- coloc.abf(dataset1=pq_data, dataset2=gw_data)
    
    return(coloc_res$summary)

}

#filter cis-pQTL summary statistics for matched proteins
pq_gwas_gene <- unique(pq_gwas_ref$Gene_name)

#Run coloc for all genes in dataset
coloc_gwas <- sapply(pq_gwas_gene[2:length(pq_gwas_gene)], run_coloc, gw=gwas, n=gw_n)

#format and annotate results
format_res <- function(coloc_res) {
    #assemble as dataframe
    coloc_df <- data.frame(t(coloc_res))

    #annotate output
    coloc_df$Gene_name <- rownames(coloc_df)
    ano <- subset(pro_ref, select=c(Assay, Gene_name, OlinkID))
    coloc_ano <- merge(ano, coloc_df, by="Gene_name")

    return(coloc_ano)

}

gwas_out <- format_res(coloc_gwas)

#export output
out_fn <- paste("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/GWAS/STARNET_Olink_20180383_20181235_cases_combined_cis-pQTLs_5FDR_", p, "_coloc_res.tsv", sep="")
write.table(gwas_out, file=out_fn, quote=FALSE, sep="\t")