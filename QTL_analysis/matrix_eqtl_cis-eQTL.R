suppressMessages({
    library(dplyr)
    library(data.table)
    library(MatrixEQTL)
    })
#set tissue paramter
tis <- commandArgs(trailingOnly=TRUE)[1]

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste("~/sean/INTRePID/STARNET/genotypes/STARNET_tissues/QTL_genotypes/", tis, "_STARNET_genotypedosageMAF05_cases_matched_QTLs_cis-SNPs_500Kb.tsv", sep="");
snps_location_file_name = paste("~/sean/INTRePID/STARNET/proteins/pQTLs/STARNET_Olink_cis-SNPs_500Kb_matqtl-snpsloc.tsv", sep="");

# Gene expression file name
expression_file_name = paste("~/sean/INTRePID/STARNET/gene_expression/", tis,".exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN_Olink.tsv", sep="");
gene_location_file_name = paste("~/sean/INTRePID/STARNET/gene_expression/STARNET_exp_matqtl_tss.tsv", sep="");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# No covariates
covariates_file_name = character();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0.05;

# Distance for local gene-SNP pairs
cisDist = 1e5;

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      
snps$fileOmitCharacters = "NA"; 
snps$fileSkipRows = 1;          
snps$fileSkipColumns = 1;       
snps$fileSliceSize = 2000;      
snps$LoadFile(SNP_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      
gene$fileOmitCharacters = "NA"; 
gene$fileSkipRows = 1;          
gene$fileSkipColumns = 1;       
gene$fileSliceSize = 2000;      
gene$LoadFile(expression_file_name);

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

## Run the analysis
me = Matrix_eQTL_main(
    snps = snps, 
    gene = gene, 
    output_file_name      = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel, 
    verbose = TRUE, 
    output_file_name.cis  = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

unlink(output_file_name_cis);


#calculate standard error
me$cis$eqtls$beta_se = me$cis$eqtls$beta / me$cis$eqtls$statistic;

# gene annotation for output
gene_ref <- read.table('~/sean/resources/gene_annotation_GRCh37_biomart.tsv', header=TRUE, sep='\t', check.names=FALSE, quote = "")
gene_ref <- gene_ref[c("Gene stable ID", "Gene name")]
gene_ref <- gene_ref[!duplicated(gene_ref), ]

df <- me$cis$eqtls
df_anno <- merge(df, gene_ref, by.x = 'gene', by.y='Gene stable ID')
df_anno <- df_anno[c("snps", "gene", "Gene name", "statistic", "pvalue", "FDR", "beta", "beta_se")]

#edit column names
setnames(df_anno, old = c('snps', 'pvalue', 'gene', 'Gene name'), new = c('rs_number','p-value', 'Gene stable ID', 'Gene_name'))

#export data
out_file_name <- paste("~/sean/INTRePID/STARNET/gene_expression/eQTLs/", tis, "_STARNET_cases_matched_QTLs_cis-eQTLs_500Kb.tsv", sep="")
write.table(df_anno, file=out_file_name, quote=FALSE, sep='\t', row.names=FALSE)