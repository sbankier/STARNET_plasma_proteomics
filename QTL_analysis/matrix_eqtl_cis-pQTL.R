suppressMessages({
    library(dplyr)
    library(data.table)
    library(MatrixEQTL)
})
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste("~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_adj-age-sex-year_samples_unique_SNPs.tsv", sep="\t");
snps_location_file_name = paste("~/sean/INTRePID/STARNET/proteins/pQTLs/STARNET_Olink_cis-SNPs_500Kb_matqtl-snpsloc.tsv", sep="\t");

# Gene expression file name
expression_file_name = paste("~/sean/INTRePID/STARNET/proteins/adjusted/20180383_20181235_cases_combined_adj-age-sex-year_genotypedosageMAF05_unique_proteins.tsv", sep="\t");
gene_location_file_name = paste("~/sean/INTRePID/STARNET/proteins/pQTLs/STARNET_Olink_proteins_matqtl-geneloc_tss.tsv", sep="\t");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# No covariates
covariates_file_name = character();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0.05;

# Distance for local gene-SNP pairs
cisDist = 5e5;

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
gene_ref <- read.table('~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv', header=TRUE, sep='\t', check.names=FALSE, quote = "")
df <- me$cis$eqtls
df_anno <- merge(df, gene_ref, by.x = 'gene', by.y='OlinkID')
df_anno <- df_anno[c("snps", "gene", "Assay", "UniProt", "Gene stable ID", "Gene_name", "statistic", "pvalue", "FDR", "beta", "beta_se")]

#edit column names
setnames(df_anno, old = c('snps', 'pvalue', 'gene', 'Assay'), new = c('rs_number','p-value', 'OlinkID', 'Protein_name'))

#export data
write.table(df_anno, file='~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_P1.tsv', quote=FALSE, sep='\t', row.names=FALSE)