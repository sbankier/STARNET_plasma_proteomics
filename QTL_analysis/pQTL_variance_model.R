suppressMessages({
  library(data.table)
  library(dplyr)
})

#import protein expression levels and reference
pro <- fread("~/sean/INTRePID/STARNET/proteins/adjusted/20180383_20181235_cases_combined_adj-age-sex-year_genotypedosageMAF05_unique_proteins.tsv", data.table=FALSE, header=TRUE)
ref <- fread("~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv", data.table=FALSE, header=TRUE) %>% select(OlinkID, Gene_name)

#annotate with gene names
pro_ref <- merge(ref, pro, by.x="OlinkID", by.y="starnet_id")

#import LD-pruned cis-pQTL genotypes
geno <- fread("~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_cis-pQTLs_5FDR_LD_pruned.tsv", data.table=FALSE, header=TRUE)

#get unique genes
gene_l <- unique(geno$Gene_name)
r2_l <- numeric(length(gene_l))
r2_adj_l <- numeric(length(gene_l))
for (i in seq_along(gene_l)) {

  #get gene name
  gene <- gene_l[i]

  #filter, format and transpose protein and genotypes
  pro_filt <- filter(pro_ref, Gene_name==gene) %>% select(-OlinkID)
  rownames(pro_filt) <- pro_filt$Gene_name
  pro_filt$Gene_name <- NULL
  pro_filt_t <- t(pro_filt) %>% as.data.frame()

  geno_filt <- filter(geno, Gene_name==gene) %>% select(-Gene_name)
  rownames(geno_filt) <- geno_filt$rs_number
  geno_filt$rs_number <- NULL
  geno_filt_t <- t(geno_filt) %>% as.data.frame()

  #combine as single dataframe
  pro_geno <- cbind(pro_filt_t, geno_filt_t)

  #fit to linear model
  formula_str <- paste(colnames(pro_geno)[1], "~", paste(colnames(pro_geno)[-1], collapse = " + "))
  fit <- lm(formula_str, data = pro_geno)

  #get variance explained from model
  r2_l[i] <- summary(fit)$r.squared
  r2_adj_l[i] <- summary(fit)$adj.r.squared
}

#format results as dataframe
var_df <- data.frame(
  Gene_name = gene_l,
  R_squared = r2_l,
  Adjusted_R_squared = r2_adj_l
)

#export results
fwrite(var_df, "~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_5FDR_LD_pruned_variance_explained.tsv", sep="\t")