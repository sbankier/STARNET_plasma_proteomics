suppressMessages({
  library(data.table)
  library(dplyr)
  library(qvalue)
})

#import STARNET protein expression data and protein-gene reference
pro <- fread("~/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year_unique_proteins.tsv", sep="\t", header=TRUE, check.names=TRUE, data.table=FALSE)
ref <- fread("~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart_STARNET_tissues.tsv", sep="\t", header=TRUE, check.names=TRUE)

#transpose protein expression dataframe
prot <- transpose(select(pro, -c('starnet_id', 'status', 'Sex', 'Age', 'year')))
colnames(prot) <- pro$starnet_id
prot$starnet_id <- colnames(pro)[6:length(pro)]
prot <- prot %>% select(starnet_id, everything())

#import STARNET gene expression data
ge = list()
for (x in c("LIV", "SKLM", "AOR", "MAM", "VAF", "SF", "Blood")){
  df <- fread(paste0("~/sean/INTRePID/STARNET/gene_expression/", x, ".exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN"), sep="\t", header=TRUE, check.names=TRUE, data.table=FALSE)
  cols <- gsub("X", "", colnames(df))
  colnames(df) <- cols
  df_filt <- df[df$Gene.stable.ID %in% ref$Gene.stable.ID, ]
  df_filt$tissue <- x
  ge[[x]] <- df_filt
} 

#assemble gene expression data as single dataframe
ge_df <- bind_rows(ge)

#obtain common samples
com <- ge
com[["pro"]] <- prot
com_sam <- Reduce(intersect, lapply(com, names))

#filter gene and protein expression for common samples
gec <- c(c("Gene.stable.ID", "tissue"), com_sam)
ge_com <- subset(ge_df, select=gec)
proc <- c(c("starnet_id"), com_sam)
pro_com <- subset(prot, select=proc)

#get proteins that are expressed in at least one tissue
ref_tis <- select(ref, c("LIV", "SKLM", "AOR", "MAM", "VAF", "SF", "Blood", "OlinkID"))
filtered_ref <- ref_tis %>% filter(if_any(LIV:Blood, ~ .x))

#iteratively regress protein expression on gene expression
fit_res <- list()
for (protein in filtered_ref$OlinkID) {

  #get corresponding ensembl ID
  gene_id <- ref$Gene.stable.ID[ref$OlinkID == protein]
  
  #filter for protein and gene expression
  protein_df <- filter(pro_com, starnet_id == protein)
  gene_df <- filter(ge_com, Gene.stable.ID == gene_id)
  
  #skip protein if there is no data available
  if (nrow(protein_df) == 0) {
    print(paste("no data for", protein, "skipping...", sep=" "))
    next
  }
  
  #convert protein values to a vector
  protein_v <- protein_df %>% slice(1) %>% select(-1) %>% as.numeric
  
  #obtain list of vectors for tissue gene expression
  lm_v <- list()
  for (t in gene_df$tissue) {
    gene_v <- gene_df %>% filter(tissue == t) %>% select(-c(1:2)) %>% as.numeric
    lm_v[[t]] <- gene_v
  }
  
  #assemble all expression values in single list
  lm_v[["pro"]] <- protein_v
  
  #construct linear model formula
  tl <- lm_v
  tl[["pro"]] <- NULL
  formula_str <- paste("pro ~", paste(names(tl), collapse = "+"))
  
  #fit model
  fit <- lm(as.formula(formula_str), data = lm_v)
  coff <- summary(fit)$coefficients
  fit_res[[protein]] <- coff
  print(protein)
}

#assemble results as single dataframe and remove intercept
res_mat <- do.call(rbind, Map(cbind, OlinkID = names(fit_res), fit_res))
res_df <- cbind(tissue = rownames(res_mat), res_mat) %>% as.data.frame
rownames(res_df) <- NULL
res_filt <- filter(res_df, tissue != "(Intercept)")

#calculate q-values
pv <- res_filt[, 6] %>% as.numeric
qv <- qvalue(pv)
res_filt$qvalues <- qv$qvalues

#annotate results
ano <- select(ref, c("OlinkID", "Gene.stable.ID", "Gene_name"))
res_ano <- merge(ano, res_filt, by="OlinkID")

#export results
write.table(res_ano, "~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/STARNET_LM_Olink_GE_tissues.tsv", sep="\t", quote=FALSE)