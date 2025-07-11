suppressMessages({
  library(data.table)
  library(dplyr)
})

#set n paramter
n <- commandArgs(trailingOnly=TRUE)[1]

#import random eigenproteins and traits
grn_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/eigengenes/eigenproteins/random/random_n", n, "_eigenproteins.tsv", sep="")
grn <- fread(grn_fn, data.table=FALSE)

#import traits
ph <- fread("~/sean/INTRePID/STARNET/phenotype/STARNET_pheno_GRN_traits.tsv", data.table=FALSE)

#combine as single dataframe
grn_ph <- merge(grn, ph, by="starnet_id")

#get variable labels
var <- colnames(grn)[-1]

#linear model to test protein variables on phenotypes
var_lm <- function(trait){
  fit_res <- list()
  for (x in var) {

    #remove Nans
    grn_ph_na <- grn_ph %>% filter(!is.na(.data[[trait]]))

    #fit linear model
    formula_str <- paste(trait, "~", x, collapse=" ")
    fit <- lm(as.formula(formula_str), data = grn_ph_na)
    coff <- summary(fit)$coefficients

    #directly compute -log10 p-values
    df <- df.residual(fit)
    tv <- coff[2, "t value"]
    log_pv <- -(log10(2) + pt(abs(tv), df, lower.tail = FALSE, log.p = TRUE) / log(10))
    coff_out <- cbind(coff, log_pv = c(NA, log_pv))
    fit_res[[x]] <- coff_out
  }

  #assemble results as single dataframe
  res_mat <- do.call(rbind, Map(cbind, variable = names(fit_res), fit_res))
  res_df <- cbind(intr = rownames(res_mat), res_mat) %>% as.data.frame

  #remove intercept and add trait label
  rownames(res_df) <- NULL
  res_filt <- filter(res_df, intr != "(Intercept)")
  res_filt$intr <- NULL
  res_filt$trait <- trait
  res_filt$n <- nrow(grn_ph_na)
  
  return(res_filt)
}

#test all traits
bmi_df <- var_lm("BMI")
hb_df <- var_lm("HbA1c")
tg_df <- var_lm("TG")
ldl_df <- var_lm("LDL")
hdl_df <- var_lm("HDL")
chol_df <- var_lm("Chol")
crp_df <- var_lm("CRP")
wp_df <- var_lm("Waist_hip")

#assemble as single dataframe
all_df <- rbind(bmi_df, hb_df, tg_df, ldl_df, hdl_df, chol_df, crp_df, wp_df)

#export results
out_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/eigengenes/eigenproteins/random/random_n", n, "_eigenproteins_trait_associations.tsv", sep="")
fwrite(all_df, out_fn, sep="\t")