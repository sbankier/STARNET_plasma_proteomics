suppressMessages({
  library(data.table)
  library(dplyr)
  library(qvalue)
})

#import protein data
pro_all <- fread("~/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year_unique_proteins.tsv", data.table=FALSE)
ref <- fread("~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart_STARNET_tissues.tsv", sep="\t", header=TRUE, check.names=TRUE)

#remove redundent columns
pro <- select(pro_all, -c('status', 'Sex', 'Age', 'year'))

#import traits
ph <- fread("~/sean/INTRePID/STARNET/phenotype/STARNET_pheno_GRN_traits.tsv", data.table=FALSE)

#combine as single dataframe
pro_ph <- merge(pro, ph, by="starnet_id")

#get variable labels
var <- colnames(pro)[-1]

#linear model to test protein variables on phenotypes
var_lm <- function(trait){
  fit_res <- list()
  for (x in var) {

    #remove Nans
    pro_ph_na <- pro_ph %>% filter(!is.na(.data[[trait]]))

    #fit linear model
    formula_str <- paste(trait, "~", x, collapse=" ")
    fit <- lm(as.formula(formula_str), data = pro_ph_na)
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
  res_filt$n <- nrow(pro_ph_na)

  return(res_filt)
}

#test all traits
bmi_df <- var_lm("BMI")
wh_df <- var_lm("Waist_hip")
hb_df <- var_lm("HbA1c")
tg_df <- var_lm("TG")
ldl_df <- var_lm("LDL")
hdl_df <- var_lm("HDL")
chol_df <- var_lm("Chol")
crp_df <- var_lm("CRP")

#assemble as single dataframe
all_df <- rbind(bmi_df, wh_df, hb_df, tg_df, ldl_df, hdl_df, chol_df, crp_df)

#obtain q-values
pv <- all_df[, 5] %>% as.numeric
qv <- qvalue(pv)
all_df$qvalues <- qv$qvalues

#annotate results
ano <- select(ref, c("OlinkID", "Gene_name"))
colnames(ano) <- c("variable", 'Gene_name')
all_ano <- merge(ano, all_df, by="variable")

#export results
fwrite(all_ano, "~/sean/INTRePID/STARNET/proteins/traits/protein_trait_associations.tsv", sep="\t")