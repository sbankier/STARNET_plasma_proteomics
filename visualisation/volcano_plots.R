suppressMessages({
       library(dplyr)
       library(data.table)
       library(EnhancedVolcano)
       library(patchwork)
})

#import tissue-protein LM
lm_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/STARNET_LM_Olink_GE_tissues.tsv", sep="")
lm <- fread(lm_fn, data.table=FALSE)

#import tissue-protein correlations
cor_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/total_STARNET_pro-exp_correlations.tsv", sep="")
cor <- fread(cor_fn, data.table=FALSE)

#volcano plots for individual tissues
plot_volcano <- function(df, x, xlab, y, ylab, title, col){
       vp <- EnhancedVolcano(df,
                             lab = df$Gene_name,
                             x = x,
                             xlab = xlab,
                             xlim = c(-1,1),
                             y = y,
                             ylab= ylab,
                             title = title,
                             subtitle = NULL,
                             pCutoff = 0.05,
                             FCcutoff = 0.2,
                             pointSize = 3.0,
                             axisLabSize = 28,
                             labSize = 7,
                             titleLabSize = 32,
                             col=c("gray", "gray", "gray", col),
                             legendPosition = 'none',
                             caption = NULL,
                             gridlines.major = FALSE,
                             gridlines.minor = FALSE,
                             widthConnectors = 1.0,
                             colConnectors = 'black')
}

#get colour code for tissues
tis <- c("LIV", "SKLM", "AOR", "MAM", "VAF", "SF", "Blood")
tis_col <- c("brown", "green", "red", "orange", "purple", "pink", "blue")

#plot volcano plots for individual tissues
lm_l <- list()
cor_l <- list()
for (i in seq_along(tis))
{
       #for linear model results
       lm_filt <- filter(lm, tissue==tis[i])
       lm_filt_plot <- plot_volcano(lm_filt, "Estimate", "Beta", "qvalues", "-log10 q-values", tis[i], tis_col[i])
       lm_l[[tis[i]]] <- lm_filt_plot
       lm_filt_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/",tis[i] , "_LM_protein_assoc_volcano_plot.png")
       ggsave(filename=lm_filt_fn, plot=lm_filt_plot, units="cm", width=24, height=25, dpi=300)

       #for correlation results
       cor_filt <- filter(cor, tissue==tis[i])
       cor_filt_plot <- plot_volcano(cor_filt, "SpearmanR", "Spearman R", "Spearmanqv", "-log10 q-values", tis[i], tis_col[i])
       cor_l[[tis[i]]] <- cor_filt_plot
       cor_filt_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/",tis[i] , "_STARNET_pro_exp-match_correlations_volcano.png")
       ggsave(filename=cor_filt_fn, plot=cor_filt_plot, units="cm", width=24, height=25, dpi=300)
}

plot_volcano_all <- function(df, x, xlab, y, ylab, title){
       #colour code for tissues
       keyvals <- ifelse(df$tissue == "LIV", "brown",
                         ifelse(df$tissue == "SKLM", "green",
                                ifelse(df$tissue == "AOR", "red",
                                       ifelse(df$tissue == "MAM", "orange",
                                              ifelse(df$tissue == "VAF", "purple",
                                                     ifelse(df$tissue == "SF", "pink",
                                                            ifelse(df$tissue == "Blood", "blue", NA)
                                                     )
                                              )
                                       )
                                )
                         )
                 )
       names(keyvals) <- df$tissue

       #generate volcano plot
       vp <- EnhancedVolcano(df,
                             lab = df$Gene_name,
                             x = x,
                             xlab = xlab,
                             xlim = c(-1,1),
                             y = y,
                             ylab= ylab,
                             title = title,
                             subtitle = NULL,
                             colCustom = keyvals,
                             pCutoff = 0.05,
                             FCcutoff = 0.2,
                             pointSize = 3.0,
                             axisLabSize = 28,
                             labSize = 7,
                             titleLabSize = 32,
                             legendPosition = 'none',
                             caption = NULL,
                             gridlines.major = FALSE,
                             gridlines.minor = FALSE)
}

#construct plots
lm_plot <- plot_volcano_all(lm, "Estimate", "Beta", "qvalues", "-log10 q-values", "All")
cor_plot <- plot_volcano_all(cor, "SpearmanR", "Spearman R", "Spearmanqv", "-log10 q-values", "All")

#export plots
lm_plot_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/LM_protein_assoc_volcano_plot.png", sep="")
ggsave(filename=lm_plot_fn, plot=lm_plot, units="cm", width=24, height=25, dpi=300)
cor_plot_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/total_STARNET_pro_exp-match_correlations_volcano.png", sep="")
ggsave(filename=cor_plot_fn, plot=cor_plot, units="cm", width=24, height=25, dpi=300)

#combined plots
lm_patch <- (lm_plot | lm_l[["LIV"]] | lm_l[["SKLM"]] | lm_l[["AOR"]]) / (lm_l[["MAM"]] | lm_l[["VAF"]] | lm_l[["SF"]] | lm_l[["Blood"]])
lm_patch_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/LM_protein_assoc_volcano_plot_patch.png", sep="")
ggsave(filename=lm_patch_fn, plot=lm_patch, units="cm", width=80, height=40, dpi=300)
cor_patch <- (cor_plot | cor_l[["LIV"]] | cor_l[["SKLM"]] | cor_l[["AOR"]]) / (cor_l[["MAM"]] | cor_l[["VAF"]] | cor_l[["SF"]] | cor_l[["Blood"]])
cor_patch_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/total_STARNET_pro_exp-match_correlations_volcano_patch.png", sep="")
ggsave(filename=cor_patch_fn, plot=cor_patch, units="cm", width=80, height=40, dpi=300)