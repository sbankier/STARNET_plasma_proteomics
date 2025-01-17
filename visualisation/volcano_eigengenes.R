 suppressMessages({
       library(dplyr)
       library(data.table)
       library(EnhancedVolcano)
       library(patchwork)
})

#import protein-eigengene correlations
cor_fn <- paste("~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations.tsv", sep="")
cor <- fread(cor_fn, data.table=FALSE)
cor <- cor %>% mutate(network_tissue = ifelse(network_tissue == "BLOOD", "Blood", network_tissue))

#add network labels 
labs=cor[, c("Gene_name", "network_id")]
cor$labs <- apply(labs, 1, function(row) {
  paste(row, collapse = " -")
})

#filter for network type
cor_ts <- filter(cor, network_type == "specific")
cor_ct <- filter(cor, network_type == "cross")

#volcano plots for individual tissues
plot_volcano <- function(df, x, xlab, y, ylab, title, col, labels){
       vp <- EnhancedVolcano(df,
                             lab = df$Gene_name,
                             x = x,
                             xlab = xlab,
                             xlim = c(-1,1),
                             y = y,
                             ylab= ylab,
                             selectLab = labels,
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
                             drawConnectors = TRUE,
                             widthConnectors = 1.0,
                             colConnectors = 'black',
                             max.overlaps=30)
}

#get information about lead correlations
lead_cor <- fread("~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_tissue_transcript_lead_correlations.tsv", data.table=FALSE)

#remove duplicate proteins keeping strongest association
cor_ct_un <- cor_ct %>% arrange(spearman_qv) %>% distinct(Gene_name, .keep_all = TRUE) 

#plot volcano plot for cross-tissue
ct_lab <- lead_cor %>%
  rowwise() %>%
  filter(cross == max(c(transcript, cross, specific), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(cross)) %>%
  head(30) %>%
  as.data.frame() 

ct_plot <- plot_volcano(cor_ct_un, "spearman_r", "Spearman R", "spearman_qv", "-log10 q-values", "Cross-tissue", "magenta", ct_lab$Gene_name)
ct_fn <- paste("/Home/ii/seanb/sean/INTRePID/STARNET/gene_expression/eigengenes/figures/volcano/GRN_protein_correlations_CT_volcano.png")
ggsave(filename=ct_fn, plot=ct_plot, units="cm", width=24, height=25, dpi=300)

#get colour code for tissues
tis <- c("LIV", "SKLM", "AOR", "MAM", "VAF", "SF", "Blood")
tis_col <- c("brown", "green", "red", "orange", "purple", "pink", "blue")

#plot volcano plots for individual tissues
ts_l <- list()
for (i in seq_along(tis))
{

       #filter tissue specific correlations for tissue and get absolute R
       ts_filt <- filter(cor_ts, network_tissue==tis[i])
       ts_filt_un <- ts_filt %>% arrange(spearman_qv) %>% distinct(Gene_name, .keep_all = TRUE)
       ts_filt_un$r_abs <- abs(ts_filt_un$spearman_r)

       #get lead gene that is highest out of all correlations
       ts_lab <- lead_cor %>%
         rowwise() %>%
         filter(spec_tissue == tis[i] & specific == max(c(transcript, cross, specific), na.rm = TRUE)) %>%
         ungroup() %>%
         arrange(desc(specific)) %>%
         head(10) %>%
         merge(ts_filt_un, by="Gene_name") %>%
         filter(r_abs >= 0.2) %>%
         as.data.frame()

       #generate volcano plots and store as list
       ts_filt_plot <- plot_volcano(ts_filt_un, "spearman_r", "Spearman R", "spearman_qv", "-log10 q-values", tis[i], tis_col[i], ts_lab$Gene_name)
       ts_l[[tis[i]]] <- ts_filt_plot
       ts_filt_fn <- paste("/Home/ii/seanb/sean/INTRePID/STARNET/gene_expression/eigengenes/figures/volcano/", tis[i], "GRN_protein_correlations_TS_volcano.png")
       ggsave(filename=ts_filt_fn, plot=ts_filt_plot, units="cm", width=24, height=25, dpi=300)
}

#construct combined plots
patch <- (
  ct_plot | 
  (ts_l[["LIV"]] + ts_l[["SKLM"]] + ts_l[["AOR"]] + ts_l[["MAM"]] + 
   ts_l[["VAF"]] + ts_l[["SF"]] + ts_l[["Blood"]]) +
  plot_layout(ncol = 2, widths = c(1, 1))
)

#export patch plot
patch_fn <- paste("/Home/ii/seanb/sean/INTRePID/STARNET/gene_expression/eigengenes/figures/volcano/GRN_protein_correlations_volcano_patch.png", sep="")
ggsave(filename=patch_fn, plot=patch, units="cm", width=60, height=40, dpi=300)