suppressMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

#import colocalisation results and filter for significance
coloc <- fread('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/total_STARNET_Olink_cis-SNPs_500Kb_cases_combined_cis-pQTLs_FDR5_eQTL_coloc.tsv', sep="\t", data.table=FALSE)
coloc_sig <- filter(coloc, PP.H4.abf >= 0.5)

#import gene position reference
ref <- fread("~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_geneposition.tsv", sep="\t", data.table=FALSE)
colnames(ref)[colnames(ref) %in% c("Chromosome/scaffold name", "Transcription start site (TSS)")] <- c("chr", "TSS")

#append gene position to coloc results
ref_lab <- select(ref, c("Gene_name", "chr", "TSS"))
coloc_ref <- merge(coloc_sig, ref_lab, by="Gene_name")
coloc_ref$chr <- as.numeric(coloc_ref$chr)

#obtain chromosome lengths  
chrlen <- data.frame(chr = 1:22, 
                    chr_len = c(249225077, 
                                243185679, 
                                197842617, 
                                190922257, 
                                180714439, 
                                170919735, 
                                159125187, 
                                146300622, 
                                141066490, 
                                135436023, 
                                134945765, 
                                133820598, 
                                115108598, 
                                107285437, 
                                102429049,  
                                90170495,  
                                81065254,  
                                78000441,  
                                59098134,  
                                62960229,  
                                48099610,  
                                51208568))

#prepare x-axis as a cumulative sum of chromosome lengths
chrlen$tot = cumsum(chrlen$chr_len)-chrlen$chr_len
plotdata = merge(coloc_ref, chrlen, by="chr", sort = F)
plotdata$BPcum = plotdata$TSS + plotdata$tot

#define grey bands to divide chromosomes
chrbreaks = unique(chrlen[order(chrlen$chr),c('chr', 'tot', 'chr_len')])
colnames(chrbreaks)[2] = 'start'
chrbreaks$end = chrbreaks$start + chrbreaks$chr_len
rep_col <- rep(c('white', 'grey'), 11)
chrbreaks$col <- rep_col
chrbreaks$mid = chrbreaks$start + ((chrbreaks$end - chrbreaks$start) / 2)

#import GWAS colocalisation results and filter based on PP threshold
gw_l <- list()
for (x in c("CAD", "LDL", "HbA1C", "BMI", "CRP")) {
  gw_fn <- paste("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/GWAS/STARNET_Olink_20180383_20181235_cases_combined_cis-pQTLs_5FDR_", x, "_coloc_res.tsv", sep="")
  gw_df <- read.table(gw_fn, sep="\t")
  gw_df <- mutate(gw_df, x = PP.H4.abf >= 0.5)
  colnames(gw_df)[10] <- x
  gw_df <- subset(gw_df, select=c("Gene_name", x))
  gw_l[[x]] <- gw_df
}

#label QTL colocalisation results with GWAS results
plotdatagw <- left_join(plotdata, gw_l[["CAD"]], by="Gene_name") %>% 
              left_join(y=gw_l[["LDL"]], by="Gene_name") %>%
              left_join(y=gw_l[["CRP"]], by="Gene_name") %>%
              left_join(y=gw_l[["HbA1C"]], by="Gene_name") %>%
              left_join(y=gw_l[["BMI"]], by="Gene_name") %>% 
              replace(is.na(.), FALSE) %>%
              mutate(
                Conditions = rowSums(select(., CAD:BMI)),
                Status = case_when(
                  Conditions > 1 ~ "multiple",
                  CAD ~ "CAD",
                  LDL ~ "LDL",
                  CRP ~ "CRP",
                  HbA1C ~ "HbA1C",
                  BMI ~ "BMI",
                  TRUE ~ "FALSE"
                )
              ) %>%
              select(-Conditions)

#select rows where Status is "multiple"
multiple_rows <- plotdatagw[plotdatagw$Status == "multiple", ]

#create a function to generate the combined label
generate_combined_label <- function(row) {
  conditions <- c("CAD", "LDL", "CRP", "HbA1C", "BMI")
  combined_labels <- conditions[row]
  return(paste(combined_labels, collapse = " + "))
}

#apply the function to create the combined labels
multiple_rows$Status <- apply(multiple_rows[, c("CAD", "LDL", "CRP", "HbA1C", "BMI")], 1, generate_combined_label)
plotdatagw[plotdatagw$Status == "multiple", ] <- multiple_rows

#selected protein annotations for plot
annotations_lower <- data.frame(
  x = c(1880853289, 697447638, 2761408808, 1016644164, 2397325449),
  label = c("B3GNT1", "SORCS2", "CD40", "SMAD5", "FES")
)
annotations_upper <- data.frame(
  x = c(542603234, 2435652562, 1826768303),
  label = c("SEMA3F", "IL4R", "DKK3")
)

#plot colocalisation results using genomic position
cplot = ggplot() +
  geom_rect(data = chrbreaks, mapping = aes(ymin = 0, ymax = 7, xmin = start, xmax = end, fill = col)) +
  scale_fill_manual(values = c('lightgrey', 'white')) +
  geom_point(data = plotdatagw[plotdatagw$Status == "None" | plotdatagw$Status == "FALSE", ], mapping = aes(x = BPcum, y = tissue, size = PP.H4.abf), color = 'lavenderblush4') +
  geom_point(data = plotdatagw[plotdatagw$Status != "None" & plotdatagw$Status != "FALSE", ], mapping = aes(x = BPcum, y = tissue, color = Status, size = PP.H4.abf)) +
  xlab("Chromosome position") +
  ylab('Colocalised tissue') +
  scale_x_continuous(label = chrbreaks$chr, breaks = chrbreaks$mid) +
  scale_y_discrete(breaks = unique(plotdatagw$tissue), label = unique(plotdatagw$tissue)) +
  scale_color_manual(values = c('CAD' = 'red', 'LDL' = 'blue3', 'HbA1C' = 'yellow', 'BMI' = 'darkorchid4', 'CRP' = 'mediumseagreen'),
                     breaks = c('CAD', 'LDL', 'HbA1C', 'BMI', 'CRP'),
                     labels = c('CAD', 'LDL', 'HbA1C', 'BMI', 'CRP'),
                     name = "Status") +
  guides(size = guide_legend(title = "Tissue PP.4"),
         color = guide_legend(title = "GWAS colocalisation"),
         fill = 'none') +
  theme_cowplot() +
  theme( 
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "bottom",
    legend.direction = "horizontal") +
  geom_text(data = annotations_lower, aes(label = label, x=x, y = 7.4)) +
  geom_segment(data = annotations_lower, aes(x = x, xend = x, y=7.3, yend = 7.1),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.5) +
  geom_text(data = annotations_upper, aes(label = label, x=x, y = 7.54)) +
  geom_segment(data = annotations_upper, aes(x = x, xend = x, y=7.45, yend = 7.1),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.5)

#export plot
cpfn <- paste("~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/total_STARNET_QTL_coloc_H40.5_plot.png")
ggsave(filename=cpfn, plot=cplot, units="cm", width=25, height=22, dpi=300)