suppressMessages({
    library(dplyr)
    library(data.table)
    library(gprofiler2, lib="~/temp/Rlib/")
})

#import correlation results
cor <- fread("/Home/ii/seanb/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_labelled.tsv", data.table=FALSE)

#filter for network type
cor_ts <- filter(cor, spearman_qv <= 0.05, network_type == "specific", present_in_GRN == TRUE)
cor_ct <- filter(cor, spearman_qv <= 0.05, network_type == "cross", present_in_GRN == TRUE) 

#get custom background gene set
gene_ref <- fread("~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv", data.table=FALSE)
bg <- unique(gene_ref$Gene_name)

#function to perform GO enrichment using gProfiler
target_go <- function(X, cor_sig) {

	#get test set from tissue set
	test_set <- filter(cor_sig, network_id == X)$Gene_name %>% unique()
   
    #run gprofiler2
    gostres <- gost(query = test_set, 
        organism = "hsapiens", 
        custom_bg = bg, 
        correction_method="fdr", 
        evcodes=TRUE, 
        sources=c("GO:MF", "GO:BP", "GO:CC", "KEGG", "REAC", "WP", "HP"))
    
    return(gostres$result)
}

#get network labels
ts_lab <- unique(cor_ts$network_id)
ct_lab <- unique(cor_ct$network_id)

#perform GO enrichment across all tissue sets
go_ts_l <- lapply(X=ts_lab, FUN=target_go, cor_sig=cor_ts)
names(go_ts_l) <- ts_lab
go_ts_df <- bind_rows(go_ts_l, .id = "network_id")

go_ct_l <- lapply(X=ct_lab, FUN=target_go, cor_sig=cor_ct)
names(go_ct_l) <- ct_lab
go_ct_df <- bind_rows(go_ct_l, .id = "network_id")

#get the most common terms across all tissue sets
count_terms <- function(df) {

    #count terms
    term_counts <- data.frame(table(df$term_id))
    colnames(term_counts) <- c('term_id', 'number_of_networks')
    df_count <- merge(df, term_counts, by='term_id')

    #remove columns, annotate and reorder
    df_count_filt <- subset(df_count, select = -c(parents, query, significant))
    df_anno <- df_count_filt[, c("network_id", "term_id", "term_name", "number_of_networks", "p_value", "term_size", "query_size", "intersection_size", "precision", "recall", "source", "effective_domain_size", "source_order", "intersection")]

    return(df_anno)
}

ts_anno <- count_terms(go_ts_df)
ct_anno <- count_terms(go_ct_df)

#export results
fwrite(ts_anno, "~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_tissue_specific_5FDR_GO.tsv", sep="\t", row.names=FALSE, quote=FALSE)
fwrite(ct_anno, "~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_cross_tissue_5FDR_GO.tsv", sep="\t", row.names=FALSE, quote=FALSE)