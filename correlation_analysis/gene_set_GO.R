suppressMessages({
    library(dplyr)
    library(data.table)
    library(gprofiler2, lib="~/temp/Rlib/")
})

#import gene sets
cross_unique <- fread("~/sean/INTRePID/STARNET/gene_expression/eigengenes/protein_cross_t_unique_genes_5FDR.tsv", data.table=FALSE)
spec_unique <- fread("~/sean/INTRePID/STARNET/gene_expression/eigengenes/protein_spec_t_unique_genes_5FDR.tsv", data.table=FALSE)
pg_unique <- fread("~/sean/INTRePID/STARNET/gene_expression/eigengenes/protein_pg_t_unique_genes_5FDR.tsv", data.table=FALSE)

#get custom background gene set
gene_ref <- fread("~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv", data.table=FALSE)
bg <- unique(gene_ref$Gene_name)

#function to perform GO enrichment using gProfiler
target_go <- function(x, lab) {
   
    #run gprofiler2
    gostres <- gost(query = x$Gene_name, 
        organism = "hsapiens", 
        custom_bg = bg, 
        correction_method="fdr", 
        evcodes=TRUE, 
        sources=c("GO:MF", "GO:BP", "GO:CC", "KEGG", "REAC", "WP", "HP"))
    res <- gostres$result

    #remove columns, annotate and reorder
    go_ano <- subset(res, select = -c(parents, query, significant))
    go_ano$gene_set <- lab
    go_out <- go_ano[, c("gene_set", "term_id", "term_name", "p_value", "term_size", "query_size", "intersection_size", "precision", "recall", "source", "effective_domain_size", "source_order", "intersection")]
    return(go_out)

}

#run goprofiler2 for each gene set
cross_go <- target_go(cross_unique, "cross")
spec_go <- target_go(spec_unique, "specific")
pg_go <- target_go(pg_unique, "transcript")

#export results
fwrite(cross_go, "~/sean/INTRePID/STARNET/gene_expression/eigengenes/protein_cross_t_unique_genes_5FDR_GO.tsv", sep="\t", row.names=FALSE, quote=FALSE)
fwrite(spec_go, "~/sean/INTRePID/STARNET/gene_expression/eigengenes/protein_spec_t_unique_genes_5FDR_GO.tsv", sep="\t", row.names=FALSE, quote=FALSE)
fwrite(pg_go, "~/sean/INTRePID/STARNET/gene_expression/eigengenes/protein_pg_t_unique_genes_5FDR_GO.tsv", sep="\t", row.names=FALSE, quote=FALSE)