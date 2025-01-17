import pandas as pd

#import 5% FDR filtered pQTLs and protein reference
pq=pd.read_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_5FDR.tsv', delimiter='\t')
ref=pd.read_csv('~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv', delimiter='\t')

#annotate proteins with Ensembl ID
pq.columns=['rs_number', 'OlinkID', 'Protein_name', 'UniProt', 'Gene stable ID', 'Gene_name', 'pQTL_statistic', 'pQTL_p-value', 'pQTL_FDR', 'pQTL_beta', 'pQTL_beta_se']

#match pQTLs and eQTLs
tis=['LIV', 'SKLM', 'AOR', 'MAM', 'VAF', 'SF', 'Blood']
qtl=[]
for t in tis:

	#import eQTLs and filter for 5% FDR
	eq=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/gene_expression/eQTLs/'+t+'_STARNET_cases_matched_QTLs_cis-eQTLs_500Kb.tsv', delimiter='\t')
	eq_filt=eq.loc[eq['FDR'] <= 0.05]
	eq_filt.columns=['rs_number', 'Gene stable ID', 'Gene_name', 'eQTL_statistic', 'eQTL_p-value', 'eQTL_FDR', 'eQTL_beta', 'eQTL_beta_se']

	#match common QTLs by tissue 
	qt_mer=pd.merge(pq, eq_filt, on=['rs_number', 'Gene stable ID', 'Gene_name'], how='inner')
	
	#format columns
	qt_mer['tissue']=t
	qt_out=qt_mer[['tissue', 'rs_number', 'OlinkID', 'UniProt', 'Protein_name', 'Gene_name', 'Gene stable ID', 
	'pQTL_statistic', 'pQTL_p-value', 'pQTL_FDR', 'pQTL_beta', 'pQTL_beta_se', 
	'eQTL_statistic', 'eQTL_p-value', 'eQTL_FDR', 'eQTL_beta', 'eQTL_beta_se']]
	qtl.append(qt_out)

#concatenate to complete datframe and export
qt_full=pd.concat(qtl)
qt_full.to_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_FDR5_eQTL_overlap.tsv', sep='\t', index=None)