import pandas as pd

#import common QTLs and genotypes 
geno=pd.read_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb.tsv', delimiter='\t')

#Obtain tissue sample filtered genotype files 
tis=['LIV', 'SKLM', 'MAM', 'AOR', 'VAF', 'SF', 'Blood']
for t in tis:
	#import STARNET tissue gene expression 
	exp=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/STARNET.'+t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.tsv', delimiter='\t')

	#obtain tissues and remove tissue prefix
	pref=t+'_'
	new_ids = exp.drop('id', axis=1).columns.str.replace(pref, '')

	#filter genotype file with tissue samples
	new_cols=['rs_number']
	new_cols.extend(new_ids)
	new_geno=geno.filter(new_cols)
	geno_out=new_geno.drop_duplicates('rs_number')

	#export filtered genotypes
	if len(exp.columns) == len(geno_out.columns):
		geno_out.to_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_tissues/QTL_genotypes/'+t+'_STARNET_genotypedosageMAF05_cases_matched_QTLs_cis-SNPs_500Kb.tsv', sep='\t', index=None)
		print(t+' done')
	else:
		print("WARNING: column lengths do not match")