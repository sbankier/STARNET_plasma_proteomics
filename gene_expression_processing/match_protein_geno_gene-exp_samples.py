import pandas as pd
import numpy as np
import sys
sys.path.append("~/sean/resources/")
from common_functions import flip_df

"""Take filtered gene expression files and filter so that sample IDs correspond to 
protein expression  and genotype files"""

#import protein data and obtain samples
pro=pd.read_csv('~/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year.tsv', delimiter='\t')
pro_id=pro[['starnet_id']]
pro_id.columns=['marker_id']
pro_id['marker_id'] = pro_id['marker_id'].astype(str)

#import genotype data and obtain samples
geno=pd.read_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_adj-age-sex-year_samples_unique_SNPs.tsv', delimiter='\t')
geno_flip=flip_df(geno)
geno_id=geno_flip[['rs_number']]
geno_id.columns=['marker_id']
geno_id['marker_id'] = geno_id['marker_id'].astype(str)

#get common samples between genotypes and proteins
pg_sam=pd.merge(pro_id, geno_id, on='marker_id', how='inner')

tis=['LIV', 'SKLM', 'AOR', 'MAM', 'VAF', 'SF', 'Blood']

for t in tis:
	#import gene expression data and obtain samples
	exp=pd.read_csv('~/STARNET/sean/eQTL_discovery/expression/'+t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN', delimiter='\t')
	exp_flip=flip_df(exp)
	exp_id=exp_flip[['Gene stable ID']]
	exp_id.columns=['marker_id']
	exp_id['marker_id'] = exp_id['marker_id'].astype(str)

	#get common samples between genotypes, proteins and gene expression
	pge_sam=pd.merge(pg_sam, exp_id, on='marker_id', how='inner')

	#filter protein, genotype and gene expression with common samples
	pro['starnet_id'] = pro['starnet_id'].astype(str)
	pro_pge=pd.merge(pge_sam, pro, left_on='marker_id', right_on='starnet_id', 
	how='inner').drop(['marker_id', 'status', 'Sex', 'Age', 'year'], axis=1)
	pro_pge_flip=flip_df(pro_pge)

	geno_flip['rs_number']=geno_flip['rs_number'].astype(str)
	geno_pge=pd.merge(pge_sam, geno_flip, left_on='marker_id', right_on='rs_number', 
	how='inner').drop('marker_id', axis=1)
	geno_pge_flip=flip_df(geno_pge)

	exp_pge=pd.merge(pge_sam, exp_flip, left_on='marker_id', right_on='Gene stable ID',
	how='inner').drop('marker_id', axis=1)
	exp_pge_flip=flip_df(exp_pge)

	#export filtered datasets
	pro_pge_flip.to_csv('~/sean/INTRePID/STARNET/proteins/adjusted/geno_pro_exp_match/'+t+'_Olink_20180383_20181235_full_combined_adj-age-sex-year_geno-pro-exp_match.tsv', sep='\t', index=None)
	geno_pge_flip.to_csv('~/sean/INTRePID/STARNET/genotypes/geno_pro_exp_match/'+t+'_genotypedosageMAF05_Olink_cis-SNPs_500Kb_geno-pro-exp_match.tsv', sep='\t', index=None)
	exp_pge_flip.to_csv('~/sean/INTRePID/STARNET/gene_expression/geno_pro_exp_match/'+t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN_geno-pro-exp_match.tsv', sep='\t', index=None)

	#mark complete
	print(t+' complete')