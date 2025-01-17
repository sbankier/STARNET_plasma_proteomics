import pandas as pd
import numpy as np
import sys
sys.path.append("~/sean/resources/")
from common_functions import flip_df

"""Take filtered gene expression files and filter so that sample IDs correspond to 
protein expression files"""

#import protein data and obtain samples
pro=pd.read_csv('~/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year.tsv', delimiter='\t')
pro['starnet_id']=pro['starnet_id'].astype(str)

pro_id=pro[['starnet_id']]
pro_id.columns=['marker_id']
pro_id['marker_id'] = pro_id['marker_id'].astype(str)

tis=['LIV', 'SKLM', 'AOR', 'MAM', 'VAF', 'SF', 'Blood']

for t in tis:
	#import gene expression data and obtain samples
	exp=pd.read_csv('~/STARNET/sean/eQTL_discovery/expression/'+t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN', delimiter='\t')
	exp_flip=flip_df(exp)
	exp_id=exp_flip[['Gene stable ID']]
	exp_id.columns=['marker_id']
	exp_id['marker_id'] = exp_id['marker_id'].astype(str)

	#get common samples between genotypes, proteins and gene expression
	pe_sam=pd.merge(pro_id, exp_id, on='marker_id', how='inner')

	#filter protein, genotype and gene expression with common samples
	pro_pe=pd.merge(pe_sam, pro, left_on='marker_id', right_on='starnet_id', 
		how='inner').drop(['marker_id', 'status', 'Sex', 'Age', 'year'], axis=1)
	pro_pe_flip=flip_df(pro_pe)

	exp_pe=pd.merge(pe_sam, exp_flip, left_on='marker_id', right_on='Gene stable ID',
		how='inner').drop('marker_id', axis=1)
	exp_pe_flip=flip_df(exp_pe)

	#export filtered datasets
	pro_pe_flip.to_csv('~/sean/INTRePID/STARNET/proteins/adjusted/pro_exp_match/'+t+'_Olink_20180383_20181235_cases_combined_adj-age-sex-year_pro-exp_match.tsv', sep='\t', index=None)
	exp_pe_flip.to_csv('~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/'+t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN_pro-exp_match.tsv', sep='\t', index=None)

	#mark complete
	print(t+' complete')