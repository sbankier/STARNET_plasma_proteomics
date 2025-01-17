import pandas as pd
import numpy as np
import sys
from scipy import stats
from itertools import product
from tqdm import tqdm

#set parameters
t=sys.argv[1]

#import protein and gene expression files with matched samples
pro=pd.read_csv('~/sean/INTRePID/STARNET/proteins/adjusted/pro_exp_match/'+t+'_Olink_20180383_20181235_cases_combined_adj-age-sex-year_pro-exp_match.tsv', delimiter='\t')
exp=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/'+t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN_pro-exp_match.tsv', delimiter='\t')

#import protein expression reference for unique genes
ref=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv', delimiter='\t')
ref_un=ref.drop_duplicates('Gene_name')
ref_filt=ref_un[['OlinkID', 'Gene stable ID']]

#annotate protein and gene expression features
pro_ref=pd.merge(ref_filt, pro, left_on='OlinkID', right_on='starnet_id', how='inner').drop('starnet_id', axis=1)
exp_ref=pd.merge(ref_filt, exp, on='Gene stable ID', how='inner')

#format to compute pairwise correlations
pro_cor=pro_ref.set_index('OlinkID').drop('Gene stable ID', axis=1)
exp_cor=exp_ref.set_index('Gene stable ID').drop('OlinkID', axis=1)

#compute pairwise correlations
def pairwise_correlation(g_df, p_df):
    results=[]
    for gene, protein in tqdm(product(g_df.index, p_df.index)):
        r, p = stats.spearmanr(g_df.loc[gene], p_df.loc[protein])
        results.append({'Gene stable ID': gene, 'OlinkID': protein, 'SpearmanR': r, 'Spearmanpv': p})
    return pd.DataFrame(results)

pairwise_corr=pairwise_correlation(exp_cor, pro_cor)

#identify proteins with a shared transcript
ref_ano=ref[['Gene stable ID', 'OlinkID']].drop_duplicates()
ref_ano['mother_gene']=True
pairwise_corr_ref=pd.merge(pairwise_corr, ref_ano, on=['Gene stable ID', 'OlinkID'], how='left').fillna(False)
pairwise_corr_ref['tissue']=t

#annotate results
gene_ref=ref_un[['Gene stable ID', 'Gene_name']]
pro_ref=ref_un[['OlinkID', 'Gene_name']]
pro_ref.columns=['OlinkID', 'Protein_name']
res_gene=pd.merge(pairwise_corr_ref, gene_ref, on='Gene stable ID', how='inner')
res_pro=pd.merge(res_gene, pro_ref, on='OlinkID', how='inner')
res_out=res_pro[['OlinkID', 'Protein_name', 'Gene stable ID', 'Gene_name', 'tissue', 'mother_gene','SpearmanR', 'Spearmanpv']]

#export results
res_out.to_csv('~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/full_correlations/'+t+'_STARNET_pro-exp_correlations_pairwise.tsv', sep='\t', index=None)