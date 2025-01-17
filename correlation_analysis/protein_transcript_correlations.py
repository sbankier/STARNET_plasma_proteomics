import pandas as pd
import numpy as np
import math
from scipy import stats
from multipy import fdr
import sys

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

#identify proteins with a shared transcript
genes=list(set(pro_ref['Gene stable ID']).intersection(exp_ref['Gene stable ID']))

#obtain correlations for all protein with a corresponding gene
pr_l=[]
pp_l=[]
for g in genes:
	exp_g=exp_ref.loc[exp_ref['Gene stable ID'] == g].drop(['Gene stable ID', 'OlinkID'], axis=1).T
	exp_g.columns=['exp']
	pro_g=pro_ref.loc[pro_ref['Gene stable ID'] == g].drop(['Gene stable ID', 'OlinkID'], axis=1).T
	pro_g.columns=['pro']
	sp_stat=stats.spearmanr(exp_g['exp'], pro_g['pro'])
	pr=sp_stat[0]
	pp=sp_stat[1]
	pr_l.append(pr)
	pp_l.append(pp)

full_cor=pd.DataFrame({'Gene stable ID': genes, 'SpearmanR': pr_l, 'Spearmanpv': pp_l})

#obtain -log10 p-values for Spearman correlations
lpvl=[]
for x in full_cor['Spearmanpv']:
    lpv=-math.log10(x)
    lpvl.append(lpv)
full_cor['Spearmanlogpv']=lpvl

#adjust for multiple testing using all matched proteins
qv=fdr.qvalue(full_cor['Spearmanpv'].values)
full_cor['Spearmanqv']=qv[1]

#annotate output
ref_cor=pd.merge(ref_un, full_cor, on='Gene stable ID', how='inner')
ref_cor['tissue']=t

#export output
ref_cor.to_csv('STARNET/gene_expression/pro_exp_match/correlations/'+t+'_STARNET_pro-exp_correlations.tsv', sep='\t', index=None)