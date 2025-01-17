import pandas as pd
import numpy as np
import math
from scipy import stats
from multipy import fdr
from tqdm import tqdm
import sys
sys.path.append("~/sean/resources/")
from common_functions import flip_df

#import GRN eigengenes
eig=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/eigengene_mat.csv', delimiter=',')

#obtain randomly permutated eigengenes
#eig_rand_mat=np.random.permutation(eig.values.flatten()).reshape(eig.shape)
#eig_rand=pd.DataFrame(eig_rand_mat, index=eig.index, columns=eig.columns).reset_index()
eig_rand=eig.drop('starnet_id', axis=1).sample(frac=1, random_state=0).reset_index(drop=True)
eig_rand.insert(0, 'starnet_id', eig['starnet_id'])

#import protein expression data
pro=pd.read_csv('~/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year_unique_proteins.tsv', delimiter='\t')

#transpose dataframes
pro['starnet_id']=pro['starnet_id'].astype(str)
pro_t=flip_df(pro.drop(['status', 'Sex', 'Age', 'year'], axis=1))
eig_t=flip_df(eig_rand)

#identify common samples
pro_sam=set(pro_t.columns[1:])
eig_sam=set(eig_t.columns[1:])
com_sam=['starnet_id']+list(pro_sam.intersection(eig_sam))

#filter protein and eigengene dataframes for common samples
pro_com=pro_t[com_sam]
eig_com=eig_t[com_sam]

pro_id=pro_com['starnet_id'].to_list()
eig_id=eig_com['starnet_id'].to_list()

#obtain pairwise correlations between proteins and eigenproteins
cor_dict={'OlinkID':[], 'network_id':[], 'spearman_r':[], 'spearman_pv':[]}
for pid in tqdm(pro_id):
	pro_v=pro_com.loc[pro_com['starnet_id'] == pid].values[0][1:]
	for eid in eig_id:
		cor_dict['OlinkID'].append(pid)
		cor_dict['network_id'].append(eid)
		eig_v=eig_com.loc[eig_com['starnet_id'] == eid].values[0][1:]
		cor=stats.spearmanr(pro_v, eig_v)
		cor_dict['spearman_r'].append(cor[0])
		cor_dict['spearman_pv'].append(cor[1])

cor_df=pd.DataFrame(cor_dict)

#annotate from protein reference
ref=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv', delimiter='\t')
ref_un=ref.drop_duplicates('Gene_name')
cor_ref=pd.merge(ref_un, cor_df, on='OlinkID', how='inner')

#obtain -log10 p-values for Spearman correlations
lpvl=[]
for x in cor_ref['spearman_pv']:
    lpv=-math.log10(x)
    lpvl.append(lpv)
cor_ref['spearman_logpv']=lpvl

#correct for multiple testing
qv=fdr.qvalue(cor_ref['spearman_pv'].values)
cor_ref['spearman_qv']=qv[1]

#import module information
mod_full=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/module_tab.csv', delimiter=',')
mod=mod_full.iloc[:,0:10]
mod.rename(columns={'Unnamed: 0': 'network_id'}, inplace=True)

#get tissue specific networks vs cross-tissue
spc=mod.loc[mod['purity'] >= 0.95]
crs=mod.loc[mod['purity'] < 0.95]

#label network tissues
spc_lab=spc.copy().drop(['network_id', 'mod_size', 'purity'], axis=1).idxmax(axis=1).values
spc['network_type']='specific'
spc['network_tissue']=spc_lab

#label cross-tissue
crs_lab_index=crs.copy().drop(['network_id', 'mod_size', 'purity'], axis=1).apply(lambda row: row[row != 0].index, axis=1)
crs_lab=string_values_list = [', '.join(index) for index in crs_lab_index]
crs['network_type']='cross'
crs['network_tissue']=crs_lab

#combine labeled modules and label correlations
mod_lab=pd.concat([spc, crs])[['network_id', 'mod_size', 'purity', 'network_type', 'network_tissue']]
mod_lab['network_id']=mod_lab['network_id'].astype(str)
cor_lab=pd.merge(cor_ref, mod_lab, on='network_id', how='inner')

#output results
cor_lab.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_random.tsv', sep='\t', index=None)