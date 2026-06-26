import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from functools import reduce
from scipy import stats
import math
from multipy import fdr
from tqdm import tqdm
import sys
sys.path.append("/Home/siv32/qac018/sean/resources/")
from common_functions import flip_df

#import eigengene modules
mod=pd.read_csv('/Home/siv32/qac018/sean/INTRePID/STARNET/gene_expression/eigengenes/modules.csv', delimiter=',')
mod['tissue']=mod['tissue'].str.replace('BLOOD', 'Blood')
mod['Gene stable ID']=mod['ensembl'].str.split('.', n=1).str[0]

#import all gene expression files
exp=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/expr_recast_imp.tsv', delimiter='\t')
exp['tissue']=exp['tissue'].str.replace('BLOOD', 'Blood')

#function to get eigengenes
def get_eigen(mod, net):
	
    #filter for module transcripts and transpose
    grn=mod.loc[mod['clust'] == net][['Gene stable ID', 'tissue']]
    grn_exp=pd.merge(grn, exp, on=['Gene stable ID', 'tissue'], how='inner').dropna(axis=1).drop('tissue', axis=1)
    grn_exp=grn_exp.T.reset_index()
    grn_exp.columns=grn_exp.iloc[0]
    grn_exp=grn_exp.iloc[1:]
    
    #apply standard scaler to dataset features
    ct_x=grn_exp.drop('Gene stable ID', axis=1).values
    ct_xt=StandardScaler().fit_transform(ct_x)
    
    #Obtain first prinicipal component
    ct_pca=PCA(n_components=1)
    ct_principalComponents=ct_pca.fit_transform(ct_x)
    
    #get variance explained
    var_exp=ct_pca.explained_variance_ratio_
    
    #format results as dataframe
    ct_df=pd.DataFrame({'starnet_id': grn_exp['Gene stable ID'], 'GRN_'+str(net): ct_principalComponents.flatten()})
    var_df=pd.DataFrame({'GRN': ['GRN_'+str(net)], 'variance_explained': var_exp})

    return ct_df, var_df

#get list of modules
mod_ids=mod[['clust']].drop_duplicates()['clust'].tolist()

#get eigengenes for selected modules
eig_l=[]
var_l=[]
for net in mod_ids:
    res=get_eigen(mod, net)
    eig_l.append(res[0])
    var_l.append(res[1])

#combine as single dataframes
eig=reduce(lambda left, right: pd.merge(left, right, on='starnet_id', how='outer'), eig_l)
var=pd.concat(var_l)

#import proteins
pro=pd.read_csv('~/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year_unique_proteins.tsv', delimiter='\t')

def get_cor(pro, eig):

	#transpose dataframes
	pro['starnet_id']=pro['starnet_id'].astype(str)
	pro_t=flip_df(pro.drop(['status', 'Sex', 'Age', 'year'], axis=1))
	eig_t=flip_df(eig)

	#get labels
	pro_id=pro_t['starnet_id'].to_list()
	eig_id=eig_t['starnet_id'].to_list()

	#obtain pairwise correlations between proteins and eigenproteins
	cor_dict={'OlinkID':[], 'network_id':[], 'spearman_r':[], 'spearman_pv':[]}
	for pid in tqdm(pro_id):
		for eid in eig_id:
			cor_dict['OlinkID'].append(pid)
			cor_dict['network_id'].append(eid)

			#filter for protein/ eigengene
			pro_f=pro_t.loc[pro_t['starnet_id'] == pid]
			eig_f=eig_t.loc[eig_t['starnet_id'] == eid].dropna(axis=1)

			#identify common samples
			pro_sam=set(pro_f.columns[1:])
			eig_sam=set(eig_f.columns[1:])
			com_sam=['starnet_id']+list(pro_sam.intersection(eig_sam))

			#filter protein and eigengene dataframes for common samples
			pro_v=pro_f[com_sam].values[0][1:]
			eig_v=eig_f[com_sam].values[0][1:]
			
			sp=stats.spearmanr(pro_v, eig_v)
			cor_dict['spearman_r'].append(sp[0])
			cor_dict['spearman_pv'].append(sp[1])

	cor_df=pd.DataFrame(cor_dict)
	cor_df['network_id']=cor_df['network_id'].str.split('_').str[1]

	#annotate from protein reference
	ref=pd.read_csv('~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv', delimiter='\t')
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
	mod_tab_full=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/module_tab.csv', delimiter=',')
	mod_tab=mod_tab_full.iloc[:,0:10]
	mod_tab.rename(columns={'Unnamed: 0': 'network_id'}, inplace=True)

	#get tissue specific networks vs cross-tissue
	spc=mod_tab.loc[mod_tab['purity'] >= 0.95]
	crs=mod_tab.loc[mod_tab['purity'] < 0.95]

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

	return cor_lab

#get correlations from eigengenes with seed genes
cor_full=get_cor(pro, eig)

#identify seed genes in GRN modules
pro_lab=cor_full[['Gene stable ID']].drop_duplicates()
pro_lab['seed_gene']=True
mod_seed=pd.merge(mod, pro_lab, on='Gene stable ID', how='left').fillna(False)

#remove seed genes from GRNs
mod_seed_filt=mod_seed.loc[mod_seed['seed_gene'] == False]

#get list of modules
mod_seed_filt_ids=mod_seed_filt[['clust']].drop_duplicates()['clust'].tolist()

#get eigengenes for modules without seed genes
eig_sd_l=[]
var_sd_l=[]
for net in mod_seed_filt_ids:
    res_sd=get_eigen(mod_seed_filt, net)
    eig_sd_l.append(res_sd[0])
    var_sd_l.append(res_sd[1])

#combine as single dataframes
eig_sd=reduce(lambda left, right: pd.merge(left, right, on='starnet_id', how='outer'), eig_sd_l)
var_sd=pd.concat(var_sd_l)

#get correlations from eigengenes without seed genes
cor_sd=get_cor(pro, eig_sd)

#export eigengenes
eig.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/seed_gene_removal/GRN_eigengenes_expr_recast_imp_with_sg.tsv', sep='\t', index=None)
eig_sd.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/seed_gene_removal/GRN_eigengenes_expr_recast_imp_no_sg.tsv', sep='\t', index=None)

#export variance explained
var.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/seed_gene_removal/GRN_eigengenes_variance_expr_recast_imp_with_sg.tsv', sep='\t', index=None)
var_sd.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/seed_gene_removal/GRN_eigengenes_variance_expr_recast_imp_no_sg.tsv', sep='\t', index=None)

#output results
cor_full.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/seed_gene_removal/Olink_correlations_expr_recast_imp_with_sg.tsv', sep='\t', index=None)
cor_sd.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/seed_gene_removal/Olink_correlations_expr_recast_imp_no_sg.tsv', sep='\t', index=None)
