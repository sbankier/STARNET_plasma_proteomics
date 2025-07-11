import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#import protein expression data
pro=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year_unique_proteins.tsv', delimiter='\t')

#import GRN78 proteins
grn=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/GRN78/Olink_cases_combined_adjusted_eigengene_correlations_GRN78_5FDR.tsv', delimiter='\t')
grn_filt=grn.loc[grn['present_in_GRN'] == True].drop_duplicates()

#filter for GRN78 proteins
pro_grn=pro[['starnet_id']+grn_filt['OlinkID'].tolist()]
pro_grn_nolep=pro_grn.drop('OID00463', axis=1)

def get_eign(exp, labct, labvar):

	#apply standard scaler to dataset features
	ct_x=exp.drop('starnet_id', axis=1).values
	ct_xt=StandardScaler().fit_transform(ct_x)

	#Obtain first 5 prinicipal components
	ct_pca=PCA(n_components=1)
	ct_principalComponents=ct_pca.fit_transform(ct_x)

	#get variance explained
	var_exp=ct_pca.explained_variance_ratio_

	#format results as dataframe
	ct_df=pd.DataFrame({'starnet_id': exp['starnet_id'], labct: ct_principalComponents.flatten()})
	var_df=pd.DataFrame({'GRN': [labvar], 'variance_explained': var_exp})

	return ct_df, var_df

#get eigenproteins
grn_eign=get_eign(pro_grn, 'PC1', '78')
nolep_eign=get_eign(pro_grn_nolep, 'PC1_nolep', '78_nolep')

#combine results
eign=pd.merge(grn_eign[0], nolep_eign[0], on='starnet_id', how='inner')
var=pd.concat([grn_eign[1], nolep_eign[1]])

#add protein expression values to dataframe
ct_full=pd.merge(eign, pro_grn, on='starnet_id', how='inner')

#update annotations
ids=pd.DataFrame({'OlinkID': pro_grn.columns[1:]})
ref=pd.merge(ids, grn[['OlinkID', 'Gene_name']], on='OlinkID', how='inner').drop_duplicates()
ct_full.columns=['starnet_id', 'PC1', 'PC1_nolep']+ref['Gene_name'].tolist()

#export results
ct_full.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/GRN78/GRN78_correlated_eigenprotein_inGRN.tsv', sep='\t', index=None)
var.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/GRN78/GRN78_correlated_eigenprotein_inGRN_varexp.tsv', sep='\t', index=None)