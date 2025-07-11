import pandas as pd
from functools import reduce
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#import protein expression data
pro=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year_unique_proteins.tsv', delimiter='\t')

#import GRN linked proteins
grn=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_labelled.tsv', delimiter='\t')

#count number of significant associated proteins
grn_sig=grn.loc[(grn['spearman_qv'] <= 0.05) & (grn['present_in_GRN'] == True)][['OlinkID', 'Gene_name', 'network_id']].drop_duplicates()
freq=grn_sig['network_id'].value_counts().reset_index()
freq.columns=['network_id', 'counts']

#filter based on number of associated proteins
grn_freq=pd.merge(grn_sig, freq, on='network_id', how='inner')
grn_filt=grn_freq.loc[grn_freq['counts'] >= 4]

#get eigenproteins from GRN linked proteins
def get_eign(net):

	#filter protein expression for GRN linked proteins
	net_grn=grn_filt.loc[grn_filt['network_id'] == net]
	pro_grn=pro[['starnet_id']+net_grn['OlinkID'].tolist()]

	#apply standard scaler to dataset features
	ct_x=pro_grn.drop('starnet_id', axis=1).values
	ct_xt=StandardScaler().fit_transform(ct_x)

	#Obtain first prinicipal component
	ct_pca=PCA(n_components=1)
	ct_principalComponents=ct_pca.fit_transform(ct_x)

	#get variance explained
	var_exp=ct_pca.explained_variance_ratio_

	#format results as dataframe
	ct_df=pd.DataFrame({'starnet_id': pro_grn['starnet_id'], 'GRN_'+str(net): ct_principalComponents.flatten()})
	var_df=pd.DataFrame({'GRN': ['GRN_'+str(net)], 'variance_explained': var_exp})

	return ct_df, var_df

#get network IDs
net_ids=grn_filt.sort_values('counts', ascending=False).drop_duplicates('network_id')['network_id'].tolist()

#get eigenproteins for networks
eign_l=[]
var_l=[]
for x in net_ids:
	eign=get_eign(x)
	eign_l.append(eign[0])
	var_l.append(eign[1])

#combine results
grn_eign=reduce(lambda left, right: pd.merge(left, right, on='starnet_id', how='inner'), eign_l)
grn_var=pd.concat(var_l)

#export results
grn_eign.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/GRN78/GRN_correlated_eigenprotein_inGRN_n4.tsv', sep='\t', index=None)
grn_var.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/GRN78/GRN_correlated_eigenprotein_inGRN_varexp_n4.tsv', sep='\t', index=None)