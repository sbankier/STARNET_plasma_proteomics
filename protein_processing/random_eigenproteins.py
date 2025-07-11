import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#import protein expression data
pro=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year_unique_proteins.tsv', delimiter='\t')
pro_filt=pro.drop(['starnet_id', 'status', 'Sex', 'Age', 'year'], axis=1)

def get_rand(n):
	ct_l=[]
	var_l=[]
	for x in range(0,1000):

		#randomly sample proteins
		rand=pro_filt.sample(n=n, axis=1, random_state=x)

		#apply standard scaler to dataset features
		ct_x=rand.values
		ct_xt=StandardScaler().fit_transform(ct_x)

		#Obtain first prinicipal component
		ct_pca=PCA(n_components=1)
		ct_principalComponents=ct_pca.fit_transform(ct_x)

		#get variance explained
		var_exp=ct_pca.explained_variance_ratio_

		#format results as dataframe
		ct_df=pd.DataFrame({'rand_'+str(x): ct_principalComponents.flatten()})
		var_df=pd.DataFrame({'random_PC1': ['rand_'+str(x)], 'variance_explained': var_exp})

		#add to list
		ct_l.append(ct_df)
		var_l.append(var_df)

	#assemble as dataframes
	ct_out=pd.concat(ct_l, axis=1)
	var_out=pd.concat(var_l, axis=0)

	#format output
	ct_out = pd.concat([pro[['starnet_id']], ct_out], axis=1)

	#export results
	ct_out.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/eigenproteins/random/random_n'+str(n)+'_eigenproteins.tsv', sep='\t', index=None)
	var_out.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/eigenproteins/random/random_n'+str(n)+'_eigenproteins_varexp.tsv', sep='\t', index=None)

#get random eigenproteins for differnt numbers of proteins
get_rand(5)
get_rand(10)
get_rand(15)
get_rand(20)
get_rand(25)
get_rand(30)