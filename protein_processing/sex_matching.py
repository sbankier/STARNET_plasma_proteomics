import pandas as pd
from sklearn.neighbors import NearestNeighbors

#import eigenproteins for males and females
m=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/eigenproteins/sex/GRN_correlated_eigenprotein_inGRN_n4_male.tsv', delimiter='\t')
f=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/eigenproteins/sex/GRN_correlated_eigenprotein_inGRN_n4_female.tsv', delimiter='\t')

#format eigneproteins
m['starnet_id']=m['starnet_id'].astype(str)
f['starnet_id']=f['starnet_id'].astype(str)

#import phenotypes and select for BMI
ph=pd.read_csv('~/sean/INTRePID/STARNET/phenotype/STARNET_pheno_GRN_traits.tsv', delimiter='\t')
ph_bmi=ph[['starnet_id', 'BMI']]

#add BMI to eigenproteins
m_ph=pd.merge(m, ph_bmi, on='starnet_id', how='inner').dropna()
f_ph=pd.merge(f, ph_bmi, on='starnet_id', how='inner').dropna()

#use nearest neighbors to select closest male for each female
nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(m_ph[['BMI']])
distances, indices = nbrs.kneighbors(f_ph[['BMI']])
matched_m = m_ph.iloc[indices.flatten()]

#export matched males
matched_m.drop('BMI', axis=1).to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/eigenproteins/sex/GRN_correlated_eigenprotein_inGRN_n4_male_matched.tsv', sep='\t', index=None)