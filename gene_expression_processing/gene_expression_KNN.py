import pandas as pd
from sklearn.impute import KNNImputer

#import data
exp_df=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/expr_recast.csv', delimiter=',')

#remove additional columns and transpose so that samples become rows
exp_mat=exp_df.drop(['Unnamed: 0', 'tissue', 'transcript_id'], axis=1).T

#impute
imputer = KNNImputer(n_neighbors=20)
imputed = imputer.fit_transform(exp_mat)

#restore dataframe structure
gene_lab = exp_df['transcript_id'].str.extract(r'(ENSG\d+)\.')[0]
exp_imp = pd.DataFrame(
    imputed,
    index=exp_mat.index,
    columns=gene_lab
).T.reset_index()
exp_imp['tissue']=exp_df['tissue']
exp_imp.rename(columns={0: 'Gene stable ID'}, inplace=True)

#export
exp_imp.to_csv('~/sean/INTRePID/STARNET/gene_expression/expr_recast_imp.tsv', sep='\t', index=None)
