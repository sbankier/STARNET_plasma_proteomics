import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

#import Olink proteomic and phenotype data 
star=pd.read_csv('~/sean/INTRePID/STARNET/proteins/STARNET_Olink_20180383_20181235_full_combined_pheno.tsv', delimiter='\t').dropna()
star_adj=pd.read_csv('~/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_full_combined_adj-age-sex-year.tsv', delimiter='\t').dropna()

def pca_plt(ct, pca1, pca2, save, title):

    #apply standard scaler to dataset features
    ct_x = ct.drop(['starnet_id', 'status', 'Sex', 'Age', 'year'], 1).dropna().values
    ct_xt = StandardScaler().fit_transform(ct_x)

    #reduce data across 2 principal components
    ct_pca = PCA(n_components=2)
    ct_principalComponents = ct_pca.fit_transform(ct_x)
    ct_principalDf = pd.DataFrame(data = ct_principalComponents
                 , columns = ['principal component 1', 'principal component 2'])

    #remove outliers from PCA
    if pca1 == 0:
        pass
    else:
        ct_principalDf=ct_principalDf.sort_values('principal component 1', ascending=False).iloc[pca1:]

    if pca1 == 0:
        pass
    else:
        ct_principalDf=ct_principalDf.sort_values('principal component 2', ascending=False).iloc[pca2:]

    #plot PCA and save output
    sns.lmplot(x='principal component 1', y='principal component 2', fit_reg=False,
                 palette='mako_r', data=ct_principalDf.dropna(), height=7)
    plt.tight_layout()
    if save == True:
        plt.savefig('/Home/ii/seanb/sean//INTRePID/STARNET/proteins/'+title+'.png')
    plt.show()

pca_plt(star, 0, 0, True, 'STARNET_Olink_20180383_20181235_full_combined_PCA')
pca_plt(star_adj, 0, 0, True, 'STARNET_Olink_20180383_20181235_full_combined_adj-age-sex-year_PCA')