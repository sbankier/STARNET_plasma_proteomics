import pandas as pd
import numpy as np
from functools import reduce
import sys
sys.path.append("~/sean/resources/")
from common_functions import flip_df

"""Filter genotype file so that sample IDs correspond to protein expression file"""

#define tissues and genotype file
pro=pd.read_csv('~/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year.tsv', delimiter='\t')
pro_id=pro[['starnet_id']]
pro_id.columns=['marker_id']
pro_id['marker_id'] = pro_id['marker_id'].astype(str)

#import genotype sample IDs
sam=pd.read_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_genotype_sample_IDs.tsv', delimiter='\t')

#import defined genotype file and process in chunks
geno_pro=[]
for geno in pd.read_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb.tsv',
        sep='\t', iterator=True, chunksize=10000):

        geno_info=geno[['rs_number', 'Gene_name', 'chromosome', 'position', 'ref', 'alt']].reset_index(drop=True)

        #flip dataframe to make samples columns
        geno=geno.drop(['Gene_name', 'chromosome', 'position', 'ref', 'alt'], axis=1)
        genot=geno.T
        genot=genot.reset_index()
        genot.columns=genot.iloc[0]
        genot=genot.iloc[1:]

        #filter for cases or controls and merge filtered IDs with genotypes
        geno_id=pd.merge(genot, sam[['marker_id']], left_on='rs_number',right_on='marker_id', how='inner').drop('rs_number', axis=1)

        #get genotype info and add to filtered dataframe
        mer=pd.merge(pro_id,geno_id,on='marker_id',how='inner')
        mer_t=flip_df(mer)
        mer_t.columns.name = None 
        mer_in=mer_t.reset_index(drop=True)
        geno_pro.append(pd.concat([geno_info, mer_in], axis=1))

mer_anno=pd.concat(geno_pro).drop('marker_id', axis=1)
mer_anno.to_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_adj-age-sex-year_samples.tsv', sep='\t', index=None)

#output unique SNPs
mer_un=mer_anno.drop_duplicates('rs_number').drop(['Gene_name', 'chromosome', 'position', 'ref', 'alt'], axis=1)
mer_un.to_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_adj-age-sex-year_samples_unique_SNPs.tsv', sep='\t', index=None)

#orientate protein dataframe
geno_mk=flip_df(mer_anno)[['rs_number']]
geno_mk.columns=['starnet_id']
geno_id=geno_mk.iloc[5:].astype(str)
pro_df=pro.drop(['status', 'Sex', 'Age', 'year'], axis=1)
pro_df['starnet_id']=pro_df['starnet_id'].astype(str)
pro_filt=pd.merge(pro_df, geno_id, on='starnet_id', how='inner')
pro_ori=flip_df(pro_filt)

#export orientated protein dataframe
pro_ori.to_csv('~/sean/INTRePID/STARNET/proteins/adjusted/20180383_20181235_cases_combined_adj-age-sex-year_genotypedosageMAF05.tsv', sep='\t', index=None)

#filter for unique proteins
ref=pd.read_csv('~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv', delimiter='\t')
ref_filt=ref.drop_duplicates('Gene_name')
pro_ref_filt=pd.merge(pro_ori, ref_filt[['OlinkID']], left_on='starnet_id', right_on='OlinkID', how='inner').drop('OlinkID', axis=1)
pro_ref_filt.to_csv('~/sean/INTRePID/STARNET/proteins/adjusted/20180383_20181235_cases_combined_adj-age-sex-year_genotypedosageMAF05_unique_proteins.tsv', sep='\t', index=None)
