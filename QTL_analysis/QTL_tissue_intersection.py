import pandas as pd

#import common QTLs
qtl=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_FDR5_eQTL_overlap.tsv', delimiter='\t')
qtl_un=qtl.drop_duplicates(['tissue', 'Gene_name'])

#import LD QTLs
qtl_ld=pd.read_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_eQTL_overlap_LD.tsv', delimiter='\t')
qtl_ld_un=qtl_ld.drop_duplicates(['tissue', 'Gene_name'])

#import colocalisation results
coloc=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/total_STARNET_Olink_cis-SNPs_500Kb_cases_combined_cis-pQTLs_FDR5_eQTL_coloc.tsv', delimiter='\t')
coloc_sig=coloc.loc[coloc['PP.H4.abf'] > 0.5]

#get tissue specific counts of common genes
def get_count(df, label):

    count=df[label].value_counts().to_dict()
    countdf=pd.DataFrame(count.items(), columns=[label, 'frequency'])

    liv_l=[]
    sklm_l=[]
    aor_l=[]
    mam_l=[]
    vaf_l=[]
    sf_l=[]
    bl_l=[]
    for x in countdf[label]:

        if df.loc[(df[label] == x) & (df['tissue'] == 'LIV')].empty:
            liv_l.append(0)
        else:
            liv_l.append(1)

        if df.loc[(df[label] == x) & (df['tissue'] == 'SKLM')].empty:
            sklm_l.append(0)
        else:
            sklm_l.append(1)

        if df.loc[(df[label] == x) & (df['tissue'] == 'AOR')].empty:
            aor_l.append(0)
        else:
            aor_l.append(1)

        if df.loc[(df[label] == x) & (df['tissue'] == 'MAM')].empty:
            mam_l.append(0)
        else:
            mam_l.append(1)

        if df.loc[(df[label] == x) & (df['tissue'] == 'VAF')].empty:
            vaf_l.append(0)
        else:
            vaf_l.append(1)

        if df.loc[(df[label] == x) & (df['tissue'] == 'SF')].empty:
            sf_l.append(0)
        else:
            sf_l.append(1)

        if df.loc[(df[label] == x) & (df['tissue'] == 'Blood')].empty:
            bl_l.append(0)
        else:
            bl_l.append(1)

    countdf['LIV']=liv_l
    countdf['SKLM']=sklm_l
    countdf['AOR']=aor_l
    countdf['MAM']=mam_l
    countdf['VAF']=vaf_l
    countdf['SF']=sf_l
    countdf['Blood']=bl_l

    return countdf

qtl_count=get_count(qtl_un, 'Gene_name')
qtl_ld_count=get_count(qtl_ld_un, 'Gene_name')
coloc_count=get_count(coloc_sig, 'Gene_name')

#export tissue intersections
qtl_count.to_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_cis-pQTL_tissue-intersection_count_20180383_20181235_cases_combined_unique.tsv', sep='\t', index=None)
qtl_ld_count.to_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_cis-pQTL_tissue-intersection_count_20180383_20181235_cases_combined_unique_LD.tsv', sep='\t', index=None)
coloc_count.to_csv('/Home/ii/seanb/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/total_STARNET_Olink_cis-SNPs_500Kb_cases_combined_cis-pQTLs_FDR5_eQTL_coloc_PP0.5_intersection.tsv', sep='\t', index=None)
