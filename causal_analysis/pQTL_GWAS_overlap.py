import pandas as pd

#import cis-pQTLs and format
pq=pd.read_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_P1.tsv', delimiter='\t')

#import SNP reference
snp_ref=pd.read_csv('~/sean/resources/dbSNP_00-common_all_human_9606_b151_GRCh37p13_rsids.tsv', 
    delimiter='\t', dtype={'chromosome': str})

#import GWAS summary statistics
cad=pd.read_csv('~/sean/INTRePID/GWAS/CARDIoGRAM/GCST90132314_buildGRCh37.tsv', delimiter='\t', dtype={'chromosome': str})
ldl=pd.read_csv('~/sean/INTRePID/GWAS/LDL/GCST90239658_buildGRCh37.tsv', delimiter='\t', dtype={'chromosome': str})
dbp=pd.read_csv('~/sean/INTRePID/GWAS/DBP/GCST90000059_buildGRCh37.tsv', delimiter='\t', dtype={'chromosome': str})
hba=pd.read_csv('~/sean/INTRePID/GWAS/HbA1C/GCST90019509_buildGRCh37.tsv', delimiter='\t', dtype={'chromosome': str, 'variant_id': str})
crp=pd.read_csv('~/sean/INTRePID/GWAS/CRP/GCST90019499_buildGRCh37.tsv', delimiter='\t', dtype={'chromosome': str})
bmi_rs=pd.read_csv('~/sean/INTRePID/GWAS/GIANT/fat-distn.giant.ukbb.meta-analysis.bmi.combined_rsID.tsv', delimiter='\t', dtype={'chromosome': str})

#annotate GWAS stats with rsIDS
cad_rs=pd.merge(cad, snp_ref[['rs_number', 'chromosome', 'position']], left_on=['chromosome', 'base_pair_location'], right_on=['chromosome', 'position'], how='inner')
cad_rs=cad_rs.loc[cad_rs['meta_analysis'] == 'Cardiogram']
ldl_rs=pd.merge(ldl, snp_ref[['rs_number', 'chromosome', 'position']], left_on=['chromosome', 'base_pair_location'], right_on=['chromosome', 'position'], how='inner')
dbp_rs=pd.merge(dbp, snp_ref[['rs_number', 'chromosome', 'position']], left_on=['chromosome', 'base_pair_location'], right_on=['chromosome', 'position'], how='inner')
crp_rs=pd.merge(crp, snp_ref[['rs_number', 'chromosome', 'position']], left_on=['chromosome', 'base_pair_location'], right_on=['chromosome', 'position'], how='inner')
hba_rs=pd.merge(hba, snp_ref[['rs_number', 'chromosome', 'position']], left_on=['chromosome', 'base_pair_location'], right_on=['chromosome', 'position'], how='inner')

#export rsID annotated GWAS results
cad_rs.to_csv('~/sean/INTRePID/GWAS/CARDIoGRAM/GCST90132314_buildGRCh37_rsID.tsv', sep='\t', index=None)
ldl_rs.to_csv('~/sean/INTRePID/GWAS/LDL/GCST90239658_buildGRCh37_rsID.tsv', sep='\t', index=None)
dbp_rs.to_csv('~/sean/INTRePID/GWAS/DBP/GCST90000059_buildGRCh37_rsID.tsv', sep='\t', index=None)
crp_rs.to_csv('~/sean/INTRePID/GWAS/CRP/GCST90019499_buildGRCh37_rsID.tsv', sep='\t', index=None)
hba_rs.to_csv('~/sean/INTRePID/GWAS/HbA1C/GCST90019509_buildGRCh37_rsID.tsv', sep='\t', index=None)

#format GWAS results
cad_filt=cad_rs[['rs_number', 'beta', 'standard_error', 'p_value']].drop_duplicates()
ldl_filt=ldl_rs[['rs_number', 'beta', 'standard_error', 'p_value']].drop_duplicates()
dbp_filt=dbp_rs[['rs_number', 'beta', 'standard_error', 'p_value']].drop_duplicates()
hba_filt=hba_rs[['rs_number', 'beta', 'standard_error', 'p_value']].drop_duplicates()
crp_filt=crp_rs[['rs_number', 'beta', 'standard_error', 'p_value']].drop_duplicates()
bmi_filt=bmi_rs[['rs_number', 'Effect', 'StdErr', 'P-value']].drop_duplicates()
bmi_filt.columns=['rs_number', 'beta', 'standard_error', 'p_value']

#filter combined results for significance
def com_res(gw):
    gw.columns=['rs_number', 'GWAS_beta', 'GWAS_se', 'GWAS_p-value']
    pq_gw=pd.merge(pq, gw, on='rs_number', how='inner')
    pq_gw_sig=pq_gw.loc[(pq_gw['FDR'] <= 0.05) & (pq_gw['GWAS_p-value'] <= 5e-8)]
    
    return pq_gw_sig

#obtain combined pQTL/ GWAS results
pq_cad_sig=com_res(cad_filt)
pq_ldl_sig=com_res(ldl_filt)
pq_dbp_sig=com_res(dbp_filt)
pq_hba_sig=com_res(hba_filt)
pq_crp_sig=com_res(crp_filt)
pq_bmi_sig=com_res(bmi_filt)

#export combined filtered results
pq_cad_sig.to_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/GWAS/STARNET_Olink_20180383_20181235_cases_combined_cis-pQTLs_5FDR_CAD_overlap.tsv', sep='\t', index=None)
pq_ldl_sig.to_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/GWAS/STARNET_Olink_20180383_20181235_cases_combined_cis-pQTLs_5FDR_LDL_overlap.tsv', sep='\t', index=None)
pq_dbp_sig.to_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/GWAS/STARNET_Olink_20180383_20181235_cases_combined_cis-pQTLs_5FDR_DBP_overlap.tsv', sep='\t', index=None)
pq_hba_sig.to_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/GWAS/STARNET_Olink_20180383_20181235_cases_combined_cis-pQTLs_5FDR_HbA1C_overlap.tsv', sep='\t', index=None)
pq_crp_sig.to_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/GWAS/STARNET_Olink_20180383_20181235_cases_combined_cis-pQTLs_5FDR_CRP_overlap.tsv', sep='\t', index=None)
pq_bmi_sig.to_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/coloc/GWAS/STARNET_Olink_20180383_20181235_cases_combined_cis-pQTLs_5FDR_BMI_overlap.tsv', sep='\t', index=None)