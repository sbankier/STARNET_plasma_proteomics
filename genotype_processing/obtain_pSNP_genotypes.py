import pandas as pd
from tqdm import tqdm

#import protein and SNP reference files
pro_tss=pd.read_csv('~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_geneposition.tsv', delimiter='\t')
snp_ref=pd.read_csv('~/sean/resources/dbSNP_00-common_all_human_9606_b151_GRCh37p13_rsids.tsv', delimiter='\t', dtype={'chromosome': str})

#obtain SNPs within threshold of TSS
th=500000
gene_l=pro_tss.drop_duplicates('Gene_name')['Gene_name'].values

cis_df=[]
for g in tqdm(gene_l):
	g_tss=pro_tss.loc[pro_tss['Gene_name'] == g]['Transcription start site (TSS)'].values[0]
	g_ch=pro_tss.loc[pro_tss['Gene_name'] == g]['Chromosome/scaffold name'].values[0]
	up_lim=g_tss+th
	dn_lim=g_tss-th
	snp_filt=snp_ref.loc[(snp_ref['chromosome'] == g_ch) & (snp_ref['position'] >= dn_lim) & (snp_ref['position'] <= up_lim)].copy()
	snp_filt['Gene_name'] = g
	cis_df.append(snp_filt[['Gene_name', 'rs_number', 'chromosome', 'position', 'alt']])

cis_snps=pd.concat(cis_df)

#output dataframe of ALL SNPs for EACH gene that are within threshold of TSS
cis_snps.to_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/STARNET_Olink_cis-SNPs_500Kb.tsv', sep='\t', index=None)

#read in STARNET genotypes in chunks and filter on rsid
pos=cis_snps[['rs_number', 'chromosome', 'position']].drop_duplicates()

df_list=[]
for cnk in pd.read_csv('~/STARNET/genotype/genotype_all/STARNET_genotypedosageMAF05_annotate.txt',
        sep='\t', iterator=True, chunksize=10000, dtype={'chromosome': str}):
        mer=pd.merge(pos, cnk, on=['chromosome', 'position'], how='inner')
        df_list.append(mer)
df=pd.concat(df_list)

#round imputed genotypes to whole number
df_imp=df.round()

#export genotypes
cis_geno=pd.merge(cis_snps[['Gene_name', 'rs_number']], df_imp, on='rs_number', how='inner')
cis_geno.to_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb.tsv', sep='\t', index=None)