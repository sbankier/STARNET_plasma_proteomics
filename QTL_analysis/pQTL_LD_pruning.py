import pandas as pd

#import significant cis-pQTLs
pq=pd.read_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_5FDR.tsv', delimiter='\t')

#import cis-pQTL genotypes
cnk_l=[]
gfn='~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_adj-age-sex-year_samples_unique_SNPs.tsv'
psnp=pq['rs_number'].tolist()
for cnk in pd.read_csv(gfn, delimiter='\t',chunksize=10000):
	geno_cnk=cnk.loc[cnk['rs_number'].isin(psnp)]
	cnk_l.append(geno_cnk)
geno=pd.concat(cnk_l)
pq_geno=pd.merge(pq[['rs_number', 'Gene_name']], geno, on='rs_number', how='inner')

#get unique gene names
gene_un=pq.drop_duplicates('Gene_name')['Gene_name'].tolist()

#get independent pQTLs for all proteins
ind_pq_l=[]
for gene in gene_un:

    #filter for protein and obtain correlation matrix of all pQTLs 
    pro_geno=pq_geno.loc[pq_geno['Gene_name'] == gene].drop('Gene_name', axis=1)
    pro_geno.set_index('rs_number', inplace=True)
    corr=pro_geno.T.corr(method='pearson')

    #get Pearson R2
    corr2=corr**2

    #get corresponding p-values
    pro_pq=pq.loc[pq['Gene_name'] == gene]
    pro_pv=pro_pq.set_index('rs_number')['p-value']

    #filter for SNPs
    independent_snps=[]
    remaining_snps=list(corr2.columns)
    while remaining_snps:

    	#select SNP with lowest p-value and update independent SNP list
        current_snp = min(remaining_snps, key=lambda rs_number: pro_pv[rs_number])
        independent_snps.append(current_snp)

        #remove all correlated SNPs from list
        high_corr_snps = corr2.loc[current_snp, corr2.loc[current_snp] >= 0.8].index.tolist()
        remaining_snps = [snp for snp in remaining_snps if snp not in high_corr_snps]

    #filter pQTL summary statistics for independent SNPs
    ind_pq=pq.loc[(pq['Gene_name'] == gene) & (pq['rs_number'].isin(independent_snps))]
    ind_pq_l.append(ind_pq)

#combine as full dataframe
ind_pq_full=pd.concat(ind_pq_l)

#get LD pruned genotypes
ind_pq_geno=pd.merge(ind_pq_full[['rs_number', 'Gene_name']], pq_geno, on=['rs_number', 'Gene_name'], how='inner')

#export LD pruned pQTLs and genotypes
ind_pq_full.to_csv('~/sean/INTRePID/STARNET/proteins/pQTLs/Olink2018-19_combined/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_20180383_20181235_cases_combined_cis-pQTLs_5FDR_LD_pruned.tsv', sep='\t', index=None)
ind_pq_geno.to_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_cis-pQTLs_5FDR_LD_pruned.tsv', sep='\t', index=None)

