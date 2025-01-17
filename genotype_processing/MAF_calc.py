import pandas as pd
import numpy as np

#import genotypes data
geno=pd.read_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb.tsv', delimiter='\t')

#Remove missing data and round values
geno=geno.drop(['Gene_name', 'chromosome', 'position', 'ref', 'alt'], axis=1)
geno=geno.set_index('rs_number') 
geno=geno.apply(np.round)
geno=geno.reset_index() 

#count genotype frequency per SNP
c0=(geno == 0).astype(int).sum(axis=1)
c1=(geno == 1).astype(int).sum(axis=1)
c2=(geno == 2).astype(int).sum(axis=1)
count=pd.concat([c0, c1, c2], axis=1)

#calculate minor allele frequency
MAF_calc=[]
for x in count.values:
    a, b, c = x
    af = (b+(c*2))/((a+b+c)*2)
    if af < 0.5:
        maf = af
    else:
        maf = 1-af
    MAF_calc.append(maf)

#compile output and export
out=pd.DataFrame({'rs_number': geno['rs_number'], 'MAF': MAF_calc})
out.to_csv('~/sean/INTRePID/STARNET/genotypes/STARNET_proteins/STARNET_genotypedosageMAF05_Olink_cis-SNPs_500Kb_MAFs.tsv', sep='\t', index=None)
