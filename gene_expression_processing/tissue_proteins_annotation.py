import pandas as pd
from functools import reduce

#import protein-gene reference
ref=pd.read_csv('~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv', delimiter='\t')

#import tissue genes
liv=pd.read_csv("~/sean/INTRePID/STARNET/gene_expression/LIV.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN", 
	delimiter='\t', usecols=[0]).assign(LIV=True)
sklm=pd.read_csv("~/sean/INTRePID/STARNET/gene_expression/SKLM.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN", 
	delimiter='\t', usecols=[0]).assign(SKLM=True)
aor=pd.read_csv("~/sean/INTRePID/STARNET/gene_expression/AOR.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN", 
	delimiter='\t', usecols=[0]).assign(AOR=True)
mam=pd.read_csv("~/sean/INTRePID/STARNET/gene_expression/MAM.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN", 
	delimiter='\t', usecols=[0]).assign(MAM=True)
vaf=pd.read_csv("~/sean/INTRePID/STARNET/gene_expression/VAF.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN", 
	delimiter='\t', usecols=[0]).assign(VAF=True)
sf=pd.read_csv("~/sean/INTRePID/STARNET/gene_expression/SF.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN", 
	delimiter='\t', usecols=[0]).assign(SF=True)
blood=pd.read_csv("~/sean/INTRePID/STARNET/gene_expression/Blood.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN", 
	delimiter='\t', usecols=[0]).assign(Blood=True)

#assemble tissue genes as single dataframe
tis_df=[liv, sklm, aor, mam, vaf, sf, blood]
tis_full=reduce(lambda left,right: pd.merge(left,right,on='Gene stable ID',how='outer'), tis_df)

#filter based on protein matched genes
ref_tis=pd.merge(ref, tis_full, on='Gene stable ID', how='left')

#export updated reference
ref_tis.to_csv('~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart_STARNET_tissues.tsv', sep='\t', index=None)