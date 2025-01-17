import pandas as pd

#import protein reference
ref=pd.read_csv('~/sean/INTRePID/STARNET/proteins/Olink_ref_target96_biomart.tsv', delimiter='\t')
ref_un=ref[['Gene stable ID']].drop_duplicates()

#filter gene expression transcripts for corresponding proteins
tis=['LIV', 'SKLM', 'AOR', 'MAM', 'VAF', 'SF', 'Blood']
for t in tis:
	ge=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/'+t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN', delimiter='\t')
	ge_filt=pd.merge(ge, ref_un, on='Gene stable ID', how='inner')
	ge_filt.to_csv('~/sean/INTRePID/STARNET/gene_expression/'+t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN_Olink.tsv', sep='\t', index=None)
