import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

#import results from protein-gene expression correlations/ multiple linear regression
corr=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/total_STARNET_pro-exp_correlations.tsv', delimiter='\t')
lm=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/STARNET_LM_Olink_GE_tissues.tsv', delimiter='\t')

#obtain -log p-values
def get_logpv(df, pvlab):
    pvs=df[pvlab].tolist()
    l=[]
    for pv in pvs:
        if pv == 0:
            l.append(1e-300)
        else:
            lpv=-math.log10(pv)
            l.append(lpv)
    df[pvlab+'_log_pv']=l
    
    return(df)

corr_log=corr[['Gene stable ID', 'Gene_name', 'tissue', 'Spearmanlogpv', 'Spearmanqv']]
lm_log=get_logpv(lm, 'Pr(>|t|)')[['Gene.stable.ID', 'tissue', 'Pr(>|t|)_log_pv', 'qvalues']]

#merge to single dataframe and export
corr_lm=pd.merge(corr_log, lm_log, 
	left_on=['Gene stable ID', 'tissue'], 
	right_on=['Gene.stable.ID', 'tissue'], 
	how='inner')

corr_lm.to_csv('/export/michoelfs/users/sean/INTRePID/STARNET/gene_expression/pro_exp_match/STARNET_LM_Spearman_Olink_GE_tissues.tsv', sep='\t', index=None)

#get -log10 threshold that corresponds to 5% FDR
corr_thd=corr_log.loc[corr_log['Spearmanqv'] <= 0.05].sort_values('Spearmanlogpv').head(1)['Spearmanlogpv'].values[0]
lm_thd=lm_log.loc[lm_log['qvalues'] <= 0.05].sort_values('Pr(>|t|)_log_pv').head(1)['Pr(>|t|)_log_pv'].values[0]

#set hue order
hue_order=['LIV', 'SKLM', 'AOR', 'MAM', 'VAF', 'SF', 'Blood']
hue_palette = {'LIV': 'brown', 'SKLM': 'green', 'AOR': 'red', 'MAM': 'orange', 'VAF': 'purple', 'SF': 'pink', 'Blood': 'blue'}

#plot scatterplot of -log10 p-values for correlation vs LM results
plt.figure(figsize=(8,6))
sns.scatterplot(data=corr_lm, x='Spearmanlogpv', y='Pr(>|t|)_log_pv', hue='tissue', hue_order=hue_order, palette=hue_palette)
plt.xlabel('Spearman correlations (-log10 p-value)', fontsize=12)
plt.ylabel('multiple linear regression (-log10 p-value)', fontsize=12)
plt.legend(title='Tissue')
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12)
plt.axvline(x = corr_thd, color='r', linestyle='dashed')
plt.axhline(y = lm_thd, color='r', linestyle='dashed')
texts=[]
for (xi, yi, si) in zip(corr_lm['Spearmanlogpv'], corr_lm['Pr(>|t|)_log_pv'], corr_lm['Gene_name']):
    if xi > 11 and yi > 5:
        texts.append(plt.annotate(
            si,
            (xi, yi),
            textcoords="offset points",
            xytext=(10, 10),
            ha='center',
            arrowprops=dict(arrowstyle="->", color='gray', lw=1)
        ))
adjust_text(texts)
plt.savefig('/export/michoelfs/users/sean/INTRePID/STARNET/gene_expression/pro_exp_match/STARNET_LM_vs_Spearman_Olink_GE_tissues.png')