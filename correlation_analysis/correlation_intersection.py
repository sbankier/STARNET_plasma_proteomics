import pandas as pd

#import protein gene correlations
corr=pd.read_csv('/Home/ii/seanb/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/total_STARNET_pro-exp_correlations.tsv', delimiter='\t')
lm=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/STARNET_LM_Olink_GE_tissues.tsv', delimiter='\t')

#import protein-eigengene correlations
ei=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_labelled.tsv', delimiter='\t')

#filter Spearman correlations by q-value
corr_filt=corr.loc[corr['Spearmanqv'] <= 0.05]
lm_filt=lm.loc[lm['qvalues'] <= 0.05]

#filter for significant cross-tissue networks
spec_filt=ei.loc[(ei['spearman_qv'] <= 0.05) & (ei['network_type'] == 'specific')]
cross_filt=ei.loc[(ei['spearman_qv'] <= 0.05) & (ei['network_type'] == 'cross')]

#filter to count tissues per network
cross_net=cross_filt[['network_id', 'tissue']].drop_duplicates()
cross_net['tissue'] = cross_net['tissue'].str.split(', ')
cross_net=cross_net.explode('tissue').reset_index(drop=True)

#import GRN modules
mod=pd.read_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/modules.csv', delimiter=',')

#filter cross tissue modules for those containing protein in GRN
cross_t=cross_filt.loc[cross_filt['present_in_GRN'] == True]

#get tissue for cross tissue linked proteins
cross_t_mer=cross_t[['Gene_name', 'network_id']]
mod_mer=mod[['gene_symbol', 'clust', 'tissue']]
cross_t_mod=pd.merge(cross_t_mer, mod_mer, left_on=['Gene_name', 'network_id'], right_on=['gene_symbol', 'clust'], how='inner').drop(['clust', 'gene_symbol'], axis=1)

#count overlapping tissues for genes
def get_count(df, clab, tlab):

    #count gene frequency in dataframe
    freq=df[clab].value_counts().to_dict()
    freqdf=pd.DataFrame(freq.items(), columns=[clab, 'frequency'])

    #create dictionary for tissues
    tissues=df[[tlab]].drop_duplicates()[tlab].tolist()
    tissue_lists = {tissue: [] for tissue in tissues}

    #loop through the dataframe and populate the tissue lists
    for index, row in df.iterrows():
        for tissue in tissues:
            if row[tlab] == tissue:
                tissue_lists[tissue].append(1)
            else:
                tissue_lists[tissue].append(0)

    #assemble tissue counts as a single dataframe
    count_df = pd.DataFrame({clab: df[clab], **tissue_lists})
    sum_df = count_df.groupby(clab, as_index=False).max()
    inc_df=pd.merge(freqdf, sum_df, on=clab, how='inner')

    return inc_df

corr_count=get_count(corr_filt, 'Gene_name', 'tissue')
lm_count=get_count(lm_filt, 'Gene_name', 'tissue')
spec_count=get_count(spec_filt, 'Gene_name', 'tissue')
cross_count=get_count(cross_net, 'network_id', 'tissue')
cross_t_count=get_count(cross_t_mod, 'Gene_name', 'tissue')

#export tissue intersections
corr_count.to_csv('/Home/ii/seanb/sean/INTRePID/STARNET/gene_expression/pro_exp_match/correlations/total_STARNET_pro_exp-match_correlations_intersection_q0.05.tsv', sep='\t', index=None)
lm_count.to_csv('~/sean/INTRePID/STARNET/gene_expression/pro_exp_match/linear_reg/STARNET_LM_Olink_GE_tissues_intersection_5FDR.tsv', sep='\t', index=None)
spec_count.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_tissue_specific_intersection_5FDR.tsv', sep='\t', index=None)
cross_count.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_cross_tissue_intersection_networks_5FDR.tsv', sep='\t', index=None)
cross_t_count.to_csv('~/sean/INTRePID/STARNET/gene_expression/eigengenes/Olink_cases_combined_adjusted_eigengene_correlations_cross_tissue_intersection_5FDR_GRN_present.tsv', sep='\t', index=None)