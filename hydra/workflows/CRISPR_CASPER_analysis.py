from DataModel import *
import pandas as pd
import re 

field_names = ['cell_line', 'type' , 'replica']

# if zeros more percent drop sample
zero_threshold = 30/100.0
# drop guides if low counts
cutoff = 10

control_tag = 'NOCAS'

filename = '/lustre/scratch117/casm/team82/np13/CRISPR_CASPER_171030/counts_manos_CRISPR_CASPER_171030_hsCASPERlibAll170329.csv'

df = pd.read_csv(filename, index_col = 'libID')

seq_data_df = df[['geneID','guide']]
df.drop(['geneID','guide'], axis=1, inplace = True)



df.drop(['Auto_library','Cas_Library','D_Library'], inplace = True, axis = 1)

print('Add NOCAS label')
for sample in df.columns:
    attrs = re.split('_',sample)
    if len(attrs) < len(field_names):
      #No Type attribute so add type 
      attrs.insert(1,'NOCAS')
      df.rename(columns = { sample : '_'.join(attrs) }, inplace = True)



# Fill missing lines
df['A35_NOCAS_A'] = df['644_NOCAS_A']
df['A35_NOCAS_B'] = df['644_NOCAS_B'] 
df['A35_NOCAS_C'] = df['644_NOCAS_C'] 


medianRatioNormalization(df)

dm = DataModel(df, field_names)

pvalue_cutoff = 0.1

for key, index in dm.iterateData('cell_line'):
    for group, samples in dm.splitNestedData(index.values(),['type']):
        if control_tag == group:
            control_group = df[list(samples)]
        else:
            treat_group = df[list(samples)]
        
    print('Control:', list(control_group))
    print('Treat:', list(treat_group))
    
    res_df = analyzeData(control_group, treat_group)
    res_df['listID'] = 'dummy_list'
    res_df['geneID'] = seq_data_df['geneID']
    

    print('Cell line:', key, 'Pvalues', res_df[res_df['pvalue_low'] < 0.05].shape, 'for cutoff', cutoff)
    res_df[res_df['pvalue_low'] < pvalue_cutoff][['geneID', 'listID', 'pvalue_low']].to_csv(key + '_low.gene_example.txt', sep = '\t')
    

            
    



