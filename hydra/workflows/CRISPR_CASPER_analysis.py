from hydra.workflow import Workflow
from hydra.analysis import *

import pandas as pd
import re 
import ntpath
import os

field_names = ['cell_line', 'type' , 'replica']

control_tag = 'NOCAS'

filename = '/lustre/scratch117/casm/team82/np13/CRISPR_CASPER_171030/counts_manos_CRISPR_CASPER_171030_hsCASPERlibAll170329.csv'
project_name = 'CRISPR_CASPER_171030'
results_root = os.path.dirname(filename)  + '/' + project_name


df = pd.read_csv(filename, index_col = 'libID')

seq_data_df = df[['geneID','guide']]
df.drop(['geneID','guide'], axis=1, inplace = True)


# Experiment specific please remove
df.drop(['Auto_library','Cas_Library','D_Library'], inplace = True, axis = 1)
print('Add NOCAS label')
for sample in df.columns:
    attrs = re.split('_',sample)
    if len(attrs) < len(field_names):
        #No Type attribute so add type 
        attrs.insert(1,'NOCAS')
        df.rename(columns = { sample : '_'.join(attrs) }, inplace = True)

# Fill missing ctrl for A35
df['A35_NOCAS_A'] = df['644_NOCAS_A']
df['A35_NOCAS_B'] = df['644_NOCAS_B'] 
df['A35_NOCAS_C'] = df['644_NOCAS_C'] 



# Normal pipeline continues

dm = Workflow(df, field_names,control_tag)

#dm.analyzeExperiment(DESeq2(results_root), df, seq_data_df)

dm.analyzeExperiment(MeanVarianceKNN(results_root), df, seq_data_df)
#dm.analyzeExperiment(MeanVarianceEmpirical(results_root), df, seq_data_df)
