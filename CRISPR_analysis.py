from DataModel import *
import pandas as pd
import re 

field_names = ['cell_type', 'replicant', 'protein', 'type', 'time']

# if zeros more percent drop sample
zero_threshold = 30/100.0
# drop guides if low counts
cutoff = 10

filename = '/lustre/scratch117/casm/team82/np13/CRISPR_170523_human/counts_olly_CRISPR_170523_human_hsLibV1Kosuke19nt.csv'

df = pd.read_csv(filename, index_col = 'libID')

# Remember to remove this when data is fixed
# Workaround for the corrupt data
# df = df.select(lambda col: (not col.endswith('7')), axis=1)

seq_data_df = df[['geneID','guide']]
df = df.drop(['geneID','guide'], axis=1)
dm = DataModel(df, field_names)

for column in list(df):
    col_zeros = (df[column] ==  0).sum()
    if col_zeros / df.shape[0] > zero_threshold:
        df = df.drop(column, axis = 1)
        dm.remove(column)
        print('Dropped column', column, ': Too many zeros,',col_zeros,'out of',df.shape[0])


# Normalize Data
medianRatioNormalization(df, dm)


# Exploratory analysis on time sequence data
avg_results,best_groups, best_values = timeSequenceGroups(df,dm)
exp_col = {}


for column in list(best_groups):
    exp_col[column] = best_groups[column].value_counts().index[0]
    # Inspect Data
    print(best_groups[column].value_counts())

# print(exp_col)

for k, v in avg_results.items():
    print(k,': best time point combination is:',min(v, key = lambda x:x[0]))


control_df = extractBestData(df,exp_col,'OCI2_MT')
treat_df = extractBestData(df,exp_col,'OCI2_WT')

model = empiricalRegression(control_df)

#test = np.log(control_df.T.mean())




