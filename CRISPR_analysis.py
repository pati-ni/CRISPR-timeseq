from DataModel import *
from operations import *
import pandas as pd
import numpy as np
import itertools
import copy

field_names = ['cell_type', 'replicant', 'protein', 'type', 'time']
zero_threshold = 0.3
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
best_groups, best_values = timeSequenceGroups(df,dm)



exp_col = {}
for column in list(best_groups):
    exp_col[column] = best_groups[column].value_counts().index[0]
    # Inspect Data
    print(best_groups[column].value_counts())



    
#time_df = dm.compare_time_sequence(df)

#ts_error = dm.exploreTimesequence(time_df)


#timestamps = calculate_timestamps(df, data_attr)


#Perform reduction in time (only pairs)
# result_set = {}
# for t1,t2 in itertools.product(timestamps, repeat=2):
#     if t2 >= t1:
#         continue
#     tmp_attr = copy.deepcopy(data_attr)
#     for group_name, _ in data_attr['time'].iteritems():
#         if str(t1) != str(group_name) or str(t2) != str(group_name):
#             del tmp_attr['time'][group_name]
#     dataGroupReduction(old_df, reduced_attr)
#     result_set[(t1,t2)] = dataGroupReduction(df, data)
#     time_df,reduceDictionary(df,)



