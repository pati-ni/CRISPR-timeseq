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


# # Exploratory analysis on time sequence data
# avg_results,best_groups, best_values = timeSequenceGroups(df,dm)

# # My results on election
# first_col = {}
# for column in list(best_groups):
#     first_col[column] = best_groups[column].value_counts().index[0]

# # Arthur's results on averaging
# avg_col = {}
# for k, v in avg_results.items():
#     #get the argmin of the min tuple
#     avg_col[k] = min(v, key = lambda x:x[0])[1]


# You gotta choose motherfucker
avg_col = {'OCI2_MT':"['OCI2_A_IDH1_MT_25'] -> ['OCI2_B_IDH1_MT_25'] ", 'OCI2_WT':"['OCI2_A_IDH1_WT_25'] -> ['OCI2_B_IDH1_WT_25']"}
exp_col = avg_col



# LinearRegression for the control group
control_df = extractBestData(df,exp_col,'OCI2_MT')
model = empiricalRegression(control_df)

# Predict adj_var and pval for the treatment
results_df = pd.DataFrame(index = control_df.index)

treat_df = extractBestData(df,exp_col,'OCI2_WT')


results_df['control_mean'] = control_df.T.mean()

model_df = results_df[results_df['control_mean'] > 0].copy()
print('Ignoring zero counters... New size:', model_df.shape)

model_df['control_var'] = control_df[results_df['control_mean'] > 0].T.var()
model_df['treat_mean'] = treat_df[results_df['control_mean'] > 0].T.mean()

model_df['adj_var'] = np.exp(model.predict(np.log(model_df['control_mean']).values.reshape(-1,1))[:,0]) + model_df['control_mean']


pmf = calculatePValues(model_df.head(), 'treat_mean', 'control_mean', 'adj_var')


#test = np.log(control_df.T.mean())
