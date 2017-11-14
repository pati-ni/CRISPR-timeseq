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
medianRatioNormalization(df)



# Comparisons in the same cell line

# Exploratory analysis on time sequence data
min_var, avg_results, best_groups, best_values = timeSequenceGroups(df,dm)

# My results on election
first_col = {}
for column in list(best_groups):
    first_col[column] = best_groups[column].value_counts().index[0]

# Arthur's results on averaging
avg_col = {}
for k, v in avg_results.items():
    #get the argmin of the min tuple
    avg_col[k] = min(v, key = lambda x:x[0])[1]

min_var_col = {}
for k, v in min_var.items():
    min_var_col[k] = min(v, key = lambda x:x[0])[1]

# You gotta choose motherfucker
exp_col = min_var_col

# Use control as t0
for group in dm.getT0samples():
    # filter if there are not enough replica data
    if list(filter(lambda x : len(x) < 2, list(group[-1]))):
        print(group,'not enough replicas')
        continue
    print('Comparing,', group)
    
    
    


control_df = extractBestData(df,exp_col,'OCI2_MT')
treat_df = extractBestData(df,exp_col,'OCI2_WT')

# LinearRegression for the control group
model = empiricalRegression(control_df)

# Predict adj_var and pval for the treatment
# results_df = pd.DataFrame(index = control_df.index)

# results_df['control_mean'] = control_df.T.mean()

# model_df = results_df[results_df['control_mean'] > 0].copy()
# print('Ignoring zero counters... New size:', model_df.shape)

# model_df['control_var'] = control_df[results_df['control_mean'] > 0].T.var()
# model_df['treat_mean'] = treat_df[results_df['control_mean'] > 0].T.mean()

# model_df['adj_var'] = np.exp(model.predict(np.log(model_df['control_mean']).values.reshape(-1,1))[:,0]) + model_df['control_mean']

# #error_lr = varianceErrorRegression(model_df['control_var'], model_df['adj_var'])
# #model_df['corrected_variance'] = error_lr.predict(model_df['control_var'].values.reshape(-1,1))[:,0] + model_df['control_var']


# model_df['pvalue_low'] = calculatePValues(model_df, 'treat_mean', 'control_mean', 'adj_var')
# #Just subtract to make the generate the values
# model_df['pvalue_high'] = 1 - model_df['pvalue_low']


# model_df['phigh_rank'] = model_df['pvalue_high'].rank(ascending = True)
# model_df['plow_rank'] = model_df['pvalue_low'].rank(ascending = True)

# model_df['geneID'] = seq_data_df.loc[model_df.index]['geneID']

# model_df['listID'] = 'dummy_list'


# model_df[model_df['pvalue_low'] < 0.1][['geneID', 'listID', 'pvalue_low']].to_csv('low.gene_example.txt', sep = '\t')

# model_df[model_df['pvalue_high'] < 0.1][['geneID', 'listID', 'pvalue_high']].to_csv('high.gene_example.txt', sep = '\t')


# # model_df[['geneID','dummy','pvalue_high']].to_csv('high.gene_example.txt', header = '\t')


# #low_df = model_df['geneID'].isin(model_df[model_df['pvalue_low'] < 0.05].groupby('geneID').count().index)
# #high_df = model_df['geneID'].isin(model_df[model_df['pvalue_high'] < 0.05].groupby('geneID').count().index)






# #test = np.log(control_df.T.mean())
