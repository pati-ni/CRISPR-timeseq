from sklearn.neighbors import NearestNeighbors
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from subprocess import call
from hydra.operations import _calculate_pvalues

import os
import pathlib

# Contains all possible analysis methods for a set of samples



class Analysis:

    def __init__(self,  path, dropout = True, enrichment = True):

        self.path = path.rstrip('/')
        self.path += '/' + self.name + '/'
        self.rra_path = os.path.expanduser("~np13/src_repos/rra/RRA")
        pathlib.Path(self.path).mkdir(parents=True, exist_ok=True)
        self.p_tag_low = 'pvalue_low'
        self.p_tag_high = 'pvalue_high'
        self.enrichment = enrichment
        self.dropout = dropout
    
    def exportGuideRNAData(self, df):
        if self.dropout:
            pass
        if self.enrichment:
            pass

    def rra(self, df, gene_df, group_id):
        dfs = []
        if self.dropout:
            dfs.append(self.rra_analysis(df, gene_df, group_id, self.p_tag_low, dropout = True))
        if self.enrichment:
            dfs.append(self.rra_analysis(df, gene_df, group_id, self.p_tag_high, dropout = False))
        return dfs
        

    def rra_analysis(self, df, gene_df, group_id, p_tag, dropout = True, keep_tmp_files = True ):
        kw = 'low'
        kw2 = 'neg'

        if not dropout:
            kw = 'high'
            kw2 = 'pos'

        input_file = self.path + group_id + '-' + kw + '.rra_input.txt'
        output_file = self.path + group_id + '-'+ kw + '.rra_output.txt'

        df['geneID'] = gene_df.loc[df.index]['geneID']
        df['listID'] = 'dummy_list'
        rra_df = df[['geneID', 'listID', p_tag]]
        rra_df.to_csv(input_file, sep = '\t')
        call([self.rra_path,'-i', input_file, '-o', output_file])
        
        results_df = pd.read_csv(output_file, index_col = 'group_id', sep='\t')

        if not keep_tmp_files:
            os.remove(output_file)
            os.remove(input_file)

        results_df.rename(columns = {'items_in_group' : 'num'}, inplace = True)
        results_df.rename(columns = {'p' : kw2 + '|pvalue'}, inplace = True)
        results_df.rename(columns = {'lo_value' : kw2 + '|score'}, inplace = True)
        results_df.rename(columns = {'FDR' : kw2 + '|fdr'}, inplace = True)
        results_df.to_csv(self.path + group_id + '.gene_summary.txt')
        return results_df


class DESeq2(Analysis):

    def __init__(self, path):
        # Load rpy extension
        get_ipython().run_line_magic('load_ext', 'rpy2.ipython')

        # Import snippet
        get_ipython().run_cell_magic('R', '', 
                """
                library('DESeq2')\n
                rpy2_deseq2 <- function(df,sample_types,ctrl_treat) {\n    
                    sampleCondition <- unlist(sample_type)\n    
                    df <- data.frame(experiment_df)\n    
                    ddsHTSeq <- DESeqDataSetFromMatrix(df, DataFrame(sampleCondition), ~ sampleCondition)\n    
                    colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$sampleCondition, levels=ctrl_treat)\n    
                    dds<-DESeq(ddsHTSeq)\n    
                    res<-results(dds)\n    
                    res<-res[order(res$padj),]\n
                    return(data.frame(res))\n
                }\n
                """)

        self.name = 'deseq'
        self.requires_normalization = False
        super().__init__(path)


    def analyzeData(self, control_df, treat_df):
        experiment_df = pd.concat([control_df,treat_df], axis = 1)
        sample_type = [len(list(control_df)) * ["control"] ] + [len(list(treat_df)) * ["treat"] ]
        types = ["control", "treat"]

        get_ipython().run_line_magic('R', '-i experiment_df,sample_type,types')
        get_ipython().run_line_magic('R', 'deseq_df <- rpy2_deseq2(experiment_df, sample_type, types)')
        get_ipython().run_line_magic('R', '-o deseq_df')
        ret_df = get_ipython().user_ns['deseq_df']
        ret_df.rename(columns={'padj': 'pvalue_low' },inplace = True)
        ret_df.index.name = 'libID'
        return ret_df


class MeanVarianceEmpirical(Analysis):

    def __init__(self, *args):
        
        # Check if overriden by subclass
        if not hasattr(self, 'name'):
            self.name = 'meanvar'

        self.requires_normalization = True
        super().__init__(*args)

    def meanVarianceRegression(self, df):

        print('Regression analysis at:', list(df))
        regression_df = pd.DataFrame(index = df.index)

        regression_df['var_dif'] = np.log(df.T.var())
        regression_df['mean'] = np.log(control_df.T.mean())

        lr = LinearRegression(n_jobs = -1)
        lr.fit(regression_df['mean'].values.reshape(-1,1), regression_df['var_dif'].values.reshape(-1,1))
        error = lr.score(regression_df['mean'].values.reshape(-1,1), regression_df['var_dif'].values.reshape(-1,1))
        print('Mean Variance regression converged, regression score error:',error)

        return lr
    

    def empiricalRegression(self, control_df, cutoff = 10):

        print('Regression analysis at the control:',list(control_df))

        regression_df = pd.DataFrame(index = control_df.index)
        regression_df['mean'] = control_df.T.mean()

        #Keep only overdispersed data, model limitation (Mageck)
        print('Keeping only overdispersed data...', (regression_df['mean'] < control_df.T.var()).sum())
        regression_df = regression_df.loc[regression_df['mean'] < control_df.T.var()]

        regression_df['var_dif'] = np.log(control_df.T.var() - regression_df['mean'])
        regression_df['mean'] = np.log(regression_df['mean'])

        lr = LinearRegression(n_jobs = -1)
        lr.fit(regression_df['mean'].values.reshape(-1,1), regression_df['var_dif'].values.reshape(-1,1))
        error = lr.score(regression_df['mean'].values.reshape(-1,1), regression_df['var_dif'].values.reshape(-1,1))
        print('Regression run, regression score error:',error)
        return lr

    def analyzeData(self, control_df, treat_df):
        # LinearRegression for the control group
        model = self.empiricalRegression(control_df)
        
        results_df = pd.DataFrame(index = control_df.index)
        # Deal with zero counters, yeah right
        control_df += 1
        treat_df += 1
        results_df['control_mean'] = control_df.T.mean()
        results_df['treat_mean'] = treat_df.T.mean()
        x = results_df['control_mean'].values.reshape(-1,1)
        results_df['adj_var'] = np.exp(model.predict(np.log(x))[:, 0]) + results_df['control_mean']
        results_df['pvalue_low'] = self.calculatePValues(results_df, 'treat_mean', 'control_mean', 'adj_var')
        return results_df

    def calculatePValues(self, df, test_mean_label, model_mean_label, adj_var_label, low_tail = True):
        length = df[test_mean_label].shape[0]
        if not  length == df[model_mean_label].shape[0] == df[adj_var_label].shape[0]:
            raise ValueError ('Vector length do not match')
        x = df[test_mean_label].values
        p = (df[model_mean_label] / df[adj_var_label]).values
        r = (df[model_mean_label] ** 2 / (df[adj_var_label] - df[model_mean_label])).values

        return _calculate_pvalues(x, r, p, df[model_mean_label].values, length)



class MeanVarianceKNN(MeanVarianceEmpirical):

    def __init__(self, *args):
        self.name = 'meanvarknn'
        self.requires_normalization = True
        super().__init__(*args)
        

    def analyzeData(self, control_df, treat_df):

        if control_df.shape[0] != treat_df.shape[0]:
            raise ValueError ('Dataframes have different index values')

        print(list(control_df), list(treat_df))
        meta_df = pd.DataFrame(index = control_df.index)
        meta_df['control_mean'] = control_df.T.mean()
        meta_df['control_var'] = control_df.T.var()
        meta_df['treat_mean'] = treat_df.T.mean()
        meta_df['treat_var'] = treat_df.T.var()

        samples_dimension = control_df.shape[1]

        test = control_df.values.reshape(-1, samples_dimension)

        # Deal with underdispersed data
        overdispersed_data = control_df[meta_df['control_mean'] < meta_df['control_var']].copy()
        overdispersed_df = pd.DataFrame(index = overdispersed_data.index)
        overdispersed_df['control_mean'] = overdispersed_data.T.mean()
        overdispersed_df['control_var'] = overdispersed_data.T.var()

        train = overdispersed_data.values.reshape(-1, samples_dimension)
        print('Found', len(control_df) - len(overdispersed_data), 'underdispersed samples')
        print('Train shape', train.shape,'Test shape',test.shape)

        nbrs = NearestNeighbors(n_neighbors = 10, algorithm = 'ball_tree').fit(train)

        distances, neighbors = nbrs.kneighbors(test)

        adj_var = np.empty(shape=(treat_df.shape[0],))

        # TODO: optimize
        for i, indices in enumerate(neighbors):
            #print(i, indices[i])
            lr = LinearRegression(n_jobs = -1)
            X = overdispersed_df.iloc[indices]['control_mean'].values.reshape(-1,1)
            y = overdispersed_df.iloc[indices]['control_var'].values.reshape(-1,1)
            lr.fit(np.log(X), np.log(y - X))
            adj_var[i] = np.exp(lr.predict(np.log(meta_df.iloc[i]['control_mean']))) + meta_df.iloc[i]['control_mean']
            if adj_var[i]  < meta_df.iloc[i]['control_mean']:
                # This should not happen
                print(list(zip(indices, distances[i])))

        meta_df['adj_var'] = adj_var
        meta_df['pvalue_low'] = self.calculatePValues(meta_df, 'treat_mean', 'control_mean', 'adj_var')
        return meta_df



