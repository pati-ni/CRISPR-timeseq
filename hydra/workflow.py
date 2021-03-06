import re
import pandas as pd
import numpy as np
from operator import itemgetter
from hydra.operations import medianRatioNormalization
from itertools import compress, combinations, chain
from ast import literal_eval

# TODO: use logging module
"""
workflow.py 

Different types of experiment require different metadata data model structures.
For example, when we are dealing with timepoint sample data we are dealing with another dimension of analysis,
which somehow needs to be reduced.

This file contains different workflows for different cases
"""

# Workflow with multiple replicas separated in control and treatment
class Workflow:
    data_model = {}
    exp_results = {}
    def __init__(self, df, fields, control_tag, delim = '_'):
        df_columns = list(df)
        self.fields = fields[:]
        self.initial_fields = fields[:]
        self.delim = delim
        self.control_tag = control_tag

        # Initialize with empty dict for attributes our model
        for attribute in self.fields:
            self.data_model[attribute] = {}
            # Generate the indexes
        for column_name in df_columns:
            col_meta = re.split(delim, column_name)
            if len(col_meta) != len(self.fields):
                print('Attribute', col_meta[0], 'skipped: Not containing required metadata')
                continue

            for index, (attribute, meta_val) in enumerate(zip(self.fields, col_meta)):
                if not (meta_val in self.data_model[attribute]):
                    self.data_model[attribute][meta_val] = {}
                self.data_model[attribute][meta_val][delim.join(col_meta[:index] + col_meta[index + 1:])] = column_name

    def analyzeExperiment(self, method, df, gene_df, cell_group = 'cell_line', ctrl_treat_id = 'type'):
        
        if method.requires_normalization:
                self.normalize(df)

        for key, index in self.iterateData(cell_group):
            # Control treat id
            for group, samples in self.splitNestedData(index.values(), [ctrl_treat_id]):
                if self.control_tag == group:
                    control_group = df[list(samples)]
                else:
                    treat_group = df[list(samples)]

            # Copy dataframes to keep original counters in tact
            self.coreAnalysis(method, control_group.copy() + 1, treat_group.copy() + 1, gene_df, key)

                
    def coreAnalysis(self, method, control_df, treat_df, gene_df, key):
        print('Control:', list(control_df))
        print('Treat:', list(treat_df))

        results_df = method.analyzeData(control_df, treat_df)
        results_df['pvalue_high'] = 1 - results_df['pvalue_low']

        method.exportGuideRNAData(results_df)
        gene_info = method.rra(results_df, gene_df, key)
        # Meta analysis
        # if not dfs[key]:
            # dfs[key] = {}
        # self.exp_results[key][method.name] = {'gRNA': results_df , 'gene': gene_info}

    def iterateData(self, attribute):
        return self.data_model[attribute].items()
    
    
    def getKeys(self, group_keys, attribute):
        return list(group_keys)

    # Uses the data model and a set of desired fields
    # to produce groups of samples
    def splitNestedData(self, samples, fields, original_fields = None, null_yield = True):
        if not samples:
            if null_yield:
                yield '', samples
            raise StopIteration
            #raise ValueError('No samples available')

        if original_fields is None:
            original_fields = list(fields)

        if not fields:
            # Hooray! We reached the final level of nesting
            # Yield results and stop iteration
            yield self.maskID(next(iter(samples)),original_fields), samples
            raise StopIteration

        # Iterate over possible groups
        for group_id, group in self.iterateData(fields[0]):
            group_members = set(group.values()).intersection(samples)
            # Not in the level we want so evaluate recursive calls
            for split_group in self.splitNestedData(group_members, fields[1:], original_fields = original_fields):
                yield split_group

    # Removes specific field from sample_id
    def maskID(self, sample_id, active_fields, discard = False, fields = None):
        # Disable all fields in starting mask
        if fields is None:
            fields = self.fields

        pre_value = discard
        mask = [pre_value]*len(fields)
        for field in active_fields:
            try:
                field_index = fields.index(field)
                # Enable field in mask
                mask[field_index] = not pre_value
            except ValueError:
                print('Field',field ,'not found')

        split_attrs = re.split(self.delim,sample_id)
        return self.delim.join(compress(split_attrs,mask))
    
    # Removes sample id from data model
    def remove(self, rem_id):
        for field, groups in self.data_model.items():
            for key, group in groups.items():
                self.data_model[field][key] = { k:v for k, v in group.items() if v != rem_id }

    # Generator that iterates groups
    def iterate(self, attribute):
        try:
            group = self.data_model[attribute]
            data_index = self.fields.index(attribute)
        except KeyError:
            print('Key error:',attribute, ' not available in data model')
            

        reduced_mask = [1] * len(self.fields)
        reduced_mask[data_index] = 0
        
        indexes = getKeys(group.keys(),attribute)
        
        base_keys = group[indexes[0]].keys()
        for key in base_keys:
            name_t = re.split(self.delim, key)
            group_name = self.delim.join(compress(name_t, reduced_mask))
            yield (group_name, indexes, [group[index][key] for index in indexes])

    # TODO: Not used anywhere, maybe cleanup?
    def reduceDataModel(self, data_name):   
        try:
            reduced_attr = self.data_model[data_name]
            data_index = self.fields.index(data_name)
            
        except KeyError:
            print('Data field',data_type,'not found')

        #Build a new data model
        new_dm = {}
        reduced_mask = [1] * len(self.fields)
        reduced_mask[data_index] = 0

        for field, groups in self.data_model.items():

            #Do not add the reduced field to the new data model
            if field  == data_name:
                continue
            try:
                current_field_index = self.fields.index(field)
            except ValueError:
                print('Data field',data_name,'not available in the list')


            # Generate bit mask to remove irrelevant fields for the alias
            alias_mask = reduced_mask[:]
            alias_mask[current_field_index] = 0

            new_dm[field] = {}
            for group, alias_dict in groups.items():
                new_group = {}
                for alias, ref_name in alias_dict.items():
                    old_col_name = re.split(self.delim, ref_name)
                    col_name = self.delim.join(compress(old_col_name, reduced_mask))
                    alias_name = self.delim.join(compress(old_col_name, alias_mask))
                    new_group[alias_name] = col_name
                    new_dm[field][group] = new_group

        # Modify the current data model
        self.fields.remove(data_name)
        return new_dm

    def normalize(self, df, *args):
        for key, index in self.iterateData('cell_line'):
            samples = list(index.values())
            print('Median Ratio Normalization', samples)
            medianRatioNormalization(df, samples)
    

# Use this to generate pairs
def all_combinations(any_list):
    return chain.from_iterable(
        combinations(any_list, i)
        for i in range(1, len(any_list) + 1))

def extractBestData(df, exp_col, key):
    exp_df = pd.DataFrame(index = df.index)
    print('Extracting mean for:', exp_col[key])
    for index, samples in enumerate(re.split('->', exp_col[key])):
        # Will it blow ?? 3 , 2 , 1 ...
        exp_df['_'.join([key, str(index)])] = df[literal_eval(samples)].T.mean()
        return exp_df


class TimeSequenceWorkflow(Workflow):

    def __init__(self, df, fields, control_tag, *args):
        
        if not ('time' in fields):
            raise('This does not look like time sequence meta data')
        super().__init__(df, fields, control_tag)
    
    # generate pairs of timesamples on day 0
    def getT0samples(self):
        t0 = str(sorted(map(int, self.data_model['time'].keys()))[0])
        print('Found minimum timepoint,', t0)

        unique_set0 = set(self.data_model['time'][t0].values())

        for key in self.data_model['time'].keys():
            if t0 == key:
                continue
            ts_set = unique_set0 | set(self.data_model['time'][key].values())
            for group_id, unique_set in self.splitNestedData(ts_set, ['cell_type','type']):
                a = self.splitNestedData(unique_set, ['time'])
                yield group_id + t0 + '-' + key , list(filter(None, map(itemgetter(1), a )))

    def getKeys(self, group_keys, attribute):
        if attribute == 'time':
            return map(str, sorted(map(int, group_keys)))
        else:
            return Workflow.getKeys(group_keys, attribute)

    # Returns a argmin histogram of the samples combinations with minimized mean deviation
    # Each count represent a guide RNA
    def timeSequenceGroups(self, df, fields = ['cell_type', 'type']):
        print('Rating Samples')
        test_avg = {}
        min_var = {}
        # print('Array shape:', df.shape)
        df = df.replace([np.inf, -np.inf], np.nan).dropna()
        #print('Finding best sequence based on:',df.shape,'values')
        compare_df = pd.DataFrame(index = df.index)
        arg_compare_df = pd.DataFrame(index = df.index)

        for time_set in all_combinations(dict(self.iterateData('time')).values()):
            time_samples = list(chain.from_iterable(ts.values() for ts in time_set))
            for group_id, unique_set in self.splitNestedData(time_samples, fields):
                # If we do not have more than 2 replicants, skip the sample
                if len(unique_set) < 2:
                    print('warn: Skipping', unique_set,', no replicant available')
                    #print(list(dm.splitNestedData(unique_set,['replicant'])))
                    continue
                mean_df = pd.DataFrame(index = df.index)

                # the id of argmin
                gid = []
                for i, (_, replicant) in enumerate(self.splitNestedData(unique_set, ['replicant'])):
                    # if replicant not available, continue
                    if _ == '':
                        continue

                    # Compute the mean accross replicant's timepoints
                    mean_df[str(i) + _] = df[list(replicant)].T.mean()

                    gid.append(str(sorted(list(replicant))))

                # Eye candy for the argmin
                gid = '->'.join(gid)

                # mad: mean absolute deviation
                mad = mean_df.T.mad()

                # Arthur's case, Much simpler huh?
                if not group_id in test_avg:
                    test_avg[group_id] = []
                    min_var[group_id] = []

                test_avg[group_id].append((mad.mean(),gid))

                # Custom metric to minimize impact for the mean variance modeling
                variance = np.var(np.log(mean_df.T.mean()) - np.log(mean_df.T.var()))
                min_var[group_id].append((variance, gid))

                if group_id in compare_df.columns:
                    # key already present compare existing data
                    better_mean = compare_df[group_id] > mad
                    #print(unique_set,':',better_mean.sum())
                    arg_compare_df[group_id].loc[better_mean] = gid # alias for str(unique_set)
                    compare_df[group_id].loc[better_mean] = mad
                else:
                    # previous group not present
                    # set min and argmin
                    compare_df[group_id] = mad
                    arg_compare_df[group_id] = gid

        return min_var, test_avg, arg_compare_df, compare_df


    def calculate_timestamps(self, df):
        timestamps = []
        for ts in self.data_model['time'].keys():
            timestamps.append((int(ts), ts))
        return sorted(timestamps, key = lambda x : x[0])

    def compare_time_sequence(self, df):
        pop_df = pd.DataFrame(index = df.index)
        timestamps = self.calculate_timestamps(df)
        time_groups = self.data_model['time']
        min_ts, min_ts_key = timestamps[0]

        for alias, base_index in time_groups[min_ts_key].items():
            #filter to avoid duplication
            base_df = df[df[base_index] > 0][base_index]
            print(base_df.shape, df[base_index].shape)
            for moment, col_map in time_groups.items():
                print('Comparing', base_index, col_map[alias])
                # this line may be a source of bugs if the data model order changes
                moments_key = '_'.join([alias, moment])
                pop_df[moments_key] = 1 - (base_df - df[col_map[alias]]) / base_df
                pop_df[moments_key] = - pop_df[moments_key].apply(np.log)
        return pop_df

    # exploratory time analysis
    def exploreTimesequence(self, df):
        res_df = pd.DataFrame(index=df.index)
        for agg_name, time_stamps, group in self.iterateGroup('time'):
            valid_df = df[group].replace([-np.inf,np.inf],np.nan).dropna()
            print(group, 'Valid records:', valid_df.shape,'Meta:', agg_name, time_stamps)
            # valid_df.agg( agg_func, axis=1)
            res_df[agg_name] = valid_df.aggregate(meanSquareError, axis = 1, x = time_stamps)
            break
        return res_df
    
