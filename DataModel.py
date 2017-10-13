import re
import pandas as pd
import numpy as np
from operator import itemgetter
from itertools import compress
from scipy.optimize import leastsq
from numpy.linalg import lstsq


def meanSquareError(points, x = None):
    if x is None:
        x = np.arange(len(points))
        # Distance from plane
    A = np.vstack([x,np.ones(len(x))]).T
    y = np.array(points)
    e = lstsq(A,y)[1][0]
    return e

def dataGroupReduction(df,id1,id2,threshold = 1):    
    df1 = df[df[id1] < threshold][id1]
    df2 = df[df[id2] < threshold][id2]
        
    return df1.index.intersection(df2.index)

def dataGroupUnion():
    return



class DataModel:
    data_model = {}

    def __init__(self, df, fields, delim = '_'):
        df_columns = list(df)
        self.fields = fields[:]
        self.initial_fields = fields[:]
        self.delim = '_'

        # Initialize with empty dict for attributes our model
        for attribute in self.fields:
            self.data_model[attribute] = {}
            # Generate # TODO: he indexes
        for column_name in df_columns:
            col_meta = re.split(delim, column_name)
            if len(col_meta) != len(self.fields):
                print('Attribute', col_meta[0], 'skipped: Not containing required metadata')
                continue

            for index, (attribute, meta_val) in enumerate(zip(self.fields, col_meta)):
                if not (meta_val in self.data_model[attribute]):
                    self.data_model[attribute][meta_val] = {}
                self.data_model[attribute][meta_val][delim.join(col_meta[:index] + col_meta[index+1:])] = column_name



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


    def reduceDataModel(self, data_name, delim = '_'):   
        try:
            reduced_attr = self.data_model[data_name]
            data_index = self.fields.index(data_name)
        except KeyError:
            print('Data field',data_type,'not found')

        #Build a new data model
        new_dm ={}
        reduced_mask = [1] * len(dict_name_list)
        reduced_mask[data_index] = 0
    

        # new_df = pd.DataFrame(index = df.index)
        # attrs = re.split(delim, id1)
        # new_name = '_'.join(compress(attrs, mask))
        # new_df[new_name] = 0
        # new_df[new_name] = df[dataGroupReduction(df) ==]

        for field, groups in self.data_model.items():

            #Do not add the reduced field to the new data model
            if field  == data_name:
                continue
            try:
                current_field_index = self.fields.index(field)
            except ValueError:
                print('Data field',data_type,'not available in the list')


            # Generate bit mask to remove irrelevant fields for the alias
            alias_mask = reduced_mask[:]
            alias_mask[current_field_index] = 0

            new_dm[field] = {}
            for group, alias_dict in groups.items():
                new_group = {}
                for alias, ref_name in alias_dict.items():
                    old_col_name = re.split(delim, ref_name)
                    col_name = delim.join(compress(old_col_name, reduced_mask))
                    alias_name = delim.join(compress(old_col_name, alias_mask))
                    new_group[alias_namecol] = _name
                    new_dm[field][group] = new_group

        # Modify the current data model
        self.fields.remove(data_name)
        return new_dm

    # Generator that iterates groups
    def iterateGroup(self, attribute, delim = '_'):
        try:
            group = self.data_model[attribute]
            data_index = self.fields.index(attribute)
        except KeyError:
            print('Key error:',attribute, ' not available in data model')

        reduced_mask = [1] * len(self.fields)
        reduced_mask[data_index] = 0

        if attribute == 'time':
            indexes = sorted(map(int, group.keys()))
        else:
            indexes = group.keys()
        base_keys = group[str(indexes[0])].keys()
        for key in base_keys:
            name_t = re.split(delim, key)
            group_name = delim.join(compress(name_t, reduced_mask))
            yield (group_name,indexes,[group[str(index)][key] for index in indexes])
    def exploreTimesequence(self, df):
        # agg_func = np.vectorize(meanSquareError)
        res_df = pd.DataFrame(index=df.index)
        for agg_name, time_stamps, group in self.iterateGroup('time'):
            valid_df = df[group].replace([-np.inf,np.inf],np.nan).dropna()
            print(group, 'Valid records:', valid_df.shape,'Meta:', agg_name, time_stamps)
            # valid_df.agg( agg_func, axis=1)
            res_df[agg_name] = valid_df.aggregate(meanSquareError, axis = 1, x = time_stamps)
            break
        return res_df

