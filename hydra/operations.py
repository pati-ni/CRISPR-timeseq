import numba
import numpy as np
from scipy.stats import nbinom

@numba.jit
def x_hat(counter_values):
    result = np.empty(len(counter_values))
    N = len(counter_values[0])
    for row in range(len(counter_values)):
        result[row] = np.exp(np.sum(np.log(counter_values[row] + 1)) / N)
    return result



@numba.jit
def _calculate_pvalues(x, r, p, mean, length):
    results = np.empty(length)
    for i in range(length):
        results[i] =  nbinom.cdf(x[i], r[i], p[i])
    return results



def medianRatioNormalization(df, samples = None):
    if samples is None:
        samples = list(df)
        # work directly on dataframe
        mod_df = df
    else:
        mod_df = df[samples].copy()

    # Insert temporarily with the key
    geom_mean = x_hat(mod_df[samples].values)
    for sample in samples:
        mod_df[sample] = mod_df[sample] / (mod_df[sample] / geom_mean).median()

    if not (samples is None):
        for col in list(mod_df):
            df[col] = mod_df[col]
