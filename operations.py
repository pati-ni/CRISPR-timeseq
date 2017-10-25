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

def x_hat2_slow(counter_values):
    result = np.empty(len(counter_values))
    N = len(counter_values[0])
    for row in range(len(counter_values)):
        result[row] = np.exp(np.sum(np.log(counter_values[row] + 1)) / N)
    return result

# @numba.jit
def _calculate_pvalues(test_mean, model_mean, est_var, length):    
    results = np.empty(length)
    for i in range(length):
        results[i] = pvalue_est(test_mean[i], model_mean[i], est_var[i])
    return results

# @numba.jit
def pvalue_est(x, mean, adj_var):
    # mean and var must be shame length
    p = 1 - (mean / adj_var)
    r = (mean * mean) / (adj_var - mean)
    print(p,r,nbinom.cdf(x,r,p))
    return nbinom.cdf(x,r,p)
