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
        if x[i] < mean[i]:
            results[i] =  nbinom.cdf(x[i], r[i], p[i])
        else:
            results[i] = 1 - nbinom.cdf(x[i], r[i], p[i])
    return results
