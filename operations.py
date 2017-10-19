import numba
import numpy as np

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
