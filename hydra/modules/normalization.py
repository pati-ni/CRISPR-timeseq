import numba
import numpy as np
import argparse




@numba.jit
def x_hat(counter_values):
    result = np.empty(len(counter_values))
    N = len(counter_values[0])
    for row in range(len(counter_values)):
        # Geometric mean, add one to avoid zero values
        result[row] = np.exp(np.sum(np.log(counter_values[row] + 1)) / N)
    return result



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Median ratio normalization for ascii files")
    parser.add_argument("-i","--input", help = "Input File")
    parser.add_argument("-o","--output", help = "Output File", default=None)
    parser.add_argument("-d","--delimiter", help = "Values delimiter", default=',')
    parser.add_argument("-n","--normalize", help = "Avoid dealing with 0 values",default = True)
    args = parser.parse_args()
