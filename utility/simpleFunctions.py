__author__ = 'Haohan Wang'

import numpy as np

def centralize(x):
    m = np.mean(x)
    return x - m

def mapping2ZeroOne(x):
    maxi = np.max(x)
    mini = np.min(x)
    return (x-mini)/(maxi-mini)

def rescale(x):
    maxi = np.max(np.abs(x))
    if maxi == 0:
        return x
    return x/maxi

def normalize(x):
    m = np.mean(x)
    s = np.std(x)
    return (x - m) / s

def discretize(X):
    X[X>0.9] = 1
    X[X<0.9] = 1e-3
    return X

def reOrder(K):
    import scipy.cluster.hierarchy as hier
    Z = hier.linkage(K, 'ward')
    l = hier.leaves_list(Z)
    nK = np.zeros_like(K)
    for i in range(len(l)):
        nK[i,:] = K[l[i],:]
    return nK