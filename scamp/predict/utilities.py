"""
Utilities for predictors.
"""
import numba as nb
from numba import njit
import numpy as np

@njit
def summarize_count_distribution(counts, min_count=2):

    sub_counts = counts[counts > min_count]
    if len(sub_counts) < 2:
        return np.nan, np.nan, np.nan, [np.nan]*10, np.nan
    
    mu = np.nanmean(sub_counts)
    sig = np.nanvar(sub_counts)
    dispersion = sig / mu
    quantiles = [np.nanpercentile(sub_counts, p) for p in np.arange(0, 100, step=10)]

    iqr = np.nanpercentile(sub_counts, 75) - np.nanpercentile(sub_counts, 25)
    
    return mu, sig, dispersion, quantiles, iqr

