"""
From the wquantiles package:  
https://github.com/nudomarinero/wquantiles/weighted.py

Modified by KAF to work without weights.

"""
import numpy as np
import quantile

def weighted_median(data, weights=None):
    """
    Weighted median of an array with respect to the last axis.

    Alias for `quantile(data,0.5,weights=weights)`.
    """
    return quantile(data, 0.5,weights=weights)

