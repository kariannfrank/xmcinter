"""
Library to compute weighted quantiles, including the weighted median, of
numpy arrays.

From the wquantiles package:  
https://github.com/nudomarinero/wquantiles/weighted.py

Modified by KAF to work without weights.

"""
from __future__ import print_function
import numpy as np
import quantile_1D

__version__ = "0.3"

def quantile(data, quantile, weights = None):
    """
    Weighted quantile of an array with respect to the last axis.

    Parameters
    ----------
    data : ndarray
        Input array.
    weights : ndarray
        Array with the weights. It must have the same size of the last 
        axis of `data`.
    quantile : float
        Quantile to compute. It must have a value between 0 and 1.

    Returns
    -------
    quantile : float
        The output value.
    """


    # TODO: Allow to specify the axis
    nd = data.ndim
    #check for weights and default to all weights = 1 if none provided
    if weights == None: 
        weights = np.ones(data.shape)
    if nd == 0:
        TypeError("data must have at least one dimension")
    elif nd == 1:
        return quantile_1D(data, weights, quantile)
    elif nd > 1:
        n = data.shape
        imr = data.reshape((np.prod(n[:-1]), n[-1]))
        result = np.apply_along_axis(quantile_1D, -1, imr, weights, quantile)
        return result.reshape(n[:-1])
