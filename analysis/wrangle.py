"""
Module of functions for wrangling the data output by xmc (typically in the
form of a pandas dataframe)


Contains the following functions:

For wrangling the blob parameter dataframe:
 filterblobs()
 simplefilterblobs() [mainly a helper function for filterblobs()]

Statistics:
 (from the wquantiles package: 
  https://github.com/nudomarinero/wquantiles/weighted.py,
  modified by KAF to work without weights)
 quantile()
 quantile_1D() [mainly helper function for quantile()]
 weighted_median()


"""
#----------------------------------------------------------------
# Import Common Modules
import pandas as pd
import numpy as np
#from __future__ import print_function

#----------------------------------------------------------------
"""
 filterblobs() and simplefilterblobs()

Author: Kari A. Frank 
Date: October 20, 2015
Purpose: Filter a pandas dataframe (i.e. remove rows) based on 
         minimum and maximum values for the specified column.
         
         
         
         
Usage: filtered_df = filterblobs(inframe,colnames,minvals=minvals,maxvals=maxvals)

Input:

  inframe -- pandas dataframe
               
  colnames  -- scalar or list of (string) names of the columns containing
               the parameters be filtered by (e.g. 'blobkT' to filter based on 
               minimum or maximum temperature)            
  

  minvals,maxvals  -- scalars or lists of minimum and maximum values of
                      parameters to keep. returned value range is inclusive.
                      order must match the order of columns

Output:
  
  Returns a dataframe identical to the input dataframe but missing rows
    which do not match the criteria.

Usage Notes:
 - to return rows for which the column equals a specific value (rather
   than is within a range), set minval=maxval

 - to include only upper or lower limits for any parameter, set the 
   associated minval or maxval to None  

Example:
    filtered_df = filterblobs(blobframe,'blobkT',min=0.5,max=1.5)
    - this will return a version of blobframe which includes only 
      rows (blobs) with temperatures between 0.5 and 1.5.

    filtered_df = filterblobs(blobframe,['blobkT','blob_Si'],minvals=[0.5,1.0],
                         maxvals=[1.5,3.0])
    - this will return a version of blobframe which includes only 
      rows (blobs) with temperatures between 0.5 and 1.5 and Si between
      1.0 and 3.0.
"""

#----import modules---

#----define function to filter by 1 parameter----
#--gets called by the primary filter function below--
def simplefilterblobs(inframe,colname,minval=None,maxval=None):

    #-check if colname is valid column name-
    if colname not in inframe.columns:
        raise ValueError(colname+" is not a valid column name")

    #-filter based on value-

    #-set default minval/maxval if none provided-
    if minval == None: minval = min(inframe[colname])
    if maxval == None: maxval = max(inframe[colname])

    #-filter the dataframe-
    outframe = inframe[(inframe[colname]<=maxval) & 
                       (inframe[colname]>=minval)]

    #-warn if zero rows match-
    if len(outframe.index) == 0:
        print "Warning: Returned dataframe has zero rows."

    #-warn if all rows match-
    if len(outframe.index) == len(inframe.index):
        print "Warning: Nothing was filtered (all rows match criteria)."

    #-return filtered dataframe-
    return outframe

def filterblobs(inframe,colnames,minvals=None,maxvals=None):

    #--initialize outframe--
    outframe = inframe

    #--check if multiple parameters--
    if type(colnames) is str:
        colnames = [colnames]
        minvals = [minvals]
        maxvals = [maxvals]

    #--filter on each parameter--
    for col in range(len(colnames)):
        outframe = simplefilterblobs(outframe,colname=colnames[col],
                                minval=minvals[col],maxval=maxvals[col])

    return outframe

#----------------------------------------------------------------

def quantile_1D(data, weights, quantile):
    """
    Compute the weighted quantile of a 1D numpy array.

    Parameters
    ----------
    data : ndarray
        Input array (one dimension).
    weights : ndarray
        Array with the weights of the same size of `data`.
    quantile : float
        Quantile to compute. It must have a value between 0 and 1.

    Returns
    -------
    quantile_1D : float
        The output value.
    """
    # Check the data
    if not isinstance(data, np.matrix) :
        data = np.asarray(data)
    if not isinstance(weights, np.matrix) :
        weights = np.asarray(weights)
    nd = data.ndim
    if nd != 1:
        raise TypeError("data must be a one dimensional array")
    ndw = weights.ndim
    if ndw != 1:
        raise TypeError("weights must be a one dimensional array")
    if data.shape != weights.shape:
        raise TypeError("the length of data and weights must be the same")
    if ((quantile > 1.) or (quantile < 0.)):
        raise ValueError("quantile must have a value between 0. and 1.")
    # Sort the data
    ind_sorted = np.argsort(data)
    sorted_data = data[ind_sorted]
    sorted_weights = weights[ind_sorted]
    # Compute the auxiliary arrays
    Sn = np.cumsum(sorted_weights)
    # TODO: Check that the weights do not sum zero
    Pn = (Sn-0.5*sorted_weights)/np.sum(sorted_weights)
    # Get the value of the weighted median
    return np.interp(quantile, Pn, sorted_data)

#----------------------------------------------------------------
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

#----------------------------------------------------------------

def weighted_median(data, weights=None):
    """
    Weighted median of an array with respect to the last axis.

    Alias for `quantile(data,0.5,weights=weights)`.
    """
    return quantile(data, 0.5,weights=weights)

#----------------------------------------------------------------
