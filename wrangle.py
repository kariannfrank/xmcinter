"""
Module of functions for wrangling the data output by xmc (typically in the
form of a pandas dataframe)

Contains the following functions:

For wrangling the blob parameter dataframe:
 filterblobs()
 simplefilterblobs() [mainly a helper function for filterblobs()]
 filtercircle()
 gaussian_volume()

Statistics:
 (from the wquantiles package: 
  https://github.com/nudomarinero/wquantiles/weighted.py,
  modified by KAF to work without weights)
 quantile()
 quantile_1D() [mainly helper function for quantile()]
 weighted_median()
 weighted_posterior()
 weighted_modes()
 credible_region()

"""
#----------------------------------------------------------------
# Import Common Modules
import pandas as pd
import numpy as np
#import astropy.stats as astats #for the histogram function
#(switch to the astropy histogram and use the optimal binning (bins='knuth')
# once they implement weights)
#from __future__ import print_function

#----------------------------------------------------------------

#----define function to filter by 1 parameter----
#--gets called by the primary filter function below--
def simplefilterblobs(inframe,colname,minval=None,maxval=None,quiet=False):
    """See docstring for filterblobs()"""

    #-check if colname is valid column name-
    if colname not in inframe.columns:
        raise ValueError(colname+" is not a valid column name.")

    #-filter based on value-

    #-set default minval/maxval if none provided-
    if minval is None: minval = min(inframe[colname])
    if maxval is None: maxval = max(inframe[colname])

    #-filter the dataframe-
    outframe = inframe[(inframe[colname]<=maxval) & 
                       (inframe[colname]>=minval)]

    #--print warnings--
    if quiet == False:

        #-warn if zero rows match-
        if len(outframe.index) == 0:
            print ("simplefilterblobs: Warning: Filtered" 
                   " dataframe has zero rows.")

        #-warn if all rows match-
        if len(outframe.index) == len(inframe.index):
            print( "simplefilterblobs: Warning: Nothing was "
                   "filtered (all rows match criteria).")

    #-return filtered dataframe-
    return outframe

def filterblobs(inframe,colnames,minvals=None,maxvals=None,logic='and'):
    """
    filterblobs() and simplefilterblobs()

    Author: Kari A. Frank 
    Date: October 20, 2015

    Purpose: Filter a pandas dataframe (i.e. remove rows) based on 
             minimum and maximum values for the specified column.

    Input:

      inframe (pandas.DataFrame): dataframe to be filtered

      colnames (str or list of str): names of the columns containing 
           the parameters be filtered by (e.g. 'blobkT' to filter based on 
           minimum or maximum temperature)            

      minvals,maxvals (numerical or list of numerical): minimum and maximum
           values of parameters to keep. returned value range is inclusive.
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

        filtered_df = filterblobs(blobframe,['blobkT','blob_Si'],
                                  minvals=[0.5,1.0],maxvals=[1.5,3.0])
        - this will return a version of blobframe which includes only 
          rows (blobs) with temperatures between 0.5 and 1.5 and Si between
          1.0 and 3.0.
    """

    #--initialize outframe--
#    outframe = inframe

    #--check for valid logic--
    if logic != 'and':
        if logic != 'or':
            print ("filterblobs: Warning: Invalid logic type."
                   "  Using logic='and'")
            logic = 'and'

    #--check if multiple parameters--
    if type(colnames) is str:
        colnames = [colnames]
        minvals = [minvals]
        maxvals = [maxvals]

    #--filter on each parameter--
    filtered_frames = [simplefilterblobs(inframe,colname=colnames[c], \
                       minval=minvals[c],maxval=maxvals[c] \
                       ,quiet=True) for c in range(len(colnames))]
    
    #--combine filtered dataframes--
    if logic == 'or':
        #-keep all blobs which meet any of the criteria-
        outframe = pd.concat(filtered_frames,axis=0,join='outer')
    else: 
        #logic == 'and':
        #-keep only blobs which meet all of the criteria-
        outframe = pd.concat(filtered_frames,axis=0,join='inner')

        
#    for col in range(len(colnames)):
#        newframe = simplefilterblobs(inframe,colname=colnames[col],
#                                   minval=minvals[col],maxval=maxvals[col]
#                                     ,quiet=True) 
    #-warn if zero rows match-
    if len(outframe.index) == 0:
        print "filterblobs: Warning: Filtered dataframe has zero rows."

    #-warn if all rows match-
    if len(outframe.index) == len(inframe.index):
        print ("filterblobs: Warning: Nothing was filtered"
               " (all rows match criteria).")

    #-return filtered dataframe-
    return outframe

#----------------------------------------------------------------
def filtercircle(inframe,x='blob_phi',y='blob_psi',r='blob_sigma',
                 radius=20.0,x0=-60.0,y0=80.0,logic='exclude',
                 fraction=False):
    """
    filtercircle()

    Author: Kari A. Frank 
    Date: March 30, 2016

    Purpose: Filter a pandas dataframe (i.e. remove rows) based on 
             the given circular region.

    Input:

      inframe (pandas.DataFrame): dataframe to be filtered

      x,y (str or list of str): names of the columns containing 
           the parameters to use as x and y circle coordinates.  
           default is to use 'blob_phi' and 'blob_psi'.
       
      r (str or numerical): size of the blob in x,y dimension.
           if provided, will use radius - r to check distance of
           each blob from x0,y0. This allows non-pointlike blobs
           which contribute to the emission in the desired region
           to be kept/exluded, even if the blob center is outside 
           said region. if a string is provided, it must correspond
           to the name of column in inframe. if numerical value 
           provided, that value will be used for all blobs.
           default = 'blob_sigma'. to treat blobs like point sources,
           set r=0.

      x0,y0 (numerical): center coordinates of the circle, in
           same units as x and y

      radius (numerical): radius of the circle, in same units as x and y

      logic (string): string to specify if the circular region should be
           removed (logic='exclude', default), or if everything outside the
           circle should be removed (logic='include')

      fraction (bool): switch to add extra column that specifies fraction
           of the blob volume inside (if logic='include') or outside (if
           logic='exclude') of the specified circle, instead of dropping
           any rows. If smaller than 0.01 (<1% is in the include region),
           the fraction will be set to 0. If >0.99, will be set to one.
           this allows more flexibility in deciding what is considered
           'significant' contribution of emission to a region, and allows
           weighting of blobs by this fraction, similar to how it is 
           done in the weighted maps.
           - if fraction=True, then the r argument is ignored.
           - NOT YET IMPLEMENTED (need function for gaussian 
             integration of blob over circular region)

    Output:

      Returns a dataframe identical to the input dataframe but missing 
        rows which are inside (or outside if logic='include') the 
        defined circle.

    Usage Notes:
     - be very careful if x and y do not have the same units

    Example:
        filtered_df = filtercircle(blobframe,radius=20.0,x0=-40.0,y0=80.0)
        - this will return a version of blobframe which excludes all blobs
          within a circle (in phi, psi coordinates) of radius 20" centered
          on phi=-40.0,psi=80.0

        filtered_df = filtercircle(blobframe,x='blob_Si',y='blob_Fe',
                                   radius=0.2,r=None,
                                   x0=1.0,y0=1.0,logic='include')
        - this will return a version of blobframe which includes only
          rows (blobs) contained in a circle of radius 0.2 centered on 
          1,1 in the Si-Fe plane.

    """
    #--check for valid logic--
    if logic != 'exclude':
        if logic != 'include':
            print ("filtercircle: Warning: Invalid logic type."
                   "  Using logic='exclude'")
            logic = 'exclude'

    #--save original size--
    inblobs = len(inframe.index)

    #--add column to dataframe that is distance from center--
    inframe['R'] = ((inframe[x]-x0)**2.0+(inframe[y]-y0)**2.0)**0.5 
    #-adjust for blob size-
    if (r is not None) and (r != 0) and (fraction is False):
        if isinstance(r,str): 
            inframe['R'] = inframe['R'] - inframe[r]
        else:
            inframe['R'] = inframe['R'] - r

    #--filter--
    if fraction is False:
        # find all blobs in circle - defining blobs to drop
        if logic == 'exclude':
            circleframe = filterblobs(inframe,'R',maxvals=radius)
        if logic == 'include':
            circleframe = filterblobs(inframe,'R',minvals=radius)

        # drop blobs
        inframe.drop(circleframe.index,inplace=True)
        # remove extra column
        inframe.drop('R',1,inplace=True)

#need to make function for gaussian integral in spherical coords
    if fraction is True:
        print "ERROR: fraction argument not yet implemented."
#        y_integrals = xm.gaussian_integral()
#        inframe['circlefraction'] = 

#    if logic == 'exclude':
#        outframe = inframe[ ( (inframe[x]-x0)**2.0+(inframe[y]-y0)**2.0)**0.5 >= radius ]
#    else:
#        outframe = inframe[ ( (inframe[x]-x0)**2.0+(inframe[y]-y0)**2.0)**0.5 <= radius ]

    #-warn if zero rows match-
    if len(inframe.index) == 0:
        print "filtercircle: Warning: Filtered dataframe has zero rows."

    #-warn if all rows match-
    if len(inframe.index) == inblobs:
        print ("filtercircle: Warning: Nothing was filtered"
               " (all rows match criteria).")

    #-return filtered dataframe-
    return inframe

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
    if weights is None: 
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
#def weighted_posterior(data, weights=None, bins='knuth',normalize=False):
def weighted_posterior(data, weights=None, bins=None,normalize=False):
    """
    Construct a (weighted) posterior histogram from an array of data,
    e.g. a 1D array of blob temperatures.

    The output histogram is normalized to 1 if normalize=True.

    Returns two 1D arrays, one containing the x and y values of the 
    histogram.  The x values are the center of each bin.
    """

    #-determine number of bins if not provided-
    if bins is None:
        bins = np.ceil(np.sqrt(len(data)))

    #--Create histogram--
#    y,binedges = astats.histogram(data,weights=weights,bins=bins,
#                                  density=normalize)

    y,binedges = np.histogram(data,weights=weights,bins=bins)

    #--shift x values to center of bin--
    x = np.array([(binedges[i+1]+binedges[i])/2.0 for i in 
                  range(len(binedges)-2)])
    
    return x,y

#----------------------------------------------------------------
def weighted_modes(data, weights=None):
    """
    Calculate the mode of a weighted array.

    Returns the mode.

    To Do:
    Add capability to find all the local maxima.
    """

    #--Construct the binned posterior--
    postx,posty = weighted_posterior(data,weights=weights,
                                     normalize=True)

    #--Find the mode(s)--
#    modes = []
    
    #-find first mode-
    mi = np.argmax(posty)
#    modes[0] = postx[mi]
    modes = postx[mi]
    #-delete bin-
 #   posty = np.delete(posty,mi)
 #   postx = np.delete(postx,mi)

    #-find additional modes-

    return modes

#----------------------------------------------------------------
def credible_region(data, weights=None, frac=0.9, method='HPD'):
    """
    Calculate the Bayesian credible region from the provided posterior distribution.

    Parameters
    ----------
    data : ndarray
        Input array.
    weights : ndarray
        Array with the weights. It must have the same size of the last 
        axis of `data`.
    frac : float
        Credible region to calculate, e.g. frac=0.9 will return the 90% 
        credible region. Must be between 0 and 1.
    method : str
        Optionally specify the method for calculating the credible interval.
        Choose from 'HPD' (default) or 'ET':
        HPD = Highest Probability Density, defined as the narrowest interval
              that contains frac (e.g. 90%) of the total probability. For a 
              unimodal distribution, this will include the mode.
              This returns the relevant end point of the distribution 
              if necessary as an upper or lower limit.
        ET = Equal-tailed interval, defined as the interval where the 
             probability of being below the interval is equal to the 
             probability of being above.  This interval will include median.
             Note that using HPD is preferable, as ET cannot properly deal
             with upper or lower limits, or distributions which peak close
             to one end.


    Returns
    -------
    A tuple of floats containing the mode (HPD) or median (ET) and the 
    endpoints of the credible region.
    In the case of method='HPD' and multiple modes, it will return a list
    of tuples, containing the modes and endpoints for each mode.
    implemented]

    Usage Notes
    -------
    Multipe modes not yet implemented.
    ET method not yet implemented.
    """

    #----HPD calculation----
    
    #--find the mode(s)--
    #in future, weighted_modes will return all local maxima, so can 
    #calculate multiple credible intervals this way, one for each maxima
    modes = weighted_modes(data,weights=weights)
    
    #--get normalized posterior (probability density function)--
    postx,posty = weighted_posterior(data,weights=weights,normalize=True)

    #--step through bins around mode, alternating directions--
    lowprob = 0.0
    highprob = 0.0
    lowfound = False
    highfound = False

    #-set starting indices at the mode
    highxi = int(np.where(postx == modes)[0])
    lowxi = highxi

    while (lowfound is False) or (highfound is False):
        if lowfound is False:
            # check if at edge of histogram
            if lowxi <= 0:
                #lowest bin - complete loop for this bin, 
                # but don't repeat
                lowfound = True 
            # accumulate probability below mode
            lowprob = lowprob+posty[lowxi]
            lowxi = lowxi - 1
        if highfound is False:
            # check if at edge of histogram
            if highxi >= len(postx)-1:
                #highest bin - complete loop for this bin, 
                # but don't repeat
                highfound = True 
            # accumulate probability behigh mode
            highprob = highprob+posty[highxi]
            highxi = highxi + 1

        # sum upper and lower probabilities
        netprob = lowprob + highprob
        # check against frac-
        if netprob >= frac:
            lowfound = 1
            highfound = 1
    
    #-save mode and final interval (whether or not desired frac was reached)
    interval = (modes,postx[lowxi+1],posty[highxi-1])
            
    #-check if both upper and lower limits encountered-
    if lowxi+1 == 0: print "Lower limit reached."
    if highxi-1 == len(postx)-1: print "Upper limit reached."
    if netprob < frac: print ("WARNING: Only reached "+str(netprob),
                              " credible interval, parameter range is",
                              " too narrow.")
            

    #----ET calculation----
    
    #--find the median--
    median = weighted_median(data,weights=weights)


    return interval

#----------------------------------------------------------------
def gaussian_volume(sigma):
    """
    Calculate the volume of a spherical gaussian.

    Parameters
    ----------
    sigma : numeric or 1D array of numeric
        Gaussian width(s)

    Returns
    -------
    Volume (float) or volumes (1D array of floats).

    Usage Notes
    -------

    """

    volume = (2.0*np.pi*np.square(sigma))**1.5
    
    return volume
