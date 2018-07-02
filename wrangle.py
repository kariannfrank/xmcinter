"""
Module of functions for wrangling the data output by xmc.
Data typically should be passed as a pandas DataFrame or Series.

Contains the following functions:

Modifying or creating columns:
 filterblobs()
 simplefilterblobs() [mainly a helper function for filterblobs()]
 filtercircle()
 line_emissivities()

Mapping columns or dataframe into new array:
 make_histogram()
 normalize_histogram()
 make_spectrum()

Statistics (outputs scalars):
 (from the wquantiles package: 
  https://github.com/nudomarinero/wquantiles/weighted.py,
  modified by KAF to work without weights)
 quantile()
 quantile_1D() [mainly helper function for quantile()]
 weighted_median()
 weighted_posterior()
 weighted_modes()
 credible_region()
 weighted_std()
 iter_err()

Functions that act on a single blob (scalar) can be found in
the astro_utilities module.

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
def circlefraction(x,y,r,x0,y0,radius,shape='gauss',
                   use_ctypes=False):
    """
    Helper function for filtercircle(), to do the integration
    
    Author: Kari A. Frank
    Date: October 21, 2016

    Purpose: Calculates fraction of blob volume in specified 
             (projected) circle.
    
    Input:

      x,y (numerical): x and y position of the blob

      r (numerical): radius of the blobs

      x0,y0 (numerical): x and y position of the circle

      radius (numerical): radius of the circle

      shape (str): shape of the blob, either gauss or sphere
                   - !!sphere not yet implemented!!

      use_ctypes (bool): specify whether to use ctypes functionality in
           gaussian_integrate_dblquad (using scipy.integrate). See 
           xmcmap for more details.

    Output:

      Returns the fraction of the blob's volume within the circle.

    Usage Notes:
    
      - x,y,r,x0,y0,radius should all be in the same units, usually 
        spatial (usually arcsec)
    
    """

    #--import gaussian_integral and related functions--
    from astro_utilities import gaussian_volume
    import xmcmap as xm
        # for improved speed, use 
        #  gaussian_integral_quad instead, but it requires extra 
        #  package scipy.integrate.    

    # get blob volume in correct units
    if shape not in ['gauss','sphere']: 
        print "Warning: Unrecognized blob shape. Using shape='gauss'"
        shape = 'gauss'
    if shape == 'gauss':
        volume = gaussian_volume(r)
    elif shape == 'sphere':
        volume = (4.0/3.0)*np.pi*r**3.0

    # define functions x(y) for lower and upper x integration limits
    def bottomcircle(y,x0=x0,y0=y0,r0=radius):
        x = -1.0*((r0**2.0)-(y-y0)**2.0)**0.5 + x0
        return x
    def topcircle(y,x0=x0,y0=y0,r0=radius):
        x = ((r0**2.0)-(y-y0)**2.0)**0.5 + x0
        return x

#    print 'topcircle() = ',

    # define y integration limits
    lowery = y0-radius
    uppery = y0+radius
    
    # integrate over circle
    xy_integrals = xm.gaussian_integral_dblquad(bottomcircle,topcircle,lowery,uppery,x,y,r,use_ctypes=False)

    # z integral, -inf to inf
    z_integrals = (2.0*np.pi*r**2.0)**0.5    

    # total integral (volume in projected circle)
    involume = xy_integrals*z_integrals

    # get blob volumes in arcsec^3
    totalvolume = gaussian_volume(r)

    return involume/totalvolume
            
#----------------------------------------------------------------
def circlefraction_star(arglist):
    """Function to unpack list of arguments and pass to circlefractions"""
    return circlefraction(*arglist)
            

#----------------------------------------------------------------
def filtercircle(inframe,x='blob_phi',y='blob_psi',r='blob_sigma',
                 r0=None,x0=-60.0,y0=80.0,logic='exclude',
                 fraction=True,regname='circle',use_ctypes=True,
                 parallel=False,nproc=3):
    """
    filtercircle()

    Author: Kari A. Frank 
    Date: March 30, 2016

    Purpose: Filter a pandas dataframe (i.e. remove rows) based on 
             the given circular region.

    Input:

      inframe (pandas.DataFrame): dataframe to be filtered

      x,y (str or list of str): names of the columns containing 
           the parameters to use as x and y blob coordinates.  
           default is to use 'blob_phi' and 'blob_psi'.
       
      r (str or numerical): size of the blob in x,y dimension.
           if provided, will use r0 - r to check distance of
           each blob from x0,y0. This allows non-pointlike blobs
           which contribute to the emission in the desired region
           to be kept/excluded, even if the blob center is outside 
           said region. if a string is provided, it must correspond
           to the name of column in inframe. if numerical value 
           provided, that value will be used for all blobs.
           default = 'blob_sigma'. to treat blobs like point sources,
           set r=0.

      x0,y0 (numerical): center coordinates of the circle, in
           same units as x and y

      r0 (numerical): radius of the circle, in same units as x and y

      logic (string): string to specify if the circular region should be
           removed (logic='exclude', default), or if everything outside the
           circle should be removed (logic='include')

      fraction (bool): switch to add extra column that specifies fraction
           of the blob volume inside (if logic='include') or outside (if
           logic='exclude') of the specified circle, instead of dropping
           any rows. 
           this allows more flexibility in deciding what is considered
           'significant' contribution of emission to a region, and allows
           weighting of blobs by this fraction, similar to how it is 
           done in the weighted maps.

       regname (string): name of the circular region, to be used in new
           column names. the weight columns (if fraction=True) will be 
           named regname+'_incl_fraction' (or '_excl_', according to
           logic), and similarly the distance column regname+'_incl_dist'

        use_ctypes (bool): specify whether to use ctypes in 
           scipy.integrate.dblquad. ignored if fraction is False

        parallel (bool): specify do integration in parallel. ignored
           if fractions is False

        nproc (int): number of processors to use if parallel is True

    Output:

      Returns a dataframe identical to the input dataframe, except:
      - contains an extra column that is distance of each blob from 
        the circle center
      - if fraction is False: missing rows which are inside (or 
        outside if logic='include') the defined circle
      - if fraction is True: with an extra column that specifies the 
        fraction of the blob volume either inside or outside the 
        circular region (whichever is to be included), if fraction is True

    Usage Notes:
     - be very careful if x and y do not have the same units

    Example:
        filtered_df = filtercircle(blobframe,r0=20.0,x0=-40.0,y0=80.0)
        - this will return a version of blobframe which excludes all blobs
          within a circle (in phi, psi coordinates) of radius 20" centered
          on phi=-40.0,psi=80.0

        filtered_df = filtercircle(blobframe,x='blob_Si',y='blob_Fe',
                                   r0=0.2,r=None,
                                   x0=1.0,y0=1.0,logic='include')
        - this will return a version of blobframe which includes only
          rows (blobs) contained in a circle of radius 0.2 centered on 
          1,1 in the Si-Fe plane.

    """
    from multiprocessing import Pool
    from astro_utilities import distance

    #--check for valid logic--
    if logic != 'exclude':
        if logic != 'include':
            print ("filtercircle: Warning: Invalid logic type."
                   "  Using logic='exclude'")
            logic = 'exclude'

    #--save original size--
    inblobs = len(inframe.index)

    #--copy dataframe to avoid changing the original--
    outframe = inframe.copy(deep=False)

    #--add column to dataframe that is distance from center--
    rcol = regname+'_'+logic[0:3]+'_dist'
    outframe[rcol] = distance(outframe[x],outframe[y],x0,y0)

    #-adjust for blob size-
    if (r is not None) and (r != 0) and (fraction is False):
        if isinstance(r,str): 
            outframe[rcol] = outframe[rcol] - outframe[r]
        else:
            outframe[rcol] = outframe[rcol] - r
    else: # if fraction=True 
        if isinstance(r,str):
            outframe['blob_size'] = outframe[r]
        else:
            outframe['blob_size'] = r

    #--filter--
    if fraction is False:
        # find all blobs in circle - defining blobs to drop
        if logic == 'exclude':
            circleframe = filterblobs(outframe,rcol,maxvals=r0)
        if logic == 'include':
            circleframe = filterblobs(outframe,rcol,minvals=r0)

#        print 'number of bad blobs = ',len(circleframe.index)

        # drop blobs
        outframe.drop(circleframe.index,inplace=True)
        # remove extra column
#        outframe.drop(rcol,1,inplace=True)
#        print 'number of good blobs = ',len(outframe.index)

    if fraction is True:

        if parallel is True:
            
            # build argument list
            args=[[]]*len(outframe.index)
            i = 0
            for b in outframe.index:
                args[i] = [outframe.loc[b,x],outframe.loc[b,y],outframe.loc[b,'blob_size'],x0,y0,r0,'gauss',use_ctypes] 
                i = i+1

            # calculate fractions in parallel
            pool=Pool(nproc)
            fractions = np.array(pool.map(circlefraction_star,args))
            pool.close()
            pool.join()

        if parallel is False:
            # calculate fractions
            fractions = outframe.apply(lambda d: circlefraction(d[x],d[y],d['blob_size'],x0,y0,r0,use_ctypes=use_ctypes),axis=1)

        # weights = amount of 'good' emission
        fraccol = regname+'_'+logic[0:3]+'_fraction'
        if logic == 'exclude':
            # remove emission inside circle
            outframe[fraccol] = 1.0-fractions
        else:
            # keep only emission inside circle
            outframe[fraccol] = fractions

        # remove extra column
        outframe.drop('blob_size',1,inplace=True)

    if (fraction is False):
    #-warn if zero rows match-
        if (len(outframe.index) == 0):
            print "filtercircle: Warning: Filtered dataframe has zero rows."

    #-warn if all rows match-
        if len(outframe.index) == inblobs:
            print ("filtercircle: Warning: Nothing was filtered"
                   " (all rows match criteria).")

    #-return filtered dataframe-
    return outframe

#----------------------------------------------------------------

def quantile_1D(data,weights=None,quantile=0.5):
    """
    Compute the weighted quantile of a 1D numpy array.

    Parameters
    ----------
    data (pd.Series) : Input data values (one dimension)

    weights (pd.Series) : Weights corresponding to data

    quantile (float) : Quantile to compute. It must have a value between 
                       0 and 1.

    Returns
    -------
    quantile_1D (float): The output value.

    """

    # Coerce data into numpy array
    data = np.asarray(data)
    
    if weights is not None:
        weights = np.asarray(weights)
    else:
        weights = np.ones_like(data)

    # Verify allowed value of quantile
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
    data (pd.Series) : Input data values (one dimension)

    weights (pd.Series) : Weights corresponding to data

    quantile (float) : Quantile to compute. It must have a value between 
                       0 and 1.

    Returns
    -------
    quantile : float
        The output value.

    Usage Notes
    -----------
    If data is multi-dimensional, then the array will be flattened
    before calculating the quantile.

    """

    # Coerce data into numpy array
    data = np.asarray(data)
    
    if weights is not None:
        weights = np.asarray(weights)
    else:
        weights = np.ones_like(data)
    
    nd = data.ndim
    if nd == 1:
        return quantile_1D(data, weights, quantile)
    elif nd > 1:
        n = data.shape
        imr = data.reshape((np.prod(n[:-1]), n[-1]))
        result = np.apply_along_axis(quantile_1D, -1, imr, weights, 
                                     quantile)
        return result.reshape(n[:-1])

#----------------------------------------------------------------
def weighted_median(data, weights=None):
    """
    Weighted median of an array with respect to the last axis.

    Alias for `quantile(data,0.5,weights=weights)`.
    """
    return quantile(data,0.5,weights=weights)

#----------------------------------------------------------------
def make_histogram(dataseries,weights=None,bins=50,logbins=False,
                   datarange=None,density=False,iterations=None,
                   centers=False,normalize=False):
    """
    Create histogram, including optional errorbars.

    Author: Kari A. Frank
    Date: November 29, 2016
    Purpose: Create histogram from given series.

    Input:

     dataseries:  a pandas series of data (e.g. column of a pandas 
                  dataframe)

     weights:     optionally provided a pandas series of weights which 
                  correspond to the values in datacolumn (e.g. emission 
                  measure)

     bins:        optionally specify the number of bins (default=30)
                  or an array containing the bin edge values. bins
                  is passed directly to numpy.histogram(), so can
                  also accept any of the strings recognized by that 
                  function for special calculation of bin sizes.
                  - 'auto' (but don't use if using weights)
                  - 'fd' (Freedman Diaconis Estimator)
                  - 'doane'
                  - 'scott'
                  - 'rice' 
                  - 'sturges'
                  - 'sqrt'
                 
     logbins:     boolean to specify if bins should be in log space. if
                  set, then bins should be the number of bins.

     datarange:   minimum and maximum values to include in the histogram.
                  
     centers:     boolean to specify if the bincenters should be also be 
                  returned (in addition to the binedges)

     density:     passed to histogram. if True, then returned histogram is
                  is the probability density function. In general, will
                  not use this (call normalize_histogram() on histy 
                  instead)

     iterations:  optionally provide a series of the iterations that 
                  correspond to each element of dataseries. if provided, 
                  will be used to calculate the errorbars for each bin
                  which will be plotted on the histogram. must be the same
                  length as dataseries if provided.


    Output:
    
     histy (numpy array, 1D) : series containing the y-values of the
        histogram bins

     histx (numpy array, 1D) : series containing the values of the 
        histogram bin edges (the x-axis)

     yerrors (numpy array, 1D) : series containing the error bars 
        associated with histy

    """

    #----Set up bins----
    if datarange is None:
        datarange = (dataseries.min(),dataseries.max())

    # set up log bins
    if logbins is True: 
        bins = np.logspace(np.log10(datarange[0]),np.log10(datarange[1]),
                           num=bins)

    if isinstance(bins,float): bins = int(bins)
    
    #----Create the weighted histogram----
    histy,binedges = np.histogram(dataseries,weights=weights,bins=bins,
                               density=density,range=datarange)
    
    #----Calculate errorbars----

    if iterations is not None:

        # combine series into single dataframe for grouping
        dfr = pd.concat([dataseries,iterations],axis=1)
        dfr.columns = ['data','iteration']

        if weights is not None: dfr['weight'] = weights.values

        # initialize array to hold histogram 'stack'
        niter = len(np.unique(iterations))
        nbins = len(histy)
        histstack = np.empty((nbins,niter))

        #--create histogram for each iteration--
        i = 0 # iteration (layer) counter
        gdf = dfr.groupby('iteration',as_index=False)
        for s,g in gdf: 
            if weights is not None: 
                gweights = g['weight']
            else:
                gweights = None
            hyi,hyierrors,hyedges = make_histogram(g['data'],
                                                   weights=gweights,
                                           bins=binedges,iterations=None,
                                           #logbins=logbins,
                                           density=density,
                                           datarange=datarange)
            histstack[:,i] = hyi
            i = i+1

        #--collapse stack to standard deviation in each bin--
        yerrors = np.std(histstack,axis=1)
#        errors[:] = np.average(histy) # large errors for testing

    else:
        yerrors = None

    #--normalize--
    if normalize is True:
        histy,yerrors = normalize_histogram(histy,yerrors=yerrors)

    #--return arrays--
    if centers is False:
        return histy,yerrors,binedges
    else:
        binsizes = binedges[1:]-binedges[:-1]
        bincenters = binedges[:-1]+binsizes/2.0
        return histy,yerrors,binedges,bincenters

#----------------------------------------------------------
def normalize_histogram(histy,yerrors=None):
    """Normalize histogram to go from 0 to 1 on y-axis, including scaling the errorbars if provided."""

    if yerrors is not None:
        # normalize errors
        yerrors = yerrors/float(np.max(histy))

    # normalize histogram
    histy = histy/float(np.max(histy))
    
    return histy,yerrors

#----------------------------------------------------------------
def weighted_modes(data, weights=None,bins=75,logbins=False):
    """
    Calculate the mode of a weighted array.

    Returns the mode.

    To Do:
    Add capability to find all the local maxima.
    """

    #--Construct the binned posterior--
#    postx,posty = weighted_posterior(data,weights=weights,
#                                     normalize=True)
    posty,postyerr,edges,postx = make_histogram(data,weights=weights,
                                                normalize=True,centers=True,
                                                bins=bins,logbins=logbins)


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
    Calculate the Bayesian credible region from the provided posterior 
    distribution.

    Parameters
    ----------
    data : ndarray or pd.Series
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
#    postx,posty = weighted_posterior(data,weights=weights,normalize=True)
    posty,postyerr,edges,postx = make_histogram(data,weights=weights,
                                     normalize=True,centers=True)

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
    interval = (modes,postx[lowxi+1],postx[highxi-1])
            
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
def weighted_std(data, weights=None):
    """
    Calculate weighted standard deviation of series.

    Parameters
    ----------
    data : pd.Series or np.ndarray
        series of values to find the standard deviation of

    weights : 1D numeric array (optional)
        series of weights to apply to data. must have the same
        size as data.

    Returns
    -------
    Weighted standard deviation of data (float).

    Usage Notes
    -------

    """
    

    # coerce into numpy arrays
    data = np.asarray(data)
    
    if weights is not None:
        weights = np.asarray(weights)
    else:
        weights = np.ones_like(data)

    # calculate standard deviation
    avg = np.average(data,weights=weights)
    a = ((data-avg)**2.0)*weights
    numerator = np.sum(a)
    denominator = np.sum(weights)

    return (numerator/denominator)**0.5

#----------------------------------------------------------------
def iter_err(df,param,function,weights=None,*args,**kwargs):
    """
    Calculate statistical error on a summary statistic.


    Parameters
    ----------
    df : pd.DataFrame
       dataframe of values. must include an 'iteration' column

    param : string
       column name of the parameter to pass to function

    function : python function name
       function to calculate the error of
       to calculate error on a sum (e.g. total em), use function='total'

    weights : string
       name of dataframe column to use as weights

    Returns
    -------
    Error of the value calculated by function (standard deviation of 
    function calculated for individual iterations).

    Usage Notes
    -------

    """
        
    # calculate number of iterations
    niter = np.unique(df['iteration']).size
    vals = []

#    print 'niter = ',niter
    
    # group by iteration
    layers = df.groupby('iteration')

    # iteration over groups (i.e. iterations)
    for i,group in layers:
        if weights is None:
            w = np.ones_like(group[param])
        else:
            w = group[weights]

        if function != 'total':
            val = function(group[param],weights=w)
        else:
            val = group[param].sum()
            #        print 'i = ',i
#        print 'val = ',val
        
        # add to list
        vals = vals + [val]
        
    return np.std(vals)


#----------------------------------------------------------------
def line_emissivities(kT,tolerance=0.05,**fetchargs):
    """
    Calculate the total emissivity of all emission lines from kT gas

    Parameters
    ----------
    kT : pd.Series or np.ndarray
        temperature values, in keV

    fetch_lines() input options:

         emissivity_range (2d tuple) : range of line emissivities 
                                       to include

         wavelength_range (2d tuple) : range of line wavelengths 
                             (angstroms) to include emission lines from

         energy_range (2d tuple) : range of line energies (keV) to include
                             emission lines from

         include_lines (list of strings) : list of element names to 
                        include emission lines from. will first drop 
                        all lines from other elements, then apply 
                        remaining criteria (e.g. kT_range).

         redshift (float) : redshift of the source. if non-zero, lines 
                            will be redshifted appropriately before
                            applying any criteria or plotting.


    Returns
    -------
    Sum of emissivities from all lines with tolerance of kT.

    Usage Notes
    -------

    """
    import os
    from astro_utilities import fetch_lines
    
    #----Read in lines----
    linefile = os.path.dirname(__file__)+'/xraylines_atomdb307.txt'
    atomdb = pd.read_table(linefile,sep='\s+',comment='#',engine='python')
   
    #----Check array type, convert to dataframe if needed----
    if not isinstance(kT,pd.DataFrame):
        if not isinstance(kT,pd.Series):
            kT = pd.Series(kT,name='kT') # convert to series
        kT = kT.to_frame() # convert series to dataframe
    #-rename column-
    kT.rename(columns={kT.columns[0]:'kT'},inplace=True)

    #--define gas temperature ranges--
    kT['kTrange'] = zip(kT.kT-tolerance,kT.kT+tolerance)

    #--fetch lines from atomdb file--
    kT['emissivities'] = kT['kTrange'].map(lambda x: fetch_lines(kT_range = x,atomdb=atomdb,total=True,**fetchargs))

    return kT.emissivities.values
    
#----------------------------------------------------------------
def blob_line_photons(kT,volume,time=None,tolerance=0.01,**fetchargs):
    """
    Calculate the maximum number of photons emitted due to lines

    Parameters
    ----------
    kT (pd.Series or np.ndarray) : temperatures of blobs, in keV

    volume (pd.Series or np.ndarray) : volume of blobs, in cm^3

    time (scalar) : total exposure time. if not provided, returned
                    values will be count rates intead of total counts

    fetch_lines() input options:

         emissivity_range (2d tuple) : range of line emissivities 
                                       to include

         wavelength_range (2d tuple) : range of line wavelengths 
                             (angstroms) to include emission lines from

         energy_range (2d tuple) : range of line energies (keV) to include
                             emission lines from

         include_lines (list of strings) : list of element names to 
                        include emission lines from. will first drop 
                        all lines from other elements, then apply 
                        remaining criteria (e.g. kT_range).

         redshift (float) : redshift of the source. if non-zero, lines 
                            will be redshifted appropriately before
                            applying any criteria or plotting.


    Returns
    -------
    Max number of photon counts (or count rate, if time!=None) expected
    from each blob, given its volume and temperature. Sums photons from
    all lines within tolerance of the blob kT.

    Usage Notes
    -------

    """

    # get total line emissivities for each blob
    #  units of photons*cm^3/s
    emiss = line_emissivities(kT,tolerance=tolerance,**fetchargs)
    
    # divide by volume
    counts = emiss/volume

    # multiply by exposure time
    if time is not None:
        counts = counts * time

    return counts

#----------------------------------------------------------------
def normclean(df,nHcol = 'blob_nH',kTcol = 'blob_kT',
              normcol = 'blob_norm',itercol='iteration',
              normthreshs = None,parallel=False,nproc=4):
    """
    Remove low kT - high nH blobs that do not meet norm criteria.

    Parameters
    ----------
    
    df (pd.DataFrame) : blob data, including all relevant columns

    nHcol,kTcol,normcol,itercol (string) : names of the relevant 
        columns, if not the usual names

    normthreshs (pd.DataFrame or file name) : matrix containg all the 
        nH-kT bins and  the associated norm thresholds for (minimum 
        norm required). should have columns 'kT','kThalfbin',
        'nH','nHhalfbin','norm'
    


    Returns
    -------
    pd.DataFrame that is copy of df, but with all blobs (rows) removed 
    that do not meet the minimum norm requirement (on a per-iteration 
    basis)

    Usage Notes
    -------
    df must contain a column name 'iteration'

    """

    from multiprocessing import Pool

    #--get thresholds--
    if isinstance(normthreshs,str):
        normthreshs = pd.read_table(normthreshs,comment='#',sep=r'\s+',
                                    index_col=False)

    #--copy dataframe--
    outdf = df.copy(deep=False)

    #--group by iteration--
    gdf = df.groupby(itercol)
    fargs = []

    #--loop over iterations--
    for i,group in gdf:
        # i = iteration number
        #--loop over nH-kT bins--
        if parallel is False:
            for ri,r in normthreshs.iterrows():
                b = group[( (group[kTcol]>r[0]-r[1]) & 
                            (group[kTcol]<r[0]+r[1]) & 
                            (group[nHcol]>r[2]-r[3]) & 
                            (group[nHcol]<r[2]+r[3]))]
                net_norm = b[normcol].sum()
                if net_norm < r[4]:
                    bad = b.index
                    outdf = outdf.drop(bad)
        else: # build argument lists
            fargs.append([group,normthreshs,kTcol,nHcol,normcol])

    # run iterations in parallel and collect bad indices
    if parallel is True:
        pool=Pool(nproc)
        badind = pool.map(normfilter_star,fargs)
        print type(badind)
        pool.close()
        pool.join()

        # combine lists of bad blobs (flatten badind)
        bad = [b for l in badind for b in l]

        # filter dataframe
        outdf = outdf.drop(bad)

    return outdf
#    return badind

def normfilter(gdf,normthreshs,kTcol,nHcol,normcol):

    bad = []
    for ri,r in normthreshs.iterrows():
        b = gdf[( (gdf[kTcol]>r[0]-r[1]) & 
                (gdf[kTcol]<r[0]+r[1]) & 
                (gdf[nHcol]>r[2]-r[3]) & 
                (gdf[nHcol]<r[2]+r[3]))]
        net_norm = b[normcol].sum()
        if net_norm < r[4]:
            bad.extend(list(b.index))

        return bad

def normfilter_star(arglist):
    """Function to unpack list of arguments and pass to normfilter"""
    return normfilter(*arglist)

#----------------------------------------------------------------
def make_spectrum(df,runpath='../',suffix='99999',oversim='auto',
                  bins=0.03,xlog=False,logbins=None,datarange=None,
                  fixed_nH=None):
    """
    Create spectrum (histogram) from the blobs in given dataframe

    Author: Kari A. Frank
    Date: May 4, 2017

    Input:

     df (pd.DataFrame) : blob dataframe. must contain, at minimum,
                   all columns in parameters.txt

     suffix (str) : must be string of an integer. will be appended to end
                    of created file names (e.g. deconvolution.suffix)
     
     bins:        optionally specify the number of bins (int), binsize 
                  (float), or an array containing the bin edge values. 
                  can also accept any of the strings recognized by 
                  np.histogram() for special calculation of bin sizes.
                  - 'auto' (but don't use if using weights)
                  - 'fd' (Freedman Diaconis Estimator)
                  - 'doane'
                  - 'scott'
                  - 'rice' 
                  - 'sturges'
                  - 'sqrt'
                  default is to create equally spaced bins of width 0.03,
                  (30 eV). note that xmc typically uses a binsize of 0.015.

    logbins, datarange, xlog, ylog : see make_histogram()

    oversim (numerical scalar or 'auto') : must specify the oversim value
                   to properly scale the model spectra.
                   - 'auto' (or any string): will read in the first 
                      statistic.* file it
                     finds in runpath/ and get the oversim from there.
                   - can also explicitly pass the oversim value

     fixed_nH (scalar) : fixed absorption value in the model, usually the 
                     Galactic NH for LMC/SMC objects. required if any of
                     the absorption parameters are fixed. 


    Output:

     Returns a spectrum as a histogram, in same format as returned by
     xmcfiles.read_spectra() -- (y,yerrors,yedges)

    Usage Notes:
     - Requires input, parameters.txt, start.xmc, events, and exposure
       map files all to be located in runpath directory
     - All new files are created in the working directory, so it is best
       to NOT run this from a run directory, to avoid the possibility 
       of overwriting any files.
     - The returned spectrum will be scaled by the number of iterations,
       i.e. will be equivalent to counts expected in a single iteration,
       for correct comparison with the data spectrum. 

   
    """
    #--imports--
    import os
    from xmcfiles import fake_deconvolution,read_spectra,merge_output
    from file_utilities import parse_file_line,ls_to_list
    from astro_utilities import abs_cross_section
    
    pwd = os.getcwd()+'/'

    #--get oversim--
    if isinstance(oversim,str):
        oversim = float(parse_file_line(runpath+ls_to_list(
            runpath,'statistic.*')[0])[2])

    
    #--create 'fake' deconvolution and other xmc output files--
    # (file names will end in .99999)
    niters=len(np.unique(df.iteration)) # get number iterations
    # scale by iterations
    for col in df.columns:
        if 'norm' in col:
            df[col] = df[col]/float(niters)
    oldnorm = sum([df[col].sum() for col in df.columns if 'norm' in col])
#    oldnorm = df.blob_norm.sum() # net norm
    df = fake_deconvolution(df,runpath=runpath,suffix=suffix)

    #--calculate weighted average absorption--
    if 'blob_nH2' in df.columns:
        nH = np.average(fixed_nH+df.blob_nH2,weights=df.blob_norm)
    else:
        if fixed_nH is None:
            nH = np.average(df.blob_nH,weights=df.blob_norm)
        else:
            nH = fixed_nH
    nH = nH*1e22
    print 'nH avg = ',nH
    
    #--read fake file into xmc and create spectrum--
    
    #-copy and modify input file-
    with open(runpath+'input','r') as file:
        inputscript = file.readlines()

    # find line with '/null' and replace with '/xw'
    indices = [i for i, s in enumerate(inputscript) if '/null' in s]
    index = indices[-1]
    inputscript[index] = inputscript[index].replace('null','xw')
        
    # find line with 'run' (in case blank lines after run)
    indices = [i for i, s in enumerate(inputscript) if 'run' in s]
    index = indices[-1]

    # get start.xmc file name and modify run line
    for p in inputscript[index].split(' '): # in case extra spaces
        if '.xmc' in p:
            startfile = p.rstrip() # remove any newlines
    
    inputscript[index] = 'run spec_start.xmc'
    with open(pwd+'input_spec','w') as file:
        file.writelines(inputscript)

    #-copy and modify the start.xmc file-
    with open(runpath+startfile,'r') as file:
        runscript = file.readlines()

    # find line with 'deconvolvec' (in case blank lines after)
    indices = [i for i, s in enumerate(runscript) if 'deconvolvec' in s]
    index = indices[-1]
        
    # modify deconvolvec line
    lastlineparts = [ p.rstrip() for p in runscript[index].split(' ')]

    if len(lastlineparts)==3: # add iteration
        lastlineparts.append(suffix+'.')
    else: # replace iteration
        lastlineparts[-1] = suffix+'.'

    lastline=' '.join(lastlineparts)
        
    runscript[index] = lastline

    # add path to data/expo calls
    dindices = [i for i, s in enumerate(runscript) if 'data ' in s]
    eindices = [i for i, s in enumerate(runscript) if 'exp' == s[:3]]
    indices = dindices + eindices
    for i in indices:
        parts = runscript[i].strip().split(' ')
        line = parts[0]+' '+runpath+parts[-1]+'\n'
        runscript[i] = line

    # write modified file
    with open(pwd+'spec_start.xmc','w') as file:
        file.writelines(runscript)

    #-call xmc to create the spectrum-
    # (for now do it manually)
    raw_input("In another terminal: \n\tFrom this directory"
              " ("+pwd+"),\n"
              "\trun xmc with 'xmc < input_spec'\n"
              "\tWhen spectrum_"+suffix+".fits is produced\n "
              "\tand the *."+suffix+" files are overwritten by xmc, "
              "then kill xmc.\n"
    "\tReturn to this terminal and press Enter to continue.")

    #--read in the new deconvolution file to get xspec normalization--
    if not os.path.isfile(pwd+'parameters.txt'):
        os.link(runpath+'/parameters.txt',pwd+'parameters.txt')
    dfnew = merge_output(runpath=pwd,save=False)
    dfnew = dfnew[dfnew.iteration==int(suffix)]
    #newnorm = dfnew.blob_norm.sum()
    newnorm = sum([dfnew[col].sum() for col in dfnew.columns if 'norm' in col])
    xspecscale = oldnorm/newnorm
#    xspecscale = newnorm/oldnorm
    print 'newnorm, oldnorm, scale = ',newnorm,oldnorm, xspecscale
    print 'niters = ',niters
    
    #--Create Histogram (Read in spectrum file)--
    blobspecs = read_spectra(runpath=pwd,smin=int(suffix),
                        smax=int(suffix),logbins=logbins,
                             datarange=datarange,bins=bins,oversim=oversim)

    bspec = blobspecs[0] # tuple

    y = bspec[0]/oversim#*xspecscale/float(niters)
    #y = bspec[0]*xspecscale/oversim#/float(niters)
    #yerrors = bspec[1]/float(niters)


    #--Modify norms by Absorption--

    #-get absorption curve as function of energy-
    crosses = np.vectorize(abs_cross_section)
    sigmas = crosses(bspec[2]) # one for each spectral bin
    absorption = np.exp(-1.0*sigmas*nH)
    print 'min,max y = ',min(y),max(y)
    y = y*absorption[:-1]
    print 'min,max y = ',min(y),max(y)
    print 'min,max cross sections = ',min(sigmas),max(sigmas)
    print 'min,max absorption = ',min(absorption),max(absorption)
    
    return (y,bspec[1],bspec[2])
