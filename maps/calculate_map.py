"""
Author: Kari A. Frank
Date: October 29, 2015
Purpose: Calculate a map (2D array) from the specified columns in a dataframe.

Usage:
     calculate_map(blobparam,blobx,bloby,blobsize,blobiterations=None,blobweights=None,binsize=10,
                  iteration_type='median',ctype='median',imagesize=None,itmod=100,n_int_steps=50,
                  verbose=0,x0=None,y0=None,shape='gauss')


Input:

        Note: Every variable that begins with 'blob' should be a 1D array containing
              one element for each blob.  All such arrays must be the same length.
   
        blobparam: vector of blob parameters from which to create a map
        blobx: vector of blob x positions
        bloby: vector of blob y positions
        blobsize: vector of blob sizes (in arcsec)
        binsize: width of pixels, in same units as blobx,bloby, 
                 and blobsize (assumes square pixel, default is 60.)
        blobiterations: array of the iteration associated with each
                blob in blob array. if not provided, assumes all
                blobs are from the same iteration.
        iteration_type: type=median, average, or total.  Determines how to 
                combine the blobs in each iteration to create the
                iteration image. (default is  median, but note that it 
                should be set to 'total' if making emission measure map)
        ctype: ctype=median, average, total, or error to specify how to
                combine the different iteration images when making the 
                final image. error will produce a map of the blobeter
                error (default=median)
        blobweights: array of weights for each blob (e.g. emission measures)
        imagesize: length of one side of image, assumes square image 
              ( default=max(blobx)-min(blobx) ).
        itmod: set to >1 to use only every itmod'th
                iteration (default = 10)
        x0,y0: cluster position, in same coordinate system as blobx and bloby
                (typically xmc coordinates, default is
                x0=median(blobx),y0=median(bloby)).
        shape: shape of the blobs, 'gauss' (default) or 'sphere'
        n_int_steps: number of steps in each integration loop (default = 50)
        verbose: switch to turn on printing of checkpoints (default = 0)

Output:

     Returns a 2D array containing the map.

Usage Notes
     
Example:


"""
#----Import Global Modules----
import pandas as pd
import numpy as np


#----Main Map Function----
def calculate_map(blobparam,blobx,bloby,blobsize,blobiterations=None,blobweights=None,binsize=10,
                  iteration_type='median',ctype='median',imagesize=None,itmod=100,n_int_steps=50,
                  verbose=0,x0=None,y0=None,shape='gauss'):

    #----Import Modules----
    from ..analysis.filter import filter

    #----Set Defaults----
    types = ['median','average','total','error']
    if ctype not in types:
        print "Warning: Unrecognized ctype. Using ctype='median'"
        ctype = 'median'
    if iteration_type not in types:
        print "Warning: Unrecognized iteration_type. Using iteration_type='median'"
        iteration_type = 'median'
    if blobiterations == None:
        blobiterations = np.zeros_like(blobparam)
    if (shape != 'gauss') and (shape != 'sphere' ):
        print "Warning: Unrecognized blob shape. Using shape='gauss'"
        shape = 'gauss'
    if imagesize == None:
        imagesize = 1.2*(max(blobx) - min(blobx))
    if x0 == None:
        x0 = np.median(blobx)
    if y0 == None:
        y0 = np.median(bloby)
    if blobweights == None: # default is equally weight all blobs
        blobweights = np.ones_like(blobparam)

    #----Set Up Map Parameters----
    imageradius = imagesize/2
    xmin = x0 - imageradius
    xmax = x0 + imageradius
    ymin = y0 - imageradius
    ymax = y0 + imageradius
    imageradius = (xmax-xmin)/2.0
    imagesize = imageradius*2.0

    #-number of map layers (one per iteration) and number of pixels-
    nlayers = np.unique(blobiterations).size/itmod
    nbins_x = np.floor((xmax - xmin)/binsize)
    nbins_y = np.floor((ymax - ymin)/binsize)

    #----Concatenate into a single dataframe and filter out iterations----
    df = pd.concat([blobparam,blobx,bloby,blobsize,blobweights,blobiterations])
    df.columns = ['param','x','y','size','weight','iteration']
    
    #--Remove iterations according to itmod--

    #-make list of iterations to use-
    its = []
    for i in range(nlayers):
        its = its + [i*itmod]
    itstr = ['iteration']*len(its)

    #-remove non-matching iterations-
    df = filter(df,istr,minvals=its,maxvals=its)
       
    #----Calculate Blob Volumes----
    if shape == 'gauss':
        df['volume'] = (2.0*np.pi*df['size']**2.0)**1.5
    if shape == 'sphere':
        df['volume'] = (4.0/3.0)*np.pi*df['size']**3.0

    #----Initialize Image Stack----
    image_stack = np.zeros((nbins_x,nbins_y,nlayers))


    #----Group by Iteration----
    layers = df.groupby('iteration')

    #----Iterate over groups----
    for i, group in layers: #i=iteration number (string or int?), group = subset of dataframe
        #--calculate fraction of each blob in each pixel and store as 2D array--
        #image_stack[:][:][layer]=gauss_fraction_image(nbinsx,nbinsy,binsize,n_int_steps,group['x'],group['y'],group['size'],iteration_type)
        #-loop over x-
        for x in range(nbins_x):
            #get x integral
            lowerx = xmin + x*binsize
            upperx = xmin + x*binsize + binsize
            dx = (upperx - lowerx)/float(n_int_steps)
            x_blob_integrals = gaussian_integral(lowerx,upperx,dx,group['x'],group['size'])
            for y in range(nbins_y):
                #get y integral
                lowery = ymin + y*binsize
                uppery = ymin + y*binsize + binsize
                dy = (uppery - lowery)/float(n_int_steps)
                y_blob_integrals = gaussian_integral(lowery,uppery,dy,group['y'],group['size'])
                fractions = x_blob_integrals*y_blob_integrals#times dz integral
                                
    #--Collapse Image Stack (combine iterations)----

        
    #----Return map----
    return themap

