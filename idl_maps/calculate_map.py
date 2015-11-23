"""
Author: Kari A. Frank
Date: October 29, 2015
Purpose: Calculate a map (2D array) from the specified columns in a dataframe.

Usage:
     calculate_map(blobparam,blobx,bloby,blobsize,blobiterations=None,
                  blobweights=None,binsize=10,iteration_type='median',
                  ctype='median',imagesize=None,itmod=100,n_int_steps=50,
                  verbose=0,x0=None,y0=None,shape='gauss')

Input:

        Note: Every variable that begins with 'blob' should be a 1D array 
              containing one element for each blob.  All such arrays must 
              be the same length.
   
        blobparam : vector of blob parameters from which to create a map
        blobx : vector of blob x positions
        bloby : vector of blob y positions
        blobsize : vector of blob sizes (in arcsec)
        binsize : width of pixels, in same units as blobx,bloby, 
                  and blobsize (assumes square pixel, default is 60.)
        blobiterations : array of the iteration associated with each
                         blob in blob array. if not provided, assumes all
                         blobs are from the same iteration.
        iteration_type : type=median, average, or total.  Determines how to 
                combine the blobs in each iteration to create the
                iteration image. (default is  median, but note that it 
                should be set to 'total' if making emission measure map)
        ctype : ctype=median, average, total, or error to specify how to
                combine the different iteration images when making the 
                final image. error will produce a map of the blobeter
                error (default=median)
        blobweights : array of weights for each blob (e.g. emission measures)
        imagesize : length of one side of image, assumes square image 
              ( default=max(blobx)-min(blobx) ).
        itmod : set to >1 to use only every itmod'th
                iteration (default = 10)
        x0,y0 : cluster position, in same coordinate system as blobx and bloby
                (typically xmc coordinates, default is
                x0=median(blobx),y0=median(bloby)).
        shape : shape of the blobs, 'gauss' (default) or 'sphere'
        n_int_steps : number of steps in each integration loop (default = 50)
        verbose: switch to turn on printing of checkpoints (default = 0)

Output:

     Returns a 2D array containing the map.

Usage Notes
     
Example:


"""
#----Import Global Modules----
import pandas as pd
import numpy as np
from ..analysis.wrangle import filterblobs,weighted_median

#----Functions to calculate integrals----

def gaussian_integral(lowerx,upperx,dx,x,sigma):

    for xp in xrange(lowerx,upperx,dx):
        integ = integ + np.exp(-1.0/2.0*(xp-x)**2.0/sigma**2.0)*dx

    return integ

#----Function to combine blobs from single iteration into 1 image----
def iteration_image(data,nbins_x,nbins_y,binsize,n_int_steps,iteration_type,
                    shape):

    for x in xrange(nbins_x):
            #get x integral
            lowerx = xmin + x*binsize
            upperx = xmin + x*binsize + binsize
            dx = (upperx - lowerx)/float(n_int_steps)
            if shape == 'gauss':
                x_blob_integrals = gaussian_integral(lowerx,upperx,dx,data['x'],
                                                     data['size'])
            elif shape == 'sphere':
                x_blob_integrals = spherical_integral(lowerx,upperx,dx,
                                                      data['x'],data['size'])
            for y in xrange(nbins_y):
                #get y integral
                lowery = ymin + y*binsize
                uppery = ymin + y*binsize + binsize
                dy = (uppery - lowery)/float(n_int_steps)
                if shape == 'gauss':
                    y_blob_integrals = gaussian_integral(lowery,uppery,dy,
                                                         data['y'],data['size'])
                elif shape == 'sphere':
                    y_blob_integrals = spherical_integral(lowery,uppery,dy, \
                                                         data['y'],data['size'])
                #calculate fraction of blob volume in this pixel

                if shape != 'points':
                # !! for now this assumes gaussian volume !!
                    fractions = (x_blob_integrals*y_blob_integrals*
                                 (2.0*np.pi*data['size']**2.0)**.5 / 
                                 data['volume'])
                #times dz integral to get total volume in pixel, 
                #then divided by total volume
                else:
                    print "points is not yet implemented"
                    # if assuming points, then fraction=1 or 0
                    # !! in progress !!

                #-combine blobs in this pixel-
                if iteration_type == 'median':
                    iterimage=weighted_median(data['param'],
                                              weights=data['weight']*fractions)
                elif iteration_type == 'average':
                    iterimage=np.average(data['param'],
                                         weights=data['weight']*fractions)
                elif iteration_type == 'total':
                    iterimage=np.sum(data['param']*data['weight']*fractions)
                else:
                    print "ERROR: unrecognized iteration_type"

    return iterimage

#----Function to collapse a stack of images into a single image----
def collapse_stack(img_stack,nbins_x,nbins_y,ctype):

    if ctype == 'average':
        collapsed_img = np.average(img_stack,axis=2)
    elif ctype == 'median':
        collapsed_img = np.median(img_stack,axis=2)
    elif ctype == 'total':
        collapsed_img = np.sum(img_stack,axis=2)
    elif ctype == 'error':
        collapsed_img = np.std(img_stack,axis=2)
    else: 
        print "ERROR: unrecognized iteration_type"

    return collapsed_img

#----Main Map Function----
def calculate_map(blobparam,blobx,bloby,blobsize,blobiterations=None,
                  blobweights=None,binsize=10,iteration_type='median',
                  ctype='median',imagesize=None,itmod=100,n_int_steps=50,
                  x0=None,y0=None,shape='gauss'):

    #----Import Modules----

    #----Set Defaults----
    types = ['median','average','total','error']
    if ctype not in types:
        print "Warning: Unrecognized ctype. Using ctype='median'"
        ctype = 'median'
    if iteration_type not in types:
        print ("Warning: Unrecognized iteration_type. "
               "Using iteration_type='median'")
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

    #-initialize image stack-
    image_stack = np.zeros((nbins_x,nbins_y,nlayers))

    #----Concatenate into a single dataframe and filter out iterations----
    df = pd.concat([blobparam,blobx,bloby,blobsize,blobweights,blobiterations])
    df.columns = ['param','x','y','size','weight','iteration']
    
    #--Remove iterations according to itmod--

    #-make list of iterations to use-
    its = []
    for i in range(nlayers):
        its = its + [i*itmod]
    itstr = ['iteration']*len(its)

    #-keep only matching iterations-
    df = filterblobs(df,itstr,minvals=its,maxvals=its)
       
    #----Calculate Blob Volumes----
    if shape == 'gauss':
        df['volume'] = (2.0*np.pi*df['size']**2.0)**1.5
    if shape == 'sphere':
        df['volume'] = (4.0/3.0)*np.pi*df['size']**3.0

    #----Group by Iteration----
    layers = df.groupby('iteration')

    #----Iterate over groups (i.e. iterations)----
    for i, group in layers: 
    #i=iteration number (string or int?), group = subset of dataframe
        image_stack[:,:,layer] = iteration_image(group,nbins_x,nbins_y,binsize,
                                                 n_int_steps,iteration_type,
                                                 shape)

    #--Collapse Image Stack (combine iterations)----
    themap = collapse_stack(image_stack,nbins_x,nbins_y,ctype)
        
    #----Return map----
    return themap

