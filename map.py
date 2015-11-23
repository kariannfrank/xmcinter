"""
Module of functions needed to create maps of xmc blob parameters

Includes:

 make_map
 calculate_map
 gaussian_integral
 iteration_image
 collapse stack

The main function that should be called is make_map.  The others are 
essentially just helper functions for make_map. The function that does most
of the work is calculate_map.

"""
#----------------------------------------------------------------------------

#----Import Global Modules----
import os
import pandas as pd
import numpy as np
import astropy.io.fits as fits

#----------------------------------------------------------------------------
def make_map(indata,outfile=None,paramname='blob_kT',paramweights=None,
             binsize=10.0,itmod=100,paramshape='gauss',ctype='median',
             x0=None,y0=None,imagesize=None,sigthresh=0.0,paramx='blob_phi',
             paramy='blob_psi',paramsize='blob_sigma',
             iteration_type='median',clobber=False):
    """
    Author: Kari A. Frank
    Date: November 19, 2015
    Purpose: Read a file containing blob parameters and create an 
             associated map.

    Input:

      indata (string or DataFrame):   name of file containing the blob 
                        parameters (string), formatted as output from 
                        the python function
                        xmcinter.files.xmcrun.merge_output() OR
                        a pandas dataframe of blob parameters in same format.

      outfile (string):  name of the output fits file

      paramname (string) : name of the column (parameter name) to be 
                           mapped (default='blob_kT')

      paramx,paramy (strings) : names of the column of the x and y
                                blob positions (default='blob_phi',
                                'blob_psi')

      paramshape (string) : shape of the blobs, 'gauss','sphere', or 
                            'points'. 'points' simply assumes every blob
                            is a point, even if the model was not 
                            (default='gauss')

      paramsize : string name of the column containing the size of 
                 the blobs (default='blob_sigma' for shape='gauss'
                 or 'blob_radius' for shape='sphere')

      paramweights : string name of the column containing the weights
                     (default='', no weighting)

      itmod :   set to an integer > 1 to use only every
               it_mod'th iteration (defaul=100)

      binsize : size of a pixel, in same units as paramx and paramy,
               default=60 (fast)

      iteration_type (string) : 'median', 'average', or 'total'.  Determines
                how to combine the blobs in each iteration to create the
                iteration image. (default is  'median', but note that it 
                should be set to 'total' if making emission measure map)

      ctype (string) : 'median', 'average', 'total', or 'error' to specify 
                how to combine the different iteration images when making 
                the final image. error will produce a map of the parameter
                error (default='median')

      imagesize (float) : optionally specify the size of output image 
                          (length of one side of square image) in same u
                          nits as paramx,y

      x0,y0 (floats):  center coordinates of the desired output image, in
                       same units as paramx,paramy

      sigthresh (float): optionally specify a significance threshold (in
                number of sigma) for the final map.  all pixels
                that do not meet this significance threshold are set
                to zero.  note this takes much longer, as it
                requires calculating an error map in addition to the desired
                output map
               
       clobber (bool) : specify whether any existing fits file of the same
                        name as outfile should be overwritten. 

    Output:

         Saves a fits file in the same directory as infile, containing the
         calculated map.

    Usage Notes:

     The implementation for shape='sphere' is not yet functional.
    """
    
    #----Import Modules----
    import time

    #----Store blob information in DataFrame and set output file----
    if type(indata) is str:
        df = pd.read_table(indata,sep='\t',index_col = 0)
        if outfile is None:
            (fname,ext) = os.path.splitext(indata)
            outfile = fname + '_'+paramname+'_'+ctype+'.fits'
        indatastr = indata
    else:
        df = indata
        if outfile is None:
            outfile = paramname+'_'+ctype+'.fits'
        indatastr = 'DataFrame'

    #----Check for weights----
    if paramweights is None:
        weights = None
    else:
        weights = df[paramweights]

    #----Calculate the map----
    img = calculate_map(df[paramname],df[paramx],df[paramy],df[paramsize],
                        blobiterations=df['iteration'],
                        blobweights=weights,binsize=binsize,
                        iteration_type=iteration_type,ctype=ctype,
                        imagesize=imagesize,itmod=itmod,
                        x0=x0,y0=y0,shape=paramshape)

    print "max, min img = ",np.max(img), np.min(img)

    #----Save map to fits file----

    #--write history--
    history1 = ('make_map,'+indatastr+',outfile='+nstr(outfile)
                +',paramname='+nstr(paramname)
                +',paramweights='+nstr(paramweights)
                +',paramx='+nstr(paramx)+',paramy='+nstr(paramy)
                +',paramsize='+nstr(paramsize)+',binsize='
                +nstr(binsize)+',itmod='+nstr(itmod)+',paramshape='
                +nstr(paramshape)+',ctype='+nstr(ctype)+',iteration_type='
                +nstr(iteration_type)+',x0='+nstr(x0)+',y0='+nstr(y0)
                +',imagesize='+nstr(imagesize)+',sigthresh='
                +nstr(sigthresh))
    
    history2 = 'Created '+str(time.strftime("%x %H:%M:%S"))

    #--write file--
    hdr = fits.Header()
    hdr['HISTORY']=history2
    hdr['HISTORY']=history1
    hdu = fits.PrimaryHDU(img,header=hdr)
    
    hdu.writeto(outfile,clobber=clobber)
    
    return img

#--------------------------------------------------------------------------
def nstr(x):
    """Wrapper for built-in str() that returns the string 'None' if 
       passed None"""
    if x is None:
        return 'None'
    else:
        return str(x)

#--------------------------------------------------------------------------
def gaussian_integral(lowerx,upperx,nsteps,x,sigma):
    """Function to calculate gaussoam integral."""
    integ = 0.0
    (steps, dx) = np.linspace(lowerx,upperx,nsteps,retstep=True)
    for xp in steps:
        integ = integ + np.exp(-1.0/2.0*(xp-x)**2.0/sigma**2.0)*dx

    return integ

#--------------------------------------------------------------------------
def iteration_image(data,nbins_x,nbins_y,binsize,xmin,ymin,n_int_steps,
                    iteration_type,shape):
    """Function to combine blobs from single iteration into 1 image."""
    from ..analysis.wrangle import weighted_median

    #--initialize 2D image array--
    iterimage = np.zeros((nbins_x,nbins_y))

    #--loop over image--
    for x in xrange(nbins_x):
            #get x integral
            lowerx = int(xmin + x*binsize)
            upperx = int(xmin + x*binsize + binsize)
#            dx = (upperx - lowerx)/float(n_int_steps)
            if shape == 'gauss':
                x_blob_integrals = gaussian_integral(lowerx,upperx,
                                                     n_int_steps,
                                                     data['x'],data['size'])
            elif shape == 'sphere':
                x_blob_integrals = spherical_integral(lowerx,upperx,\
                                                      n_int_steps,\
                                                     data['x'],data['size'])
            for y in xrange(nbins_y):
                #get y integral
                lowery = int(ymin + y*binsize)
                uppery = int(ymin + y*binsize + binsize)
#                dy = int((uppery - lowery)/float(n_int_steps))
                if shape == 'gauss':
                    y_blob_integrals = gaussian_integral(lowery,uppery, \
                                                     n_int_steps,\
                                                     data['y'],data['size'])
                elif shape == 'sphere':
                    y_blob_integrals = spherical_integral(lowery,uppery,\
                                                     n_int_steps,\
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
                    iterimage[x,y]=weighted_median(data['param'],
                                              weights=data['weight']
                                              *fractions)
                elif iteration_type == 'average':
                    iterimage[x,y]=np.average(data['param'],
                                         weights=data['weight']*fractions)
                elif iteration_type == 'total':
                    iterimage[x,y]=np.sum(data['param']*data['weight']
                                          *fractions)
                else:
                    print "ERROR: unrecognized iteration_type"

    return iterimage

#--------------------------------------------------------------------------
def collapse_stack(img_stack,nbins_x,nbins_y,ctype):
    """Function to collapse a stack of images into a single image."""

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

#-------------------------------------------------------------------------
def calculate_map(blobparam,blobx,bloby,blobsize,blobiterations=None,
                  blobweights=None,binsize=10,iteration_type='median',
                  ctype='median',imagesize=None,itmod=100,n_int_steps=50,
                  x0=None,y0=None,shape='gauss'):
    """
    The main mapping function.

    Author: Kari A. Frank
    Date: October 29, 2015
    Purpose: Calculate a map (2D array) from the provided columns a 
             dataframe.

    Usage:
         calculate_map(blobparam,blobx,bloby,blobsize,blobiterations=None,
                      blobweights=None,binsize=10,iteration_type='median',
                      ctype='median',imagesize=None,itmod=100,
                      n_int_steps=50,verbose=0,x0=None,y0=None,
                      shape='gauss')

    Input:

            Note: Every variable that begins with 'blob' should be a 1D 
                  array containing one element for each blob.  All such 
                  arrays must be the same length.

            blobparam : vector of blob parameters from which to create a map
            blobx : vector of blob x positions
            bloby : vector of blob y positions
            blobsize : vector of blob sizes (in arcsec)
            binsize : width of pixels, in same units as blobx,bloby, 
                      and blobsize (assumes square pixel, default is 60.)
            blobiterations : array of the iteration associated with each
                             blob in blob array. if not provided, assumes 
                             all blobs are from the same iteration.
            iteration_type : type=median, average, or total.  Determines 
                    how to combine the blobs in each iteration to create the
                    iteration image. (default is  median, but note that it 
                    should be set to 'total' if making emission measure map)
            ctype : ctype=median, average, total, or error to specify how to
                    combine the different iteration images when making the 
                    final image. error will produce a map of the blobeter
                    error (default=median)
            blobweights : array of weights for each blob (e.g. emission 
                          measures)
            imagesize : length of one side of image, assumes square image 
                        ( default=max(blobx)-min(blobx) ).
            itmod : set to >1 to use only every itmod'th
                    iteration (default = 10)
            x0,y0 : cluster position, in same coordinate system as blobx 
                    and blob (typically xmc coordinates, default is
                    x0=median(blobx),y0=median(bloby)).
            shape : shape of the blobs, 'gauss' (default) or 'sphere'
            n_int_steps : number of steps in each integration loop 
                          (default = 50)

    Output:

         Returns a 2D array containing the map.

    Usage Notes

    Example:


    """
    #----Import Modules----
    from ..analysis.wrangle import filterblobs

    #----Set Defaults----
    types = ['median','average','total','error']
    if ctype not in types:
        print "Warning: Unrecognized ctype. Using ctype='median'"
        ctype = 'median'
    if iteration_type not in types:
        print ("Warning: Unrecognized iteration_type. "
               "Using iteration_type='median'")
        iteration_type = 'median'
    if blobiterations is None:
        blobiterations = np.zeros_like(blobparam)
    if (shape != 'gauss') and (shape != 'sphere' ):
        print "Warning: Unrecognized blob shape. Using shape='gauss'"
        shape = 'gauss'
    if imagesize is None:
        imagesize = 1.2*(max(blobx) - min(blobx))
    if x0 is None:
        x0 = np.median(blobx)
    if y0 is None:
        y0 = np.median(bloby)
    if blobweights is None: # default is equally weight all blobs
        blobweights = pd.Series(np.ones_like(blobparam)
                                ,index=blobparam.index)

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
    if nlayers == 0: 
        nlayers = 1
    nbins_x = int(np.floor((xmax - xmin)/binsize))
    nbins_y = int(np.floor((ymax - ymin)/binsize))
    print 'nbins_x, nbins_y = ',nbins_x,nbins_y

    #-initialize image stack-
    image_stack = np.zeros((nbins_x,nbins_y,nlayers))

    #----Concatenate into a single dataframe and filter out iterations----
    df = blobparam.to_frame(name='param')
    df['x'] = blobx
    df['y'] = bloby
    df['size'] = blobsize
    df['weight'] = blobweights
    df['iteration'] = blobiterations

    #--Remove iterations according to itmod--

    #-make list of iterations to use-
    its = []
    itmin = np.min(df['iteration'])
    for i in range(nlayers):
        its = its + [i*itmod + itmin]
    itstr = ['iteration']*len(its)

    #-keep only matching iterations-
    df = filterblobs(df,itstr,minvals=its,maxvals=its,logic='or')
   
    #----Calculate Blob Volumes----
    if shape == 'gauss':
        df['volume'] = (2.0*np.pi*np.square(df['size']))**1.5
    if shape == 'sphere':
        df['volume'] = (4.0/3.0)*np.pi*df['size']**3.0

    #----Group by Iteration----
    layers = df.groupby('iteration')

    #----Iterate over groups (i.e. iterations)----
    layer = 0
    for i, group in layers: 
        print 'layer = ',layer
        print 'min, max group ',np.min(group),np.max(group)
        #i=iteration number (string or int?), group = subset of dataframe
        image_stack[:,:,layer] = iteration_image(group,nbins_x,nbins_y,
                                                 binsize,xmin,ymin,
                                                 n_int_steps,
                                                 iteration_type,shape)
        print 'min, max layer = ',np.min(image_stack[:,:,layer]),np.max(image_stack[:,:,layer])
        layer = layer + 1

    #--Collapse Image Stack (combine iterations)----
    themap = collapse_stack(image_stack,nbins_x,nbins_y,ctype)
        
    #----Return map----
    return themap



