"""
Module of functions needed to create maps of xmc blob parameters

Includes:

 make_map
 calculate_map
 gaussian_integral
 point_integral
 iteration_image
 collapse stack
 mask_circle

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
             paramy='blob_psi',paramsize='blob_sigma',exclude_region=None,
             iteration_type='median',clobber=False,nlayers=None):
    """
    Author: Kari A. Frank
    Date: November 19, 2015
    Purpose: Read a file containing blob parameters and create an 
             associated map.

    Input:

      indata (string or DataFrame):   name of file containing the blob 
                        parameters (string), formatted as output from 
                        the python function
                        xmcinter.xmcfiles.merge_output() OR
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

      nlayers: optionally set number of layers (number of iterations to use
               rather than itmod.  if set, will override itmod. default=None

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

      exclude_region (3D float tuple): tuple of form (x0,y0,radius)
                specifying a circular region to mask (set image
                values to zero). x0,y0 are the center coordinates
                of the circle. all units are in the same
                units as paramx and paramy.

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
                        x0=x0,y0=y0,shape=paramshape,nlayers=nlayers)

#    print "max, min img = ",np.max(img), np.min(img)

    #----Mask Region----
    if exclude_region is not None:
        msk = circle_mask(df,paramx,paramy,exclude_region,binsize,
                          imagesize,x0,y0)
        img = img*msk

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
    """Function to calculate gaussian integral."""
    integ = 0.0
    (steps, dx) = np.linspace(lowerx,upperx,nsteps,retstep=True)
    for xp in steps:
        integ = integ + np.exp(-1.0/2.0*(xp-x)**2.0/sigma**2.0)*dx

    return integ

#--------------------------------------------------------------------------
def point_integral(lowerx,upperx,lowery,uppery,x,y):
    """Function to determine if blob center lies in pixel."""
    
    print "ERROR: point_integral is not yet functional."
    f = x.to_frame('x')
    f['y'] = y
    f.ycheck = np.where((lowery < f['y'] & f['y'] < uppery),1.0,0.0)
    f.xcheck = np.where(lowerx < f['x'],1.0,0.0)
    f['integ'] = f.ycheck*f.xcheck

#    for i in xrange(len(x)):
#        if (lowerx < x[i] < upperx) and (lowery < y[i] < uppery):
#            integ[i] = 1.0
#        else: 
#            integ[i] = 0.0

#    print 'len integ = ',len(integ)
#    print 'len x = ',len(x)
    return f['integ'].values

#--------------------------------------------------------------------------
def circle_mask(df,paramx,paramy,exclude_region,binsize,imagesize,x0,y0):
    """Function to create mask (set to zero) of a circular region 
    in an image array. To apply mask, multiply the image array by mask."""

    # dummy parameter (all ones)
    dummy = pd.Series(np.ones_like(df[paramx])
                      ,index=df[paramx].index)
    dummy=dummy.to_frame('mask')
    print 'len dummy = ',len(dummy.index)
    print 'len df = ',len(df.index)

    dummy[paramx] = df[paramx].values
    dummy[paramy] = df[paramy].values

    # set mask within circle to zero
    dummy['mask'][( (dummy[paramx]-exclude_region[0])**2.0+(dummy['blob_psi']-exclude_region[1])**2.0 < exclude_region[2]**2.0 )] = 0.0

    # -since shape='points', paramsize (2nd instance of df[paramx]) is
    #  ignored, so just need a dummy column
    mask = calculate_map(dummy['mask'],dummy[paramx],dummy[paramy],
                         dummy[paramx],
                         blobiterations=df['iteration'],
                         binsize=binsize,iteration_type='max',
                         imagesize=imagesize,nlayers=1,x0=x0,y0=y0,
                         shape='points') 
                
    return mask
    


#--------------------------------------------------------------------------
def iteration_image(data,nbins_x,nbins_y,binsize,xmin,ymin,n_int_steps,
                    iteration_type,shape):
    """Function to combine blobs from single iteration into 1 image."""
    from wrangle import weighted_median

    #--initialize 2D image array--
    iterimage = np.zeros((nbins_x,nbins_y))

    #--loop over image--
    for x in xrange(nbins_x):
            #get x integral
            lowerx = int(xmin + x*binsize)
            upperx = int(xmin + x*binsize + binsize)
#            dx = (upperx - lowerx)/float(n_int_steps)
            if shape == 'gauss' or shape == 'points':
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
                if shape == 'gauss' or shape == 'points':
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
                    # for now, points is implemented by setting the volumes 
                    #   to be much smaller than a pixel size
                    fractions = (x_blob_integrals*y_blob_integrals*
                                 (2.0*np.pi*data['size']**2.0)**.5 / 
                                 data['volume'])
#                    print "points is not yet implemented"
                    # if assuming points, then fraction=1 or 0
#                    fractions = point_integral(lowerx,upperx,lowery,uppery,
#                                               data['x'],data['y'])

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
                elif iteration_type == 'max':
                    iterimage[x,y]=np.max(data['param']*data['weight']
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
                  x0=None,y0=None,shape='gauss',nlayers=None):
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
                      and blobsize (assumes square pixel, default is 10.)
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
            nlayers : optionally can set number of layers instead of itmod.
                    if set, will override itmod. default=None
            x0,y0 : cluster position, in same coordinate system as blobx 
                    and blob (typically xmc coordinates, default is
                    x0=median(blobx),y0=median(bloby)).
            shape : shape of the blobs, 'gauss' (default),'sphere', 
                    or 'points'
            n_int_steps : number of steps in each integration loop 
                          (default = 50)

    Output:

         Returns a 2D array containing the map.

    Usage Notes

    Example:


    """
    #----Import Modules----
    from wrangle import filterblobs

    #----Set Defaults----
    types = ['median','average','total','error','max']
    if ctype not in types:
        print "Warning: Unrecognized ctype. Using ctype='median'"
        ctype = 'median'
    if iteration_type not in types:
        print ("Warning: Unrecognized iteration_type. "
               "Using iteration_type='median'")
        iteration_type = 'median'
    if blobiterations is None:
        blobiterations = np.zeros_like(blobparam)
    if (shape != 'gauss') and (shape != 'sphere' ) and (shape != 'points'):
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
    niter = np.unique(blobiterations).size
    if nlayers is None:
        nlayers = niter/itmod
        if nlayers == 0: 
            nlayers = 1
    else:
        if nlayers > niter: #max nlayers = number iterations
            nlayers = niter
        itmod = niter/nlayers
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
    if shape == 'points':
        df['volume'] = (0.1*binsize)**3.0 # set to much smaller than pixel

    #----Group by Iteration----
    layers = df.groupby('iteration')

    #----Iterate over groups (i.e. iterations)----
    layer = 0
    for i, group in layers: 
        print 'layer = ',layer
#        print 'min, max group ',np.min(group),np.max(group)
        #i=iteration number (string or int?), group = subset of dataframe
        image_stack[:,:,layer] = iteration_image(group,nbins_x,nbins_y,
                                                 binsize,xmin,ymin,
                                                 n_int_steps,
                                                 iteration_type,shape)
        layer = layer + 1

    #--Collapse Image Stack (combine iterations)----
    themap = collapse_stack(image_stack,nbins_x,nbins_y,ctype)
        
    #----Return map----
    return themap



