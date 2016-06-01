"""
Module of functions needed to create maps of xmc blob parameters

Includes:

 make_map
 gaussian_integral
 point_integral
 iteration_image
 collapse stack
 mask_circle

The main function that should be called is make_map.  The others are 
essentially just helper functions for make_map. The function that does most
of the work is iteration_image.

IMPORTANT NOTE ABOUT DENSITY MAPS
make_map() does not properly handle density maps (i.e. number density or
mass density) since it is not linearly dependent on EM.  Intead, use
em_to_map_density() to calculate the correct parameter columns to pass to 
make_map(), with the weights set to 'densityspecial'.

"""
#----------------------------------------------------------------------------

#----Import Global Modules----
import os,inspect
import pandas as pd
import numpy as np
import astropy.io.fits as fits
from scipy.integrate import quad,nquad
from scipy.ndimage.interpolation import rotate
from multiprocessing import Pool
import ctypes

#----------------------------------------------------------------------------
def make_map(indata,outfile=None,paramname='blob_kT',paramweights=None,
             binsize=10.0,itmod=100,paramshape='gauss',ctype='median',
             x0=None,y0=None,imagesize=None,witherror=True,sigthresh=0.0,
             sigthreshparam=None,imgstdthresh=None,imgstdthreshparam=None,
             paramx='blob_phi',
             paramy='blob_psi',paramsize='blob_sigma',exclude_region=None,
             iteration_type='median',clobber=False,nlayers=None,
             parallel=True,nproc=3,cint=True,movie=False,moviedir=None,
             cumulativemovie=False,withsignificance=False,rotation=0.0):
    """
    Author: Kari A. Frank
    Date: November 19, 2015
    Purpose: Read a file containing blob parameters and create an 
             associated map.

    Input:

      indata (string or DataFrame):   name of file containing the blob 
                      parameters (string), formatted as output from 
                      the python function xmcinter.xmcfiles.merge_output() 
                      OR a pandas dataframe of blob parameters in same 
                      format.

      outfile (string):  name of the output fits file

      paramname (string or list of str) : name of the column (parameter 
                          name) to be mapped (default='blob_kT')

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
                (default='', no weighting). if paramname is a list, 
                then can also give paramweights as a list of the same
                length, specifying a different weights column for each
                map -- this is important if, e.g. one of the paramnames 
                is 'blob_em', which should is typically also used as the 
                weights for the other parameters.  passing a value of None
                will result in an unweighted map.

      itmod :   set to an integer > 1 to use only every
               it_mod'th iteration (defaul=100)

      nlayers: optionally set number of layers (number of iterations to use
               rather than itmod.  if set, will override itmod. default=None

      binsize : size of a pixel, in same units as paramx and paramy,
               default=60 (fast)

      iteration_type (string) : 'median', 'average', or 'total'.  Determines
                how to combine the blobs within each iteration to create the
                iteration image. (default is  'median', but note that it 
                should be set to 'total' if making emission measure map)
                if paramname is a list, then can also give iteration_type 
                as a list of the same length, specifying a different 
                iteration_type for each map -- this is important if, e.g.
                one of the paramnames is 'blob_em', which should generally
                use iteration_type='total'.

      ctype (string) : 'median', 'average', 'total', or 'error' to specify 
                how to combine the different iteration images when making 
                the final image. error will produce a map of the parameter
                error (default='median'). if paramname is a list, then can
                also give ctype as a list of the same length, specifying
                a differenty ctype for each map.

      witherror (bool) : switch to also return a map of the error in 
                         each pixel (standard deviation)

      withsignificance (bool) : switch to also return a map of the 
                significance (in #sigma, i.e. img/errimg) in 
                each pixel. If True, then will set witherror=True, 
                regardless of whether the witherror argument was explicitly
                set.

      imagesize (float) : optionally specify the size of output image 
                          (length of one side of square image) in same
                          units as paramx,y. if paramnames is list, then 
                          all maps will have the same image size as the
                          first one in the list.

      x0,y0 (floats):  center coordinates of the desired output image, in
                       same units as paramx,paramy. if paramname is a list,
                       then all maps will have the same x0,y0 coordinates
                       as the first map.

      exclude_region (3D float tuple): tuple of form (x0,y0,radius)
                specifying a circular region to mask (set image
                values to zero). x0,y0 are the center coordinates
                of the circle. all units are in the same
                units as paramx and paramy.

      sigthresh (float): optionally specify a significance threshold (in
                number of sigma) for the final map.  all pixels
                that do not meet this significance threshold are set
                to nan.
               
      sigthreshparam (string): specify which parameter should be used to 
                calculate the significance for the significance threshold.
                Ignored if sigthresh=0.0.  Most commonly, 
                sigthreshparam=None (default) or sigthreshparam='blob_em'. 
                The latter will then only map the 
                regions (on a per pixel basis) for which the emission 
                measure significance was greater than sigthresh. If not 
                None, then sigthreshparam must be an element of paramname 
                list.

      imgstdthresh (float) : similar to sigthresh, except the threshold is
                set as number of standard deviations below the maximum
                value in the image, where the stdev is calculated from the 
                image array itself.  imgstdthresh is the number of standard
                deviations to subtract from the maximum to obtain the 
                minimum allowed pixel value.  If this minimum is < 0, then 
                it will be ignored.

      imgstdthreshparam (string) : same as sigthreshparam, but associated
                with the imgstdthresh argument. typically, this should be
                either None (do all parameters independently, default), or 
                'blob_em'. 

      clobber (bool) : specify whether any existing fits file of the same
                       name as outfile should be overwritten. 

      parallel (bool) : switch to specify if the iteration images
                    should be computed in serial or in parallel using
                    multiprocessing (default=True)

      nproc (int) :  if parallel=True, then nproc sets the number of 
              processors to use (default=3). ignored if parallel=False

      cint (bool) :  turn on/off the use of ctypes for integration 
              (default=True). set cint=False if gaussian.c is not 
              compiled on your machine.

      movie (bool) :  save each layer image individually in order to 
              create a movie from the images. Number of frames=nlayers. 
              (default=False).  If paramname is a list, then a movie will
              created for each parameter map, or can pass a list of bool 
              to movie specifying which maps in paramname should get
              an associated movie.

      moviedir (str) : optionally specify the folder in which to save
                   the frames for the movie. ignored if 
                   movie=False (default=outfile_base_parname_movie/).
                   Be careful - if moviedir already exists, any frames
                   in it will be overwritten!

      cumulativemovie (bool) : create the movie using cumulative images,
                i.e. recreate the image using all available iterations each
                time.  ignored if movie=False (default=False)

      rotation (numeric) : number of degrees to rotate the final images. If
               not a multiple of 90, then the output image size will be
               greater than the imagesize parameter (but with empty 
               corners), to avoid dropping any pixels.

    Output:

         Saves a fits file in the same directory as infile, containing the
         calculated map.

    Usage Notes:

     - The implementation for shape='sphere' is not yet functional.
     - If given multiple parameters to map, then all will mapped on the
       same x,y grid (imagesize, binsizes, and x0,y0 will be the same)
     - The image does not have to be square, but each pixel is always 
       square.
     - If the input dataframe has no columns 'iteration', then all blobs 
       will be assumed to come from a single iteration.
     - If both sigthresh and imgstdthresh are used, then sigthresh 
       will be applied first.

    """
    
    #----Import Modules----
    from wrangle import filterblobs,gaussian_volume
    import time

    #----Set any defaults----
    if withsignificance is True: witherror = True

    #----Check if lists----
    if not isinstance(paramname,list):
        paramname=[paramname]
    if not isinstance(paramweights,list):
        paramweights = [paramweights]*len(paramname)
    if not isinstance(iteration_type,list):
        iteration_type = [iteration_type]*len(paramname)
    if not isinstance(ctype,list):
        ctype = [ctype]*len(paramname)
    if not isinstance(movie,list):
        movie = [movie]*len(paramname)

    #----Verify inputs----
    types = ['median','average','total','error','max']
    for i in xrange(len(paramname)):
        if ctype[i] not in types:
            print "Warning: Unrecognized ctype. Using ctype='median'"
            ctype[i] = 'median'
        if iteration_type[i] not in types:
            print ("Warning: Unrecognized iteration_type. "
                   "Using iteration_type='median'")
            iteration_type[i] = 'median'

    if (paramshape != 'gauss') and (paramshape != 'sphere' ) and \
    (paramshape != 'points'):
        print "Warning: Unrecognized paramshape. Using paramshape='gauss'"
        paramshape = 'gauss'

    if (sigthreshparam is not None) and (sigthreshparam not in paramname):
        print ("Warning: "+sigthreshparam+" is not in paramname. "
               "Resetting sigthreshparam=None.")
        sigthreshparam=None

    if (imgstdthreshparam is not None) and \
    (imgstdthreshparam not in paramname):
        print ("Warning: "+imgstdthreshparam+" is not in paramname. "
               "Resetting imgstdthreshparam=None.")
        imgstdthreshparam=None

    #----Store blob information in DataFrame and set output file----
    if outfile is not None:
        outfile_base,ext = os.path.splitext(outfile)
        if outfile_base[-1] != '_': outfile_base = outfile_base+'_'
    if isinstance(indata,str):
        df = pd.read_table(indata,sep='\t',index_col = 0)
        if outfile is None:
            (fname,ext) = os.path.splitext(indata)
            outfile_base = fname+'_bin'+str(int(binsize))+'_'
        indatastr = indata
    else:
        df = indata
        if outfile is None:
            outfile_base = 'bin'+str(int(binsize))+'_'
        indatastr = 'DataFrame'

    if 'iteration' not in df.columns:
        df['iteration'] = np.zeros_like(df[paramname[0]])

    #--set output file names and moviedirs
    outfiles = [outfile]*len(paramname)
    moviedirs = [None]*len(paramname)
    badparams = []
    for p in xrange(len(paramname)):
        outfiles[p] = outfile_base+ctype[p]+'_'+paramname[p]+'.fits'
        moviedirs[p] = outfile_base+ctype[p]+'_'+paramname[p]+'_movie/'

        #--check if output file already exists--
        if os.path.isfile(outfiles[p]) and clobber is not True:
            print ("Warning: "+outfile+" exists and clobber=False. "
                   "Not mapping "+paramname[p]+".")
            badparams = badparams + [paramname[p]]
            #-check if sigthreshparam is being removed-
            if paramname[p] == sigthreshparam:
                print ("Warning: sigthreshparam is not being mapped. "
                       "Resetting sigthresh=0.0")
                sigthreshparam = None
                sigthresh = 0.0

            #-check if imgstdthreshparam is being removed-
            if paramname[p] == imgstdthreshparam:
                print ("Warning: imgstdthreshparam is not being mapped. "
                       "Resetting imgstdthresh=None")
                imgstdthreshparam = None
                imgstdthresh = None

    #--remove parameters that would be clobbered if clobber=False--
    for b in badparams:
        bi = paramname.index(b)
        outfiles.remove(outfiles[bi])
        moviedirs.remove(moviedirs[bi])
        paramname.remove(paramname[bi])
        iteration_type.remove(iteration_type[bi])
        ctype.remove(ctype[bi])
        paramweights.remove(paramweights[bi])
        movie.remove(movie[bi])

    #----Set default image size and center----
    if imagesize is None:
        ximagesize = 1.1*(max(df[paramx] - min(df[paramx])))
        yimagesize = 1.1*(max(df[paramx] - min(df[paramy])))
    elif isinstance(imagesize,tuple) or isinstance(imagesize,list):
        if len(imagesize)>2: 
            print ("calculate_map: Warning: imagesize has too many"+ 
                   " elements, using first two only")
        if len(imagesize)>=2:
            ximagesize=imagesize[0]
            yimagesize=imagesize[1]
        if len(imagesize)==1:
            ximagesize=imagesize[0]
            yimagesize=imagesize[0]
    else:
        ximagesize=imagesize
        yimagesize=imagesize
    if x0 is None:
        x0 = (max(df[paramx])-min(df[paramx]))/2.0+min(df[paramx])
    if y0 is None:
        y0 = (max(df[paramy])-min(df[paramy]))/2.0+min(df[paramy])

    ximageradius = ximagesize/2
    yimageradius = yimagesize/2
    xmin = x0 - ximageradius
    xmax = x0 + ximageradius
    ymin = y0 - yimageradius
    ymax = y0 + yimageradius
    ximageradius = (xmax-xmin)/2.0
    yimageradius = (ymax-ymin)/2.0
    ximagesize = ximageradius*2.0
    yimagesize = yimageradius*2.0

    print 'x,yimagesize,x0,y0,xmin,ymin = ',ximagesize,yimagesize,x0,y0,xmin,ymin

    #-number of map layers (one per iteration) and number of pixels-
    niter = np.unique(df['iteration']).size
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
    print 'nbins_x, nbins_y, nlayers = ',nbins_x,nbins_y,nlayers

    imgs = [] #empty list of image arrays (one per parameter)
    errimgs = [] #empty list of image arrays (one per parameter)

    #-initialize image stack or arguments list-
    nparams = len(paramname)
    if parallel is False:
        image_stacks = np.zeros((nbins_x,nbins_y,nparams,nlayers))
    else:
        imgargs = [[]]*nlayers

    #--Remove iterations according to itmod--

    #-make list of iterations to use-
    # randomly chooses the required number of iterations
    #  from iterations which exist in the dataframe
    its = np.random.choice(df['iteration'].unique(),size=nlayers,
                           replace=False)
    itstr = ['iteration']*len(its)

    #-keep only matching iterations-
    df = filterblobs(df,itstr,minvals=its,maxvals=its,logic='or')

    #----Calculate Blob Volumes----
    if 'blob_volume' not in df.columns:
        if paramshape == 'gauss':
            df['blob_volume'] = (2.0*np.pi*np.square(df[paramsize]))**1.5
        if shape == 'sphere':
            df['blob_volume'] = (4.0/3.0)*np.pi*df[paramsize]**3.0
        if shape == 'points':
            df['blob_volume'] = (0.1*binsize)**3.0 # set to much smaller 
                                                   # than pixel

    #----Group by Iteration----
    layers = df.groupby('iteration')

    #----Iterate over groups (i.e. iterations)----
    layer = 0
    for i, group in layers: 

        if parallel is False: # create iteration images in serial
            print 'layer = ',layer
            #i=iteration number, group = subset of dataframe
            image_stacks[:,:,:,layer] = iteration_image(group,paramname,
                                 paramweights,
                                 nbins_x,nbins_y,binsize,xmin,ymin,
                                 iteration_type,paramshape,paramx,paramy,
                                 paramsize,cint,fast=True,
                                 n_int_steps=n_int_steps)
        else: # construct argument lists for multiprocessing
            imgargs[layer] = [group,paramname,paramweights,nbins_x,nbins_y,
                              binsize,xmin,ymin,iteration_type,paramshape,
                              paramx,paramy,paramsize,cint]
        layer = layer + 1

    # using multiprocessing package
    if parallel is True:
        pool=Pool(nproc)
        image_stacks = np.array(pool.map(iteration_image_star,
                                                 imgargs))
        pool.close()
        pool.join()
        image_stacks = image_stacks.swapaxes(0,2).swapaxes(0,1).swapaxes(2,3)        
        
    #----Loop through parameters to create and manipulate final images----

    #--Collapse Image Stack (combine iterations)--
    collapsed_images = np.zeros((nbins_x,nbins_y,nparams))
    err_images = np.zeros((nbins_x,nbins_y,nparams))
    for p in xrange(len(paramname)):
        themap = collapse_stack(image_stacks[:,:,p,:],
                                                 ctype=ctype[p])

        #--Apply significance threshold and Create Error Maps--
        if (sigthresh != 0.0) or (witherror is True):
            # - compute error (standard deviation) map -
            errmap = collapse_stack(image_stacks[:,:,p,:],
                                               ctype='error')
        else:
            errmap = None

        collapsed_images[:,:,p] = themap
        err_images[:,:,p] = errmap

        #--Mask Region--
        # not yet functional
        if exclude_region is not None:
            msk = circle_mask(df,paramx,paramy,exclude_region,binsize,
                              imagesize,x0,y0)
            themap = themap*msk
            if errmap is not None: errmap = errmap*msk

        #--Rotate Image--
        if rotation != 0.0:
            themap = rotate(themap,rotation,axes=(0,1))
            errmap = rotate(errmap,rotation,axes=(0,1))
        #after testing, add cval=np.nan argument                

        #--Save images to list--
        imgs = imgs+[themap]
        errimgs = errimgs+[errmap]

        #--Make movie--
        if movie[p] is True: movie_from_stack(themap,moviedirs[p],
                                       cumulativemovie=cumulativemovie,
                                       parallel=parallel)
        

    #--Loop through images, apply thresholds, and save to fits--
    # (must be separate loop to allow for sigthreshparam!=param[p])
    if sigthreshparam is not None:
        #-sigthreshparam has just been mapped-
        if sigthreshparam in paramname:
            sigp = paramname.index(sigthreshparam)
            sigmap = abs(imgs[sigp])/errimgs[sigp]        

    if imgstdthreshparam is not None:
        #-sigthreshparam has just been mapped-
        if imgstdthreshparam in paramname:
            imgstdp = paramname.index(imgstdthreshparam)
            imgstdmap = imgs[imgstdp]       

    for p in xrange(len(paramname)):

        themap=imgs[p]
        errmap=errimgs[p]

        #--Apply significance threshold--
        if sigthreshparam is None:
            sigmap = abs(imgs[p])/errimgs[p]

            # - set pixels with significance < threshold to Nan - 
        if sigthresh != 0.0:
            themap[sigmap < sigthresh] = np.nan

        #--Apply img std threshold--
        if imgstdthresh != None:
            if imgstdthreshparam is None:
                imgstdmap = imgs[p]

            # - set pixels with value < threshold to Nan - 
            imgmin = np.nanmax(imgstdmap)-imgstdthresh*np.nanstd(imgstdmap)
            if imgmin < 0.0: imgmin = 0.0
            themap[imgstdmap < imgmin] = np.nan

        #--Save map to fits file--

        #--write history--
        history1 = ('make_map,'+indatastr+',outfile='+nstr(outfiles[p])
                    +',paramname='+nstr(paramname[p])
                    +',paramweights='+nstr(paramweights[p])
                    +',paramx='+nstr(paramx)+',paramy='+nstr(paramy)
                    +',paramsize='+nstr(paramsize)+',binsize='
                    +nstr(binsize)+',itmod='+nstr(itmod)+',paramshape='
                    +nstr(paramshape)+',ctype='+nstr(ctype[p])
                    +',iteration_type='
                    +nstr(iteration_type[p])+',x0='+nstr(x0)+',y0='
                    +nstr(y0)
                    +',imagesize='+nstr(imagesize)+',sigthresh='
                    +nstr(sigthresh)+',sigthreshparam='
                    +nstr(sigthreshparam)+',imgstdthresh='
                    +nstr(imgstdthresh)+',imgstdthreshparam='
                    +nstr(imgstdthreshparam)+',movie='
                    +str(movie[p])
                    +',moviedir='+moviedirs[p])
        history3 = 'ximagesize = '+nstr(ximagesize)
        history4 = 'yimagesize = '+nstr(yimagesize)
        history2 = 'Created '+str(time.strftime("%x %H:%M:%S"))

        #--write file--
        hdr = fits.Header()
        hdr['HISTORY']=history2
        hdr['HISTORY']=history1
        hdr['HISTORY']=history3
        hdr['HISTORY']=history4
        hdu = fits.PrimaryHDU(themap,header=hdr)

        hdu.writeto(outfiles[p],clobber=clobber)

        if witherror is True:
            hdr=fits.Header()
            hdr['HISTORY']=history2
            hdr['HISTORY']=history1
            hdr['HISTORY']=history3
            hdr['HISTORY']=history4
            hdr['HISTORY']='error map'
            fits.append(outfiles[p],errmap,hdr)

        if withsignificance is True:
            hdr=fits.Header()
            hdr['HISTORY']=history2
            hdr['HISTORY']=history1
            hdr['HISTORY']=history3
            hdr['HISTORY']=history4
            hdr['HISTORY']='significance (img/errimg) map'
            fits.append(outfiles[p],sigmap,hdr)

#        imgs = imgs+[themap]
#        errimgs = errimgs+[errmap]
    return imgs

#--------------------------------------------------------------------------
def movie_from_stack(stack,moviedir,cumulativemovie=False,ctype='median',
                     delay=20,cmap='CMRmap',parallel=True):
    """
    Save each image in a stack to file for creating a movie

    See http://matplotlib.org/examples/color/colormaps_reference.html for
    list of colormaps.

    """

    # - import plotting functions - 
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
#    from scipy.misc import imsave

    # - create directory if it doesn't exist -
    if not os.path.exists(moviedir):
        os.makedirs(moviedir)

    # - set up figure object - 
    fig = plt.figure()

    # - loop over layers and save images to file -
    nlayers = stack.shape[2]
    for layer in xrange(nlayers):
        framenum="%03d" % (layer,)
        print 'frame '+framenum

        fig.set_size_inches(5,5)
        ax = plt.Axes(fig,[0.,0.,1.,1.])
        ax.set_axis_off()
        fig.add_axes(ax)

        if (cumulativemovie is False) or (layer == 0):
            im = stack[:,:,layer]
        else:
            if ctype == 'average':
                collapsed_img = np.average(stack[:,:,:layer],axis=2)
            elif ctype == 'median':
                collapsed_img = np.median(stack[:,:,:layer],axis=2)
            elif ctype == 'total':
                collapsed_img = np.sum(stack[:,:,:layer],axis=2)
            elif ctype == 'error':
                collapsed_img = np.std(stack[:,:,:layer],axis=2)
            else: 
                print "movie_from_stack: ERROR: unrecognized ctype"
            im = collapsed_img
            print 'min,max im = ',np.min(im),np.max(im)

        # - convert to color image and plot - 
#        fig,ax=plt.subplots()
        ax.imshow(im,cmap=cmap,aspect='auto')

        # - save to file -
        fig.savefig(moviedir+'frame'+framenum+'.png',dpi=300)
#                    bbox_inches=0)

        fig.clf() # clear figure

    plt.close(fig) # close figure (release memory)
        
    # - combine frames into movie - 
    cmd = 'convert -set delay '+str(int(delay))+' '+moviedir+'/frame*.png '+moviedir+'/movie.gif'
    os.system(cmd)

    return None


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
def gaussian_integral_quad(lowerx,upperx,blobx,blobsize,use_ctypes=True):
    """Function to calculate gaussian integral using scipy.integrate.quad()."""

    if use_ctypes is False:
        result = quad(lambda x: gaussian1D(x,blobx,blobsize),lowerx,upperx )
    else:
        # - get gaussian function information -
        gausspath = os.path.dirname(os.path.abspath(inspect.getfile(
                    inspect.currentframe())))
        lib = ctypes.CDLL(gausspath+'/gaussian.so')
        cgauss = lib.gaussian # get function name from c library
        cgauss.restype = ctypes.c_double
        cgauss.argtypes = (ctypes.c_int,ctypes.c_double)
    
        # - integrate -
        result = quad(cgauss,lowerx,upperx,(blobx,blobsize) )
        #result = nquad(cgauss,[[lowerx,upperx]],args=[blobx,blobsize])

    return result[0]
    
#--------------------------------------------------------------------------
def gaussian1D(x,mux,sigma):
    """1-D Gaussian function"""
    return np.exp(-1.0/2.0*(x-mux)**2.0/sigma**2.0)

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
    dummy = pd.Series(np.ones_like(df[paramx].values)
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
def iteration_image(data,params,weights,nbins_x,nbins_y,binsize,xmin,ymin,
                    iteration_type,shape,blobx,bloby,blobsize,use_ctypes,
                    fast=True,
                    n_int_steps=1000):
    """Function to combine blobs from single iteration into 1 image."""
    from wrangle import weighted_median,gaussian_volume

    #--initialize stack of 2D images, one for each parameter--
    iterimages = np.zeros((nbins_x,nbins_y,len(params)))

   #----Calculate blob volumes in correct units----
    if shape == 'gauss':
        volumes = gaussian_volume(data[blobsize]) 
    elif shape == 'sphere':
        volumes = (4.0/3.0)*np.pi*data[blobsize]**3.0
    else: # points
        volumes = (0.1*binsize)**3.0 # set to much smaller than pixel

    #--loop over image--
    for x in xrange(nbins_x):
        #get x integral
        lowerx = int(xmin + x*binsize)
        upperx = int(xmin + x*binsize + binsize)
        if shape == 'gauss' or shape == 'points':
            if fast is False: 
                # only use fast=False if scipy.integrate is
                #   not available
                x_blob_integrals = gaussian_integral(lowerx,upperx,
                                                     n_int_steps,
                                                     data[blobx],
                                                     data[blobsize])
            else:
                x_blob_integrals = data.apply(lambda d: \
                                gaussian_integral_quad(lowerx,\
                                upperx,d[blobx],d[blobsize],\
                                use_ctypes=use_ctypes),\
                                axis=1)
        elif shape == 'sphere':
            print "ERROR: spherical_integral() not yet implemented"
            x_blob_integrals = spherical_integral(lowerx,upperx,\
                                                      n_int_steps,\
                                                     data[blobx],
                                                  data[blobsize])
        for y in xrange(nbins_y):
            #get y integral
            lowery = int(ymin + y*binsize)
            uppery = int(ymin + y*binsize + binsize)
            if shape == 'gauss' or shape == 'points':
                if fast is False:
                    y_blob_integrals = gaussian_integral(lowery,uppery,\
                                                     n_int_steps,\
                                                     data[bloby],
                                                         data[blobsize])
                else:
                    y_blob_integrals = data.apply(lambda d: \
                                gaussian_integral_quad(lowery,\
                                uppery,d[bloby],d[blobsize],\
                                use_ctypes=use_ctypes),\
                                axis=1)

            elif shape == 'sphere':
                y_blob_integrals = spherical_integral(lowery,uppery,\
                                                     n_int_steps,\
                                                     data[bloby],
                                                      data[blobsize])
                #calculate fraction of blob volume in this pixel

            if shape != 'points':
                # !! for now this assumes gaussian volume !!
                fractions = (x_blob_integrals*y_blob_integrals*
                             (2.0*np.pi*data[blobsize]**2.0)**.5 / volumes)
                #times dz integral to get total volume in pixel, 
                #then divided by total volume

            else:
                # for now, points is implemented by setting the volumes 
                #   to be much smaller than a pixel size
                fractions = (x_blob_integrals*y_blob_integrals*
                             (2.0*np.pi*data[blobsize]**2.0)**.5 / 
                             volumes)
#                    print "points is not yet implemented"
                    # if assuming points, then fraction=1 or 0
#                    fractions = point_integral(lowerx,upperx,lowery,uppery,
#                                               data['x'],data['y'])

            #-combine blobs in this pixel (loop over parameters)-
            for p in xrange(len(params)):
                if weights[p] is None: # default is equal weights
                    w = pd.Series(np.ones_like(data[params[p]]),
                                  index=data[params[p]].index)
#                elif weights[p] == 'densityspecial': 
                    # for plotting density from EM - assumes the column passed
                    # was sqrt(EM*Volume/1.21), so weights=1/Vblobpix
#                    w = 1.0/(fractions*data['blob_volume'])
                else: 
                    w = data[weights[p]]
                if iteration_type[p] == 'median':
                    iterimages[x,y,p]=weighted_median(data[params[p]],
                                                  weights=w*fractions)
                elif iteration_type[p] == 'average':
                    iterimages[x,y,p]=np.average(data[params[p]],
                                                 weights=w*fractions)
                elif iteration_type[p] == 'total':
                    iterimages[x,y,p]=np.sum(data[params[p]]*w*fractions)
                elif iteration_type[p] == 'max':
                    iterimages[x,y,p]=np.max(data[params[p]]*w*fractions)
                else:
                    print "ERROR: unrecognized iteration_type"

    return iterimages

#--------------------------------------------------------------------------
def iteration_image_star(arglist):
    """Function to unpack list of arguments and pass to iteration_image()"""
    # for use with multiprocessing package
#    print 'iteration = ',arglist[0].iteration[0]
    return iteration_image(*arglist)

#--------------------------------------------------------------------------
def collapse_stack(img_stack,ctype):
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
