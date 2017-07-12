"""
Scripts that call xmcinter functions to do common diagnostics tasks.

Contents:
 clean
 check
"""
#----------------------------------------------------------

#-import common modules-
import pandas as pd
import numpy as np
import plots as xplt
import wrangle as xw
from xmcfiles import merge_output

#----------------------------------------------------------
def clean(runpath='./',itmin=0,itmax=None,distance=8.0):
    """
    Perform basic processing of blobs and save to merged file.

    Name: clean
    Author: Kari A. Frank
    Date: November 1, 2015
    Purpose: Create a dataframe and associated saved file that includes
             an iteration and emission measure column, and only includes
             iterations after convergence 

    Usage: 
      import xmcinter.diagnostics as xd
      xd.clean(runpath='./',itmin=0,itmax=None,distance=8.0)

    Input:

     runpath: string of path to the deconvolution files, or a dataframe
              resulting from a previous call to xw.merge_output()

     itmin: minimum iteration to keep

     distance: distance to the object in kpc (default=8.0), used to 
               calculate the emission measure


    Output:

     Returns the dataframe of deconvolution parameters, filtered by iteration
      add with the emission measure and iteration columns included, plus a 
      column with the blob sizes in arcsec (if blob shape = gaussian)


    Usage Notes:
     - typically this is run after xplt.chi2, to determine the minimum 
       iteration 
     - assumes that relevant column names begin with 'blob'. if not found,
       will skip adding the new column.

    Example:
    """

    # -- import modules --
    import astro_utilities as astro

    # -- read deconvolution files --
    if isinstance(runpath,str):
        df = merge_output(runpath,save=False)
    else:
        df = runpath

    # -- remove iterations before convergence --
    if itmax == None:
        itmax = np.max(df['iteration'])
    if itmin == None:
        itmin = np.min(df['iteration'])
    df = xw.filterblobs(df,'iteration',minvals=itmin,maxvals=itmax)
#    df = df[(df['iteration'] >= itmin) & (df['iteration']<=itmax)]
        
    # -- add blob size in arcsec column --
    if 'blob_lnsigma' in df.columns:
        df['blob_sigma'] = np.exp(df['blob_lnsigma'])

    # -- add tau column, if used lvpshock --
    if 'blob_logtau' in df.columns:
        df['blob_tau'] = 10.0**(df['blob_logtau'])

    # -- add emission measure column --
    if 'blob_norm' in df.columns:
        df['blob_em'] = astro.norm_to_em(df['blob_norm'],
                                         astro.convert_distance(distance,
                                                                'kpc',
                                                                'cm'))
    print 'saved EM = ',df.blob_em.iloc[0]

    # -- add hydrogen number densities of blobs in cm^-3, hydrogen mass --
    if 'blob_sigma' in df.columns:
        df['blob_volume'] = astro.gaussian_volume(astro.convert_arcsec(\
                df['blob_sigma'],distance,'kpc','cm'))
        df['blob_numberdensity'] = astro.em_to_density(df['blob_em'],\
                                   df['blob_volume'],density_type='number')

        df['blob_mass'] =astro.em_to_mass(df['blob_em'],df['blob_volume'],
                                          tounit='sol')        

    # -- save as file --
    outfile = ('deconvolution_merged_iter'
               +str(int(itmin))+'-'+str(int(itmax))+'.txt')
    df.to_csv(outfile,sep='\t')

    # -- make traceplots --
#    tracefigs = xplt.traceplots(df)

    # -- check EM calculation --
    testem = astro.norm_to_em(df.iloc[0]['blob_norm'],
                              astro.convert_distance(distance,'kpc','cm'))
    print 'test EM = ',testem

    return df

#----------------------------------------------------------
def check(runpath='./',outpath='./',itmin=0,itmax=None,kTthresh=None,
          cint=False,display=False,init_file=None,skipspectrum=False,
          legacy=False):
    """
    Name: check
    Author: Kari A. Frank
    Date: November 23, 2015
    Purpose: Make and display some diagnostic plots to check
             ongoing xmc runs.  Plots chi2, norm image, spectrum, and
             histograms.
            
    Usage: 
      import xmcinter.diagnostics as xd
      xd.check(runpath='./',itmin=0,itmax=None)

    Input:

     runpath: string of path to the deconvolution files

     outpath: string of path to store output files

     display (bool) : if False, then will not display the figures 

     itmin/itmax: minimum iteration to use in check

     kTthresh: optionally specify a minimum kT. if set to 
               numeric value, then all blobs with kT<kTthresh will 
               not be included in the map or the kTthresh histograms
               (but will be included in all other plots). 
               if both kTthresh and init_file are provided, then 
               nHkTthresh will be applied first, and then kTthresh.


    init_file: if set to a string file name, that file must exist 
               in xmcinter/scripts/ directory, and must contain a 
               definition of nHkTthresh() function, distance, pixelsize, 
               mapsize, and x0,y0. 
               - typically this file should be the the init file,
                 i.e. <object>_<obsid>_init.py

    skipspectrum: skip plotting of the spectrum

    legacy: assume that spectrum file names are iteration/100. ignored
            if skipspectrum=True

    Output:

     - Displays plots.
     - Returns DataFrames of blob parameters and statistics as tuple (df,sf)

    Usage Notes:


    Example:
    """

    # -- import modules --
    import os
    from file_utilities import ls_to_list
    import xmcmap as xm
    from xmcfiles import read_spectra
   
    # -- define defaults -- 
    # - these should be overwritten if init_file is provided -
    distkpc = 3.0 
    x0 = None # automatically calculate in make_map()
    y0 = None # automatically calculate in make_map()
    rotation = 0.0
    pixelsize = None # default reset later to match data size
    mapsize = None # automatically calculate in make_map()
    nbins = 75
    w = 500
    h = 200
    z = 0.0

    # -- read init file if provided --
    if init_file is not None:
        f,ext = os.path.splitext(init_file) # remove extension
        import importlib
        il = importlib.import_module('scripts.'+f)
        #        __import__('scripts.'+f,fromlist=[''])
        distkpc = il.distkpc
        x0 = il.x0
        y0 = il.y0
        rotation = il.rotation 
        pixelsize = il.pixelsize
        mapsize = il.mapsize
        nbins = il.nbins
        w = il.w
        h = il.h
        z = il.z
        
    if kTthresh is None:
        kTthresh = 0.0

    # -- plot chi2 --
    print "\nPlotting chi2 ...\n"
    sf = xplt.chi2(runpath,itmax=itmax,display=display,
                   outfile=outpath+'/chi2_vs_iteration.html')
    
    # -- calculate median chi2 --
    if itmax is None:
        itmax = np.max(sf.iteration)
    sf = xw.filterblobs(sf,'iteration',minvals=itmin,maxvals=itmax)
    print 'Filtered sf'
    medchi2 = xw.weighted_median(sf['redchi2'])

    # -- read deconvolution files --
    print "\nReading deconvolution files ...\n"
    dfall=merge_output(runpath,save=False,itmin=itmin,itmax=itmax)

    
    # -- set default pixelsize --
    if pixelsize is None:
        pixelsize = (dfall.blob_phi.max()-dfall.blob_phi.min())/50.0
    
    # -- add derivative columns --
    dfall = clean(dfall,itmin=itmin,itmax=itmax,distance=distkpc)
    if itmax == None:
        itmax = np.max(dfall['iteration'])
    if itmin == None:
        itmin = np.min(dfall['iteration'])

    print '\nIterations '+str(itmin)+' - '+str(itmax)+' : '
    print "Total Number of Blobs = ",len(dfall.index)
    print 'Median chi2/dof = '+str(medchi2)+'\n'

    # -- plot model and data spectra --
    if skipspectrum is False:
        print "\nPlotting spectrum ...\n"

        sfig = xplt.standard_spectra(runpath=runpath,display=display,
                                     itmin=itmin,legacy=legacy,
                                 itmax=itmax,
                                 outfile='spectra_all.html',ylog=True,
                                 xlog=False,logbins=None,bins=0.03,
                                 lines=True,emissivity_range=(1e-17,1.0),
                                 nlines=100,energy_range=(0.87,10.0))
    
    # -- make median traceplots --
    print "\nPlotting traces ...\n"
    efig = xplt.trace(dfall,weights=None,display=display,
                      outfile=outpath+'/trace_plots.html')

    # -- make histograms --
    print "\nPlotting posteriors ...\n"
    
    hfigs = xplt.histogram_grid([dfall,dfall],weights=[None,'blob_mass'],
                                bins=nbins,ncols=3,norm=True,
                                display=display,
                                outfile=outpath+'/histogram_grid_allblobs.html',
                                legends=['Unweighted','Mass weighted'],
                                width=w,height=h,iterations='iteration')

    print "\nPlotting filtered posteriors ...\n"
    if init_file is not None:
        dfgood = il.nHkTthresh(dfall)
        dfgood = dfgood[dfgood.blob_kT>=kTthresh]
    else:
        dfgood = dfall[dfall.blob_kT>=kTthresh]
    print "Total Number of Filtered Blobs = ",len(dfgood.index)
        
    hfigs = xplt.histogram_grid([dfgood,dfgood],
                                weights=[None,'blob_mass'],
                                display=display,
                                bins=nbins,ncols=3,norm=True,
                            outfile=outpath+'/histogram_grid_filtered.html',
                                legends=['Unweighted','Mass weighted'],
                                width=w,height=h,iterations='iteration')

    # -- scatter plots--
    print "\nPlotting scatter plots ...\n"
    blobcols = [c for c in dfall.columns if 'blob' in c]
    sfigs2 = xplt.scatter_grid(dfall[blobcols],agg=None,sampling=1000,
                               display=display)

    # -- make norm map from most recent iteration --
    print "\nMaking blob norm map ...\n"
        
    img1file = (outpath+'/bin'+str(int(pixelsize))+
                '_iter'+str(itmax))

    img = xm.make_map(dfgood,x0=x0,y0=y0,
                      paramname='blob_norm',cint=cint,
                      paramweights=None,iteration_type='total',
                      binsize=pixelsize,nlayers=1,
                      rotation=rotation,imagesize=mapsize,
                      withsignificance=True,nproc=2,
                      outfile=img1file,clobber=True)

    return (dfall,sf)

