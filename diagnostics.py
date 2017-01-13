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

    # -- add blob size in arcsec column --
    if 'blob_lnsigma' in df.columns:
        df['blob_sigma'] = np.exp(df['blob_lnsigma'])

    # -- add linear tau column if used lvpshock --
    if 'blob_logtau' in df.columns:
        df['blob_tau'] = 10.0**df.blob_logtau

    # -- add tau column, if used lvpshock --
    if 'blob_logtau' in df.columns:
        df['blob_tau'] = 10.0**(df['blob_tau'])

    # -- add emission measure column --
    if 'blob_norm' in df.columns:
        df['blob_em'] = astro.norm_to_em(df['blob_norm'],
                                         astro.convert_distance(distance,
                                                                'kpc',
                                                                'cm'))
    
    # -- add hydrogen number densities of blobs in cm^-3, hydrogen mass --
    if 'blob_sigma' in df.columns:
        df['blob_volume'] = xw.gaussian_volume(astro.convert_arcsec(\
                df['blob_sigma'],distance,'kpc','cm'))
        df['blob_numberdensity'] = astro.em_to_density(df['blob_em'],\
                                   df['blob_volume'],density_type='number')

        df['blob_mass'] =astro.em_to_mass(df['blob_em'],df['blob_volume'],
                                          tounit='sol')

    # -- remove iterations before convergence --
    if itmax == None:
        itmax = np.max(df['iteration'])
    df = xw.filterblobs(df,'iteration',minvals=itmin,maxvals=itmax)

    # -- save as file --
    outfile = ('deconvolution_merged_iter'
               +str(int(itmin))+'-'+str(int(itmax))+'.txt')
    df.to_csv(outfile,sep='\t')

    # -- make traceplots --
#    tracefigs = xplt.traceplots(df)

    return df

#----------------------------------------------------------
def check(runpath='./',outpath='./',itmin=0,itmax=None,kTthresh=0.17,
          cint=False):
    """
    Name: clean
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

     itmin/itmax: minimum iteration to use in check

     kTthresh: float to optionally specify a minimum kT. if set,
               then all blobs with kT<kTthresh will not be included
               in the map (but will be included in all other plots)

    Output:

     - Displays plots.
     - Returns DataFrame of blob parameters

    Usage Notes:


    Example:
    """

    # -- import modules --
    import os
    from file_utilities import ls_to_list
    import xmcmap as xm

    if kTthresh is None:
        kTthresh = 0.0

    # -- plot chi2 --
    print "\nPlotting chi2 ...\n"
    sf = xplt.chi2(runpath,itmax=itmax,
                   outfile=outpath+'/chi2_vs_iteration.html')
    
    # -- calculate median chi2 --
    if itmax is None:
        itmax = np.max(sf.iteration)
    sf = xw.filterblobs(sf,'iteration',minvals=itmin,maxvals=itmax)
    print 'Filtered sf'
    medchi2 = xw.weighted_median(sf['redchi2'])

    # -- read deconvolution files --
    print "\nReading deconvolution files ...\n"
    dfall=merge_output(runpath,save=False)

    # -- add derivative columns --
    dfall = clean(dfall,itmin=itmin,itmax=itmax)
    print '\nIterations '+str(itmin)+' - '+str(itmax)+' : '
    print "Total Number of Blobs = ",len(dfall.index)
    print 'Median chi2/dof = '+str(medchi2)+'\n'

    # -- plot model and data spectra --
    print "\nPlotting spectrum ...\n"
    smin = itmin/100
    if itmax is None: 
        smax = None
    else:
        smax = itmax/100
    sfig = xplt.spectrum(runpath=runpath,smin=smin,smax=smax,bins=0.03,
                         ylog=True,xlog=True,
                         outfile=outpath+'/spectrum.html',
                         lines=True,nlines=100,energy_range=(0.5,10.0))

    # -- make median traceplots --
    print "\nPlotting traces ...\n"
    efig = xplt.trace(dfall,weights=None,
                      outfile=outpath+'/trace_plots.html')

    # -- make histograms --
    print "\nPlotting posteriors ...\n"
    nbins = 75
    w = 500
    h = 200
    hfigs = xplt.histogram_grid([dfall,dfall],weights=[None,'blob_em'],
                                bins=nbins,ncols=2,norm=True,
                                outfile=outpath+'/histogram_grid.html',
                                legends=['Unweighted','EM weighted'],
                                width=w,height=h,iterations='iteration')

    print "\nPlotting posteriors with kT threshold ...\n"
    hfigs = xplt.histogram_grid([xw.filterblobs(dfall,'blob_kT',
                                                minvals=kTthresh),
                                 xw.filterblobs(dfall,'blob_kT',
                                                minvals=kTthresh)],
                                weights=[None,'blob_em'],
                                bins=nbins,ncols=2,norm=True,
                                outfile=outpath+'/histogram_grid_kTthresh.html',
                                legends=['Unweighted','EM weighted'],
                                width=w,height=h,iterations='iteration')

    # -- scatter plots--
    print "\nPlotting scatter plots ...\n"
    blobcols = [c for c in dfall.columns if 'blob' in c]
    sfigs2 = xplt.scatter_grid(dfall[blobcols],agg=None,sampling=2000)

    # -- make norm map from most recent iteration --
    print "\nMaking blob em map ...\n"
    pixelsize = (dfall.blob_phi.max()-dfall.blob_phi.min())/50.0
#    img1file = (outpath+'/bin'+str(int(pixelsize))+
#                '_iter'+str(itmin)+'-'+str(itmax))
    img1file = (outpath+'/bin'+str(int(pixelsize))+
                '_iter'+str(itmax))
    img = xm.make_map(xw.filterblobs(dfall,['blob_kT'],
                                     minvals=[kTthresh,itmax],
                                     maxvals=[None,itmax]),
                      paramname='blob_em',cint=cint,
                      paramweights=None,iteration_type='total',
                      binsize=pixelsize,nlayers=1,
                      withsignificance=True,nproc=2,
                      outfile=img1file,clobber=True)

    return (dfall,sf)

