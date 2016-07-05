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

     runpath: string of path to the deconvolution files

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

    Example:
    """

    # -- import modules --
    import astro_utilities as astro

    # -- read deconvolution files --
    df = merge_output(runpath,save=False)

    # -- add log10(arcsec) blob size column --
    if 'blob_lnsigma' in df.columns:
        df['blob_sigma'] = np.exp(df['blob_lnsigma'])

    # -- add emission measure column --
    df['blob_em'] = astro.norm_to_em(df['blob_norm'],
                                     astro.convert_distance(distance,'kpc',
                                                            'cm'))
    
    # -- add hydrogen number densities of blobs in cm^-3, hydrogen mass --
    df['blob_volume'] = xw.gaussian_volume(astro.convert_arcsec(\
            df['blob_sigma'],distance,'kpc','cm'))
    df['blob_numberdensity'] = astro.em_to_density(df['blob_em'],\
                               df['blob_volume'],density_type='number')

    df['blob_mass'] =astro.em_to_mass(df['blob_em'],df['blob_volume'],tounit='sol')
                               
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
def check(runpath='./',itmin=0,itmax=None):
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

     itmin: minimum iteration to use in check

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

    # -- read deconvolution files --
    df = merge_output(runpath,save=False)

    # -- add log10(arcsec) blob size column --
    if 'blob_lnsigma' in df.columns:
        df['blob_sigma'] = np.exp(df['blob_lnsigma'])

    # -- remove iterations before itmin --
    if itmax == None:
        itmax = np.max(df['iteration'])
    df = xw.filterblobs(df,'iteration',minvals=itmin)

    # -- plot chi2 --
    sf = xplt.chi2(runpath)
    
    # -- print median chi2 --
    medchi2 = xw.weighted_median(sf['redchi2'])
    print ('Median chi2/dof for iteration > '+str(itmin)+': '
           +str(medchi2))

    # -- plot model and data spectra --
    data_wave = xplt.spectra(runpath,smin=itmin/100,smax=itmax/100)

    # -- make norm maps --
#    normmap =  'itmin'+str(itmin)+'_median_norm.fits'
#    median_norm_img = xm.make_map(df,outfile=normmap,paramname='blob_norm'
#                                  ,binsize=10.0,itmod=(itmax-itmin)/500
#                                  ,iteration_type='total')

    normmaplatest = 'iter'+str(itmax)+'_norm.fits'
    df_latest = xw.filterblobs(df,'iteration',minvals=itmax,maxvals=itmax)
    latest_norm_img = xm.make_map(df_latest,outfile=normmaplatest,
                                  paramname='blob_norm',binsize=10.0,
                                  itmod=1,iteration_type='total',nproc=1)
    

    # -- plot parameter histograms --
    histfigs = xplt.histogram_grid(df)

    # -- plot norm maps and event file --

    #-get event file names-
    evfiles = ls_to_list('./',ls_args='*events*.fits')
    cmd = 'ds9'
    for evf in evfiles:
        cmd = cmd + ' '+ evf + ' -bin factor 32'
#    cmd = cmd + ' ' + normmap + ' ' + normmaplatest
    cmd = cmd + ' ' + normmaplatest
    os.system(cmd)


    # -- make traceplots --
#    tracefigs = xplt.traceplots(df)

    return (df,sf)

