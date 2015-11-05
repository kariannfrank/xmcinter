"""
Scripts that call xmcinter functions to do common diagnostics tasks.

A typical workflow might be as follows:

# start python
ipython
# import scripts and common modules
import xmcinter.scripts.diagnostics as xd 
# check convergence (determine minimum iteration to use)
statframe = xplt.chi2('./')
# filter by iteration and add emission measure column
# (change itmin and distance as appropriate)
df = xd.clean(itmin=1000,distance=3.3)
# check traceplots
tracefigs = xplt.traceplots(df)
# check weighted histograms
histfigs = xplt.histogram_grid(df,weights=df['blob_em'])

"""
#----------------------------------------------------------

#-import common modules-
#import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xmcinter.plots as xplt
from xmcinter.files.xmcrun import merge_output
from xmcinter.analysis.wrangle import filterblobs
import xmcinter.astro_utilities as astro

#----------------------------------------------------------
"""
Name: clean
Author: Kari A. Frank
Date: November 1, 2015
Purpose: Create a dataframe and associated saved file that includes
         an iteration and emission measure column, and only includes
         iterations after convergence 

Usage: 
  import xmcinter.diagnostics as xd
  xd.clean(runpath='./',itmin=0,distance=8.0)

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
 - typically this is run after xplt.chi2, to determine the minimum iteration 

Example:

"""
def clean(runpath='./',itmin=0,distance=8.0):

    # -- read deconvolution files --
    df = merge_output(runpath,save=False)

    # -- add log10(arcsec) blob size column --
    if 'blob_lnsigma' in df.columns:
        df['blob_sigma'] = np.exp(df['blob_lnsigma'])

    # -- add emission measure column
    df['blob_em'] = astro.norm_to_em(df['blob_norm'],
                                     astro.convert_distance(distance,'kpc',
                                                            'cm'))/(10.0**55.0)

    # -- remove iterations before convergence --
    df = filterblobs(df,'iteration',minvals=itmin)

    # -- save as file --
    outfile = 'deconvolution_merged_itmin'+str(int(itmin))+'.txt'
    df.to_csv(outfile,sep='\t')

    # -- make traceplots --
    tracefigs = xplt.traceplots(df)

    return df
