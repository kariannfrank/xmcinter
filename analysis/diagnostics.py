"""
Scripts that call xmcinter functions to do common diagnostics tasks.

To run any of these scripts, first:
>> import xmcinter.scripts.diagnostics as xd
This should import all the common modules.

Then, call any of the following by prepending with xd

A typical workflow might be as follows:

# start python
ipython
# import scripts and common modules
import xmcinter.scripts.diagnostics as xd 
# check convergence (determine minimum iteration to use)
statframe = xplt.chi2('./')
# filter by iteration and add emission measure column
# (change itmin and distance as appropriate)
df = xd.diagnostics(itmin=1000,distance=3.3)
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
from xmcinter.analysis.filter import filter
import xmcinter.astro_utilities as astro


#----------------------------------------------------------
#Name: iteration_filter
#Author: Kari A. Frank
#Date: November 1, 2015
#Purpose: Create a dataframe and associated saved file that includes
#         an iteration and emission measure column, and only includes
#         iterations after convergence 
#
#Usage: 
#  import xmcinter.diagnostics as xd
#  xd.iteration_filter(runpath='./',itmin=0,distance=8.0)
#
#Input:
# 
# runpath: string of path to the deconvolution files
# 
# itmin: minimum iteration to keep
#
# distance: distance to the object in kpc (default=8.0), used to 
#           calculate the emission measure
#
#Output:
# 
# Returns the dataframe of deconvolution parameters, filtered by iteration
#  add with the emission measure and iteration columns included, plus a 
#  column with the blob sizes in log10(arcsec)
#   
#
#Usage Notes:
# - typically this is run after xplt.chi2, to determine the minimum iteration 
#
#Example:
# 
#
def diagnostics(runpath='./',itmin=0,distance=8.0):

    # -- read deconvolution files --
    df = merge_output(runpath,save=False)

    # -- add log10(arcsec) blob size column --
    df['blob_log10sigma'] = np.log10(np.exp(df['blob_lnsigma']))

    # -- add emission measure column
    df['blob_em'] = astro.norm_to_em(df['blob_norm'],astro.convert_distance(distance,'kpc','cm'))
    
    # -- remove iterations before convergence --
    df = filter(df,'iteration',minvals=itmin)

    # -- save as file --
    outfile = 'deconvolution_merged_itmin'+str(int(itmin))+'.txt'
    df.to_csv(outfile,sep='\t')

    # -- make traceplots --
    tracefigs = xplt.traceplots(df)

    return df
