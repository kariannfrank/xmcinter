# xmcinter
Code to manipulate and interpret output from xmc.  This includes some IDL 
code (mainly just for creating maps), along with the python package xmcinter.

The python package requires the following third-party packages:
- pandas
- numpy
- astropy
- bokeh

At the beginning of analysis python session, typical imports are:
import pandas as pd
import numpy as np
import xmcinter.plots as xplt
from xmcinter.files.xmcrun import merge_output
#import xmcinter.analysis.wrangle as wrangle
from xmcinter.analysis.wrangle import filterblobs
import xmcinter.astro_utilities as astro
import xmcinter.analysis.diagnostics as xd

A typical workflow for early diagnostics is:

# check convergence
xplt.chi2('./')

# filter by iteration and add emission measure (change itmin and distance)
df = xd.clean(runpath='./',itmin=1000,distance=3.3)

# check traceplots
tracefigs = xplt.traceplots(df)

# check weighted histograms
histfigs = xplt.histogram_grid(df[df.columns[:-1]],weights=df['blob_em'])

# filter by some parameter (e.g. temperature), and then re-check plots
df_temp = filterblobs(df,'blob_kT',minvals=0.3,maxvals=3.0)
temphistfigs = xplt.histogram_grid(df_temp[df_temp.columns[:-1]],weights=df_temp['blob_em'])