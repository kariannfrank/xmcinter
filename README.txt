# xmcinter
Code to manipulate and interpret output from xmc.  This includes some IDL 
code (mainly just for creating maps), along with the python package xmcinter.

The xmcinter python package requires the following third-party packages:
- pandas
- numpy
- astropy
- bokeh

###################################################################
At the beginning of analysis python session, typical imports are:
(can copy the following import statements and paste into ipython
 with %paste)

import pandas as pd
import numpy as np
import xmcinter.plots as xplt
from xmcinter.files.xmcrun import merge_output
import xmcinter.analysis.wrangle as xw
from xmcinter.analysis.wrangle import filterblobs
import xmcinter.astro_utilities as astro
import xmcinter.analysis.diagnostics as xd
import xmcinter.map as xm

###################################################################
A typical workflow for early diagnostics
(requires access to the original deconvolution.* and statistic.* files)

# check convergence
sf = xplt.chi2('./')

# determine median chi2 after convergence
sf = filterblobs(sf,'iteration',minvals=1000)
final_chi2 = xw.weighted_median(sf['redchi2'])

# filter by iteration and add emission measure (change itmin and distance)
df = xd.clean(runpath='./',itmin=1000,distance=3.3)

# check traceplots
tracefigs = xplt.traceplots(df)

# check weighted histograms
histfigs = xplt.histogram_grid(df[df.columns[:-1]],weights=df['blob_em'])

###################################################################
Examples of some common tasks
(all executed from ipython unless specified otherwise with IDL> prompt)

# read in a typical merged deconvolution file
df = pd.read_table('deconvolution_merged.txt',sep='\t',index_col=0)

# filter by some parameter (e.g. temperature)
df_filt = filterblobs(df,'blob_kT',minvals=0.3,maxvals=3.0)

# or filter by multiple parameters at once
# (df_filt will include only blobs with 0.5 <= kT <= 1.0 AND 1.0 <= Si 
#  abundance <= 5.0)
df_filt = filterblobs(df,['blob_kT','blob_Si'],minvals=[0.5,1.0],
#		      maxvals=[1.0,5.0])

# write filtered dataframe to file
df.to_csv('deconvolution_merged_filtered.txt',sep='\t')

# plot weighted histograms of parameters
histfigs = xplt.histogram_grid(df[df.columns[:-1]],weights=df['blob_em'])

# make a map of one of the parameters
IDL> make_map,infile='deconvolution_merged.txt',paramname='blob_kT',weights='blob_em',binsize=4.0
