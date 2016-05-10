# xmcinter
Code to manipulate and interpret output from xmc.

The xmcinter python package requires the following third-party packages:
- pandas
- numpy
- astropy
- bokeh
- scipy

You may need to add xmcinter's parent directory to your PYTHONPATH, e.g.,
if you placed xmcinter inside ~/python_programs/,

setenv PYTHONPATH ~/python_programs/:$PYTHONPATH

###################################################################
At the beginning of analysis python session, some common imports are:
(can copy any of the following import statements and paste into ipython
 with %paste)

import pandas as pd
import numpy as np
import xmcinter.plots as xplt
import xmcinter.xmcfiles as xf
import xmcinter.wrangle as xw
import xmcinter.astro_utilities as astro
import xmcinter.diagnostics as xd
import xmcinter.xmcmap as xm

###################################################################
Usage Notes:

- To use the fastest (and default) method for making maps 
  (via xmcinter.xmcmap.make_map), it is necessary to compile the gaussian.c
  function on your local machine, e.g. on linux:

  gcc -shared -o gaussian.so -fPIC gaussian.c

  If you choose not to use this, make sure to always pass the cint=False
  argument to make_map(). It should still run relatively fast.  

###################################################################
A typical workflow for early diagnostics
(requires access to the original deconvolution.* and statistic.* files)

# check progress
df,sf = xd.check(itmin=500)

# check median chi2 after convergence
# (change 1000 to first converged iteration)
sf = xplt.chi2()
sf = xw.filterblobs(sf,'iteration',minvals=1000)
medchi2 = xw.weighted_median(sf['redchi2'])
print medchi2

# filter by iteration and add emission measure, blob_sigma 
# (change itmin and distance)
df = xd.clean(itmin=1000,distance=3.3)

# check traceplots
tracefigs = xplt.traceplots(df)

# check weighted histograms
histfigs = xplt.histogram_grid(df[df.columns[:-1]],weights=df['blob_em'])

###################################################################
Examples of some other common tasks
(all executed from ipython unless specified otherwise with IDL> prompt)

# read in a typical merged deconvolution file
df = pd.read_table('deconvolution_merged.txt',sep='\t',index_col=0)

# filter by some parameter (e.g. temperature)
df_filt = xw.filterblobs(df,'blob_kT',minvals=0.3,maxvals=3.0)

# or filter by multiple parameters at once
# (df_filt will include only blobs with 0.5 <= kT <= 1.0 AND 1.0 <= Si 
#  abundance <= 5.0)
df_filt = xw.filterblobs(df,['blob_kT','blob_Si'],minvals=[0.5,1.0],
#		      maxvals=[1.0,5.0])

# write filtered dataframe to file
df.to_csv('deconvolution_merged_filtered.txt',sep='\t')

# plot weighted histograms of parameters
histfigs = xplt.histogram_grid(df[df.columns[:-1]],weights=df['blob_em'])

# make scatter plots of specified parameters
scatterfigs = xplt.traceplots(df,columns=['blob_kT','blob_tau','blob_norm','blob_lnsigma'])

# make a map of one of the parameters (saved as fits file)
image = xm.make_map(df,paramname='blob_kT',paramweights='blob_em',binsize=10.0)
