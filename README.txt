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

import pandas as pd
import numpy as np
import xmcinter.plots as xplt
import xmcinter.xmcfiles as xf
import xmcinter.wrangle as xw
import xmcinter.astro_utilities as astro
import xmcinter.diagnostics as xd
import xmcinter.xmcmap as xm

You can copy any of the following import statements and paste into ipython
 with %paste, or just run the provided initialization script which will 
 import the above:
 
%run diagnostic_init.py

###################################################################
Usage Notes:

- To use the fastest (and default) method for making maps 
  (via xmcinter.xmcmap.make_map), it is necessary to compile the gaussian.c
  function on your local machine, e.g. on linux:

  gcc -shared -o gaussian.so -fPIC gaussian.c

  If you choose not to use this, make sure to always pass the cint=False
  argument to make_map(). It should only run slightly slower.

- If you encounter an error that looks like the following, the cause is 
  likely an incorrect parameters.txt file

  ValueError: Length mismatch: Expected axis has 23 elements, new values have 22 elements

###################################################################
A typical workflow for early diagnostics
- requires access to the original deconvolution.* and statistic.* files
- assumes running from the <run_directory>/analysis/ subfolder, one level
  below the directory containing the xmc output files
- for more examples, including using with jupyter notebooks, see the 
  scripts in the scripts/ directory

# check progress (makes some basic diagnostic figures, including a 
#  current EM map, and prints basic info)
df,sf = xd.check(runpath='../',itmin=500)

# or check manually:
# check median chi2 after convergence
# (change 1000 to first converged iteration)
sf = xplt.chi2(runpath='../')
sf = xw.filterblobs(sf,'iteration',minvals=1000)
medchi2 = xw.weighted_median(sf['redchi2'])
print medchi2

# filter by iteration and add emission measure, blob_sigma 
# (change itmin and distance)
df = xd.clean(itmin=1000,distance=3.3)

# check traceplots -- the plotting can take awhile
tracefigs = xplt.traceplots(df)

# check weighted histograms
histfigs = xplt.histogram_grid(df,columns=df.columns[:-1],weights='blob_em')

###################################################################
Examples of some other common tasks
(all executed from ipython)

# read in a typical merged deconvolution file
df = pd.read_table('deconvolution_merged.txt',sep='\t',index_col=0)

# filter by some parameter (e.g. temperature)
dfnew = df[df['blob_kT'] < 1.0]
dfnew = df[(df['blob_kT'] < 1.0) & (df['blob_tau']>5e11)]
dfnew = df[~((df['blob_kT'] < 1.0) & (df['blob_tau']>5e11))]

# write filtered dataframe to file
df.to_csv('deconvolution_merged_filtered.txt',sep='\t')

# plot weighted histograms of parameters
histfigs = xplt.histogram_grid(df[df.columns[:-1]],weights=df['blob_em'])

# make scatter plots of specified parameters
scatterfigs = xplt.scatter_grid(df,columns=['blob_kT','blob_tau','blob_norm','blob_lnsigma'])

# make a map of one of the parameters (saved as fits file)
image = xm.make_map(df,paramname='blob_kT',paramweights='blob_em',binsize=10.0,withsignificance=True)
image = xm.make_map(df,paramname='blob_em',paramweights=None,binsize=10.0,withsignificance=True,iteration_type='total')

# make a map of multiple parameters, suppress pixels with faintest emission, rotate output images by -270deg
images = xm.make_map(df,paramname=['blob_em''blob_kT','blob_nH'],paramweights=[None,'blob_em','blob_em'],binsize=10.0,withsignificance=True,iteration_type=['total','median','median'],imagesize=8.0*60.0,nproc=5,imgstdthresh=4.2,imgstdtreshparam='blob_em',rotation=-270,clobber=True)

# calculate the [X/Y] abundance ratio (where Y can also be Si, O, Mg, and 
# X can be any free abundance)
df['blob_SFe']=astro.xspec_abund_to_nomoto_dataframe(df,'blob_S','blob_Fe')
