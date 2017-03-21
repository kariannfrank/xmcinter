#### Common Parameters
distkpc = 8.5 # distance to object in kpc
outroot = './' # output location for new files

## Histogram Parameters
nbins=75
w = 500
h = 200

## Map Parameters
# start with center coords from start.xmc file and rot=0 if don't know
x0=-60.
y0=80.
rotation=-90.
r0=145. # size of (significant) object
pixelsize = 3.0 # should be size of typical blob (see histograms)
mapsize = 300.0 # should include the entire phi/psi range

import xmcinter.astro_utilities as astro
holex = astro.wcs2xmc(280.33033,-4.9366564)[0]
holey = astro.wcs2xmc(280.33033,-4.9366564)[1]
holer = 25.0

def nHkTthresh(df):
    import numpy as np
    return df[df.blob_kT >= 0.46*np.log(df.blob_nH+0.95)-0.02]
