#### Common Parameters
distkpc = 8.0 # distance to object in kpc
outroot = './' # output location for new files

## Histogram Parameters
nbins=75
w = 500
h = 200

## Map Parameters
# start with center coords from start.xmc file and rot=0 if don't know
x0=-60.
y0=80.
z = 0.0
rotation=270.
r0=100. # size of (significant) object
pixelsize = 5.0 # should be size of typical blob (see histograms)
mapsize = 200.0 # should include the entire phi/psi range

def nHkTthresh(df):
    import numpy as np
    return df[df.blob_kT >= 0.3*np.log(df.blob_nH)+0.11]
