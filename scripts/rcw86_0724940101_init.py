#### Common Parameters
distkpc = 2.5 # distance to object in kpc
outroot = './' # output location for new files

## Histogram Parameters
nbins=75
w = 500
h = 200

## Map Parameters
# start with center coords from start.xmc file and rot=0 if don't know
x0=0. # this number is not accurate yet
y0=80.
r0=750. # radius of significant object
rotation=0.
z = 0.0
pixelsize = 10.0 # should be size of typical blob (see histograms)
mapsize = 750.0 # should include the entire phi/psi range

def nHkTthresh(df):
#    import numpy as np
    return df[df.blob_kT >= 0.1*np.log(df.blob_nH+0.2)+0.11+0.025]
    return df
