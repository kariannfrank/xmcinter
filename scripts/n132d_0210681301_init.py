#### Common Parameters
distkpc = 50.0 # distance to object in kpc
outroot = './' # output location for new files

## Histogram Parameters
nbins=75
w = 500
h = 200

## Map Parameters
# start with center coords from start.xmc file and rot=0 if don't know
x0=-70.
y0=-10.
z =  0.000927
rotation=0.
r0=200. # size of (significant) object
pixelsize = 3.0 # should be size of typical blob (see histograms)
mapsize = 200.0 # should include the entire phi/psi range

def nHkTthresh(df,nHgal=0.062):
    # nH needs to be total NH = Galactic + LMC,
    # as was used to calculate the NH-kT curve
    import numpy as np
#    return df[df.blob_kT >= 0.28*np.log(df.blob_nH+0.03)+0.23]
    return df[df.blob_kT >= 0.28*np.log(nHgal+df.blob_nH2+0.03)+0.23]
