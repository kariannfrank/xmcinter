#### Common Parameters
distkpc = 60.0 # distance to object in kpc
outroot = './' # output location for new files

## Histogram Parameters
nbins=75
w = 500
h = 200

## Map Parameters
# start with center coords from start.xmc file and rot=0 if don't know
x0=-75
y0=-5.
z = 0.000527
rotation=0.
pixelsize = 1.0 # should be size of typical blob (see histograms)
mapsize = 100.0 # should include the entire phi/psi range

def nHkTthresh(df,nHgal=0.046):
    import numpy as np
#    print "Emission lines detectable from all kT-nH in range."
#    print "Not filtering any blobs."
    return df[df.blob_kT >= 0.04*np.log(nHgal+df.blob_nH2+0.01)+0.155]

