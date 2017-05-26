#### Common Parameters
distkpc = 50.0 # distance to object in kpc
outroot = './' # output location for new files

## Histogram Parameters
nbins=75
w = 500
h = 200

## Map Parameters
# start with center coords from start.xmc file and rot=0 if don't know
x0=-140.
y0=70.
r0=42.0 # radius of significant object
rotation=0.
z = 0.000927
pixelsize = 5.0 # should be size of typical blob (see histograms)
mapsize = 150.0 # should include the entire phi/psi range

def nHkTthresh(df):
#    import numpy as np
    #return df[df.blob_kT >= 0.28*np.log(df.blob_nH+0.03)+0.23]
    print "Emission lines detectable from all kT-nH in range."
    print "Not filtering any blobs."
    return df
