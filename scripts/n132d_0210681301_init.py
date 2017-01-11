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
pixelsize = 5.0 # should be size of typical blob (see histograms)
mapsize = 200.0 # should include the entire phi/psi range
