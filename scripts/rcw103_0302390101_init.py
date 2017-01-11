#### Common Parameters
distkpc = 3.3 # distance to object in kpc
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
r0=340. # size of (significant) object
pixelsize = 10.0 # should be size of typical blob (see histograms)
mapsize = 700.0 # should include the entire phi/psi range

holex = astro.wcs2xmc(244.40073,-51.04045)[0]
holey = astro.wcs2xmc(244.40073,-51.04045)[1]
holer = 8.0
