#### Common Parameters
distkpc = 4.4 # distance to object in kpc
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
rotation=0.0 #270.
r0=140. # size of (significant) object
pixelsize = 5.0 # should be size of typical blob (see histograms)
mapsize = 400.0 # should include the entire phi/psi range

# Central Pulsar
# must copy the atthk.fits file into the analysis directory
import xmcinter.astro_utilities as astro
holex = astro.wcs2xmc(272.8716327,-19.42388673)[0] # RA, DEC in degrees
holey = astro.wcs2xmc(272.8716327,-19.42388673)[1]
holer = 10.972 # arcsec

# PWN Region, approximate (used for powerlaw spatial prior)
# must copy the atthk.fits file into the analysis directory
import xmcinter.astro_utilities as astro
pwnx = astro.wcs2xmc(272.8716327,-19.42388673)[0] # RA, DEC in degrees
pwny = astro.wcs2xmc(272.8716327,-19.42388673)[1]
pwnr = 45.0 # arcsec

def nHkTthresh(df):
    import numpy as np
    return df[df.blob_kT >= 0.3*np.log(df.blob_nH+0.2)+0.21]
