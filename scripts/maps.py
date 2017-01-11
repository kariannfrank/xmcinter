### Initialize Session ###

# copy paste this at beginning of EVERY session
%run ../../diagnostic_init.py
%load_ext autoreload
%autoreload 2

# read in cleaned blobs
df = pd.read_table('deconvolution_merged_iter8000-16112_cleaned.txt',index_col=0,sep=r'\s+',engine='python',comment='#')

# set map parameters
distkpc = 50.0 # distance to object in kpc
outroot = './' # output location for new files

## Histogram Parameters
nbins=75
w = 500
h = 200

## Map Parameters
# start with center coords from start.xmc file and rot=0 if don't know
z = 0.000927
x0=-140.
y0=70.
rotation=0.
r0=150. # size of (significant) object
pixelsize = 1.0 # should be size of typical blob (see histograms)
mapsize = 150.0 # should include the entire phi/psi range

print "Total Number of Blobs = ",len(df.index)

# (re)define blobcols to include new columns
blobcols = [c for c in df.columns if 'blob' in c]

itermin = np.min(df.iteration)
itermax = np.max(df.iteration)
niter = itermax-itermin
smin = itermin/100
smax = itermax/100

# define columns to map
mapcols=list(blobcols)
#mapcols.remove('blob_norm')
#mapcols.remove('blob_sigma')
mapcols.remove('blob_phi')
mapcols.remove('blob_psi')
mapcols.remove('blob_volume')
mapcols.remove('blob_lnsigma')
mapcols.remove('blob_frac')

# set weights for each column
pweights = ['blob_em']*len(mapcols)
# unweight em and density
pweights[mapcols.index('blob_em')]=None
pweights[mapcols.index('blob_mass')]=None

# set iteration combination types
itypes=['median']*len(mapcols)
itypes[mapcols.index('blob_em')]='total'
itypes[mapcols.index('blob_mass')]='total'

#### Basic Emission Measure Map (to set threshold and check map size)
img1file = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_allblobs'
imgs = xm.make_map(df,paramname='blob_em',
                   paramweights=None,iteration_type='total',
                   binsize=pixelsize,nlayers=20,imagesize=mapsize,
                   withsignificance=True,nproc=4,
                   outfile=img1file,x0=x0,y0=y0,clobber=True,
                   rotation=rotation)

#### Map all parameters
imgfile = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_cleaned'
imgs=xm.make_map(df,paramname=mapcols,paramweights=pweights,iteration_type=itypes,binsize=pixelsize,nlayers=20,imagesize=mapsize,withsignificance=True,nproc=4,outfile=imgfile,x0=x0,y0=y0,clobber=True)
