############################################################
#  Initialization                                          #
#  (run at beginning of every session)                     #
############################################################

## Assumes cleaning.py has been completed

#### Imports
%run ~/Python_Programs/xmcinter/scripts/diagnostic_init.py
%load_ext autoreload
%autoreload 2

#### Import common settings and object specific settings (e.g. distance)
# change object name and obsid as appropriate
# if the following file does exist in xmcinter/scripts/ then create it
# by copying, renaming, and modifying another init file.
%run ~/Python_Programs/xmcinter/scripts/rcw103_0302390101_init.py

############################################################
#  Read Merged, Cleaned Deconvolution File                 #
############################################################

#### Read File
# edit filename as appropriate
df = pd.read_table('deconvolution_merged_iter3500-9292.txt',index_col=0,sep=r'\s+',engine='python',comment='#')

print "Total Number of Blobs = ",len(df.index)

# (re)define blobcols to include new columns
blobcols = [c for c in df.columns if 'blob' in c]

itermin = np.min(df.iteration)
itermax = np.max(df.iteration)
niter = itermax-itermin
smin = itermin/100
smax = itermax/100

############################################################
#  Set up map parameters                                   #
############################################################

#### define columns to map
mapcols=list(blobcols)
#mapcols.remove('blob_norm')
#mapcols.remove('blob_sigma')
mapcols.remove('blob_phi')
mapcols.remove('blob_psi')
mapcols.remove('blob_volume')
mapcols.remove('blob_lnsigma')
mapcols.remove('blob_frac')

#### set weights for each column
pweights = ['blob_em']*len(mapcols)
# unweight em and density
pweights[mapcols.index('blob_em')]=None
pweights[mapcols.index('blob_mass')]=None

#### set iteration combination types
itypes=['median']*len(mapcols)
itypes[mapcols.index('blob_em')]='total'
itypes[mapcols.index('blob_mass')]='total'

############################################################
#  Basic Emission Measure Map                              #
#  (to set threshold and check map size)                   #
############################################################

#### Set file location and name prefix
img1file = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_allblobs'

#### Make map
imgs = xm.make_map(df,paramname='blob_em',
                   paramweights=None,iteration_type='total',
                   binsize=pixelsize,nlayers=20,imagesize=mapsize,
                   withsignificance=True,nproc=4,
                   outfile=img1file,x0=x0,y0=y0,clobber=True,
                   rotation=rotation)

#### Use EM map to determine EM threshold
emthresh = 


############################################################
#  Map all parameters                                      #
############################################################

#### Set file location and name prefix
imgfile = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_cleaned'

#### Make maps
imgs=xm.make_map(df,paramname=mapcols,paramweights=pweights,iteration_type=itypes,binsize=pixelsize,nlayers=20,imagesize=mapsize,withsignificance=True,nproc=4,outfile=imgfile,x0=x0,y0=y0,clobber=True,imgthresh=emthresh,imgthreshparam='blob_em')
