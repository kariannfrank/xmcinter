############################################################
#  Initialization                                          #
#  (run at beginning of every session)                     #
############################################################

#### Imports
%run diagnostic_init.py
#automatically reload functions before execution
# if the reloading doesn't happen automatically, just run %autoreload
# to do it manually 
%load_ext autoreload
%autoreload 2 

#### Import common settings and object specific settings (e.g. distance)
# change object name and obsid as appropriate
# if the following file does exist in xmcinter/scripts/ then create it
# by copying, renaming, and modifying another init file.
%run rcw103_0302390101_init.py

############################################################
#  Check Run Progress                                      #
#  (skip if already have merged deconvolution file)        #
############################################################

#### Create basic figures with script
dfall,sf = xd.check(runpath='../')

#### OR Create figures manually

#### Check Chi2
sf=xplt.chi2(runpath='../')

#### Check traces
dfall=xf.merge_output('../',save=False)
efig = xplt.trace(dfall,weights=None)

#### Decide minimum converged iteration (optionally maximum)
itermin = 3500
itermax = None
# determine converged chi2/dof
sf = xw.filterblobs(sf,'iteration',minvals=itermin,maxvals=itermax)
print xw.weighted_median(sf['redchi2']) 
print max(sf.iteration)

#### Check spectrum
smin = itermin/100
if itermax is None: 
    smax = None
else:
    smax = itermax/100
sfig = xplt.spectrum(runpath='../',smin=smin,smax=smax,bins=0.03,
                     ylog=True,xlog=True,
                     lines=True,nlines=50)

############################################################
#  Filter Data                     #
#  (skip if already have merged deconvolution file)        #
############################################################

#### Initial Filtering and Derived Columns

itermin =  # set minimum converged iteration
itermax =  None # set maximum iteration to use

# remove iterations before convergence, add new derived columns, 
#  and save merged deconvolution file
dfall = xd.clean(dfall,itmin=itermin,distance=distkpc,itmax=itermax)

############################################################
#  Read Merged Deconvolution File                          #
############################################################

#### Read File
# edit filename as appropriate
dfall = pd.read_table('deconvolution_merged_iter3500-9292.txt',index_col=0,sep=r'\s+',engine='python',comment='#')


############################################################
#  Basic Inspection                                        #
############################################################

#### Check number of blobs
print "Total Number of Blobs = ",len(dfall.index)

# (re)define blobcols to include new columns
blobcols = [c for c in dfall.columns if 'blob' in c]

#### Check spectrum
itermin = dfall.iteration.min()
itermax = dfall.iteration.max()
if itermax is None: 
    smax = None
else:
    smax = itermax/100
sfig = xplt.spectrum(runpath='../',smin=smin,smax=smax,
                     ylog=True,xlog=True,
                     lines=True,nlines=50,
                     kT_range=(dfall.blob_kT.min(),dfall.blob_kT.max()))

#### Overall Histograms
hfigs = xplt.histogram_grid([dfall,dfall],weights=[None,'blob_em'],
                            bins=nbins,ncols=2,norm=True,
                            legends=['Unweighted','EM weighted'],
                            width=w,height=h,iterations='iteration')

#### Trace Plots
efig = xplt.trace(dfall,weights=None)

#### Interactive scatter plots of all parameters
# (will take awhile to plot)
tfigs2 = xplt.scatter_grid(dfall[blobcols],agg=None,sampling=2000)

#### Define Map Parameters

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

#### Basic Emission Measure Map
img1file = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_allblobs'
imgs = xm.make_map(dfall,paramname='blob_em',
                   paramweights=None,iteration_type='total',
                   binsize=pixelsize,nlayers=70,imagesize=mapsize,
                   withsignificance=True,nproc=4,
                   outfile=img1file,x0=x0,y0=y0,clobber=True,
                   rotation=rotation)


############################################################
#  Remove Insignificant High EM Blobs                      #
############################################################
# Use the histograms and scatter plots to determine good guesses
# for an EM threshold. Then create EM significance maps
# for the 'bad' blobs to see if the emission from those blobs was 
# significant.

# Repeat this section, lowering the EM until significant emission
# starts to show up (significance>1sigma). Set that EM as the final
# final threshold below.

#### Stretched out EM histogram to guess threshold
hfigs = xplt.histogram(dfall['blob_em'],weights=None,bins=4*nbins,
                       norm=False,width=w*2,height=2*h)

#### Define Shortcut Functions to Filter by EM
emthresh = 5e59
def dflowem(df,emthresh):
    return xw.filterblobs(df,'blob_em',maxvals=emthresh)
def dfhighem(df,emthresh):
    return xw.filterblobs(df,'blob_em',minvals=emthresh)

#### Create EM Significance Maps of High EM Blobs
imgbadfile = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_badem'
imgs = xm.make_map(dfhighem(dfall,emthresh),paramname='blob_em',
                   paramweights=None,iteration_type='total',
                   binsize=pixelsize,nlayers=70,imagesize=mapsize,
                   withsignificance=True,nproc=4,
                   outfile=imgbadfile,x0=x0,y0=y0,clobber=True,
                   rotation=rotation)

# open in ds9
# ds9 -multiframe bin10_700arcsec_*_median_blob_em.fits -cmap rainbow -match colorbars &

#### Set Final EM Threshold and Filter
emthresh = 5e59
dfgood = dflowem(dfall,emthresh)

############################################################
#  Remove Insignificant Low kT Blobs                       #
############################################################
# Use the histograms and scatter plots to determine good guesses
# for a temperature threshold. Then create EM significance maps
# for the 'bad' blobs to see if the emission from those blobs was 
# significant.

# Repeat this section, increasing the threshold until significant emission
# starts to show up (significance>1sigma). Set that kT as the final
# final threshold below.

#### Stretched out kT histogram to guess threshold
hfigs = xplt.histogram_grid([dfgood,dfgood],columns=['blob_kT'],
                            weights=[None,'blob_em'],
                       bins=4*nbins,legends=['unweighted','EM weighted'],
                       norm=False,width=w*2,height=2*h,ncols=1)

#### Define Shortcut Functions to Filter by kT
kTthresh = 0.17
def dflowkT(df,kTthresh):
    return xw.filterblobs(df,'blob_kT',maxvals=kTthresh)
def dfhighkT(df,kTthresh):
    return xw.filterblobs(df,'blob_kT',minvals=kTthresh)

#### Create EM Significance Maps of Low kT Blobs
imgbadfile = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_badkT'
imgs = xm.make_map(dflowkT(dfgood,kTthresh),paramname='blob_em',
                   paramweights=None,iteration_type='total',
                   binsize=pixelsize,nlayers=70,imagesize=mapsize,
                   withsignificance=True,nproc=4,
                   outfile=imgbadfile,x0=x0,y0=y0,clobber=True,
                   rotation=rotation)

# open in ds9
# ds9 -multiframe bin10_700arcsec_*_median_blob_em.fits -cmap rainbow -match colorbars &

#### Set Final kT Threshold and Filter df
kTthresh1 = 0.17
dfgood = dfhighkT(dfgood,kTthresh1)

############################################################
#  Remove Insignificant High kT Blobs                       #
############################################################
# Use the histograms and scatter plots to determine good guesses
# for a temperature threshold. Then create EM significance maps
# for the 'bad' blobs to see if the emission from those blobs was 
# significant.

# Repeat this section, decreasing the threshold until significant emission
# starts to show up (significance>1sigma). Set that kT as the final
# final threshold below.

#### Stretched out kT histogram to guess threshold
hfigs = xplt.histogram(dfall['blob_em'],weights=None,bins=4*nbins,
                       norm=False,width=w*2,height=2*h)

kTthresh=4.0

#### Create EM Significance Maps of High kT Blobs
imgbadfile = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_badkT'
imgs = xm.make_map(dfhighkT(dfgood,kTthresh),paramname='blob_em',
                   paramweights=None,iteration_type='total',
                   binsize=pixelsize,nlayers=70,imagesize=mapsize,
                   withsignificance=True,nproc=4,
                   outfile=imgbadfile,x0=x0,y0=y0,clobber=True,
                   rotation=rotation)

# open in ds9
# ds9 -multiframe bin10_700arcsec_*_median_blob_em.fits -cmap rainbow -match colorbars &

#### Set Final kT Threshold and Filter df
kTthresh2 = 4.0
dfgood = dflowkT(dfgood,kTthresh2)


############################################################
#  Create New Cleaned EM Map                               #
############################################################
imgcleanfile = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_clean'
imgs = xm.make_map(dfgood,paramname='blob_em',
                   paramweights=None,iteration_type='total',
                   binsize=pixelsize,nlayers=70,imagesize=mapsize,
                   withsignificance=True,nproc=4,
                   outfile=imgcleanfile,x0=x0,y0=y0,clobber=True,
                   rotation=rotation)


############################################################
#  Optionally Create Spatial Weights                       #
#  (to de-weight emission outside of the remnant in        #
#   future histograms)                                     #
############################################################
# Make sure x0,y0,r0 are accurate. Check the EM map above (and the 
# corresponding significance map) to verify. The circle defined by
# x0,y0,r0 should be such that it encompasses all of the >1sigma
# significant emission, and as little as possible of everything else.
# Make sure to choose carefully because calculating the weights 
# takes a very long time - you only want to do it once.

#### Calculate Fractional Weight and Add Column to dataframe
# (this will take a very long time! probably several hours)
dfgood = xw.filtercircle(dfgood,r0=r0,x0=x0,y0=y0,logic='include',
                         fraction=True,regname='blob_spatweight',
                         use_ctypes=True,parallel=True,nproc=5)

#### Compare Histograms With and Without the New Weights
dfgood['blob_emspatweight'] = dfgood.blob_em*dfgood.blob_spatweight
hfigs = xplt.histogram_grid([dfgood,dfgood,dfgood,dfgood],columns=blobcols,weights=[None,'blob_em','blob_spatweight','blob_emspatweight'],bins=nbins,ncols=2,norm=False,legends=['Unweighted','EM weighted','Spatial weighted','EM*Spatial weighted'],width=w,height=h)

############################################################
#  Create Spatial Weights to de-weight emission in masked  #
#  regions of the remnant .  Repeat for each masked region #
############################################################

#### If using region from event file, get the wcs coordinates
#### in decimal degrees from the ds9 region, and the radius in 
#### in arcsec. Must copy the atthk.fits file to analysis folder.

## Convert to xmc coords:
holex = astro.wcs2xmc(244.4311,-51.082313)[0]
holey = astro.wcs2xmc(244.4311,-51.082313)[1]
holer = 123.

#### Calculate Fractional Weight and Add Column to dataframe
# (this will take a very long time! probably several hours)
dfgood = xw.filtercircle(dfgood,r0=holer,x0=holex,y0=holey,logic='exclude',
                         fraction=True,regname='blob_hole',
                         use_ctypes=True,parallel=True,nproc=5)

#### Update weights
dfgood['blob_emspatweight'] = dfgood.blob_emspatweight*dfgood.blob_hole_exc_fraction

#### Compare Histograms With and Without the New Weights
hfigs = xplt.histogram_grid([dfgood,dfgood,dfgood,dfgood],columns=blobcols,weights=[None,'blob_em','blob_emspatweight'],bins=nbins,ncols=2,norm=False,legends=['Unweighted','EM weighted','EM*Spatial weighted'],width=w,height=h)


############################################################
#  Save Cleaned and Augmented DataFrame to File            #
############################################################
outfile = ('deconvolution_merged_iter'+str(int(min(dfgood.iteration)))+'-'+str(int(max(dfgood.iteration)))+'_cleaned.txt')

# edit header information as appropriate
f = open(outfile,'w+')
f.write('# Cleaning criteria: \n')
f.write('# EM<'+str(emthresh)+' \n')
f.write('# '+str(kTthresh1)+'<kT<'+str(kTthresh2)+' keV \n')
f.write('# spatweights: x0= '+str(x0)+', y0='+str(y0)+', r0='+str(r0)+'\n')
dfgood.to_csv(f,sep='\t')
f.close()
