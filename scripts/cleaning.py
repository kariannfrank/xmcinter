############################################################
#  Initialization                                          #
#  (run at beginning of every session)                     #
############################################################

#### Imports
%run ~/Python_Programs/xmcinter/scripts/diagnostic_init.py
#automatically reload functions before execution
# if the reloading doesn't happen automatically, just run %autoreload
# to do it manually 
%load_ext autoreload
%autoreload 2 

#### Import common settings and object specific settings (e.g. distance)
# change object name and obsid as appropriate
# if the following file does exist in xmcinter/scripts/ then create it
# by copying, renaming, and modifying another init file.
%run ~/Python_Programs/xmcinter/scripts/kes73_0013340201_init.py 

############################################################
#  Check Run Progress                                      #
#  (skip if already have merged deconvolution file)        #
############################################################

#### Create basic figures with script
dfall,sf = xd.check(runpath='../',itmin=None,display=True)

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
sfig = xplt.standard_spectra(runpath=runpath,display=display,itmin=itermin,
                             itmax=itermax,
                             outfile='spectra_all.html',ylog=True,
                             xlog=False,logbins=None,bins=0.03,
                             lines=True,emissivity_range=(1e-17,1.0),
                             nlines=100,energy_range=(0.87,10.0))

############################################################
#  Filter Data                     #
#  (skip if already have merged deconvolution file)        #
############################################################

#### Initial Filtering and Derived Columns

itermin =  4000 # set minimum converged iteration
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
nblobs = len(dfall.index)
print "Total Number of Blobs = ",len(dfall.index)

# (re)define blobcols to include new columns
blobcols = [c for c in dfall.columns if 'blob' in c]

#### Check spectrum
itermin = dfall.iteration.min()
itermax = dfall.iteration.max()
niter = itermax-itermin
sfig = xplt.standard_spectra(runpath='../',display=True,itmin=itermin,
                             itmax=itermax,
                             outfile='spectra_all.html',ylog=True,
                             xlog=False,logbins=None,bins=0.03,
                             lines=True,emissivity_range=(1e-18,1.0),
                             nlines=1000,energy_range=(0.9,10.0),
                             kT_range=(dfall.blob_kT.min(),dfall.blob_kT.max()))

    
#### Clean by nH-kT
dfgood = nHkTthresh(dfall)
ngoodblobs = len(dfgood.index)
print "Number all blobs = ",nblobs
print "Number good blobs = ",ngoodblobs
print "Fraction good blobs = ",float(ngoodblobs)/float(nblobs)
print "Total EM = ",dfall.blob_em.sum()/niter
print "Total Good EM = ",dfgood.blob_em.sum()/niter

#### Compare Spectra of Good and Bad Blobs

# - read spectra if already made -
goodspec = xf.read_spectra(runpath='../',smin=99999,smax=99999,
                           average=False,
                           bins=0.03,logbins=None)
badspec = xf.read_spectra(runpath='../',smin=99998,smax=89999,
                          average=False,
                          bins=0.03,logbins=None)
# - OR make spectra -
goodspec = xw.make_spectrum(dfgood,runpath='../',suffix='99999',
                            bins=0.03,
                            xlog=False,logbins=None)

badspec = xw.make_spectrum(dfall[~dfall.isin(dfgood)].dropna(),
                           runpath='../',suffix='89999',bins=0.03,
                           xlog=False,logbins=None)

# - read data spectrum -
dataspec = xf.read_spectra(runpath='../',smin=0,smax=0,average=False,
                          bins=0.03,logbins=None)[0]
# - read average (total) model spectrum -
allspec = read_spectra(runpath='../',smin=smin,smax=smax,
                          average=True,bins=0.03)[-1]


# - plot spectra -
sfig = xplt.spectra([dataspec,allspec,goodspec,badspec],
                colors=['black','steelblue','darkolivegreen','firebrick'],
                    labels=['Data','All','Good','Bad'],
                    outfile='spectra_compare.html',
                    ylog=True,xlog=False,lines=True,nlines=100,
                    emissivity_range=(1e-18,1.0),
                    kT_range=(dfgood.blob_kT.min(),dfgood.blob_kT.max()))

#### Compare all and cleaned blob histograms
hfigs = xplt.histogram_grid([dfall,dfgood],weights=[None,None],
                            bins=nbins,ncols=3,norm=True,
                            legends=['All Blobs','Good Blobs'],
                            outfile='histogram_grid_allvsgood_unweighted.html',
                            width=w,height=h,iterations='iteration')

hfigs = xplt.histogram_grid([dfall,dfgood],weights='blob_mass',
                            bins=nbins,ncols=3,norm=True,
                            legends=['All Blobs','Good Blobs'],
                            outfile='histogram_grid_allvsgood_massweighted.html',
                            width=w,height=h,iterations='iteration')


#### Overall Histograms
hfigs = xplt.histogram_grid([dfgood,dfgood],weights=[None,'blob_em'],
                            bins=nbins,ncols=3,norm=True,
                            legends=['Unweighted','EM weighted'],
                            width=w,height=h,iterations='iteration')

#### Trace Plots
efig = xplt.trace(dfall,weights=None)

#### Interactive scatter plots of all parameters
# (will take awhile to plot)
tfigs2 = xplt.scatter_grid(dfall[blobcols],agg=None,sampling=2000,outfile='scatter_grid_allblobs.html')
tfigs2 = xplt.scatter_grid(dfgood[blobcols],agg=None,sampling=2000,outfile='scatter_grid_goodblobs.html')

#### Basic Norm Map
img1file = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_goodblobs'
imgs = xm.make_map(dfgood,paramname='blob_norm',
                   paramweights=None,iteration_type='total',
                   binsize=pixelsize,nlayers=70,imagesize=mapsize,
                   withsignificance=True,nproc=4,
                   outfile=img1file,x0=x0,y0=y0,clobber=True,
                   rotation=rotation)


############################################################
#  Optionally Check Significance of Blob Subsets           #
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
imgbadfile = outroot+'bin'+str(int(pixelsize))+'_'+str(int(mapsize))+'arcsec_badem_allblobs'
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
dfgood['blob_massspatweight'] = dfgood.blob_mass*dfgood.blob_spatweight
hfigs = xplt.histogram_grid([dfgood,dfgood,dfgood,dfgood],columns=blobcols,weights=[None,'blob_mass','blob_spatweight','blob_massspatweight'],bins=nbins,ncols=3,norm=False,legends=['Unweighted','Mass weighted','Spatial weighted','Mass*Spatial weighted'],width=w,height=h)

############################################################
#  Create Spatial Weights to de-weight emission in masked  #
#  regions of the remnant .  Repeat for each masked region #
############################################################

#### If using region from event file, get the wcs coordinates
#### in decimal degrees from the ds9 region, and the radius in 
#### in arcsec. Must copy the atthk.fits file to analysis folder.

## These should be defined in xmcinter/scripts/ object initialization script
## Convert to xmc coords:
#holex = astro.wcs2xmc(244.4311,-51.082313)[0]
#holey = astro.wcs2xmc(244.4311,-51.082313)[1]
#holer = 123.

#### Calculate Fractional Weight and Add Column to dataframe
# (this will take a very long time! probably several hours)
dfgood = xw.filtercircle(dfgood,r0=holer,x0=holex,y0=holey,logic='exclude',
                         fraction=True,regname='blob_hole',
                         use_ctypes=True,parallel=True,nproc=5)

#### Update weights
dfgood['blob_massspatweight'] = dfgood.blob_massspatweight*dfgood.blob_hole_exc_fraction

#### Compare Histograms With and Without the New Weights
hfigs = xplt.histogram_grid([dfgood,dfgood,dfgood,dfgood],columns=blobcols,weights=[None,'blob_mass','blob_massspatweight'],bins=nbins,ncols=3,norm=False,legends=['Unweighted','Mass weighted','Mass*Spatial weighted'],width=w,height=h)


############################################################
#  Save Cleaned and Augmented DataFrame to File            #
############################################################
outfile = ('deconvolution_merged_iter'+str(int(min(dfgood.iteration)))+'-'+str(int(max(dfgood.iteration)))+'_cleaned.txt')

# edit header information as appropriate
f = open(outfile,'w+')
f.write('# Cleaning criteria: \n')
#f.write('# EM<'+str(emthresh)+' \n')
#f.write('# '+str(kTthresh1)+'<kT<'+str(kTthresh2)+' keV \n')
f.write('# spatweights: x0= '+str(x0)+', y0='+str(y0)+', r0='+str(r0)+'\n')
dfgood.to_csv(f,sep='\t')
f.close()
