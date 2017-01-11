### Initialize Session ###

# copy paste this at beginning of EVERY session
%run ../diagnostic_init.py

# if level 1 cleaning already has been, read in cleaned blobs
dfall = pd.read_table('deconvolution_merged_iter3000-5336.txt',index_col=0,sep=r'\s+',engine='python')

### Checking Progress of Ongoing Run and Level 1 Cleaning ###

# check chi2 
# - determine minimum converged iteration from plot and set below
sf = xplt.chi2()
itermin = 2000

# remove early iterations and record median chi2
sf = xw.filterblobs(sf,'iteration',minvals=itermin)
xw.weighted_median(sf['redchi2']) 
# median chi2/dof = 

# perform level 1 cleaning (also writes cleaned blobs to file)
# - change distance for each object
dfall = xd.clean(itmin=itermin,distance=50.0)

# check number of blobs remaining
len(dfall.index)


### Level 2 Cleaning - Remove high EM/low kT blobs ###
# (repeat as needed with different threshold values until the removed blobs
#  are not significant and there are no questionable blobs dominated the 
#  weighted histograms)

# plot histograms and scatter plot of blob_em, blob_kT
hfigs = xplt.histogram_grid(dfall[['blob_kT','blob_em']],bins=1000)
whfigs = xplt.histogram_grid(dfall[['blob_kT','blob_em']],bins=1000,weights=dfall['blob_em'])
sfig = xplt.scatter(dfall,'blob_kT','blob_em',size=1,npoints=10000)

# from plots, guess threshold values for kT and EM and separate out 
# the 'bad' blobs. change x0,y0, and imagesize (and possibly binsize for
# large objects) as appropriate for that object and xmc run.
dfbad = xw.filterblobs(dfall,'blob_kT',maxvals=0.17)
dfbad = xw.filterblobs(dfbad,'blob_em',minvals=0.8e59)
imgs=xm.make_map(dfbad,paramname='blob_em',paramweights=None,iteration_type='total',binsize=2.0,nlayers=20,imagesize=160.,withsignificance=True,nproc=4,outfile=outroot+'bin2_160arcsec_bad',x0=-60.,y0=-10.,clobber=True)
