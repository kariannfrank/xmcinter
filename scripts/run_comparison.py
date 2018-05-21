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
%run ~/Python_Programs/xmcinter/scripts/ctb109_0038140101_init.py 

############################################################
#  Read Merged Deconvolution File                          #
############################################################

#### Read Files
# edit filename as appropriate
dfa = pd.read_table('deconvolution_merged_iter2000-17063.txt',index_col=0,sep=r'\s+',engine='python',comment='#')
dfb = pd.read_table('/data/drive1/staging/CTB109_0038140101_b_b100_EPIC/analysis/deconvolution_merged_iter3000-17691.txt',index_col=0,sep=r'\s+',engine='python',comment='#')


############################################################
#  Basic Inspection                                        #
############################################################

#### Check number of blobs
nblobs = len(dfa.index)
print "Total Number of Blobs in A = ",len(dfa.index)

nblobs = len(dfb.index)
print "Total Number of Blobs in B = ",len(dfb.index)

# (re)define blobcols to include new columns
blobcolsa = [c for c in dfa.columns if 'blob' in c]

#### Compare histograms
w = 500
h = 200

hfigs = xplt.histogram_grid([dfa,dfb],weights=[None,None],
                            bins=nbins,ncols=2,norm=True,
                            legends=['A','B'],
                            outfile='histogram_grid_unweighted_ab100_b100.html',
                            width=w,height=h,iterations='iteration')

hfigs = xplt.histogram_grid([dfa,dfb],weights='blob_em',
                            bins=nbins,ncols=2,norm=True,
                            legends=['A','B'],
                            outfile='histogram_grid_AB_emweighted.html',
                            width=w,height=h,iterations='iteration')


#### Trace Plots
efig = xplt.trace(dfall,weights=None)

#### Interactive scatter plots of all parameters
# (will take awhile to plot)
tfigs2 = xplt.scatter_grid(dfall[blobcols],agg=None,sampling=2000,outfile='scatter_grid_allblobs.html')
tfigs2 = xplt.scatter_grid(dfgood[blobcols],agg=None,sampling=2000,outfile='scatter_grid_goodblobs.html')
