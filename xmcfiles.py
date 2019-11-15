"""
Module of functions for dealing with files output from xmc

Contains the following functions:

 parse_file_line
 read_parnames
 merge_output 
 remove_nans
 add_header
 calculate_prior
 fake_devonvolution
 read_spectra
"""
#----------------------------------------------------------
#Import common modules
import os
import numpy as np
import pandas as pd
import file_utilities as fu


#----------------------------------------------------------
def remove_nans(datatable,filename=None,verbose=1):
    """Remove rows containing nans from dataframe and print warning."""
    
    nb1 = len(datatable.index)
    datatable.dropna(inplace=True)
    nb2 = len(datatable.index)
    ndropped = nb1-nb2
    if ndropped <= 1: 
        row=' row'
        was=' was'
    else:
        row=' rows'
        was=' were'
    if verbose > 1: print "number of dropped rows = ",ndropped
    if (ndropped > 0) and (verbose > 0): 
        if filename == None: 
            print "Warning: "+str(nb1-nb2)+row+" contained nan "+\
                   "values and"+was+" dropped."
        else:
            print "Warning: "+str(nb1-nb2)+row+" in "+filename+\
                   " contained nan values and"+was+" dropped."
    return datatable

#----------------------------------------------------------
def read_parnames(runpath):
    """
    Read in parameters.txt file from an xmc run to get names of columns in deconvolution files

    Author: Kari A. Frank 
    Date: March 20, 2014

    Usage: parnames = read_parnames(runpath)

    Input:

      runpath -- Path to xmc run directory, containing the file parameters.txt.
                 Only the first line of the file is read.  Parameter names
                 should contain no spaces or punctuation and are assumed to be 
                 separated by commas. They should correspond to the columns
                 of the deconvolution files.



    Output:

      returns string list of column (parameter) names

    Usage Notes:

    """

    #--Set file paths--
    parfile = runpath+'/parameters.txt'

    outpars = fu.parse_file_line(parfile)

    #--return string list--
    return outpars

#----------------------------------------------------------
def merge_output(runpath='./',filetype='deconvolution',save=True,sep='\t',
                 itmin=0,itmax=None):
    """
    Author: Kari A. Frank 
    Date: October 16, 2015
    Purpose: Read in deconvolution, statistic, or other files from 
             an xmc run and merge into a single dataframe, with an extra
             column to specify iteration. optionally can also save as 
             a text file.

    Usage: merge_output(runpath='./',filetype='deconvolution',save=True)

    Input:

      runpath     -- Path to xmc run directory, containing the
                      deconvolution, statistic, etc. files

      filetype     -- string containing the type of xmc file to merge
                      deconvolution, statistic, sigma, mean, changed

      save         -- boolean switch to prevent saving of the merged 
                      file. use save=False if all that is need is the merged
                      dataframe (default=True).

      sep          -- optionally specify the delimiter in the output file
                      (default='\t', tab-separated)

      itmin/itmax  -- optionally specify the minimum and maximum iterations
                      to read

    Output:
      - if save=True, writes a text file <filetype>_merged.txt into the 
        runpath directory
      - returns a pandas dataframe containing data from all the input files

    Usage Notes:
      - overwrites the file if it already exists
      - automatically includes all iterations present if no 
        itmin/itmax given
      - iterations are not written in order, but an extra column is added
        to specify which iteration each row is associated with
      - the created file can be read into a dataframe using
        datatable = pd.read_table(mergedfile,sep='\t',index_col=0)

    """

    # -- defaults --
    if itmax is None:
        itmax = 100000. # dummy itmax value
    
    # - get list of files - 
#    filelist_raw = fu.ls_to_list(runpath,ls_args = filetype+'.*')
    filelist_raw = fu.files_to_list(runpath,search_str = filetype+'.*')

    # - remove less than itmin, greater than itmax -
    filelist = [f for f in filelist_raw if (int(f.split('.')[-1]) >= itmin and int(f.split('.')[-1]) <= itmax)]

    # --Initialize dataframe with first file --
    
    # - read file - 
    datatable = pd.read_table(runpath+'/'+filelist[0],sep='\s+',header=None)

    # - check for nans and drop them -
    datatable = remove_nans(datatable,filename=filelist[0])

    # - add headers -
    parnames = []
    # read header information
    if filetype == 'statistic':
        parnames =  ['stat1','dof','oversim','nblobs','alpha','chi2']
        datatable.columns = parnames
    if filetype == 'deconvolution':
        if os.path.isfile(runpath+'/parameters.txt'):
            parnames = fu.parse_file_line(runpath+'/parameters.txt') 
            datatable.columns = parnames

    # - get iteration number and add as column - 
    iternum = int(filelist[0].split('.')[-1]) # returns integer
    datatable['iteration'] = iternum

    # -- Loop through files and concatenate into single dataframe --
    for f in filelist[1:]:        
        if os.stat(runpath+'/'+f).st_size > 0:
            newframe = pd.read_table(runpath+'/'+f,sep='\s+',
                                     header=None)
            iternum = int(f.split('.')[-1]) # returns integer
            if len(parnames) > 0: newframe.columns=parnames
            newframe['iteration'] = iternum
            # - check for nans and drop them -
            newframe = remove_nans(newframe,filename=f)
            datatable = pd.concat([datatable,newframe],
                                  ignore_index=True)
            if iternum%500 == 0:
                print 'Read iteration '+str(iternum)
        else:
            print 'Warning: '+f+' is missing or empty. Skipping.'
    # -- Write to file --
    if save == True: datatable.to_csv(runpath+'/'+filetype+'_merged.txt',
                                      sep=sep)

    return datatable

#----------------------------------------------------------
def calculate_prior(minval,maxval):
    """Given endpoints of a range, returns the center and stepsize."""

    width = maxval-minval
    stepsize = width/2.0
    center = minval+stepsize

    return center,stepsize

#----------------------------------------------------------
def fake_deconvolution(df,suffix='99999',runpath='../'):
    """
    Write a deconvolution file from a dataframe

    Author: Kari A. Frank 
    Date: April 27, 2016

    Input:
    
      df           -- dataframe
    
      suffix       -- string to attach to end of file names

      runpath      -- Path to xmc run directory, containing the
                      parameters.txt file


    Output:
      - writes a text file with blob parameters formatted as a 
        deconvolution file, and corresponding empty sigma, statistic,
        changed, and mean files.
      - returns the new dataframe

    Usage Notes:
      - overwrites the files if they already exist
      - requires a parameters.txt file containing the deconvolution file
        column names in the correct order

    """

    #--read parameters.txt--
    pnames = read_parnames(runpath)
    
    #--create dataframe with only parameters.txt columns, in correct order--
    dfout = df[pnames].copy()

    #--check if file exists--
    clobber = 'y'
    if os.path.isfile('deconvolution.'+suffix):
        clobber = raw_input("File deconvolution."+suffix+" already exists. Overwrite file? (y/n)  ")
    if clobber in ['y','yes','Yes','Y']:
        clobber = True
    else:
        clobber = False
        
    #--write dataframe to file--
    if clobber:
        dfout.to_csv('deconvolution.'+suffix,sep=' ',header=False,
                     index=False)

    #--create auxiliary (empty) files--
    open('statistic.'+suffix,'w').close()
    open('sigma.'+suffix,'w').close()
    open('mean.'+suffix,'w').close()
    open('changed.'+suffix,'w').close()
    
    return dfout
#----------------------------------------------------------
def read_spectra(runpath='../',itmin=1,itmax=None,logbins=False,
                 legacy=False,
                 average=False,bins=0.03,datarange=None,oversim='auto'):
    """
    Read xmc spectrum files into histograms.

    Author: Kari A. Frank 
    Date: May 3, 2017

    Input:
    
     runpath (string) : the relative path to xmc run folder, which
                        contains the spectrum.* files. alternatively,
                        can pass a dataframe with the x and y values, e.g.
                        as would be read from a text file written by
                        xspec with the wd command.

     itmin/itmax (int) : minimum and/or maximum iterations included in the
                       averaging of the model spectra. corresponds to the 
                       spectrum* file names, e.g. itmax=3 will average the 
                       files spectrum_1.fits, spectrum_2.fits, and 
                       spectrum_3.fits. default is all available spectra.

     bins:        optionally specify the number of bins (int), binsize 
                  (float), or an array containing the bin edge values. 
                  can also accept any of the strings recognized by 
                  np.histogram() for special calculation of bin sizes.
                  - 'auto' (but don't use if using weights)
                  - 'fd' (Freedman Diaconis Estimator)
                  - 'doane'
                  - 'scott'
                  - 'rice' 
                  - 'sturges'
                  - 'sqrt'
                  default is to create equally spaced bins of width 0.03,
                  (30 eV). note that xmc typically uses a binsize of 0.015.

     average (bool) : if True, will return the average of the spectra 
                      between itmin and itmax, as the last tuple in the 
                      returned list.

     datarange () : passed directly to make_histogram()

     oversim (numerical scalar or 'auto') : must specify the oversim value
                   to properly scale the model spectra.
                   - 'auto' (or any string): will read in the first 
                      statistic.* file it
                     finds in runpath/ and get the oversim from there.
                   - can also explicitly pass the oversim value

     legacy (bool) : assume old-style spectrum file names, i.e.
                     spectrum_iter/100.fits

    Output:
      - returns a tuple containing the histogram information, 
        (y,yerrors,yedges) for each iteration read, in a list of tuples

    Usage Notes:
     - to read in a single spectrum set itmin=itmax=iteration number
     - to read in the data spectrum, set itmin=0 and itmax=0

    """

    #----Import Modules----
    import os,re
    import astropy.io.fits as fits
    from file_utilities import ls_to_list,parse_file_line
    from wrangle import make_histogram

    #----Set defaults----
        
    # check for MPI file names (if xmc was run with mpi)
    if os.path.isfile(runpath+'/spectrum0_0.fits'):
        specname = 'spectrum0_'
    else:
        specname = 'spectrum_'

    # get oversim
    if isinstance(oversim,str):
        oversim = float(parse_file_line(runpath+ls_to_list(runpath,'statistic.*')[0])[2])

    # set file name type
    if legacy is True:
        itmin = itmin/100
        if itmax is not None: itmax = itmax/100

    if itmax is None:
#        itmax = len(ls_to_list(runpath,'spectrum*')) - 1
        itmax = int(re.search(r'_(.*)\.fits',
                             ls_to_list(runpath,'-tr spectrum*')[-1]).group(1))
                
    #----Read in first model spectrum----
    foundspec = False
    sm = itmin
    hists = []
    print itmin,itmax
    while (sm<=itmax and foundspec is False):
        specfile = runpath+'/'+specname+str(sm)+'.fits'
        print specfile
        if os.path.isfile(specfile):
            table = fits.getdata(specfile,0)
            wave = table.field('wave')
            wave_avg = table.field('wave')
            iters_avg = np.ones_like(wave)
            iters_avg.fill(sm)
            foundspec = True

            #--convert pd.Series--
            if np.max(wave) < 15.0:
                xlabel = 'Energy (keV)'
            else:
                xlabel = 'Wavelength (Angstroms)'
 
            if sm == 0:
                wave = pd.Series(wave.byteswap().newbyteorder(),
                                 name=xlabel)
            else:
                wave = pd.Series(wave+0,name=xlabel)

            #--Create and Save Histogram--
            # y units are counts per bin
            if isinstance(bins,float): # assume binsize, convert
                                       # to number bins
                nbins = np.ceil((np.max(wave.values)-
                        np.min(wave.values))/bins)    

            y,yerrors,yedges = \
                       make_histogram(wave,bins=nbins,
                       logbins=logbins,
                       datarange=datarange,
                       density=False,iterations=None)

            hists = hists + [(y,yerrors,yedges)]

        else:
#            print ("Warning: "+specfile+" not found.  Skipping +"
#                   "to next spectrum.")
            sm = sm+1

        if foundspec is False:
            #print sm
            print "ERROR: no spectrum files found in range."
#            print "WARNING: No spectrum files found in range. Plotting data spectrum only."
            
    #----Loop over remaining spectra----
    for s in xrange(sm+1,itmax+1):
        specfile = runpath+'/'+specname+str(s)+'.fits'
        if os.path.isfile(specfile):
            table = fits.getdata(specfile,0)
            wave = table.field('wave')
            iters = np.ones_like(wave)
            iters.fill(s) #assign 'iteration' number
            wave_avg = np.hstack((wave_avg,wave))
            iters_avg = np.hstack((iters_avg,iters))

            #--convert pd.Series--
            if np.max(wave) < 15.0:
                xlabel = 'Energy (keV)'
            else:
                xlabel = 'Wavelength (Angstroms)'
            wave = pd.Series(wave+0,name=xlabel)

            #--Create and Save Histogram--
            # y units are counts per bin
            if isinstance(bins,float): # assume binsize, convert
                                       # to number bins
                nbins = np.ceil((np.max(wave.values)-
                        np.min(wave.values))/bins)    

            y,yerrors,yedges = \
                       make_histogram(wave,bins=nbins,
                       logbins=logbins,
                       datarange=datarange,
                       density=False,iterations=None)
            hists = hists + [(y/oversim,yerrors,yedges)]
            
#        else:
#            print "Warning: "+specfile+" does not exist. Skipping."

    #----Save average----
    if average is True and foundspec is True:
        # --convert average to pd.Series--
        wave_avg = pd.Series(wave_avg+0,name=xlabel)
        iters_avg = pd.Series(iters_avg+0,name='iteration')

        #--Create and Save Histogram--
        # y units are counts per bin
        if isinstance(bins,float): # assume binsize, convert
            # to number bins
            nbins = np.ceil((np.max(wave_avg.values)-
                            np.min(wave_avg.values))/bins)    
                
        y,yerrors,yedges = make_histogram(wave_avg,bins=nbins,
                                          logbins=logbins,
                                          datarange=datarange,
                                          density=False,
                                          iterations=iters_avg)

        #-scale by number of iterations-
        nspec = len(np.unique(iters_avg))
        print 'nspec = ',nspec
        y = y/float(nspec)/oversim
        yerrors = yerrors/float(nspec)/oversim
        
        #-save-
        hists = hists + [(y,yerrors,yedges)]
 
    return hists
