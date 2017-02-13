"""
Module of functions for dealing with files output from xmc

Contains the following functions:

 parse_file_line
 read_parnames

 functions using pandas:
  merge_output 
  remove_nans
  add_header
  calculate_prior
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
    Author: Kari A. Frank 
    Date: March 20, 2014
    Purpose: Read in parameters.txt file from an xmc run to get 
              names of columns in deconvolution files.

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
def merge_output(runpath='./',filetype='deconvolution',save=True,sep='\t'):
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

    Output:
      - if save=True, writes a text file <filetype>_merged.txt into the 
        runpath directory
      - returns a pandas dataframe containing data from all the input files

    Usage Notes:
      - overwrites the file if it already exists
      - automatically includes all iterations present
      - iterations are not written in order, but an extra column is added
        to specify which iteration each row is associated with
      - the created file can be read into a dataframe using
        datatable = pd.read_table(mergedfile,sep='\t',index_col=0)

    """

    # - get list of files - 
    filelist = fu.ls_to_list(runpath,ls_args = filetype+'.*')
    
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
