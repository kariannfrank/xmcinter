#----------------------------------------------------------
#Module of functions for dealing with files output from xmc
#
#Contains the following functions:
#
# parse_file_line
# read_parnames
# get_xmcoutput

# functions using pandas:
#  merge_output -- replaces get_deconvolution, get_statistic, and get_xmcoutput
#
#----------------------------------------------------------
#Import common modules
import os
import numpy as np
import pandas as pd
import utilities

#----------------------------------------------------------
# read_parnames()
#
#Author: Kari A. Frank 
#Date: March 20, 2014
#Purpose: Read in parameters.txt file from an xmc run to get 
#          names of columns in deconvolution files.
#
#Usage: parnames = read_parnames(runpath)
#
#Input:
#
#  runpath -- Path to xmc run directory, containing the file parameters.txt.
#             Only the first line of the file is read.  Parameter names
#             should contain no spaces or punctuation and are assumed to be 
#             separated by commas. They should correspond to the columns
#             of the deconvolution files.
#
#Output:
#  
#  returns string list of column (parameter) names
#
#Usage Notes:
#

def read_parnames(runpath):

    #--Set file paths--
    parfile = runpath+'/parameters.txt'

    outpars = utilities.parse_file_line(parfile)

    #--return string list--
    return outpars

#----------------------------------------------------------
# merge_output
#
#Author: Kari A. Frank 
#Date: October 16, 2015
#Purpose: Read in deconvolution, statistic, or other files from an xmc run
#          and merge into a single tab-separated text file, with an extra
#          column to specify iteration
#          
#Usage: merge_output(filepath,filetype='deconvolution',save=True)
#
#Input:
#
#  filepath     -- Path to xmc run directory, containing the
#                  deconvolution, statistic, etc. files
# 
#  filetype     -- string containing the type of xmc file to merge
#                  deconvolution, statistic, sigma, mean, changed
#
#  save         -- boolean switch to prevent saving of the merged 
#                  file. use save=False if all that is need is the merged 
#                  dataframe (default=True).
#
#Output:
#  - if save=True, writes a text file <filetype>_merged.txt into the filepath 
#    directory
#  - returns a pandas dataframe containing data from all the input files
#
#Usage Notes:
#  - overwrites the file if it already exists
#  - automatically includes all iterations present
#  - iterations are not written in order, but and extra column is added
#    to specify which iteration each row is associated with
#  - the created file can be read into a dataframe using
#    datatable = pd.read_table(mergedfile,sep='\t')
#

def merge_output(filepath,filetype='deconvolution',save=True):

    # - get list of files - 
    filelist = utilities.ls_to_list(filepath,ls_args = filetype+'*')
    
    # --Initialize dataframe with first file --

    # - read file - 
    datatable = pd.read_table(filepath+'/'+filelist[0],sep='\s+',header=None)

    # - add headers -
    parnames = []
    # read header information
    if filetype == 'statistic':
        parnames =  ['stat1','dof','oversim','nblobs','alpha','chi2']
        datatable.columns = parnames
    if filetype == 'deconvolution':
        if os.path.isfile(filepath+'/parameters.txt'):
            parnames = utilities.parse_file_line(filepath+'/parameters.txt') 
            datatable.columns = parnames

    # - get iteration number and add as column - 
    iternum = int(filelist[0].split('.')[-1]) # returns integer
    datatable['iteration'] = iternum

    # -- Loop through files and concatenate into single dataframe --
    for f in filelist[1:]:        
        newframe = pd.read_table(filepath+'/'+f,sep='\s+',header=None)
        iternum = int(f.split('.')[-1]) # returns integer
        if len(parnames) > 0: newframe.columns=parnames
        newframe['iteration'] = iternum
        datatable = pd.concat([datatable,newframe],ignore_index=True)
        
    # -- Write to file --
    if save == True: datatable.to_csv(filepath+'/'+filetype+'_merged.txt',
                                      sep='\t')

    return datatable

#----------------------------------------------------------
# get_xmcoutput
#
#Author: Kari A. Frank 
#Date: March 24, 2014
#Purpose: Read in deconvolution or statistic files from an xmc run
#          and return as a large array, with one parameter
#          per column and an extra column for iteration number
#          
#
#Usage: deconvolution_array = get_xmcoutput(runpath,itmin=itmin,itmax=itmax,
#                                           kind='deconvolution')
#
#Input:
#
#  runpath     -- Path to xmc run directory, containing the
#                  deconvolution files.
# 
#  itmin,itmax -- minimum and maximum number of iterations to read
#                 default is 0 to maximum available
#
#  kind        -- optional string to specify which type of output
#                 file should be read, 'deconvolution' (default)
#                 of 'statistic'.
#
#Output:
#  
#  returns a tuple of the float array of deconvolution (or statistic)
#     file content for all iterations in range itmin-itmax, and a 
#     string list of the parameter(column) names: (array,parnames)
#
#Usage Notes:
#

def get_xmcoutput(runpath,itmin=0,itmax=-1,kind='deconvolution'):

    if itmax == -1 :
        itmax = len(utilities.ls_to_list(runpath,'deconvolution.*'))-1

    #--set base filenames--
    dfile = runpath+'/'+kind+'.'
    parfile = runpath+'/parameters.txt'

    #--get parameter (column) names--
    if kind == 'deconvolution':
        parnames = read_parnames(runpath)
    else:
        parnames = ['stat1','dof','oversim','nblobs','alpha','chi2']
    npars = len(parnames)

    #--blob (row) numbers--
    if kind == 'deconvolution':
        nblobs = utilities.file_lines(dfile+'0')
#        blobnums = range(0,nblobs-1) #blob numbers start at 0
    else:
        nblobs = 1 #only one line in statistic files

    #--loop through deconvolution files--
    
    for i in range(itmin,itmax):
        
        ifile = dfile+str(i)
#        print ifile
        
        #-define new empty array-
        total_arr = np.empty((nblobs,npars+1))

        #-read file content into numpy array-
        # first dimension is blobnumber (row), second
        # dimension is parameter (column)

        ## eventually turn this into a separate function
        ## which takes as arguments the parnames (or can 
        ## optionally read them from a header line in the
        ## file itself) and a single file name, so that it
        ## can also read in files written by write_blobs()

        #-check if file exists and is not empty-
        if (path.isfile(ifile) == False) or (path.getsize(ifile) == 0):
            print "WARNING: "+ifile+" is missing or empty and will be skipped."
        else:
            total_arr[:,:-1] = np.loadtxt(ifile)

            #-set last column to iteration number-
            total_arr[:,-1] = i

            #-concatenate onto master array-
            if i == itmin:
                master_array = total_arr
            else:
                master_array = np.vstack((master_array,total_arr))

    #--append iteration to parnames--
    parnames.append('iteration')

    #--convert parnames list to an array--
#    parnames_array = numpy.array(parnames)

    #--return array and parameter names--
    return (master_array,parnames)

#----------------------------------------------------------
