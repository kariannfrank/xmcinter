#----------------------------------------------------------
#Module of functions for dealing with files output from xmc
#
#Contains the following functions:
#
# parse_file_line
# read_parnames

# functions using pandas:
#  merge_output -- replaces get_deconvolution, get_statistic, and get_xmcoutput
#
#----------------------------------------------------------
#Import common modules
#import os
#import numpy as np
#import pandas as pd
import file_utilities as fu

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
#
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

    outpars = fu.parse_file_line(parfile)

    #--return string list--
    return outpars

#----------------------------------------------------------
# merge_output
#
#Author: Kari A. Frank 
#Date: October 16, 2015
#Purpose: Read in deconvolution, statistic, or other files from an xmc run
#          and merge into a single dataframe, with an extra
#          column to specify iteration. optionally can also save as a text file.
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
#  sep          -- optionally specify the delimiter in the output file
#                  (default='\t', tab-separated)
#
#Output:
#  - if save=True, writes a text file <filetype>_merged.txt into the filepath 
#    directory
#  - returns a pandas dataframe containing data from all the input files
#
#Usage Notes:
#  - overwrites the file if it already exists
#  - automatically includes all iterations present
#  - iterations are not written in order, but an extra column is added
#    to specify which iteration each row is associated with
#  - the created file can be read into a dataframe using
#    datatable = pd.read_table(mergedfile,sep='\t',index_col=0)
#

def merge_output(filepath,filetype='deconvolution',save=True,sep='\t'):

    # - get list of files - 
    filelist = fu.ls_to_list(filepath,ls_args = filetype+'.*')
    
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
            parnames = fu.parse_file_line(filepath+'/parameters.txt') 
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
                                      sep=sep)

    return datatable

#----------------------------------------------------------
