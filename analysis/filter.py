# filter()
#
#Author: Kari A. Frank 
#Date: October 20, 2015
#Purpose: Filter a pandas dataframe (i.e. remove rows) based on 
#         minimum and maximum values for the specified column.
#         
#         
#         
#         
#Usage: filtered_df = filter(inframe,colnames,minvals,maxvals)
#
#Input:
#
#  inframe -- pandas dataframe
#               
#  colnames  -- scalar or list of (string) names of the columns containing
#               the parameters be filtered by (e.g. 'blobkT' to filter based on 
#               minimum or maximum temperature)            
#  
#
#  minvals,maxvals  -- scalars or lists of minimum and maximum values of
#                      parameters to keep. returned value range is inclusive.
#                      order must match the order of columns
#
#Output:
#  
#  Returns a dataframe identical to the input dataframe but missing rows
#    which do not match the criteria.
#
#Usage Notes:
# - to return rows for which the column equals a specific value (rather
#   than is within a range), set minval=maxval
#
# - to include only upper or lower limits for any parameter, set the 
#   associated minval or maxval to None  
#
#Example:
#    filtered_df = filter(blobframe,'blobkT',min=0.5,max=1.5)
#    - this will return a version of blobframe which includes only 
#      rows (blobs) with temperatures between 0.5 and 1.5.
#
#    filtered_df = filter(blobframe,['blobkT','blob_Si'],minvals=[0.5,1.0],
#                         maxvals=[1.5,3.0])
#    - this will return a version of blobframe which includes only 
#      rows (blobs) with temperatures between 0.5 and 1.5 and Si between
#      1.0 and 3.0.
#
#


#----define exceptions---
#class UnmatchedDimensions(Exception): pass

#----import modules---
import pandas as pd

#----define simple function to filter by 1 parameter----
#--gets called by the primary filter function below--
def simplefilter(inframe,colname,minval=None,maxval=None):

    #-check if colname is valid column name-
    if colname not in inframe.columns:
        raise ValueError(colname+" is not a valid column name")

    #-filter based on value-

    #-set default minval/maxval if none provided-
    if minval == None: minval = min(inframe[colname])
    if maxval == None: maxval = max(inframe[colname])

    #-filter the dataframe-
    outframe = inframe[(inframe[colname]<=maxval) & 
                       (inframe[colname]>=minval)]

    #-warn if zero rows match-
    if len(outframe.index) == 0:
        print "Warning: Returned dataframe has zero rows."

    #-warn if all rows match-
    if len(outframe.index) == len(inframe.index):
        print "Warning: Nothing was filtered (all rows match criteria)."

    #-return filtered dataframe-
    return outframe

def filter(inframe,colnames,minvals=None,maxvals=None):

    #--initialize outframe--
    outframe = inframe

    #--check if multiple parameters--
    if type(colnames) is str:
        colnames = [colnames]
        minvals = [minvals]
        maxvals = [maxvals]

    #--filter on each parameter--
    for col in range(len(colnames)):
        outframe = simplefilter(outframe,colname=colnames[col],
                                minval=minvals[col],maxval=maxvals[col])

    return outframe
