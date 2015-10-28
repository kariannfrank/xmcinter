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
#Usage: filtered_df = filter(inframe,colname,minval=min,maxval=max,equalval=equalval)
#
#Input:
#
#  inframe -- pandas dataframe
#               
#  colname  -- (string) name of the column containing the values to
#               be filtered by (e.g. 'blobkT' to filter based on 
#               minimum or maximum temperature)            
#  
#
#  minval,maxval   -- minimum and maximum values of the specified parameter 
#                     to keep. returned value range is inclusive.
#  
#  equalval     -- if provided, filter will return only rows for which the value
#               colname is equal to this. if provided, then min and max 
#               arguments are ignored.
#
#Output:
#  
#  Returns a dataframe identical to the input dataframe but missing rows
#    which do not match the criteria.
#
#Usage Notes:
#
#  
#
#Example:
#    filtered_df = filter(blobframe,'blobkT',min=0.5,max=1.5)
#    - this will return a version of blobframe which includes only 
#      rows (blobs) with temperatures between 0.5 and 1.5.
#


#----define exceptions---
#class UnmatchedDimensions(Exception): pass

#----import modules---
import pandas as pd

#----define function----
def filter(inframe,colname,minval=None,maxval=None,equalval=None):

    #-check if colname is valid column name-
    if colname not in inframe.columns:
        raise ValueError(colname+" is not a valid column name")

    #-filter based on value-
    if equalval == None:

        #-set default minval/maxval if none provided-
        if minval == None: minval = min(inframe[colname])
        if maxval == None: maxval = max(inframe[colname])

        #-filter the dataframe-
        outframe = inframe[(inframe[colname]<=maxval) & 
                           (inframe[colname]>=minval)]
    else:
        outframe = inframe[inframe[colname]==equalval]

    #-warn if zero rows match-
    if len(outframe.index) == 0:
        print "Warning: Returned dataframe has zero rows."

    #-warn if all rows match-
    if len(outframe.index) == len(inframe.index):
        print "Warning: Nothing was filtered (all rows match criteria)."

    #-return filtered dataframe-
    return outframe
