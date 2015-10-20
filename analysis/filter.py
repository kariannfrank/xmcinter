#----------------------------------------------------------
#Module of functions for dealing with intermediate
#  products of xmc analysis procedures (e.g. the 
#  dictionaries output by get_deconvolution() or 
#  get_statistic() ).
#
#Contains the following functions:
#
# weighted_median()
# make_mask() 
# abund_ratio()
#
#----------------------------------------------------------
#Import common modules


#----------------------------------------------------------
# weighted_median()
#
#Author: Kari A. Frank 
#Date: March 28, 2014
#Purpose: 
#         Return the weighted median of the provided array and weights.
#         Based on the xmc idl function weighted_median()        
#
#Usage: middle = weighted_median(values,weights)
#
#Input:
#
#  values  -- array of the values to find the median of 
#              
#  weights -- array of weights corresponding to values            
#              
#  
#Output:
#  
#  returns the weighted median of values
#
#Usage Notes:
#
#  - values and weights must have the same dimensions
#  - should work for lists or numpy arrays
#  - numpy arrays are converted to 1D lists
#
#Example:
#
#

#----define exceptions---
class UnmatchedDimensions(Exception): pass

#----define function----
def weighted_median(values,weights):

    try:
        #--check for matching dimensions--
        if type(values) is list:    
            if len(values) != len(weights):
                exerror = UnmatchedDimensions()
                exerror.message = ("ERROR: values and weights have different"
                                   " dimensions.")
                raise exerror
        else:
            if type(values) is np.ndarray:
                if values.shape != weights.shape:
                    exerror = UnmatchedDimensions()
                    exerror.message = ("ERROR: values and weights have "
                                       "different dimensions.")
                    raise exerror
    except UnmatchedDimensions,ex:
        print ex.message

    #--if numpy array, flatten to 1D list--
    if type(values) is numpy.ndarray:
        values.flatten()
        weights.flatten()
        values.tolist()
        weights.tolist()

##eventually add capability to ignore zero values, do it here
    trimmed_values = values
    trimmed_weights = weights

#    print values
#    print weights
        
    #--calculate median--
        
    #-zip the values and weights lists for sorting-
    zipped = zip(trimmed_values,trimmed_weights)
        
    #-sort on values-
    sorted_zipped = sorted(zipped,key=lambda val: val[0])

    #-unzip-
    sorted_values,sorted_weights = zip(*sorted_zipped)
        
    #-create cumulative weight array, quit when cumul_weight crosses 0.5-
    net_weight = sum(sorted_weights)
    i = 0
    total = 0
    while total <= 0.5:
        total = total + sorted_weights[i]/net_weight
        i = i + 1
    #-check distance from 0.5-
    dist = abs(0.5 - total)
    #-get next weight, to see if it is closer to 0.5-
    new_total = total + sorted_weights[i]/net_weight
    if abs(new_total - 0.5) > dist:
        half_sub = i
    else:
        half_sub = i-1


#    cumul_weights = [sum(sorted_weights[:i+1])/net_weight for i in range(len(sorted_weights))]

    #-find where cumulative weight = 0.5-
#    half_sub,half_cumul_weight = min(enumerate(cumul_weights),key=lambda x:abs(x[1]-0.5))

    #-return weighted median-
    return sorted_values[half_sub]

        
#----------------------------------------------------------
# make_mask()
#
#Author: Kari A. Frank 
#Date: April 11, 2014
#Purpose: 
#         Return a boolean mask array based on 
#          with blobs (rows) masked 
#          based on given parameter conditions.
#         
#Usage: masd = xmcanalyze.make_mask(data_array,col_name,minparam,maxparam)
#
#Input:
#
#  data_array -- array of blob parameters (e.g. from read_blobs.pro),
#                assumes one blob per line
#              
#  param_i    -- column index of the parameter (column in data_array) to
#                which to apply the conditions
#
#  minparam,maxparam -- minimum and maximum values of the
#                       specified parameter to keep             
#  
#Output:
#  
#  Returns a boolean mask the same shape as data_array, with 
#    all blobs masked for which col_name parameter lies outside
#    the range [minparam,maxparam]
#
#Usage Notes:
#
#  
#
#Example:
#
#

#----define exceptions---
#class UnmatchedDimensions(Exception): pass

#----define function----
def make_mask(data_array,param_i,minparam,maxparam):

    import numpy

    #--make empty masks (nothing masked)--
    max_mask = numpy.ma.make_mask_none(numpy.shape(data_array))
    min_mask = numpy.ma.make_mask_none(numpy.shape(data_array))

    #-find 'bad' rows(blobs)--
    max_rows = numpy.where(data_array[:,param_i] > maxparam)
    min_rows = numpy.where(data_array[:,param_i] < minparam)

    #-mask all 'bad' blobs-
    max_mask[max_rows,:] = True
    min_mask[min_rows,:] = True
    themask = numpy.ma.mask_or(max_mask,min_mask)

    #-return mask-
    return themask

#----------------------------------------------------------
# abund_ratio()
#
#Author: Kari A. Frank 
#Date: May 15, 2014
#Purpose: 
#         Return a boolean mask array based on 
#          with blobs (rows) masked 
#          based on given parameter conditions.
#         
#Usage: masd = xmcanalyze.abund_ratio(data_array,col_name,weights=weights,stat='average',ref_elem=ref_elem,min=min)
#
#Input:
#
#  data_array -- array of blob parameters (e.g. from read_blobs.pro),
#                assumes one blob per line
#              
#  param_i    -- column index of the parameter (column in data_array) to
#                which to apply the conditions
#
#  weights    -- optional weights for each blob, typically emission measure
#  
#  min        -- only include blobs with abundance > min
#  
#Output:
#  
#  Returns the (weighted) median or average of the specified element
#  
#  
#
#Usage Notes:
#
#  
#
#Example:
#
#

#----define exceptions---
#class UnmatchedDimensions(Exception): pass

#----define function----
def abund_ratio(data_array,col_i,weights=1,stat='average',minval=0.0):
    
    import numpy

#    weights=1
    abund = data_array[:,col_i]
#    if weights == 1: weights = numpy.ones_like(abund)
    ej = numpy.where(abund >= minval)
    avg = numpy.average(abund[ej],weights=weights[ej])
    med = weighted_median(abund[ej],weights=weights[ej])

    print 'col,avg,med = ',col_i,avg,med
    
    if stat == 'average':
        return avg
    else:
        return med
    
        
