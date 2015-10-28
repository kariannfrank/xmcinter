#----------------------------------------------------------
#Module of handy general astronomy functions
#
#Contains the following functions:
#
# norm_to_em
# transfer_header()
# normalize_image()
# stack_images()
# angstroms2keV()
# K2keV()
# convert_distance()
#
#----------------------------------------------------------
#Import common modules
import shutil
import numpy as np
from astropy.io import fits


#----------------------------------------------------------
#norm_to_em()

#Author: Kari Frank
#Date: March 28, 2014
#Purpose: given model normalization (as defined for Xspec mekal model),
#          distance in cm, and redshift, return the emission measure
#         
#Usage: norm_to_em(norm,dist,redshift=redshift)
#Input: 
#  norm -- model normalization, in units output by Xspec 
#  dist -- distance to the object in cm
#  redshift -- redshift of the object. default = 0.0
#
#Output:
#  returns emission measure corresponding to the given norm (float).
#Usage Notes:
# 
# - Assumes  
#      norm = 10^14/4pi[dist(1+z)]^2 * EM
# - Units of returned emission measure are cm^-3
# - EM = n_e*n_H*V

def norm_to_em(norm,dist,redshift=0.0):

  #-convert to emission measure-
  em = norm*10.0**(14.0)*dist*4.0*np.pi*dist*(1.0+redshift)**2.0
  
  #-return emission measure-
  return em


#----------------------------------------------------------
#transfer_header()

#Author: Kari Frank
#Date: April 4, 2014
#Purpose: copy a header from one fits file into another fits file
#         
#         
#Usage: transfer_header(sourcefile,targetfile,newfile)
#Input: 
#  sourcefile -- fits file from which to read the header
#  targetfile -- fits file which will have its header replaced
#  newfile    -- new fits file name to save the edited targetfile
#
#Output:
#  makes a copy of targetfile with the header replaced by that
#   from sourcefile.
#  returns the name of the new file (newfile)
#
#Usage Notes:
#  
# 
# 

def transfer_header(sourcefile,targetfile,newfile):

  #-get sourcefile header-
  source_hdus = fits.open(sourcefile)
  source_hdr = source_hdus[0].header
  new_hdr = source_hdr.copy()
  source_hdus.close()

  #-open targetfile-
  target_hdus = fits.open(targetfile)
  target_img = target_hdus[0].data
  new_img = target_img.copy()
  target_hdus.close()
  
  #-write to file-
  fits.writeto(newfile,new_img,new_hdr,clobber=True)

  return newfile


#----------------------------------------------------------
#normalize_image()

#Author: Kari Frank
#Date: February 2, 2015
#Purpose: Normalize a fits image such that is sums to 1
#         
#Usage: normalize_image(infile,outfile='default')
#Input: 
#  infile -- fits file to normalize
#  outfile -- name of output normalized fits image. default is infile+'.norm'
#
#Output:
#  makes a copy of infile that is normalized
#  returns normalized image array
#
#Usage Notes:
#  
# 
# 

def normalize_image(infile,outfile='default'):

  #-set output file name-
  if outfile=='default':
    outfile=infile+'.norm'

  #-copy image file-
  shutil.copy(infile,outfile)

  #-open file-
  hdus = fits.open(outfile,mode='update')
  
  #-normalize-
  out_img = hdus[0].data
  out_img = out_img.astype(float)
  out_img=out_img/np.sum(out_img)
  hdus[0].data = out_img

  #-write normalized file-
  #in_hdus.writeto(outfile,clobber=True)
 
  hdus.close()

  return out_img

#----------------------------------------------------------
#stack_images()

#Author: Kari Frank
#Date: February 3, 2015
#Purpose: Read in stack of fits images and add them together.
#         
#Usage: stack_images(infiles,outfile='default')
#Input: 
#  infile -- list of fits image filenames to stack
#  outfile -- name of output stacked fits image. default is 'stacked.fits'
#
#Output:
#  fits image file containing sum of input images
#  returns the stacked image numpy array
#
#Usage Notes:
# - assumes images are already aligned (in image coordinates) 
# - the output file, including header, is a copy of the first
#   input image with modified image data

def stack_images(infiles,outfile='default'):

  #-set output file name-
  if outfile=='default':
    outfile='stacked.fits'

  #-copy image file-
  shutil.copy(infiles[0],outfile)

  #-open output file-
  hdus = fits.open(outfile,mode='update')
  out_img = hdus[0].data


  #-loop through remaining images-
  for infile in infiles[1:]:
    inhdus = fits.open(infile)
    in_img = inhdus[0].data
    out_img = out_img+in_img
    inhdus.close()

  #-flush updated image array to file and close-
  hdus[0].data = out_img
  hdus.close()

  return out_img

#----------------------------------------------------------
#angstroms2keV()

#Author: Kari Frank
#Date: July 13, 2015
#Purpose: quickly convert angstroms to keV
#         
#Usage: angstroms2keV(inwave)
#Input: 
#  inwave -- input wavelength, in angstroms
#
#Output:
#  returns the wavelength in keV (if input was angstroms)
#   or the wavelength in angstrom (if the input was keV)
#
#Usage Notes:
# - uses E=h*c/lambda

def angstroms2keV(inwave):

  #-set conversion constant-
  C = 12.398521
  
  #-convert (same equation regardless of keV or A input)-
  outwave = C / inwave

  return outwave

#----------------------------------------------------------
#Author: Kari A. Frank
#Date: October 9, 2013
#Purpose: 
#Usage: energy_keV = K2keV(temperature,reverse=False)
#
#Input:
#
# temperature: float temperature in Kelvins or keV. 
#
# reverse: switch to convert input from keV to Kelvins
#          (default=False).  if True, then input 
#          should be the temperature in keV.
#
#Output:
# - returns the converted temperature
#
#Usage Notes:
# 

def K2keV(temperature,reverse=False):

  #-define constant conversion factor-
  k = 8.62*10**(-8) #constant in keV/K

  #-convert from K to keV-
  if not reverse:
    energy_keV = k*float(temperature)
    return energy_keV
  #-convert from keV to K-
  else:
    energy_keV = float(temperature)
    #convert to Kelvins
    kelvins = energy_keV/k
    return kelvins

#----------------------------------------------------------
# convert_units()
#Author: Kari A. Frank
#Date: October 28, 2015
#Purpose: Convert between distance units
#Usage: convert_distance(val,fromunit,tounit)
#
#Input:
#
# val      -- value to be converted
#
# fromunit -- string specifying units of the input val:
#             'kpc', 'km','cm','ly','pc'
#
# tounit   -- string specifying the units to convert val into:
#             same values as fromunit.
#
#Output:
# - returns the converted val
#
#Usage Notes:
# - 
#   

def convert_distance(val,fromunit,tounit):

    #--set constants--
    kpc_to_km = 1000.0*3.0857*10.0**13.0
    pc_to_km = kpc_to_km/1000.0
    km_to_cm = 100000.0
    km_to_ly = 1.057e-13

    #--convert values--

    if fromunit == 'pc':
      if tounit == 'km':
        return val*pc_to_km
      if tounit == 'cm':
        return val*pc_to_km*km_to_cm
      if tounit == 'kpc':
        return val*pc_to_km/kpc_to_km
      if tounit == 'ly':
        return val*pc_to_km*km_to_ly

    if fromunit == 'km':
      if tounit == 'pc':
        return val/pc_to_km
      if tounit == 'cm':
        return val*km_to_cm
      if tounit == 'kpc':
        return val/kpc_to_km
      if tounit == 'ly':
        return val*km_to_ly

    if fromunit == 'cm':
        if tounit == 'pc':
            return val/km_to_cm/pc_to_km
        if tounit == 'kpc':
            return val/km_to_cm/kpc_to_km
        if tounit == 'km':
            return val/km_to_cm
        if tounit == 'ly':
            return val/km_to_cm*km_to_ly

#----------------------------------------------------------
