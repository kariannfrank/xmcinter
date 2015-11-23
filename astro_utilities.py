"""
Module of handy general astronomy functions

Contains the following functions:

 norm_to_em
 transfer_header()
 normalize_image()
 stack_images()
 angstroms2keV()
 K2keV()
 convert_distance()
 get_xmm_attitude()
 hms2deg()
 deg2hms()
 xmc_wcs_convert()

"""
#Import common modules
import shutil
#import numpy as np
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

def norm_to_em(norm,dist_cm,redshift=0.0):

  #-convert to emission measure-
  em = norm*10.0**(14.0)*dist_cm*4.0*np.pi*dist_cm*(1.0+redshift)**2.0
  
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
# val      -- numerical value to be converted
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

    if fromunit == 'kpc':
        if tounit == 'pc':
            return val*1000.0
        if tounit == 'cm':
            return val*kpc_to_km*km_to_cm
        if tounit == 'km':
            return val*kpc_to_km
        if tounit == 'ly':
            return val*kpc_to_km*km_to_ly

#----------------------------------------------------------
def get_xmm_attitude(attfile='atthk.fits',hms=False):
  """
  get_xmm_attitude

  Author: Kari Frank
  Date: November 17, 2015
  Purpose: Return an XMM observation boresight coordinates and rotation 
           angle from an attitude (e.g. atthk.fits) file.

  Usage: obscoords = get_xmm_attitude(attfile='atthk.fits',hms=FALSE)
  Input:
      attfile : [string] filename (including path) of the attitude file

      hms : [boolean] return the coordinates as [hour,minutes,seconds] 
            instead of degrees

  Output:
      Returns list of coordinates and angle, 
      [right_ascension,declination,angle], where all are in decimal degrees

      If hms=True, then right_ascension and declination are lists of the 
      form [hh,mm,ss]

  Usage Notes:
     - the attitude file should be generated from the odf files using 
       the XMM-SAS tool atthkgen
  """
  
  #-ra, dec, and rotation from attitude file header-
  att_hdus = fits.open(attfile)
  att_hdr = att_hdus[0].header
  # this retrieves the median values in case they were not constant, e.g. 
  # due to a slew failure
  att_ra = att_hdr['MAHFRA']
  att_dec = att_hdr['MAHFDEC']
  att_ang = att_hdr['MAHFPA']
  att_hdus.close()

  #-convert to hh,mm,ss if hms=True-
  if hms:
    coords = [deg2hms(att_ra),deg2hms(att_dec),att_ang]
  else:
    coords = [att_ra,att_dec,att_ang]

  #-return coordinates and angle-
  return coords

#----------------------------------------------------------
def deg2hms(x_deg):
  """
  Author: Kari Frank
  Date: November 17, 2015
  Purpose: Convert a value from decimal degrees into hours,minutes,seconds.

  Usage: y = deg2hms(x_deg)

  Input:
      x_deg : value in decimal degrees (e.g. right ascension)

  Output:
      Returns the converted value as a list in the form [hh,mm,ss]

  Usage Notes:
  """
  
  #-set constants-
  day_length = 24.0

  #-convert to hours-
  x_hour_full = x_deg*day_length/360.0
  x_hour = np.floor(x_hour_full)
  #-convert remainder to minutes
  x_min_full = (x_hour_full - x_hour)*60.0
  x_min = np.floor(x_min)
  #-convert remainder to seconds-
  x_sec = (x_min_full - x_min)*60.0

  return [x_hour,x_min,x_sec]
#----------------------------------------------------------
def hms2deg(x_hms):
  """
  Author: Kari Frank
  Date: November 17, 2015
  Purpose: Convert a value from hours,minutes,seconds into decimal degrees.
           (This is the inverse of deg2hms)

  Usage: y = hms2deg(x_hms)

  Input:
      x_hms : list in format [hh,mm,ss] (e.g. right ascension)

  Output:
      Returns the input as a single value in decimal degrees

  Usage Notes:
  """
  
  #-set constant-
  day_length = 24.0

  #-convert hours to degrees-
  x_deg = x_hms[0]*360.0/day_length
  #-convert minutes to degrees-
  x_deg = x_deg + x_hms[1]*360.0/(60.0*day_length)
  #-convert seconds to degrees-
  x_deg = x_deg + x_hms[2]*360.0/(60.0*60.0*day_length)

  return x_deg
#----------------------------------------------------------
def xmc2wcs(phi,psi,attfile='atthk.fits'):
  """
  Author: Kari Frank
  Date: November 18, 2015
  Purpose: Convert from xmc coordinates to RA and DEC.

  Usage: y = xmc2wcs(phi,psi,attfile='atthk.fits')

  Input:
     phi (float): xmc phi (x) coordinate (arcsec)
     psi (float): xmc psi (y) coordinate (arcsec)
     attfile (string) : path to the observation's attitude file, 
                        as created by the XMM-SAS tool atthkgen

  Output:
     Returns list [ra,dec] of the input xmc coordinates in decimal degrees.

  Usage Notes:
    - See also wcs2xmc()
  """
  
  #--Set Center Coordinates--

  #-get observation boresight and rotation angle-
  obscoords = get_xmm_attitude(attfile=attfile) #=[ra,dec,angle]
  obs_ra = obscoords[0]
  obs_dec = obscoords[1]
  obs_ang = obscoords[2]*np.pi/180.0 - np.pi/2.0

  #-set xmc center coordinates-
  phi0 = -60.0
  psi0 = 80.0

  #-convert xmc center coords to ra and dec in degrees-
  ra0 = (obs_ra + (1.0/np.cos(obs_dec*np.pi/180.))*
         (-1.0*phi0*np.cos(obs_ang)-psi0*np.sin(obs_ang))/3600.)
  dec0 = obs_dec - (-1.0*phi0*np.sin(obs_ang)+psi0*np.cos(obs_ang))/3600.

  #--Convert to WCS--
  
  #-shift phi and psi-
  psi = -(psi-psi0)
  phi = phi - phi0

  #-rotate phi and psi-
  phi_rot = phi*np.cos(obs_ang) - psi*np.sin(obs_ang)
  psi_rot = phi*np.sin(obs_ang) - psi*np.cos(obs_ang)

  #-shift to ra and dec and convert units-
  ra = ra0 - (1.0/np.cos(dec0*np.pi/180.0))*phi_rot/3600.0
  dec = dec0 + psi_rot/3600.

  return [ra,dec]
#----------------------------------------------------------
def wcs2xmc(ra,dec,attfile='atthk.fits'):
  """
  Author: Kari Frank
  Date: November 18, 2015
  Purpose: Convert from RA and DEC into phi and psi coordinates used by xmc.

  Usage: y = wcs2xmc(ra,dec,attfile='atthk.fits')

  Input:
     ra (float): right ascension in decimal degrees
     dec (float): declination in decimal degrees
     attfile (string) : path to the observation's attitude file, 
                        as created by the XMM-SAS tool atthkgen

  Output:
     Returns list [phi,psi] of the input ra and dec coordinates 
     in arcsec.

  Usage Notes:
     - See also xmc2wcs()
  """
  
  #--Set Center Coordinates--

  #-get observation boresight and rotation angle-
  obscoords = get_xmm_attitude(attfile=attfile) #=[ra,dec,angle]
  obs_ra = obscoords[0]
  obs_dec = obscoords[1]
  obs_ang = obscoords[2]*np.pi/180.0 - np.pi/2.0

  #-set xmc center coordinates-
  phi0 = -60.0
  psi0 = 80.0

  #-convert xmc center coords to ra and dec in degrees-
  ra0 = (obs_ra + (1.0/np.cos(obs_dec*np.pi/180.))*
         (-1.0*phi0*np.cos(obs_ang)-psi0*np.sin(obs_ang))/3600.)
  dec0 = obs_dec - (-1.0*phi0*np.sin(obs_ang)+psi0*np.cos(obs_ang))/3600.

  #--Convert to xmc coordinates--
  
  #-distance from input coordinates to center-
  phi_rot =  -(ra - ra0)*np.cos(obs_dec*np.pi/180.)*3600. 
  psi_rot = (dec - dec0)*3600.

  #-rotated phi and psi-
  phi = phi_rot*np.cos(obs_ang)+psi_rot*np.sin(obs_ang)
  psi = -phi_rot*np.sin(obs_ang)+psi_rot*np.cos(obs_ang)

  #-shift according to xmc center coords-
  phi = phi + phi0
  psi = psi0 - psi

  return [phi,psi]
#----------------------------------------------------------
