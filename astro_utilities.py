"""
Module of handy general astronomy functions

Contains the following functions:

 norm_to_em()
 em_to_density()
 transfer_header()
 normalize_image()
 stack_images()
 angstroms2keV()
 K2keV()
 convert_distance()
 convert_arcsec()
 get_xmm_attitude()
 hms2deg()
 deg2hms()
 xmc_wcs_convert()
 chandra_psf_wing()
 xspec_abund_to_nomoto()

"""

#----------------------------------------------------------
#Import common modules
import shutil
import numpy as np
from astropy.io import fits

#----------------------------------------------------------
def norm_to_em(norm,dist_cm,redshift=0.0):
  """
  Author: Kari Frank
  Date: March 28, 2014
  Purpose: given model normalization (as defined for Xspec mekal model),
            distance in cm, and redshift, return the emission measure
          
  Input: 
    norm (numerical): model normalization, in units output by Xspec 
    dist (numerical): distance to the object in cm
    redshift (numerical): redshift of the object. default = 0.0
  
  Output:
    returns emission measure corresponding to the given norm (float).

  Usage Notes:
   
   - Assumes  
        norm = 10^-14/4pi[dist(1+z)]^2 * EM
   - Units of returned emission measure are cm^-3
   - EM = n_e*n_H*V

  """
  #-convert to emission measure-
  em = norm*10.0**(14.0)*dist_cm*4.0*np.pi*dist_cm*(1.0+redshift)**2.0
  
  #-return emission measure-
  return em

#----------------------------------------------------------
def em_to_density(em,volume,density_type='number',mu=1.0):
  """
  Author: Kari Frank
  Date: April 14, 2016
  Purpose: Given emission measure and volume, calculate the density
          
  Input: 
    em (numerical): emission measure (e.g. as output from norm_to_em)
    volume (numerical): volume of region in cm^3
    density_type (string): type of density to return (default='mass')
            - 'mass' = g/cm^3
            - 'number' = 1/cm^3

  Output:
    returns emission measure corresponding to the given norm (float).

  Usage Notes:
   
   - Assumes  
        norm = 10^14/4pi[dist(1+z)]^2 * EM
        EM = n_e*n_H*V
   - Calculates only the hydrogen densities (even if mu!=1.0; setting a 
     different mu only changes the electron density used to convert from 
     emission measure)
   - For typical cosmic abundances, set mu=1.21 (n_He = 0.1*n_H, all 
     other elements negligible)

  """
  #-constants-
  proton_mass = 1.67*10.0**(-24.0) #grams

  #-convert to density
  if density_type == 'number':
    return (em/(volume*mu))**0.5  
  else:
    return proton_mass*(em/(volume*mu))**0.5  

#----------------------------------------------------------
def em_to_mass(em,volume,mu=1.0,tounit='g'):
  """
  Author: Kari Frank
  Date: July 5, 2016
  Purpose: Given emission measure and volume, calculate the hydrogen
           gas density.
  Input: 
    em (numerical): emission measure (e.g. as output from norm_to_em)
    volume (numerical): volume of region in cm^3
    mu (numerical): mean weight (cosmic abundance mu = 1.21, pure 
                    hydrogen=1.0)
    tounit (string): specify units of output mass. options are 'g' 
                     (grams, default), or solar masses ('sol')
    
  Output:
    returns mass in grams (default) or solar masses

  Usage Notes:
  """

  # solar mass in g
  Msol_g = 1.989*10.0**(33.0)
  
  nH = em_to_density(em,volume,density_type='mass',mu=mu)

  if tounit == 'sol':
    return nH*volume/Msol_g
  else:
    return nH*volume
    
#----------------------------------------------------------
def transfer_header(sourcefile,targetfile,newfile):
  """
  Author: Kari Frank
  Date: April 4, 2014

  Purpose: copy a header from one fits file into another fits file

  Input: 
    sourcefile (str): name of fits file from which to read the header
    targetfile (str): name of fits file which will have its header replaced
    newfile (str): name of new fits file to save the edited targetfile

  Output:
    makes a copy of targetfile with the header replaced by that
     from sourcefile.
    returns the name of the new file (newfile)

  Usage Notes:
  """   

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
def normalize_image(infile,outfile='default'):
  """
  Author: Kari Frank
  Date: February 2, 2015

  Purpose: Normalize a fits image such that is sums to 1

  Input: 
    infile (str): name of fits file to normalize
    outfile (str): name of output normalized fits image file. 
         default is infile+'.norm'

  Output:
    makes a copy of infile that is normalized
    returns normalized image array

  Usage Notes:
  """  
 
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
def stack_images(infiles,outfile='default'):
  """
  Author: Kari Frank
  Date: February 3, 2015

  Purpose: Read in stack of fits images and add them together.

  Input: 
    infile (list of str): list of fits image filenames to stack
    outfile (str): name of output stacked fits image. default is 
         'stacked.fits'

  Output:
    fits image file containing sum of input images
    returns the stacked image numpy array

  Usage Notes:
   - assumes images are already aligned (in image coordinates) 
   - the output file, including header, is a copy of the first
     input image with modified image data
  """

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
def angstroms2keV(inwave):
  """
  Author: Kari Frank
  Date: July 13, 2015
  Purpose: quickly convert angstroms to keV

  Usage: angstroms2keV(inwave)
  Input: 
    inwave (numerical): input wavelength, in angstroms

  Output:
    returns the wavelength in keV (if input was angstroms),
         or the wavelength in angstrom (if the input was keV)

  Usage Notes:
   - uses E=h*c/lambda
  """

  #-set conversion constant-
  C = 12.398521
  
  #-convert (same equation regardless of keV or A input)-
  outwave = C / inwave

  return outwave

#----------------------------------------------------------
def K2keV(temperature,reverse=False):
  """
  Author: Kari A. Frank
  Date: October 9, 2013
  Purpose: Convert temperature units between Kelvin and keV.

  Input:

   temperature (numerical): temperature in Kelvins (or keV if reverse=True)

   reverse (optional bool): switch to convert input from keV to Kelvins
            (default=False).  if True, then input 
            should be the temperature in keV.

  Output:
   numerical: returns the converted temperature

  Usage Notes:

  """

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
def convert_distance(val,fromunit,tounit):
  """
  Author: Kari A. Frank
  Date: October 28, 2015

  Purpose: Convert between distance units

  Input:

       val (numerical): value to be converted

       fromunit (str): units of the input val, choose from
            'kpc', 'km','cm','ly','pc'

       tounit (str): units to convert val into, same options
            as fromunit.

  Output:
       numerical: returns the converted val

  Usage Notes:

  """   

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
def convert_arcsec(theta,distance,distanceunit,tounit):
  """
  Author: Kari A. Frank
  Date: April 14, 2016

  Purpose: Convert angular length in arcsec to a physical length

  Input:

       theta (numerical): value to be converted, in arcsec

       distance (numerical): distance to the object

       distanceunit (str): units of distance, choose from
            'kpc', 'km','cm','ly','pc'

       tounit (str): units to convert theta into, same options
            as distanceunit

  Output:
       numerical: returns theta converted to the specified units

  Usage Notes:

  """   

  #--set constants--
  kpc_to_km = 1000.0*3.0857*10.0**13.0
  pc_to_km = kpc_to_km/1000.0
  km_to_cm = 100000.0
  km_to_ly = 1.057e-13
  arcsec_to_rad = np.pi/60.0/60.0/180.0

  #--convert distance into desired unit--
  if distanceunit != tounit:
    distance = convert_distance(distance,distanceunit,tounit)

  #--convert arcsec to radians--
  theta = theta*arcsec_to_rad 
  
  #--return angular distance--
  return distance*theta

#----------------------------------------------------------
def get_xmm_attitude(attfile='atthk.fits',hms=False):
  """
  get_xmm_attitude

  Author: Kari Frank
  Date: November 17, 2015
  Purpose: Return an XMM observation boresight coordinates and rotation 
           angle from an attitude (e.g. atthk.fits) file.

  Input:
      attfile (str): file name (including path) of the attitude file

      hms (bool): if True, return the coordinates as [hour,minutes,seconds] 
            instead of degrees (default=False)

  Output:
      list of numerical: returns list of coordinates and angle, 
           [right_ascension,declination,angle], 
           where all are in decimal degrees

           If hms=True, then right_ascension and declination are lists 
           of the form [hh,mm,ss]

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

  Input:
      x_deg (numerical): value in decimal degrees (e.g. right ascension)

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

  Input:
      x_hms (length 3 list of numerical): list in format [hh,mm,ss] 
           (e.g. right ascension)

  Output:
      Returns the input as a single value converted to decimal degrees

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

  Input:
     phi (numerical): xmc phi (x) coordinate (arcsec)
     psi (numerical): xmc psi (y) coordinate (arcsec)
     attfile (str) : path to the observation's attitude file, 
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

  Input:
     ra (numerical): right ascension in decimal degrees
     dec (numerical): declination in decimal degrees
     attfile (str) : path to the observation's attitude file, 
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
def chandra_psf_wing(dist,y,a,c,counts=1):
  """
  Author: Kari Frank
  Date: July 24, 2015
  Purpose: estimate the surface brightness due to the chandra psf wings
            at a given distance from a bright point source

  Usage: sb=chandra_psf_wing(dist,counts,y,a,c)
  Input: 
    dist -- distance in arcsec from the on-axis source
    counts -- total number of counts in the point source
    y,a,c -- parameters that describe the psf wing profile
             (from cxc.harvard.edu/cal/Hrma/HerX1-Wings.html)

  Output:
    returns the surface brightness in counts/arcsec^2 expected due to the 
      psf wings

  Usage Notes:

   - Only works for dist>15"
   - Profiles are of the form psi=a * (theta/theta0)^(-y) * exp(-c*theta)
       - theta0 is nominally 10"
   - set counts=1 to return the normalized surface brightness
  """

  from math import exp

  #--set constants--
  theta0=10.0

  #--calculate normalized surface brightness--
  sb = a * ((dist/theta0)**(-1.*y)) * exp(-1.*c*dist)

  #--convert to counts/arcsec^2--
  sb = counts*sb

  #--return final value--
  return sb

#----------------------------------------------------------
def xspec_abund_to_nomoto_dataframe(df,specZcol,refZcol,specZerrcol=None,
                             refZerrcol=None):
                             
  """
  Wrapper to call primary function on a dataframe
  
  Takes as input the dataframe and string names of relevant columns.
  Returns a pd.Series object with the new column

  Assumes errors are symmetric.
  """
  import pandas as pd

  if (specZerrcol is None) and (refZerrcol is None):
    ds = df.apply(lambda x: \
                     xspec_abund_to_nomoto(x[specZcol],x[refZcol]), axis=1)
  elif (specZerrcol is None) and (refZerrcol is not None):
    ds = df.apply(lambda x: \
        xspec_abund_to_nomoto(x[specZcol],x[refZcol],refZerr=x[refZerrcol]), 
                                        axis=1)
  elif (specZerrcol is not None) and (refZerrcol is None):
    ds = df.apply(lambda x: \
        xspec_abund_to_nomoto(x[specZcol],x[refZcol],specZerr=x[specZerrcol]), 
                                        axis=1)
  else:
    ds = df.apply(lambda x: \
        xspec_abund_to_nomoto(x[specZcol],x[refZcol],
                              specZerr=x[specZerrcol],
                              refZerr=x[refZerrcol]), axis=1)

  if (specZerrcol is None) and (refZerrcol is None):
    # return series of scalars
    df2 = pd.DataFrame(ds.tolist(),columns=['ratio','err1','err2'],index=ds.index)
    ds2=df2['ratio']
#    df['scalar'] = df.apply(lambda x: x[specZcol+'/'+refZcol][0])
    return ds2 #df['scalar']
  else:
    # return series of tuples
    return df #df[specZcol+'/'+refZcol]

def xspec_abund_to_nomoto(specZ,refZ,specZerr=0.0,refZerr=0.0):
  """
  Author: Kari Frank
  Date: May 12, 2016

  Purpose: Convert an abundance as defined in XSPEC with 
           the abundance ratio (relative to the element refZ)
           as defined in Nomoto2006:
           [X/Fe] = log(N_X/N_Fe) - log(N_X/N_Fe)_solar = log(Z_X/Z_Fe)

  Input:
     specZ (float) : abundance value relative to solar, as used in XSPEC
     refZ (float) : abundance, defined the same as specZ, to use as the
                    reference abundance (denominator in the ratio).
                    This should be an abundance associated with the same
                    region as specZ. Typically refZ is the Fe, Si, or O
                    abundance.
     specZerr,refZerr (float or 2-element float tuple) : measurement 
                    errors for the abundances. If provided, will also
                    return the errors on the abundance
                    ratio. If only one of specZerr, refZerr is provided, 
                    the other will be assumed to be zero.
                    If a tuple, then assumed to be of the form 
                    (lowerr,higherr) or (symmetricerr)
                  

  Output:
     Returns float containing the calculated abundance ratio, and optionally
     the associated error bars.

  Usage Notes:
     - 
  """

  #--convert errs to tuples if not--
  if not isinstance(specZerr,tuple): specZerr=(specZerr,specZerr)
  if not isinstance(refZerr,tuple): refZerr=(refZerr,refZerr)

  #--make sure given tuples are of correct form--
  if isinstance(specZerr,tuple) and (len(specZerr) > 2): 
    print "ERROR: specZerr tuple has more than 2 elements."
    return None
  if isinstance(specZerr,tuple) and (len(specZerr) == 1): 
    specZerr = (specZerr[0],specZerr[0]) # assume symmetric errors

  if isinstance(refZerr,tuple) and (len(refZerr) > 2):
    print "ERROR: refZerr tuple has more than 2 elements."
    return None
  if isinstance(refZerr,tuple) and (len(refZerr) == 1):
    refZerr = (refZerr[0],refZerr[0]) # assume symmetric errors
  

  #--convert to float if int--
  if isinstance(specZ,int): specZ = float(specZ)
  if isinstance(refZ,int): refZ = float(refZ)

  if isinstance(specZerr[0],int): specZerr[0] = float(specZerr[0])
  if isinstance(specZerr[1],int): specZerr[1] = float(specZerr[1])

  if isinstance(refZerr[0],int): refZerr[0] = float(refZerr[0])
  if isinstance(refZerr[1],int): refZerr[1] = float(refZerr[1])

  #--Calculate ratio--
  ratio = np.log10(specZ/refZ)
  
  #--Calculate errors--

  ratioerrlow = ( (specZerr[0]/(specZ*np.log(10.0)))**2.0 + 
                  (refZerr[0]/(refZ*np.log(10.0)))**2.0 )**0.5
  ratioerrhigh = ( (specZerr[1]/(specZ*np.log(10.0)))**2.0 + 
                  (refZerr[1]/(refZ*np.log(10.0)))**2.0 )**0.5


  #--Return abundance ratio and errors--
  return ratio,ratioerrlow,ratioerrhigh


#----------------------------------------------------------
def ratio_error_bars(top,bottom,top_low,toperr=0.0,bottomerr=0.0):
  """
  Author: Kari A. Frank
  Date: May 12, 2016
  
  Purpose: Calculate the error bars for a ratio.

  Input: 
     top,bottom (float) : the numerator and denominator of the ratio
     toperr,bottomerr (float or 2-element float tuple) : measurement errors
                    for top and bottom. If provided, will also return 
                    the errors on the ratio. If only one of toperr,
                    bottomerr is provided, error for the other is 
                    assumed to be zero.
                    If a tuple, then assumed to be of the form 
                    (lowerr,higherr)
                  
  """

  #--make sure given tuples are of correct form--
  if isinstance(toperr,tuple) and (len(toperr) > 2):
    print "ERROR: toperr tuple has more than 2 elements."
    return
  if isinstance(toperr,tuple) and (len(toperr) == 1):
    toperr = (toperr[0],toperr[0]) # assume symmetric errors
  if isinstance(bottomerr,tuple) and (len(bottomerr) > 2):
    print "ERROR: bottomerr tuple has more than 2 elements."
    return
  if isinstance(bottomerr,tuple) and (len(bottomerr) == 1):
    bottomerr = (bottomerr[0],bottomerr[0]) # assume symmetric errors
  

  #--convert errs to tuples if not--
  if not isinstance(toperr,tuple): toperr=(toperr,toperr)
  if not isinstance(bottomerr,tuple): toperr=(bottomerr,bottomerr)

  #--convert to float if int--
  if isinstance(top,int): top = float(top)
  if isinstance(bottom,int): bottom = float(bottom)

  if isinstance(toperr[0],int): toperr[0] = float(toperr[0])
  if isinstance(toperr[1],int): toperr[1] = float(toperr[1])

  if isinstance(bottomerr[0],int): bottomerr[0] = float(bottomerr[0])
  if isinstance(bottomerr[1],int): bottomerr[1] = float(bottomerr[1])


  #--calculate the ratio--
  ratio = top/bottom

  #--calculate ratio errors--
  ratio_lowerr = ratio - ( (toperr[0]/bottom)**2.0 + 
                        (top/bottom**2.0*bottomerr[0])**2.0 )**0.5
  ratio_higherr = ( (toperr[1]/bottom)**2.0 + 
                    (top/(bottom)**2.0*bottomerr[1])**2.0 ) **0.5 + ratio
  
  return ratio_lowerr,ratio_higherr

#----------------------------------------------------------
def error_range_to_bars(x,xlow,xhigh):
  """Function to convert given error range to error bars"""

  xerrlow = x-xlow
  xerrhigh = xhigh-x

  return xerrlow,xerrhigh
