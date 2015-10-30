"""
Author: Kari A. Frank
Date: October 29, 2015
Purpose: A wrapper to call the IDL function make_map.pro, which creates a
         map of the given blob parameters and writes it to a fits file.

Usage:
     make_map(infile,outfile=outfile,param=param,weights=weights,binsize=binsize,
            itmod=itmod,paramshape=paramshape,type=type,points=points,x0=x0,y0=y0,
            imagesize=imagesize,sigthresh=sigthresh)


Input: (same as input to make_map.pro)
    
      infile:   string of the file name containing the 
                blob parameters, formatted as output from 
                the python function
                xmcinter.files.xmcrun.merge_output()
    
      outfile:  optional string name of the output fits file

      param:    string name of the column (parameter name) to be
                mapped (default='blob_kT')

      paramx,paramy: string name of the column of the x and y
                blob positions (default='blob_phi','blob_psi')
   
      paramshape: string specifying the shape of the blobs, 
                  'gauss' or 'sphere' (default='gauss')

      paramsize: string name of the column containing the size of 
                 the blobs (default='blob_lnsigma' for shape='gauss'
                 or 'blob_radius' for shape='sphere')

      weights: string name of the column containing the weights
               (default='', no weighting)

      itmod:   set to an integer > 1 to use only every
               it_mod'th iteration (defaul=100)

      binsize: size of a pixel, in same units as paramx and paramy,
               default=60 (fast)

      iteration_type: 'median', 'average', or 'total'.  Determines how to 
                combine the blobs in each iteration to create the
                iteration image. (default is  'median', but note that it 
                should be set to 'total' if making emission measure map)

      mtype:    'median', 'average', 'total', or 'error' to specify how to
                combine the different iteration images when making the 
                final image. error will produce a map of the parameter
                error (default='median')

      points: switch (0 or 1) to assume that all blobs are points,
              i.e. ignore paramsize and skip integration over pixels
              (default=0)

      imagesize: optionally specify the size of output image (length
                 of one side of square image) in same units as paramx,y
 
      x0,y0:  center coordinates of the desired output image, in
              same units as paramx,paramy

      sigthresh: optionally specify a significance threshold (in
                number of sigma) for the final map.  all pixels
                that do not meet this significance threshold are set
                to zero.  note this takes much longer, as it
                requires calculating an error map (default = 0)
 
      outfile:  optional string of the file name to save the
                map. default is base name of 
                infile+'_'+paramname+'_'+type+'.fits'

Output:

     Saves a fits file in the same directory as infile, containing the
     calculated map.

Usage Notes
     - Defaults are set in make_map.pro, if None is passed for the argument.
     
Example:


"""

def make_map(infile,outfile=None,param=None,weights=None,binsize=None,itmod=None,
             paramshape=None,type=None,points=None,x0=None,y0=None,imagesize=None,
             sigthresh=None,paramx=None,paramy=None,paramsize=None,iteration_type=None,mtype=None):

#----Import Modules----
    import os
    
#----Construct command line string----
#cmd = "idl -quiet -e expansion_fit -args infile='"+args.radiusfile+"',band='"+b+"',mincounts="+str(args.mincounts)

#    cmd = "idl -quiet -e make_map -args infile='"+infile+"'"
    cmd = "idl -quiet -e 'make_map,infile='"+infile+"'"
    if outfile != None:
        cmd = cmd+",outfile='"+outfile+"'"
    if param != None:
        cmd = cmd+",param='"+param+"'"
    if paramx != None:
        cmd = cmd+",paramx='"+paramx+"'"
    if paramy != None:
        cmd = cmd+",paramy='"+paramy+"'"
    if paramsize != None:
        cmd = cmd+",paramsize='"+paramsize+"'"
    if paramshape != None:
        cmd = cmd+",paramshape='"+paramshape+"'"
    if weights != None:
        cmd = cmd+",weights='"+weights+"'"
    if binsize != None:
        cmd = cmd+",binsize="+str(binsize) 
    if itmod != None:
        cmd = cmd+",itmod="+str(itmod)
    if mtype != None:
        cmd = cmd+",type='"+mtype+"'" 
    if iteration_type != None:
        cmd = cmd+",iteration_type='"+iteration_type+"'"
    if imagesize != None:
        cmd = cmd+",imagesize="+str(imagesize)
    if (points != None) and (points != False) and (points != 0):
        cmd = cmd+",points=1" 
    if x0 != None:
        cmd = cmd+",x0="+str(x0) 
    if y0 != None:
        cmd = cmd+",y0="+str(y0)
    if sigthresh != None:
        cmd = cmd+",sigthresh="+str(sigthresh)
    
    cmd = cmd+"'"
        
#----Execute command----
    print "Executing: "+cmd
    os.system(cmd)

    return True
