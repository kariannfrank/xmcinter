;; Name: make_map.pro
;;
;; Date: October 29, 2015
;; Author: Kari A. Frank
;;
;; Purpose:
;;      Read file containing blob parameters and 
;;      create associated maps as fits files.
;;      (based on blob_image.pro)
;;
;; Calling Sequence:
;;
;;      make_map,infile=infile,outfile=outfile,param=paramname,weights=weights,binsize=binsize,itmod=itmod,paramshape=paramshape,type=type,points=points,x0=x0,y0=y0,imagesize=imagesize,sigthresh=sigthresh
;;       
;;
;; Input:
;;    
;;      infile:   string of the file name containing the 
;;                blob parameters, formatted as output from 
;;                the python function
;;                xmcinter.files.xmcrun.merge_output()
;;    
;;      outfile:  optional string name of the output fits file
;;
;;      paramname: string name of the column (parameter name) to be
;;                 mapped (default='blob_kT')
;;
;;      paramx,paramy: string name of the column of the x and y
;;                blob positions (default='blob_phi','blob_psi')
;;   
;;      paramshape: string specifying the shape of the blobs, 
;;                  'gauss' or 'sphere' (default='gauss')
;;
;;      paramsize: string name of the column containing the size of 
;;                 the blobs (default='blob_lnsigma' for shape='gauss'
;;                 or 'blob_radius' for shape='sphere')
;;
;;      weights: string name of the column containing the weights
;;               (default='', no weighting)
;;
;;      itmod:   set to an integer > 1 to use only every
;;               it_mod'th iteration (defaul=100)
;;
;;      binsize: size of a pixel, in same units as paramx and paramy,
;;               default=60 (fast)
;;
;;      iteration_type: 'median', 'average', or 'total'.  Determines how to 
;;                combine the blobs in each iteration to create the
;;                iteration image. (default is  'median', but note that it 
;;                should be set to 'total' if making emission measure map)
;;
;;      type:    'median', 'average', 'total', or 'error' to specify how to
;;                combine the different iteration images when making the 
;;                final image. error will produce a map of the parameter
;;                error (default='median')
;;
;;      points: switch (0 or 1) to assume that all blobs are points,
;;              i.e. ignore paramsize and skip integration over pixels
;;              (default=0)
;;
;;      imagesize: optionally specify the size of output image (length
;;                 of one side of square image) in same units as paramx,y
;; 
;;      x0,y0:  center coordinates of the desired output image, in
;;              same units as paramx,paramy
;;
;;      sigthresh: optionally specify a significance threshold (in
;;                number of sigma) for the final map.  all pixels
;;                that do not meet this significance threshold are set
;;                to zero.  note this takes much longer, as it
;;                requires calculating an error map (default = 0)
;; 
;;      outfile:  optional string of the file name to save the
;;                map. default is base name of 
;;                infile+'_'+paramname+'_'+type+'.fits'
;;
;; Output:
;;
;;     Saves a fits file in the same directory as infile, containing the
;;     calculated map.
;;
;; Usage Notes:
;;
;;
;;
;; Example:
;;
;;

PRO make_map,infile=infile,outfile=outfile,paramname=paramname,weights=weights,binsize=binsize,itmod=itmod,paramshape=paramshape,type=type,points=points,x0=x0,y0=y0,imagesize=imagesize,sigthresh=sigthresh,paramx=paramx,paramy=paramy,paramsize=paramsize,iteration_type=iteration_type

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------

IF N_ELEMENTS(infile) EQ 0 THEN BEGIN
   MESSAGE, 'ERROR: minimum usage: make_map,infile=infile'
ENDIF
IF N_ELEMENTS(points) EQ 0 THEN points = 0
IF N_ELEMENTS(paramshape) EQ 0 THEN paramshape = 'gauss'
IF N_ELEMENTS(paramsize) EQ 0 THEN BEGIN
   CASE paramshape OF 
      'gauss': paramsize = 'blob_lnsigma'
      'sphere': paramsize = 'blob_radius'
   ENDCASE
ENDIF
IF N_ELEMENTS(paramname) EQ 0 THEN paramname = 'blob_kT'
IF N_ELEMENTS(paramx) EQ 0 THEN paramx = 'blob_phi'
IF N_ELEMENTS(paramy) EQ 0 THEN paramy = 'blob_psi'
IF N_ELEMENTS(type) EQ 0 THEN type = 'median'
IF N_ELEMENTS(iteration_type) EQ 0 THEN iteration_type = 'median'
IF N_ELEMENTS(sigthresh) EQ 0 THEN sigthresh = 0.0

IF N_ELEMENTS(itmod) EQ 0 THEN itmod = 100
IF N_ELEMENTS(binsize) EQ 0 THEN binsize = 60

basename = file_trunk(infile,/PATH)
IF N_ELEMENTS(outfile) EQ 0 THEN BEGIN
   outfile = basename+'_'+paramname+'_'+type+'.fits'
ENDIF

;-------------------------------------------
;       Read in deconvolution file
;-------------------------------------------
;as output by files.xmcrun.merge_output()

parnames = []
dframe = read_blobs(infile,colnames=parnames)

;-------------------------------------------
;              Get Columns
;-------------------------------------------

pi = WHERE(parnames EQ paramname)
xi = WHERE(parnames EQ paramx)
yi = WHERE(parnames EQ paramy)
IF N_ELEMENTS(weights) NE 0 THEN BEGIN
   weighti = WHERE(parnames EQ weights)
   weights_arr = dframe[*,weighti]
ENDIF ELSE BEGIN
   weights_arr = REPLICATE(1.0,N_ELEMENTS(dframe[*,pi]))
ENDELSE
ri = WHERE(parnames EQ paramsize)
iti = WHERE(parnames EQ 'iteration')

IF points EQ 1 THEN r = 0.001 ELSE r = dframe[*,ri]

;-------------------------------------------
;              Calculate Map
;-------------------------------------------

map = calculate_map(dframe[*,pi],dframe[*,xi],dframe[*,yi],r,PARAM_ITERATIONS=dframe[*,iti],BINSIZE=binsize,ITERATION_TYPE=iteration_type,TYPE=type,IT_MOD = itmod,WEIGHTS=weights_arr,X0=x0,Y0=y0,SIZE=imagesize)

;-------------------------------------------
;(Optionally) Suppress Low Significance Pixels
;-------------------------------------------

IF sigthresh GT 0.0 THEN BEGIN
   ;-calculate error map-
   errmap = calculate_map(dframe[*,pi],dframe[*,xi],dframe[*,yi],r,PARAM_ITERATIONS=dframe[*,iti],BINSIZE=binsize,ITERATION_TYPE=iteration_type,TYPE='error',IT_MOD = itmod,WEIGHTS=weight_arr,X0=x0,Y0=y0,SIZE=imagesize)

   ;-suppress low significance pixels-
   map[WHERE(map/errmap) LT sigthresh] = 0.0
ENDIF

;-------------------------------------------
;           Write FITS File
;-------------------------------------------

;-construct history-
history = 'make_map,'+infile+',outfile='+outfile+',paramname='+paramname+',binsize='+STRING(binsize,FORMAT='(I0)')+',itmod='+STRING(itmod,FORMAT='(I0)')+',paramshape='+paramshape+',type='+type+',points='+STRING(points,FORMAT='(I0)')+',sigthresh='+STRING(sigthresh,FORMAT='(F0)')
IF N_ELEMENTS(weights) NE 0 THEN history = history + ',weights='+weights
IF N_ELEMENTS(x0) NE 0 THEN history = history + ',x0='+STRING(x0,FORMAT='(F0)')
IF N_ELEMENTS(y0) NE 0 THEN history = history + ',y0='+STRING(y0,'(F0.0)')
IF N_ELEMENTS(imagesizse) NE 0 THEN history = history + ',imagesize='+STRING(imagesize,FORMAT='(F0)')

;-initialize-
MKHDR,head,map,/IMAGE

;-break history into multiple lines-
l = STRLEN(history)
i = 0
ll = 72
WHILE l GT 0 DO BEGIN
   SXADDHIST,STRMID(history,i,ll),head
   i = i + ll
   l = STRLEN(STRMID(history,i))
ENDWHILE

;-write file-
WRITEFITS,outfile,map,head
PRINT, "Wrote "+outfile

END
