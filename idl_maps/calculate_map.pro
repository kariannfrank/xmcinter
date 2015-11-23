;Name: create_map.pro
;
;Date: January 18, 2011
;Author: Kari Frank
;
;Purpose: Read in pipeline results and create map of blob parameter(s)
;         
;Calling Sequence: 
;         result = create_map(param,param_x,param_y,sigma[,PARAM_ITERATIONS=param_iterations][,WEIGHTS=weights][,BINSIZE=binsize][,ITERATION_TYPE=iteration_type][,TYPE=type][,SIZE=size][,IT_MOD=it_mod][,N_INT_STEPS=n_int_steps][,VERBOSE=verbose][,OUTPUT_SIZE=output_size][,X0=x0][,Y0=y0])
;      
;
;Input:
;        param: vector of blob parameters from which to create a map
;        param_x: vector of blob x positions
;        param_y: vector of blob y positions
;        param_sigma: vector of blob sizes (in arcsec)
;        binsize: width of pixels, in same units as param_x,param_y, 
;                 and param_sigma (assumes square pixel, default is 60.)
;        param_iterations: array of the iteration associated with each
;                blob in param array. if not provided, assumes all
;                blobs are from the same iteration.
;        iteration_type: type=median, average, or total.  Determines how to 
;                combine the blobs in each iteration to create the
;                iteration image. (default is  median, but note that it 
;                should be set to 'total' if making emission measure map)
;        type: type=median, average, total, or error to specify how to
;                combine the different iteration images when making the 
;                final image. error will produce a map of the parameter
;                error (default=median)
;        weights: array of weights for each blob (e.g. emission measures)
;        size: length of one side of image, assumes square image 
;              ( default=max(param_x)-min(param_x) ).
;        it_mod: set to >1 to use only every it_mod'th
;                iteration (default = 10)
;        x0,y0: cluster position, in same coordinate system as param_x and y
;                (typically xmc coordinates, default is
;                x0=median(param_x),y0=median(param_y)).
;        n_int_steps: number of steps in each integration loop (default = 50)
;        verbose: switch to turn on printing of checkpoints (default = 0)
;
;Output:
;        returns a 2D float array (image) of given parameter      
;        output_size: optional float to store the final image size for
;               later plotting, in case it has changed from the input value.
; 
;Usage Notes:
;      - param_x, param_y, param_sigma, binsize, and size must all 
;        be in the same units (usually arcsec)
;      - creates square map that is centered on x0,y0
;
;Example:


;-----------------------------------------------
;             HELPER FUNCTIONS
;-----------------------------------------------

;Name: weighted_average
;
;Date: November 16, 2010
;Author: Kari Frank
;
;Purpose: Calculates the emission-measure weighted average of a parameter
;Calling Sequence: 
;      weighted_average(parameter,weight)
;
;Input:
;      parameter: array of whichever parameter is to be average,
;                  one element for each blob.  Must correspond
;                  to the same elements in the weight array
;      weight: array of weights (e.g. 10.0^-14*emission measure) for each blob
;Output:
;     Returns the weighted averaged value of parameter.


FUNCTION weighted_average, parameter, weight

IF N_ELEMENTS(parameter) NE N_ELEMENTS(weight) THEN MESSAGE,'ERROR: parameter and weight arguments are not the same size.'

averaged=TOTAL(parameter*weight)/TOTAL(weight)

RETURN, averaged

END

;-----------------------------------------------
;Name: weighed_median
;
;Date: April 5, 2011
;Author: Kari Frank
;
;Purpose: Calculates the weighted median of a parameter, without the
;          need for binning, by calculating the 'probability' of each data point
;          (e.g. each blob).
;Calling Sequence: 
;      weighted_median(parameter,weight[,PROBABILITIES=probabilities][,TWO_SIGMA=two_sigma][,ONE_SIGMA=one_sigma][,NOZERO=nozero])
;
;Input:
;      weight:    array containing the blob weights (may be unnormalized)
;      parameter: array that the median will be found for,
;                  one element for each blob.  Must correspond
;                  to the same elements in the weight array.
;      nozero:    optional switch to ignore zeros (default=0)
;Output:
;     Returns the weighted median of parameter.
;     probabilities: optional array that stores the probabilities
;       corresponding each parameter value.
;     two_sigma: optional provided variable to store the location
;       beyond which (to the right) lie 95.44% of the blobs.   
;     one_sigma: optional provided variable to store the location
;       beyond which (to the right) lie 68.3% of the blobs.


FUNCTION weighted_median, parameter, weight,PROBABILITIES=probabilities,TWO_SIGMA=two_sigma, ONE_SIGMA=one_sigma,NOZERO=nozero

;check arguments, set defaults
IF N_PARAMS() LT 2 THEN MESSAGE,'ERROR: minimum usage is weighted_median(parameter,weight)'
IF N_ELEMENTS(parameter) NE N_ELEMENTS(weight) THEN MESSAGE,'ERROR: parameter and weight arguments are not the same size.'
IF N_ELEMENTS(nozero) EQ 0 THEN nozero = 0

;IF N_ELEMENTS(probabilities) EQ 0 THEN probabilities = FLTARR(N_ELEMENTS(parameter))

IF nozero EQ 1 THEN BEGIN
   nonzero_subs = WHERE(parameter NE 0.0)
   IF nonzero_subs[0] NE -1 THEN BEGIN
      trimmed_parameter = parameter[nonzero_subs]
      trimmed_weight = weight[nonzero_subs]
   ENDIF ELSE BEGIN ;if array zeros, keep zeros
      PRINT, 'Warning: Array all zeros.'
      trimmed_parameter = parameter
      trimmed_weight = weight
   ENDELSE 
ENDIF ELSE BEGIN
   trimmed_parameter = parameter
   trimmed_weight = weight
ENDELSE

;find probabilities the fast way (by sorting)
IF N_ELEMENTS(probabilities) EQ 0 THEN BEGIN
   sorted_indices = SORT(trimmed_parameter)
   probability = TOTAL(trimmed_weight[sorted_indices],/CUMULATIVE)/TOTAL(trimmed_weight)
   sorted_parameter=trimmed_parameter[sorted_indices]
ENDIF ELSE BEGIN

;find probabilities the slow way 
;    (preserves kT - probability vector correspondence)
   probability = FLTARR(N_ELEMENTS(trimmed_parameter))
   FOR b=0L, N_ELEMENTS(trimmed_parameter)-1 DO BEGIN

      ;;subscript of all blobs with lower value (e.g. lower temperatures)
      lower_subs = WHERE(trimmed_parameter LE trimmed_parameter[b])

      ;;probability = cumulative normalized weights
      probability[b] = TOTAL(trimmed_weight[lower_subs])/TOTAL(trimmed_weight)
   ENDFOR
ENDELSE

;PRINT, 'probability = ',probability
;find subscript of parameter with probability closest to 0.5
half_sub = WHERE( ABS(0.5-probability) EQ MIN(ABS(0.5-probability)) ) 

;find subscript of 2sigma location (95.44% lie above this location)
twosigma_sub = WHERE( ABS((1.0-0.9544)-probability) EQ MIN(ABS((1.0-0.9544)-probability)) ) 

;find subscript of 1sigma location (68.3% lie above this location)
onesigma_sub = WHERE( ABS((1.0-0.683)-probability) EQ MIN(ABS((1.0-0.683)-probability)) ) 

;median = parameter value of above parameter
IF half_sub[0] NE -1 THEN BEGIN
   IF N_ELEMENTS(probabilities) NE 0 THEN BEGIN
      median = trimmed_parameter[half_sub[0]] 
   ENDIF ELSE BEGIN
      median = sorted_parameter[half_sub[0]] 
   ENDELSE
ENDIF ELSE BEGIN 
   median = 0.0
ENDELSE

;95.44% location
IF twosigma_sub[0] NE -1 THEN BEGIN
   IF N_ELEMENTS(probabilities) NE 0 THEN BEGIN
      two_sigma = trimmed_parameter[twosigma_sub[0]] 
   ENDIF ELSE BEGIN
      two_sigma = sorted_parameter[twosigma_sub[0]] 
   ENDELSE
ENDIF ELSE BEGIN 
   two_sigma = 0.0
ENDELSE

;95.3% location
IF onesigma_sub[0] NE -1 THEN BEGIN
   IF N_ELEMENTS(probabilities) NE 0 THEN BEGIN
      one_sigma = trimmed_parameter[onesigma_sub[0]] 
   ENDIF ELSE BEGIN
      one_sigma = sorted_parameter[onesigma_sub[0]] 
   ENDELSE
ENDIF ELSE BEGIN 
   one_sigma = 0.0
ENDELSE


;save probabilities in output vector if provided
IF N_ELEMENTS(probabilities) NE 0 THEN probabilities[nonzero_subs] = probability

RETURN, median

END
;-----------------------------------------------

;-----------------------------------------------
;               MAIN FUNCTION
;-----------------------------------------------

FUNCTION calculate_map,param,param_x,param_y,param_sigma,PARAM_ITERATIONS=param_iterations,BINSIZE=binsize,ITERATION_TYPE=iteration_type,TYPE=type,FRACTIONS=fractions,SIZE=size,WEIGHTS=weights,IT_MOD=it_mod,VERBOSE=verbose,OUTPUT_SIZE=output_size,N_INT_STEPS=n_int_steps,X0=x0,Y0=y0

;check for required arguments and set defaults
IF N_PARAMS() LT 4 THEN MESSAGE, "Minimum Usage: result = calculate_map(param,param_x,paramy_,param_sigma)"
IF N_ELEMENTS(param_iterations) EQ 0 THEN param_iterations = REPLICATE(0.0,N_ELEMENTS(param))
IF (N_ELEMENTS(binsize) EQ 0) THEN binsize = 60.0
IF N_ELEMENTS(type) EQ 0 THEN type='median'
IF (type NE 'median') AND (type NE 'average') AND (type NE 'total') AND (type NE 'error') THEN BEGIN
   PRINT, 'Unrecognized TYPE. Using TYPE=median'
   type='median'
ENDIF
IF (iteration_type NE 'median') AND (iteration_type NE 'average') AND (iteration_type NE 'total') THEN BEGIN
   PRINT, 'Unrecognized ITERATION_TYPE. Using ITERATION_TYPE=median'
   type='median'
ENDIF
IF N_ELEMENTS(weights) EQ 0 THEN weights = REPLICATE(1.0,N_ELEMENTS(param))
IF N_ELEMENTS(X0) EQ 0 THEN x0 = MEDIAN(param_x)
IF N_ELEMENTS(Y0) EQ 0 THEN y0 = MEDIAN(param_y)
IF N_ELEMENTS(log_contours) EQ 0 THEN log_contours = 0
IF N_ELEMENTS(it_mod) EQ 0 THEN it_mod = 10
IF N_ELEMENTS(n_int_steps) EQ 0 THEN n_int_steps = 50
IF N_ELEMENTS(verbose) EQ 0 THEN verbose = 0
IF N_ELEMENTS(size) EQ 0 THEN  size=1.2*(MAX(param_x)-MIN(param_x))

;set up map parameters
image_radius = size/2
xmin = x0 - image_radius
xmax = x0 + image_radius
ymin = y0 - image_radius
ymax = y0 + image_radius
;trim to only include area containing blobs
;IF xmin LT MIN(param_x) THEN xmin = MIN(param_x)-MAX(param_sigma);1.2*MIN(param_x)
;IF xmax GT MAX(param_x) THEN xmax = MAX(param_x)+MAX(param_sigma);1.2*MAX(param_x)
;IF ymin LT MIN(param_y) THEN ymin = MIN(param_y)-MAX(param_sigma);1.2*MIN(param_y)
;IF ymax GT MAX(param_y) THEN ymax = MAX(param_y)+MAX(param_sigma);1.2*MAX(param_y)

PRINT,'xmin,xmax = ',xmin,xmax

image_radius = (xmax-xmin)/2.0
size = image_radius*2.0
output_size = size
;PRINT, 'param_iterations = ',param_iterations
n_layers = (MAX(param_iterations) - MIN(param_iterations))/it_mod + 1
;PRINT, 'n_layers = ',n_layers
nbins_x = FLOOR((xmax - xmin)/binsize)
nbins_y = FLOOR((ymax - ymin)/binsize)
;PRINT, 'it_mod = ',it_mod

;subtract layers for missing iterations
FOR it = MIN(param_iterations),MAX(param_iterations) DO BEGIN
   IF it MOD it_mod EQ 0 THEN BEGIN
      subs = WHERE(param_iterations EQ it)
      IF subs[0] EQ -1 THEN n_layers = n_layers -1 
      IF subs[0] EQ -1 THEN PRINT, 'iteration missing'
   ENDIF 
ENDFOR

PRINT, 'n_layers = ',n_layers

PRINT, 'max, min param = ',max(param),min(param)

;create map
blob_volumes = (2.0*!PI*param_sigma^2.0)^1.5
image_stack = MAKE_ARRAY(nbins_x,nbins_y,n_layers) ;stack of images, one per iteration (do i need only need this for emission measure?)

;calculate fraction of each blob in each bin
layer = 0
FOR i = MIN(param_iterations),MAX(param_iterations) DO BEGIN ;loop over iterations (layers)
   IF i MOD it_mod EQ 0 THEN BEGIN ;only do every it_mod'th iteration
      iteration_subs = WHERE(param_iterations EQ i)
      IF iteration_subs[0] NE -1 THEN BEGIN ;make sure iteration existed

      FOR x=0L,nbins_x-1 DO BEGIN ;loop over x
         lower_x = xmin + x*binsize
         upper_x = xmin + x*binsize + binsize
         dx = (upper_x - lower_x)/FLOAT(n_int_steps)
;         PRINT, 'dx = ',dx
         x_blob_integrals = FLTARR(N_ELEMENTS(param[iteration_subs]))
         FOR xp = FLOAT(lower_x),FLOAT(upper_x),dx DO BEGIN ;integrals over x
            x_blob_integrals = x_blob_integrals + EXP(-1.0/2.0*(xp-param_x[iteration_subs])^2.0/param_sigma[iteration_subs]^2.0)*dx
         ENDFOR 
         FOR y=0L,nbins_y-1 DO BEGIN ;loop over y
            lower_y = ymin + y*binsize
            upper_y = ymin + y*binsize + binsize
            dy = (upper_y - lower_y)/FLOAT(n_int_steps)
            y_blob_integrals = FLTARR(N_ELEMENTS(param[iteration_subs])) 
            FOR yp=FLOAT(lower_y),FLOAT(upper_y),dy DO BEGIN ;integral over y
               y_blob_integrals = y_blob_integrals + EXP(-1.0/2.0*(yp-param_y[iteration_subs])^2.0/param_sigma[iteration_subs]^2.0)*dy
            ENDFOR 
            fractions = x_blob_integrals*y_blob_integrals*(2.0*!PI*param_sigma[iteration_subs]^2.0)^0.5 / blob_volumes[iteration_subs] ;multiply by dz integral, divide by total volume
            CASE iteration_type OF
               'average': image_stack[x,y,layer] = weighted_average(param[iteration_subs],weights[iteration_subs]*fractions)
               'median': image_stack[x,y,layer] = weighted_median(param[iteration_subs],weights[iteration_subs]*fractions)
               'total': image_stack[x,y,layer] = TOTAL(param[iteration_subs]*weights[iteration_subs]*fractions)

            ENDCASE 
;            IF (verbose EQ 1) AND (image_stack[x,y,layer] EQ 0) THEN PRINT, 'x,y,layer,image_stack[x,y,layer] = ',x,y,layer,image_stack[x,y,layer]
         ENDFOR  ;endfor y
      ENDFOR  ;endfor x
      IF verbose EQ 1 THEN PRINT, 'layer ',layer,' of ',n_layers-1
      IF verbose EQ 1 THEN PRINT, 'min,max,median of image_stack[layer] = ',MIN(image_stack[*,*,layer]),MAX(image_stack[*,*,layer]),MEDIAN(image_stack[*,*,layer])
;      IF verbose EQ 1 THEN PLOTHIST,fractions,/AUTOBIN,TITLE='layer '+STRING(layer)+' fractions'
      layer = layer + 1
      ENDIF  ;iteration existed
   ENDIF ;endif it_mod
ENDFOR ;endfor i
IF verbose EQ 1 THEN PRINT, 'min, max, median of image_stack = ',MIN(image_stack),MAX(image_stack),MEDIAN(image_stack)
IF verbose EQ 1 THEN PRINT, '---image stack created---'

;collapse the image_stack array
image = MAKE_ARRAY(nbins_x,nbins_y)
iteration_weights = REPLICATE(1.0,N_ELEMENTS(image_stack[0,0,*])) ;weight each iteration uniformly
FOR x_pixel = 0L,nbins_x-1 DO BEGIN ;loop through x pixels
   FOR y_pixel = 0L,nbins_y-1 DO BEGIN
      CASE type OF
         'average': image[x_pixel,y_pixel] = weighted_average(image_stack[x_pixel,y_pixel,*],iteration_weights)
         'median': image[x_pixel,y_pixel] = weighted_median(image_stack[x_pixel,y_pixel,*],iteration_weights,/NOZERO)
         'total': image[x_pixel,y_pixel] = TOTAL(image_stack[x_pixel,y_pixel,*]*iteration_weights)
         'error': image[x_pixel,y_pixel] = weighted_stdev(image_stack[x_pixel,y_pixel,*],iteration_weights)
      ENDCASE 
      ;PRINT, 'image_stack[x_pixel,y_pixel,*]',image_stack[x_pixel,y_pixel,*]
;      IF image_stack[x_pixel,y_pixel,*] NE 0 THEN PRINT, 'image_stack[x_pixel,y_pixel,*]',image_stack[x_pixel,y_pixel,*]
;      IF (verbose EQ 1) AND (image[x_pixel,y_pixel] EQ 0) THEN PLOTHIST,image_stack[x_pixel,y_pixel,*],BIN=0.2,TITLE='histogram of layer '+STRING(layer)
   ENDFOR ;endfor y_pixel
ENDFOR ;endfor x_pixel
IF verbose EQ 1 THEN PRINT, 'min,max,median of image = ',min(image),max(image),median(image)
IF verbose EQ 1 THEN PRINT, '---image stack collapsed---'

;PRINT, 'number zeros in stack = ',N_ELEMENTS(WHERE(image_stack EQ 0.))
;PRINT, 'number elements in stack = ',N_ELEMENTS(image_stack)

;PRINT, 'image = ', image

RETURN, image

END



