;; Name: file_trunk.pro
;;
;; Date: April 8, 2014
;; Author: Kari A. Frank
;;
;; Purpose: 
;;       Read in a file name (string) and return 
;;        the file name stripped of its extension 
;;        and preceding path.
;;
;; Calling Sequence:
;;
;;       name = file_trunk(filename,/path)
;;
;; Input:
;;       filename -- string file name
;;   
;;       path     -- optional switch (0 or 1) to turn off
;;                   stripping of the path (only strips
;;                   off the extension).  default=0.
;;                   
;;                    
;;
;; Output:
;; 
;;       Returns the file name, stripped of extension
;;         path, e.g. '/path/to/myfile.fits' becomes
;;         'myfile'.
;;
;; Usage Notes:
;;
;;
;;
;; Example:
;;
;;

FUNCTION file_trunk, filename,path=path

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------

IF N_ELEMENTS(filename) EQ 0 THEN MESSAGE,'ERROR: Minimum usage: result = file_trunk(filename)\n'
IF N_ELEMENTS(path) EQ 0 THEN path = 0

;-------------------------------------------
;             Split Filename
;-------------------------------------------

;----Split Path----
IF path EQ 0 THEN fdrs = STRSPLIT(filename,'/',/EXTRACT) ELSE fdrs = [filename]

;----Separate on '.'----
file_parts = STRSPLIT(fdrs[N_ELEMENTS(fdrs)-1],'.',/EXTRACT)

;----Construct Name----
;allows for extraneous '.' in file name
name = ''
FOR p=0,N_ELEMENTS(file_parts)-2 DO BEGIN
   name = name + file_parts[p]
ENDFOR

;----Return Final Name----
RETURN, name


END
