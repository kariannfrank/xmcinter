;; Name: read_blobs.pro
;;
;; Date: April 8, 2014
;; Author: Kari A. Frank
;;
;; Purpose:
;;
;;     Read in blob parameters from file (as saved by 
;;     python function xmcinter.files.xmcrun.merge_output() ) and
;;     return as 2D array with column names.
;;
;; Calling Sequence:
;;
;;    blob_array = read_blobs(blobfile[,colnames = colnames])
;;
;; Input:
;;
;;    blobfile -- file written by write_blobs()
;;   
;; Output:
;;
;;    Returns a 2D array of blob parameters, with 
;;       blobs=rows, parameters=columns.
;;
;;    colnames  -- optionally return the parameter names
;;                  as a list of strings
;;
;; Usage Notes:
;;
;;
;;
;; Example:
;;
;;

FUNCTION  read_blobs,blobfile,colnames=colnames

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------

IF N_ELEMENTS(blobfile) EQ 0 THEN MESSAGE, 'ERROR: Minimum usage: blob_arr = read_blobs(blobfile)\n'

;-------------------------------------------
;            Read Header Line
;-------------------------------------------

;-------------------------------------------
;             Read Blobfile
;-------------------------------------------

outarr = readtext(blobfile,head)
head = STRSPLIT(head,/EXTRACT)
;colnames = head[1:N_ELEMENTS(head)-1]
colnames = ['index',head]

RETURN, TRANSPOSE(outarr)

END
