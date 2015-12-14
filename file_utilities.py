"""
file_utilities.py

Module of miscellaneous low-level useful functions involving files.

Contains the following functions:

 ls_to_list
 fetch_file
 read_list
 file_lines
 parse_file_line
"""

#----------------------------------------------------------
#-import common modules-
import re,sys,commands
import os

#----------------------------------------------------------
def ls_to_list(search_dir,ls_args=''):
  """
  Author: Kari Frank
  Date: January 22, 2013

  Purpose: read in a path and return list of all files in that 
           directory as list of strings

  Input: 
    search_dir (str): full path to desired directory
    ls_args (optional str): arguments to the ls command

  Output:
    returns a list of strings, one element per file

  Usage Notes:

  """ 
  #-set initial directory-
  pwd = os.getcwd()
    
  #-temporarily change to desired directory-
  os.chdir(search_dir)

  #-get list of files-
  longstr = commands.getoutput("ls "+ls_args) #stores as single string
  
  #-parse into separate files-
  file_list = longstr.split('\n')

  #-return to original working directory-
  os.chdir(pwd)

  return file_list


#----------------------------------------------------------
def fetch_file(search_dir,pat='evt2',prompt=True):
  """
  Author: Kari Frank
  Date: January 22, 2013

  Purpose: get a single file from given directory that 
           matches the search string pattern. if multiple 
           matching files or none, can prompt user for
           file path.

  Input: 
    search_dir (str): full path to directory
    pat (str): regular expression containing the desired
          pattern to match (e.g. 'evt2'). default='evt2'
    prompt (bool): Only matters if multiple matching files are found.
            if set to False then will not prompt user for file
            in event of multiple matching files, will instead
            just pick the first one and print out a warning.
            can also be set to a (string)number, in which case it 
            will automatically choose that numbered file from the 
            available choices (starts with 0), as if the user
            had typed a number with prompt=True. 
            default is prompt=True. 

  Output:
    string containing full path to a single file

  Usage Notes:
   - If there is no file that matches, will result in an error.
  """

  #-create list of all files in the directory
  # which match the given pattern- 
  file_list = ls_to_list(search_dir)
  good_files = []
  for f in file_list:
    match = re.search(pat,f)
    if match: #add all matching files to the good file list
      good_files.append(search_dir+'/'+f)

  #-decide on just one file-
  if len(good_files) == 1: #if only one match found, return the path
    return good_files[0] 
  elif len(good_files) > 1: #if multiple matching files
    if prompt is False:
      print 'Warning: Multiple files found matching '+pat+' ;'
      print '         choosing '+good_files[0]
      return good_files[0]
    else:
      if prompt is True: #if prompt was not set to a number
        print 'More than one file found matching '+pat+','
        print 'please choose from the following:'
        for g in range(len(good_files)):
          print good_files[g]+' ['+str(g)+']'
        ind = input('Number of file: ')
      else: #if prompt set to a number
        ind = int(prompt)
      return good_files[ind]
  else: #if no matching files found
    if prompt is False:
      sys.exit('ERROR: no file found in '+search_dir+' matching '+pat)
    else:
      inquiry = 'No file found matching '+pat+'\nEnter full path to file: '
      file = raw_input(inquiry)
      if os.path.isfile(file)==False: #make sure given file exists
        sys.exit('ERROR: Specified file does not exist.')
      else:
        return file
#----------------------------------------------------------
def read_list(listfile,comment=''):
  """
  Author: Kari Frank
  Date: July 30, 2013

  Purpose: return a list containing one element for each line
           of the specified file

  Usage: list = read_list(listfile,comment=comment)

  Input: 
    listfile (str): path to the file to be read  
    comment (optional str): comment character.
             lines beginning with comment will be ignored.

  Output:
    - Returns a list of the file contents (one element per line).
      New line characters are stripped from the end of lines.

  Usage Notes:

  """

  #-open file-
  thefile = open(listfile)
 
  #-read in lines-
  lines = thefile.readlines()
  
  #-strip newlines-
  trimmed_lines = [line.strip() for line in lines]

  #-remove commented lines-
  if comment !='':
    good_lines = []
    for line in trimmed_lines:
      if (line[0] != comment) and (line != ''): good_lines.append(line)
  else:
    good_lines = trimmed_lines

  #-close file-
  thefile.close()

  #-return list-
  return good_lines

#----------------------------------------------------------
def file_lines(infile,comment=''):
  """
  Author: Kari Frank
  Date: March 20, 2014

  Purpose: return the number of non-commented lines in a file

  Input: 
    infile (str): path and name of the file to be read  
    comment (optional str): comment character.
             lines beginning with comment will be ignored.

  Output:
    - Returns number of non-commented lines in file 

  Usage Notes:

  """

  #--get list of lines--
#  lines = read_lines(infile,comment=comment)

  counter = 0
  with open(infile,'r') as f:
    for line in f:
      counter += 1
  
  #--size of list--
  return counter

#----------------------------------------------------------
def parse_file_line(infile,delimiter=',;\s\t'):
  """
  Author: Kari A. Frank 
  Date: March 20, 2014

  Purpose: Read a first line from a file and parse it into a string list.

  Input:
    infile (str): Path to text file.
    delimiter (optional str): delimiter(s).  Default is space,
                  comma, tab, and semicolon.  Can provide more than
                  one delimiter character in the string.  Spaces and 
                  tabs should be given in regular expression syntax,
                  \s and \t.

  Output:  
    returns string list of 'words' in the line.

  Usage Notes:

  """

  #--set delimiter string--
  delstr = r'[^'+delimiter+']+'

  #--open file--
  thefile = open(infile)

  #--read line--
  rawline = thefile.readlines()

  #--close file--
  thefile.close()

  #--strip newline--
  line = rawline[0].strip()

  #--parse line--
#    outpars = line.split()
  outlist = re.findall(delstr,line)

  #--return string list--
  return outlist

#----------------------------------------------------------
