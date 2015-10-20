#----------------------------------------------------------
#Author: Kari A. Frank
#Date: October 20, 2015
#Purpose: plot simple scatter plot of chi2 vs iteration, directly from 
#         the statistic.* files output by xmc.  This is mainly intended
#         as a quick way to check the convergence of an ongoing xmc run.
#
#Usage: chi2(runpath,itmin=itmin,itmax=itmax)
#
#Input:
# 
# runpath:     string containing the relative path to xmc run folder, which
#              contains the statistic.* files.
#
# itmin/itmax: optionally provide the minimum and maximum iterations to include
#             
#Output:
# - plots chi2 vs iteration for the specified iteration to an interactive plot
#
#Usage Notes:
# - must close and save (if desired) the plot manually
#
#Example:
# 
#

#----Import Modules----

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ..files.xmcrun.merge_output

#----BEGIN FUNCTION----

def chi2(runpath,itmin=0,itmax=None):

#----Read statistic files----
    blobframe = merge_output(runpath,filetype='statistic',save=False)

#----Plot chi2 vs iteration----
    blobframe.plot(x='iteration',y='chi2')

#----Return----
    return True
