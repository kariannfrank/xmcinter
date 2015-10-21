#----------------------------------------------------------
#Module of functions for plotting xmc results.
#
#Contains the following functions:
#
# chi2
#----------------------------------------------------------

#-import common modules-
#import matplotlib.pyplot as plt
import pandas as pd

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

def chi2(runpath):

#----Import Modules----
    from .files.xmcrun import merge_output
    import bokeh
    from bokeh.plotting import figure, output_file, show

#----Read statistic files----
    statframe = merge_output(runpath,filetype='statistic',save=False)

#----Calculate reduced chi2----
    statframe['redchi2'] = statframe['chi2']/statframe['dof']

#----Set up Plot----
    output_file('chi2_vs_iteration.html')
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    fig = figure(tools=TOOLS)
    fig.xaxis.axis_label='iteration'
    fig.yaxis.axis_label='chi2/dof'

#----Plot chi2 vs iteration----
#    statframe.plot(x='iteration',y='redchi2',kind='scatter')#
#    plt.show()

    fig.circle(x=statframe['iteration'],y=statframe['redchi2'])
    show(fig,new='window')

#----Return----
    return statframe
