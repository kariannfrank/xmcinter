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
import bokeh
import bokeh.plotting as bplt

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
#    import bokeh
#    from bokeh.plotting import figure, output_file, show
    #from bokeh.mpl import to_bokeh

#----Read statistic files----
    statframe = merge_output(runpath,filetype='statistic',save=False)

#----Calculate reduced chi2----
    statframe['redchi2'] = statframe['chi2']/statframe['dof']

#----Set up Plot----
    bplt.output_file('chi2_vs_iteration.html')
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    fig = bplt.figure(tools=TOOLS)
    fig.xaxis.axis_label='iteration'
    fig.yaxis.axis_label='chi2/dof'

#----Plot chi2 vs iteration----
#    statframe.plot(x='iteration',y='redchi2',kind='scatter')
#    plt.show()

    fig.circle(x=statframe['iteration'],y=statframe['redchi2'])
    bplt.show(fig)#,new='window',browser='firefox --no-remote')

#----Return----
    return statframe

#----------------------------------------------------------
#Author: Kari A. Frank
#Date: October 26, 2015
#Purpose: plot interactive matrix of scatter plots from given dataframe
#
#Usage: traceplots(dframe,agg='sampling',npoinst=1000.0)
#
#Input:
# 
# dframe: input dataframe
#
# agg:    type of aggregration to perform before plotting, since plotting with
#         all possible rows/blobs is generally prohibitively slow options are
#         - 'none' (plot all rows)
#         - 'sampling' (take random subsample of rows, default)
#         - 'contour' (plot density contours instead of points - does not 
#            use linked brushing)
#
#  npoints: number of aggregated points to plot, ignored if agg='none' 
#           (default = 1000.0)
#           
#             
#Output:
# - plots matrix of scatter plots of all provided dataframe columns to 
#   interactive (browser-based) plot
#
#Usage Notes:
# - must close and save (if desired) the plot manually
#
#Example:
# 
#

def traceplots(dframe,agg='sampling',npoints=1000.0):

#----Import Modules----
from bokeh.io import gridplot

#----Aggregate Data----


#----Set up plot----
    output_file('traceplots.html')
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    fig = figure(tools=TOOLS)
    w = 50 # plot width 
    h = 50 # plot height
#    fig.xaxis.axis_label='iteration'
#    fig.yaxis.axis_label='chi2/dof'

#----Plot data---

#initialize plot array (half will be empty, across a diagonal)

#--loop through parameters--
#see loop in xmc_diagnostics
    for column in df.columns:

    #-create new plot-
        newfig = bplt.figure(width=w,height=h)
        newfig.circle(column1,column2,color='navy',size=1)
    #also set axis labels
    #add to array of figures

#--plot grid--
    p = gridplot(plotarray)
    bplt.show(p)

#----Return----
    return True
