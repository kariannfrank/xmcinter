#----------------------------------------------------------
#Module of functions for plotting xmc results.
#
#Contains the following functions:
#
# chi2
# traceplots
# histogram 
#
#----------------------------------------------------------

#-import common modules-
#import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import bokeh
import bokeh.plotting as bplt
import bokeh.charts as bchart

#----------------------------------------------------------
#Author: Kari A. Frank
#Date: October 20, 2015
#Purpose: plot simple scatter plot of chi2 vs iteration, directly from 
#         the statistic.* files output by xmc.  This is mainly intended
#         as a quick way to check the convergence of an ongoing xmc run.
#
#Usage: chi2(runpath)
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
#            use linked brushing) -- not yet implemented
#
#  npoints: number of aggregated points to plot, ignored if agg='none' 
#           or agg='contour' (default = 1000.0)
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
    from bokeh.models import ColumnDataSource,PrintfTickFormatter
    
#----Aggregate Data----
    if agg=='none':
        df = dframe
    if agg=='sampling':
        df = dframe.sample(npoints)
#    if agg=='contour':

#----Set up plot----
    bplt.output_file('traceplots.html')
    source = ColumnDataSource(df)
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    w = 300 # plot width 
    h = 300 # plot height

#----Plot data---

#--initialize figure array--
    dim = len(df.columns)
    figlist=[[None]*dim for i in range(dim)]

#--loop through parameters and fill in scatter matrix (empty above diagonal)--
    for col in range(dim):
        for row in range(dim):
            # increase size if on edge (provide space for axis labels)
            if col == 0:
                wi=w+20
            else:
                wi=w
            if row == dim-1:
                he=h+20
            else:
                he=h

            if col < row:
                # create scatter plot
                newfig = bplt.figure(width=wi,height=he,tools=TOOLS)
                newfig.circle(df.columns[col],df.columns[row],webgl=True,
                             color='navy',source=source,size=1)

                # add axis labels if on edge, remove tick labels if not
                if col == 0:
                    newfig.yaxis.axis_label=df.columns[row]
                    if ((max(df[df.columns[row]]) < 0.01) or 
                        (max(df[df.columns[row]]>999)) ): 
                        tformat="%4.0e"
                    else:
                        tformat="%4.2f"
                    newfig.yaxis.formatter=PrintfTickFormatter(format=tformat)
                else:
                    newfig.yaxis.major_label_text_color = None
                    newfig.yaxis.major_label_text_font_size = '0'

                if row == dim-1:
                    newfig.xaxis.axis_label=df.columns[col]
                    if ((max(df[df.columns[col]]) < 0.01) or 
                        (max(df[df.columns[col]]>999)) ): 
                        tformat="%4.0e"
                    else:
                        tformat="%4.2f"
                    newfig.xaxis.formatter=PrintfTickFormatter(format=tformat)
                else:
                    newfig.xaxis.major_label_text_color = None
                    newfig.xaxis.major_label_text_font_size = '0'

                # add to figure array
                figlist[row][col]=newfig
            if col == row:
                # plot histogram
                newfig = bchart.Histogram(df[[col]],bins=30,
                                          width=wi,height=he,tools=TOOLS,
                                          xlabel=df.columns[col])
                # add axis label if corner, remove tick labels if not
##                if col == 0:
##                    newfig.yaxis.axis_label=df.columns[row]
##                else: 
##                    newfig.yaxis.major_label_text_font_color = None
##                if row == dim-1:
##                    newfig.xaxis.axis_label=df.columns[col]
##                else:
##                    newfig.xaxis.major_label_text_font_color = None
                # add to figure array
                figlist[row][col]=newfig
            if col > row:
                # leave plot empty
                figlist[row][col]=bplt.figure(width=wi,height=he,tools=TOOLS)

#--plot grid--
    p = bplt.gridplot(figlist)
    bplt.show(p)

#----Return----
    return figlist

#----------------------------------------------------------
#----------------------------------------------------------
#Author: Kari A. Frank
#Date: October 28, 2015
#Purpose: plot a (weighted) histogram from provided series.
#
#Usage: theplot = histogram(dataseries,weights=None,bins=30)
#
#Input:
# 
# datacolumn:  a pandas series of data (e.g. column of a pandas dataframe)
#  
# weights:     optionally provided a pandas series of weights which correspond
#              to the values in datacolumn (e.g. emission measure)
#             
# bins:        optionally specify the number of bins (default=30)
#
# save:        optionally turn off opening and saving the plot as an 
#              html file - returns the figure object only (default=True)
#
# **kwargs:    pass any number of extra keyword arguments that are 
#              accepted by bokeh.plotting.quad().  some of the most
#              useful may be fill_color and line_color
#
#Output:
# - uses bokeh to open a plot of the (weighted) histogram in a browser
# - returns the figure object (allows replotting it, e.g. in a grid with
#   other figures)
#Usage Notes:
# - must close and save (if desired) the plot manually
# - axis labels will use the pandas series names (e.g. dataseries.name)
#
#Example:
#  
#

def histogram(dataseries,weights=None,bins=30,save=True,**kwargs):

#----Import Modules----


#----Set up Plot----
    if save: 
        bplt.output_file('histogram.html')
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    w = 500
    h = 400
    fig = bplt.figure(tools=TOOLS,width=w,height=h,**kwargs)
    fig.xaxis.axis_label=dataseries.name
    if weights is not None:
        fig.yaxis.axis_label=weights.name

#---Create the weighted histogram----
    histy,binedges = np.histogram(dataseries,weights=weights,bins=bins)

#---Plot the histogram----
    h = fig.quad(top=histy,bottom=0,left=binedges[:-1],right=binedges[1:],
                 **kwargs)

    if save: 
        bplt.show(fig)#,new='window',browser='firefox --no-remote')

#----Return----
    return fig

#----------------------------------------------------------
#Author: Kari A. Frank
#Date: October 28, 2015
#Purpose: plot interactive matrix of weighted histograms from given dataframe
#
#Usage: histogram_grid(dframe,weights=None,bins=30,**kwargs):
#
#Input:
# 
# dframe:      input dataframe
#
# weights:     optionally provided a pandas series of weights which correspond
#              to the values in datacolumn (e.g. emission measure)
#             
# bins:        optionally specify the number of bins (default=30)
#
# **kwargs:    pass any number of extra keyword arguments that are 
#              accepted by bokeh.plotting.quad().  some of the most
#              useful may be fill_color and line_color
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
def histogram_grid(dframe,weights=None,bins=30,**kwargs):

#----Import Modules----
    import math

#----Set up plot----
    bplt.output_file('histogram_grid.html')

#----Initialize empty figure list----
    figlist=[]

#----Fill in list of figures----
    for column in dframe:
        newfig = histogram(dframe[column],weights=weights,save=False)
        figlist=figlist+[newfig]

#----Reshape list into a 4 column array---

    #--define new shape--
    nfigs = len(figlist)
    ncols = 4
    nrows = int(math.ceil(float(nfigs)/float(ncols)))

    #--pad list with None to have nrows*ncols elements--
    figlist = figlist+[None]*(nrows*ncols-nfigs)

    #--reshape list--
    figarr = [figlist[ncol*i:ncol*(i+1)] for i in range(nrows)]

#----Plot histograms----
    p = bplt.gridplot(figarr)
    bplt.show(p)

    return figarr
