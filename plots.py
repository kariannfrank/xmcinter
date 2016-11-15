"""
Module of functions for plotting xmc results.

Contains the following functions:

 get_range
 chi2
 traceplots
 scatter
 histogram 
 histogram_grid
 format_ticks
 spectra
 evolution
"""

#-import common modules-
import pandas as pd
import numpy as np
import bokeh
import bokeh.plotting as bplt
import bokeh.charts as bchart
from bokeh.layouts import gridplot

#----------------------------------------------------------
def logaxis(minval,maxval,limit=2.0):
    """
    Given minimum and maximum values of axis range, decide if log or linear
    scale should be used.
    
    limit specifies how many orders of magnitude must separate the 
    min and max values before switching to log scale.

    Returns either True (log axis) or False (linear axis)
    """
#    rng = (dataseries.min(),dataseries.max())
    norders = np.log10(maxval) - np.log10(minval)
    if norders > float(limit):
        return True
    else:
        return False

#----------------------------------------------------------
def chi2(runpath='./',itmin=0,itmax=None,outfile='chi2_vs_iteration.html'):
    """
    Author: Kari A. Frank
    Date: October 20, 2015
    Purpose: plot simple scatter plot of chi2 vs iteration, directly from 
             the statistic.* files output by xmc.  This is mainly intended
             as a quick way to check the convergence of an ongoing xmc run.

    Usage: chi2(runpath='./')

    Input:

     runpath:     string containing the relative path to xmc run folder, 
                  which contains the statistic.* files.

     itmin/itmax: minimum/maximum iteration to plot
     
     outfile:     name of html file to save to, or set outfile='notebook'
                  to plot inline in current notebook.

    Output:
     - plots chi2 vs iteration for the specified iteration to an \
       interactive plot

    Usage Notes:
     - must close and save (if desired) the plot manually
     - if outfile='notebook', then MUST have called bplt.output_notebook()
       in current notebook session before calling this function, else it 
       will not display the plot inline.

    Example:

    """

#----Import Modules----
    from xmcfiles import merge_output
    from wrangle import filterblobs
#    import bokeh
#    from bokeh.plotting import figure, output_file, show
    #from bokeh.mpl import to_bokeh

#----Read statistic files----
    statframe = merge_output(runpath,filetype='statistic',save=False)
    if itmax is None:
        itmax = np.max(statframe['iteration'])
    statframe = filterblobs(statframe,'iteration',minvals=itmin,maxvals=itmax)

#----Calculate reduced chi2----
    statframe['redchi2'] = statframe['chi2']/statframe['dof']

#----Set up Plot----
    if outfile != 'notebook':    
        bplt.output_file(outfile)
    else:
        bplt.output_notebook()
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    fig = bplt.Figure(tools=TOOLS)
    fig.xaxis.axis_label='iteration'
    fig.yaxis.axis_label='chi2/dof'

#----Plot chi2 vs iteration----
#    statframe.plot(x='iteration',y='redchi2',kind='scatter')
#    plt.show()

    fig.circle(x=statframe['iteration'],y=statframe['redchi2'])
    bplt.show(fig)#,new='window',browser='firefox --no-remote')
    if outfile != 'notebook': bplt.curdoc().clear()

#----Return----
    return statframe

#----------------------------------------------------------
def scatter(inframe,x,y,sampling=2000.,agg=None,aggcol=None,save=True,
            width=600,height=600,source=None,tools=None,size=5,
            xlog='auto',ylog='auto',outfile=None,returnfunc=False):
    """
    scatter()
 
   Author: Kari A. Frank
    Date: March 30, 2016
    Purpose: plot simple scatter plot of the two columns

   Usage: scatter(inframe,x,y)

   Input:
 
    inframe (DataFrame):  pandas DataFrame containing the data columns to 
         be plotted   

    x,y (strings): name of the columns in inframe to plot
             
    agg (string): type of aggregration to perform before plotting, since 
         plotting all blobs individually is generally prohibitively 
         slow. options are
         - None = plot every blob as individual point (does not require
           Datashader, can output straight to html file)
         - Datashader options:
             - 'dscount' = count elements in each plot bin
             - 'dsany' = zero if bin is empty, 1 if not
             - 'dssum','dsmin','dsmax','dsmean','dsvar','dsstd' = 
               take the sum, mean, max, etc, of all elements in a 
               plot bin. must also provide the aggcol argument
               to specify which column should be used to calculate
               the min, max, etc.
             - 'dscount_cat' = count all elements in bin, grouped by 
               category. must also provide the aggcol argument
               to specify which column should be used to group by. 
               said column must have categorical data type
                 -- may not work yet

    aggcol (string): name of the column that should be passed to 
            Datashader agg function. required for several of the 
            agg options (see above)

    sampling: number of blobs to include in plot (default = 2000.0). If
            set to None will include all blobs.

    save:   optionally turn off automatic display of the plot 
            (and saving the plot as an html file if not using Datashader)
            Returns the figure object only if save=False (default=True)
            -should rename this to display or show (also in other functions)

    width/height: optionally specify the size of the figure (default=300)

    xlog/ylog: make x and y axis scales log. options are:
                  - 'auto' (default) -- try to automatically determine
                    which is best based on the data range
                  - True -- force log scale
                  - False -- force linear scale

    size: optionally specify size of data points (default=5)

    source: a ColumnDataSource object specifying the data source (to use for
            linked brushing)

    tools: optionally pass plot tools

    outfile:  name of html file to save to, or set outfile='notebook'
              to plot inline in current notebook. setting to None
              will result in file name containing the plotting parameter
              names.

     returnfunc: bool to specify returning of the image_callback function.
                 returned value will be fig,source,image_callback. If True,
                 will disable the InteractiveImage call within scatter().
                 This allows dynamic zooming by calling 
                 InteractiveImage(fig,image_callback) on the python
                 command line. Only enable this if calling 
                 scatter() directly.


   Output:
   - plots x vs y to an interactive plot 
   - if save=True, then plot is displayed.
        - if outfile = 'notebook', then displayed in current
          notebook. In this case, MUST have called 
          bplt.output_notebook() in notebook session before 
          calling this function.
        - if outfile = anything else, then displayed (and saved)
          as html file
   - returns the figure object and the ColumnDataSource object

  Usage Notes:
  - must close and save (if desired) the plot manually
  - may not work properly if trying to plot too many points, unless
    using Datashader
  - if using Datashader agg options:
    - MUST BE RUNNING FROM A JUPYTER NOTEBOOK!! bokeh plot output 
      will be set to output_notebook() and an html file will not be
      created or saved.
    - a jupyter notebook can be started by simply typing 
      'jupyter notebook' on the command line (assuming you have 
      jupyter notebooks installed)

  Example:
 
  """
#    print len(inframe[x]),len(inframe[y])
    
#----Import Modules and Parse Datashader Functions----

    from bokeh.models import ColumnDataSource
    if agg is not None:
        from datashader.bokeh_ext import InteractiveImage
        import datashader as ds 
        import datashader.transfer_functions as tf 

        # set up dictionary to fetch ds agg functions
        fdict = {'dscount':ds.count,'dsany':ds.any,'dssum':ds.sum,'dsmin':
                 ds.min,'dsmax':ds.max,'dsmean':ds.mean,'dsvar':ds.var,
                 'dsstd':ds.std,'dscount_cat':ds.count_cat,
                 'dssummary':ds.summary}
        if agg not in fdict:
            print "Warning: "+agg+" not a valid agg option. Using agg='dscount' instead"
            agg = 'dscount'

        if aggcol is not None:
            from matplotlib.cm import jet
            cm = jet
        else:
            cm = ['lightblue','darkblue']

#----Set Output Type and Location, if any----
    if save is True:
        if (agg is not None) or (outfile == 'notebook'):
            # force notebook output if Datashader
            bplt.output_notebook()
        else:
            if outfile is None: outfile = x+'_vs_'+y+'.html'
            bplt.output_file(outfile)


#----Sample the Data----
    if sampling is None: sampling = len(inframe.index)
    if len(inframe.index) > sampling:
        df = inframe.sample(sampling)
    else:
        df = inframe

#----Set default source----
    if source is None:
        source = ColumnDataSource(df)

#----Check for log scale----
    x_rng = (df[x].min(),df[x].max())
    if xlog == 'auto':
        xlog = logaxis(x_rng[0],x_rng[1])

    y_rng = (df[y].min(),df[y].max())
    if ylog == 'auto':
        ylog = logaxis(y_rng[0],y_rng[1])

    x_axis_type='linear'
    y_axis_type='linear'
    y_range = None
    x_range = None
    if xlog is True:
        x_axis_type='log'
    if ylog is True:
        y_axis_type='log'

#----Set up Plot----
    if tools is None: 
        tools = "pan,wheel_zoom,box_zoom,reset,save,box_select"
    fig = bplt.Figure(tools=tools,width=width,height=height,
                      x_axis_type=x_axis_type,y_axis_type=y_axis_type,
                      y_range=y_rng,x_range=x_rng,webgl=True)
    fig.xaxis.axis_label=x
    fig.yaxis.axis_label=y

#----Plot x vs y----
#    print 'agg = ',agg
    if agg is None:
        fig.circle(x,y,source=source,size=size)
        
        if save is True: 
            bplt.show(fig)#,new='window',browser='firefox --no-remote')
            if outfile != 'notebook': bplt.curdoc().clear()

    else:        
        def image_callback(x_range,y_range,w,h):
            # re-draw canvas
#            cvs = ds.Canvas(plot_width=w,plot_height=h,x_range=x_range,
#                           y_range=y_range)
            cvs = ds.Canvas(x_range=x_range,y_range=y_range)
            # re-aggregate
            aggr = cvs.points(df,x,y,agg=fdict[agg](aggcol))
            img = tf.interpolate(aggr,cmap=cm,
                                 how='log')
            return tf.dynspread(img,shape='circle',max_px=size,
                                threshold=0.8)#,how='saturate')
        
        if returnfunc is False:
            interimg = InteractiveImage(fig,image_callback)
#        if save is True: bplt.show(fig)
        #print "To view interactive plot, call InteractiveImage(fig,image_callback)"

#----Return----
    if returnfunc is False:
        return fig,source #this version will not allow dynamic zooming
    else:
        return fig,source,image_callback
#----------------------------------------------------------
#def image_callback(x_range,y_range,w,h):
    # re-draw canvas
    #            cvs = ds.Canvas(plot_width=w,plot_height=h,x_range=x_range,
#                           y_range=y_range)
#    cvs = ds.Canvas(x_range=x_range,y_range=y_range)
    # re-aggregate
#    aggr = cvs.points(df,x,y,agg=fdict[agg](aggcol))
#    img = tf.interpolate(aggr,cmap=['lightblue','darkblue'],
#                         how='log')
#    return tf.dynspread(img,shape='circle',max_px=size,
#                        threshold=0.9)#,how='saturate')


#----------------------------------------------------------
def traceplots(dframe,sampling=1000.0,agg=None,aggcol=None,
               outfile='traceplots.html',save=True,size=None):
# deprecated - columns
    """
    Author: Kari A. Frank
    Date: October 26, 2015
    Purpose: plot interactive matrix of scatter plots from given dataframe

    Usage: traceplots(dframe,agg='dscounts,aggcol=None,sampling=1000)

    Input:

     dframe: input dataframe

     agg:    type of aggregration to perform before plotting, since 
             plotting with all possible rows/blobs is generally 
             prohibitively slow. options are:
             - 'none' (plot all rows)
             - 'sampling' (take random subsample of rows)

      columns: DEPRECATED list of dataframe column names to include in 
               the plots. default is to use all columns

      save: optionally specify to return the figure list, without 
            saving the plot file or plotting in a browser by setting
            save=False.

      agg (string): type of aggregation to perform before plotting, since 
         plotting all blobs individually is generally prohibitively 
         slow. options are
         - None = plot every blob as individual point (does not require
           Datashader, can output straight to html file)
         - Datashader options:
             - 'dscount' = count elements in each plot bin
             - 'dsany' = zero if bin is empty, 1 if not
             - 'dssum','dsmin','dsmax','dsmean','dsvar','dsstd' = 
               take the sum, mean, max, etc, of all elements in a 
               plot bin. must also provide the aggcol argument
               to specify which column should be used to calculate
               the min, max, etc.
             - 'dscount_cat' = count all elements in bin, grouped by 
               category. must also provide the aggcol argument
               to specify which column should be used to group by. 
               said column must have categorical data type
                 -- may not work yet

    aggcol (string): name of the column that should be passed to 
            Datashader agg function. required for several of the 
            agg options (see above)

    sampling: number of blobs to include in plot (default = 2000.0). If
            set to None will include all blobs.
          
    outfile: name of output html file. set to 'notebook' to plot to 
             currently open notebook.
    
    size: size of data points to plot. if agg=None, defaults to 3, else
          defaults to 5 (for dynamic zooming)

    Output:
     - plots matrix of scatter plots of all provided dataframe columns to 
       interactive (browser-based) plot, unless save=False and agg=None, 
       then just creates the figure object
     - returns the list of figure objects and the ColumnDataSource object

    Usage Notes:
     - must close and save (if desired) the plot manually
     - See important Datashader and Jupyter Notebook notes in scatter()
    Example:

    """

#----Import Modules----
    from bokeh.models import ColumnDataSource,PrintfTickFormatter
    from bokeh.layouts import column,row    

#----Trim dataframe----
#    if columns is not None:
#        dframe = dframe[columns]
# explicitly pass columns in the dataframe, e.g. 
#df[['blob_a','blob_c','blob_d']]

#----Set Output Type and Location, if any----
    if save is True:
        if (agg is not None) or outfile == 'notebook':
            # force notebook output if Datashader
            bplt.output_notebook()
        else:
            if outfile is None: outfile = 'traceplots.html'
            bplt.output_file(outfile)

#----Sample the Data----
    if sampling is None: sampling = len(dframe.index)
    if len(dframe.index) > sampling:
        df = dframe.sample(sampling)
    else:
        df = dframe

#----Set up plot----
    source = ColumnDataSource(df)
    if (agg != None) and (size is None):
        size = 5
    else:
        if size is None: size = 3
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select"
    w = 300 # plot width 
    h = 300 # plot height

#----Plot data---

#--initialize row--
    dim = len(df.columns)
    figarr = [] #columns (each element will be a list of figs=row)
    rowlist = []

#--loop through parameters and fill in row--
    parslist = list(df.columns)
    for par in parslist:
        for otherpar in parslist:
            # plot scatter plot
            if parslist.index(par) > parslist.index(otherpar):
                newfig,newsource=scatter(
                    df,otherpar,par,size=size,
                    agg=agg,save=False,sampling=sampling,
                    width=w,height=h,source=source,
                    tools=TOOLS,aggcol=aggcol)
                rowlist.append(newfig)
            # plot histogram
            if parslist.index(par) == parslist.index(otherpar):
                newfig = histogram(
                    df[par],bins=100,width=w,
                    height=h,tools=TOOLS,save=False)
                rowlist.append(newfig)
        figarr.append(row(rowlist[:]))
        rowlist = []

#--plot grid--
    bigfig = column(figarr[:])
    if agg is None:
        if save is True: 
#            bplt.show(column(figarr[:]))
            bplt.show(bigfig)
            bplt.curdoc().clear()
    else:
        if save is True: bplt.show(bigfig)

#----Return----
#    return figarr
    return bigfig

#----------------------------------------------------------
def histogram(dataseries,weights=None,bins=100,save=True,height=600,
              width=800,tools="pan,wheel_zoom,box_zoom,reset,save",
              infig=None,color='steelblue',outfile='histogram.html',
              density=False,alpha=None,xlog='auto',legend=None,
              norm=False,xmin=None,xmax=None,**kwargs):
    """
    Author: Kari A. Frank
    Date: October 28, 2015
    Purpose: plot a (weighted) histogram from provided series.

    Usage: theplot,bintype = histogram(dataseries,weights=None,bins=30,
                                      width=400,height=300,xlog='auto',
                                      tools=tools,save=True)

    Input:

     datacolumn:  a pandas series of data (e.g. column of a pandas 
                  dataframe)

     weights:     optionally provided a pandas series of weights which 
                  correspond to the values in datacolumn (e.g. emission 
                  measure)

     bins:        optionally specify the number of bins (default=30).
                  is passed directly to numpy.histogram(), so can
                  also accept any of the strings recognized by that function
                  for special calculation of bin sizes.
                  - 'auto' (but don't use if using weights)
                  - 'fd' (Freedman Diaconis Estimator)
                  - 'doane'
                  - 'scott'
                  - 'rice' 
                  - 'sturges'
                  - 'sqrt'

     height,width: height and width of each histogram

     save:        optionally turn off opening and saving the plot,
                  returns the figure object only (default=True)

     infig:       optionally pass an initialized figure object to plot on 
                  (allows plotting multiple dataseries on the same figure) 
                  (default=None)

     density:     passed to histogram. if True, then returned histogram is
                  is the probability density function. In general, will
                  not use this.

     norm:        normalize the histogram(s) so the y-axis goes from 0 to 1
                  (default=False). useful when 
                  plotting multiple histograms for comparison. 

     xlog:        make x-axis bins uniform in log10 space. options are:
                  - 'auto' (default) -- try to automatically determine
                    which is best based on the data range
                  - True -- force log bins
                  - False -- force linear bins

     xmin,xmax:   explicitly force xaxis min and max

     outfile:     string name of html file to save figure to. set to 
                  outfile='notebook' to plot to jupyter notebook.

     legend:      string to use for legend

     **kwargs:    pass any number of extra keyword arguments that are 
                  accepted by bokeh.plotting.quad().  some of the most
                  useful may be fill_color and line_color

    Output:
     - uses bokeh to open a plot of the (weighted) histogram in a browser
     - returns the figure object (allows replotting it, e.g. in a grid with
       other figures), along with string to specify if the xbins are spaced
       linearly ('lin') or logarithmically ('log')

    Usage Notes:
     - must close and save (if desired) the plot manually
     - axis labels will use the pandas series names (e.g. dataseries.name)
     - if outfile='notebook', then MUST call bplt.output_notebook() in its
       own cell in the current notebook before calling this function, else
       it won't automatically display the figure inline.

    Example:

    """

#----Import Modules----
    from bokeh.models import PrintfTickFormatter

#----Set up opacity----
    if alpha is None:
        if (infig is not None):
            alpha = 0.7
        else: 
            alpha = 1.0

#----Set up Log Bins and Axis----
    rng = (dataseries.min(),dataseries.max())
    
    # log axis
    if xmin is None:
        xmin = rng[0]
    if xmax is None:
        xmax = rng[1]
    xaxisrng = (xmin,xmax)

    if xlog == 'auto':
        xlog = logaxis(xaxisrng[0],xaxisrng[1])

    # set up log bins
    if xlog is True: 
        logbins = np.logspace(np.log10(rng[0]),np.log10(rng[1]),bins)
        bins = logbins
        bintype = 'log'
    else:
        bintype = 'linear'

#----Create the weighted histogram----
    histy,binedges = np.histogram(dataseries,weights=weights,bins=bins,
                                  density=density,range=rng)
#    print 'histy = ',histy
#    print 'binedges = ',binedges

    if norm is True:
#        histy = histy/float(np.sum(histy))
        histy = histy/float(np.max(histy))
        yaxisrng = (0.0,1.1)
    else:
        yaxisrng = (0.0,1.1*np.max(histy))

#----Set up Plot----
    if save:
        if (outfile != 'notebook'): 
            bplt.output_file(outfile)
    #    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
        else:
            bplt.output_notebook()

    if infig is None:
        fig = bplt.Figure(tools=tools,plot_width=width,plot_height=height,
                          x_axis_type=bintype,x_range=xaxisrng,
                          y_range=yaxisrng)
        fig.xaxis.axis_label=dataseries.name

        if weights is not None:
            if norm is True:
                fig.yaxis.axis_label='normalized '+weights.name
            else:
                fig.yaxis.axis_label=weights.name
            if np.log10(max(histy))>3: 
                fig.yaxis.formatter=PrintfTickFormatter(format = "%1.1e")
        else:
            if norm is True:
                fig.yaxis.axis_label = 'normalized number of blobs'
            else:
                fig.yaxis.axis_label = 'number of blobs'
    else:
        fig = infig


#----Format Ticks----
#    fig.yaxis.formatter=PrintfTickFormatter(format="%4.1e")
#    fig.xaxis.formatter=PrintfTickFormatter(format="%4.1e")
#    if infig is None:
#        fig.yaxis.formatter=format_ticks(histy)
#        fig.xaxis.formatter=format_ticks(binedges)

#---Plot the histogram----
    h = fig.quad(top=histy,bottom=0,left=binedges[:-1],right=binedges[1:],
                 color=color,alpha=alpha,legend=legend,**kwargs)
    
    if legend is not None:
        # add some fancier treatment of legend (e.g. move outside the plot)
#        fig.legend.location='top_left'
        fig.legend.orientation='horizontal'
#        leg = fig.legend
        pass

    if save:
        bplt.show(fig)#,new='window',browser='firefox --no-remote')
        if (outfile != 'notebook'): 
            bplt.curdoc().clear()

#----Return----
    return fig

#----------------------------------------------------------
def histogram_grid(dframes,weights=None,bins=100,height=300,width=400,
                   ncols=2,outfile='histogram_grid.html',
                   colors=['steelblue','darkolivegreen',
                  'mediumpurple','darkorange','firebrick','gray'],
                   alphas=None,norm=False,legends=None,**kwargs):
    """
    Author: Kari A. Frank
    Date: October 28, 2015
    Purpose: plot interactive matrix of weighted histograms from 
             given dataframe

    Usage: histogram_grid(dframe,weights=None,bins=30,**kwargs):

    Input:

     dframes:     input dataframe, or list of dataframes

     weights:     optionally provided a pandas series of weights which 
                  correspond to the values in datacolumn (e.g. emission 
                  measure). if dframes is a list of len>1, then 
                  weights must be one of the following:
                  - a list of the same length (with
                    each element a series corresponding one of the 
                    dataframes in dframes)
                  - a series that will be applied
                    as weights to all dataframes in dframes. all 
                    dataframes must have the same number of rows
                  - a string containing a column name to use as weights
                    every dataframe in dframes must contain this column
                  - a list of strings, containing the name of a column
                    in each dframes dataframe to use a weight
                  Note that if a list is provided, None is a valid element.
                  

     bins:        numerical or list of numerical, optional. the number 
                  of bins (default=100) used to create the histogram for 
                  each dataframe in dframes. if dframes is a list of 
                  len>1, then bins must be either a scalar (the same 
                  number of bins will be used for every dataframe), or a 
                  a list the same length as dframes.

     norm:        normalize the histogram(s) so the y-axis goes from 0 to 1
                  (default=False). useful when 
                  plotting multiple histograms for comparison. 

     height,width: height and width of each histogram, passed to histogram()


     alphas:      optionally pass float or list of floats (if 
                  len(dframes>1) specifying the opacity of each plotted
                  histogram (passed to histogram()). (default=0.4 if 
                  len(dframes)>1, else default=1)
                 
     colors:      optionally pass color or list of colors (if 
                  len(dframes>1) specifying the fill color of each plotted
                  histogram (passed to histogram()). default is to choose
                  in order from ['steelblue','darksage',
                  'mediumpurple','darkorange','firebrick','gray']

     legends:     optional string or list of strings (list if 
                  len(dframes)>1) to use as legend labels.

     ncols:       number of columns in the grid of histograms (default=4)

     outfile:     string name of output html plot file. if plotting to 
                  a jupyter notebook, use outfile='notebook'.

     **kwargs:    pass any number of extra keyword arguments that are 
                  accepted by bokeh.plotting.quad().  some of the most
                  useful may be fill_color and line_color

    Output:
     - plots matrix of scatter plots of all provided dataframe columns to 
       interactive (browser-based) plot

    Usage Notes:
     - must close and save (if desired) the plot manually
     - if outfile='notebook', then MUST call bplt.output_notebook() in its
       own cell in the current notebook before calling this function, else
       it won't automatically display the figure inline.
     - if dframes is a list of len>1, then only the columns specified in 
       first dataframe will be plotted as histograms. corresponding 
       columns in the other dataframes will be overplotted if they exist. 
       additional columns in the other dataframes will be ignored.
     - if dframes is a list, and any of the other arguments which can
       be passed as lists are of longer length than dframes, any
       extra elements in those lists will ignored.

    Example:
     - Plot histograms of the parameters in blobcols list, in grid with 
       1 column, weighted by the blob_weight column, with 300 bins.

         hfigs=xplt.histogram_grid(df[blobcols],bins=300,
                     ncols=1,height=200,width=400,weights=df.blob_weight)
       
     - Plot histograms of the columns labeled 'blob_kT' and 'blob_mass'
       with the histograms from df2 overplotted on those from df1
       
         hfigs=xplt.histogram_grid([df['blob_kT','blob_mass'],df2],
                       bins=[300,100],ncols=1,height=200,width=400,
                       weights=[df.blob_weight,df2.blob_weight])


    """

#----Import Modules----
    import math
    from bokeh.models import Label,Span

#----Check for multiple dataframes and cast variables as lists----
    if isinstance(dframes,list):
        
        if not isinstance(weights,list):
            if isinstance(weights,str):
                weights = [dfr[weights] for dfr in dframes]
            else:
                weights = [weights]*len(dframes)

        if not isinstance(bins,list):
            bins = [bins]*len(dframes)

        if not isinstance(alphas,list):
            if alphas is None:
                alphas = [0.4]*len(dframes)
            else:
                alphas = [alphas]

        if not isinstance(legends,list):
            legends = [legends]*len(dframes)

    else:
        dframes = [dframes]
        if isinstance(weights,str):
            weights = [dframes[0][weights]]
        else:
            weights = [weights]
        if alphas is None:
            alphas=[1.0]
        bins = [bins]
        legends = [legends]

#----Set up plot----
    if outfile != 'notebook':
        bplt.output_file(outfile)
    else:
        bplt.output_notebook()

#----Initialize empty figure list----
    figlist=[]

#----Create empty figure with legend only----

    # create empty plot
    newfig = bplt.figure(plot_width=width,plot_height=height)

    # create labels (one per dframe)
    for i in xrange(len(dframes)):
        yi = height-30-30*i
               
        # add color bar
        legbar = Span(location=yi,dimension='width',line_color=colors[i],
                      line_alpha=alphas[i],line_width=25,
                      location_units='screen')
        newfig.add_layout(legbar)

        # add text
        leg = Label(x=70,y=yi-10,x_units='screen',y_units='screen',
                    text=legends[i],
                    text_color='black',
                    text_font_style='bold',#background_fill_color=colors[i],
                    #background_fill_alpha=alphas[i]
                    )
        newfig.add_layout(leg)        

    figlist = figlist+[newfig]

#----Fill in list of figures----
    for column in dframes[0]:

        # get min and max of xaxis
        xmin = min([dseries.min() for dseries in 
                    [dfr[column] for dfr in dframes] ])
        xmax = max([dseries.max() for dseries in 
                    [dfr[column] for dfr in dframes] ])

        # get xaxis scale
        xlog = logaxis(xmin,xmax)

        # first histogram
        if isinstance(weights[0],str):
            weights[0] = dframes[0][weights[0]]
        newfig = histogram(dframes[0][column],weights=weights[0],
                           save=False,height=height,width=width,xmin=xmin,
                           xmax=xmax,xlog=xlog,#legend=legends[0],
                           bins=bins[0],norm=norm,color=colors[0],
                           alpha=alphas[0],**kwargs)

        # loop through any remaining dataframes
        if len(dframes)>1:
            for d in xrange(len(dframes)-1):
                # proceed only if column exists
                if column in dfr.columns:
                    if isinstance(weights[d+1],str):
                        weights[d+1] = dframes[d+1][weights[d+1]]
                    # plot histogram
                    newfig = histogram(dframes[d+1][column],
                                     weights=weights[d+1],bins=bins[d+1],
                                     save=False,color=colors[d+1],
                                     alpha=alphas[d+1],
                                     norm=norm,#legend=legends[d+1],
                                     infig=newfig,xlog=xlog)

        figlist=figlist+[newfig]
        #print column,bins

#----Reshape list into a 4 column array---

    #--define new shape--
#    nfigs = len(figlist)
#    nrows = int(math.ceil(float(nfigs)/float(ncols)))

    #--pad list with None to have nrows*ncols elements--
#    figlist = figlist+[None]*(nrows*ncols-nfigs)

    #--reshape list--
#    figarr = [figlist[ncols*i:ncols*(i+1)] for i in range(nrows)]

#----Plot histograms----
#    p = gridplot(figarr)
    p = gridplot(figlist,ncols=ncols,plot_width=width,plot_height=height)
    bplt.show(p)
    if outfile != 'notebook': bplt.curdoc().clear()

    return figlist

#----------------------------------------------------------
def format_ticks(vals):
    """
    Author: Kari A. Frank
    Date: October 28, 2015
    Purpose: Adjust tick label format to look better.

    Usage: <fig>.xaxis.formatter=format_ticks(xvals)

    Input:

     vals: values used for the x or y variable

    Output:

     Returns a PrintfTickFormatter object, which the user should then 
        use to set <fig>.xaxis.formatter (or yaxis)

    Usage Notes:

    Example:

    """
#----Import Modules----
    from bokeh.models import PrintfTickFormatter

#----Determine Most Reasonable Format and Return----
    rng = (np.min(vals),np.max(vals))
    norders = np.log10(rng[1]) - np.log10(rng[0])
#    minorder = np.log10(rng[0])
#    maxorder = np.log10(rng[1])
    check1 = abs(max(vals))
    check2 = abs(np.log10(check1))
#    print 'check1,check2 = ',check1,check2
    if norders > 2.0:
#    if ( (max(vals) > 0.01) or (max(vals) > 999.0) ):
        return PrintfTickFormatter(format = "%1.1e")
    if check1 < 1.0:
#    if maxorder <= 0.0: #max value < 1
        return PrintfTickFormatter(format = "%1.2f")
    if check1 < 10.0:
#    if maxorder <= 1.0: #max value < 10
        return PrintfTickFormatter(format = "%2.1f")
    if check1 < 100.0:
#    if maxorder <= 2.0: #max value < 100
        return PrintfTickFormatter(format = "%2.0f")
    else: #max value > 100#
#        return PrintfTickFormatter(format = "%3.0f"
        return PrintfTickFormatter(format = "%1.2e")

#----------------------------------------------------------
def spectrum(runpath='./',smin=0,smax=None,datacolor='black',save=True,
            modelcolor='steelblue',lastmodelcolor='firebrick',nbins=None,
            outfile='spectrum.html'):
    """
    Author: Kari A. Frank
    Date: November 18, 2015
    Purpose: Make a plot of the data and (average) model spectra compared

    Usage: spectrum(runpath='./',smin=0,smax=None,nbins=350.0,
                   datacolor='black',modelcolor='steelblue',
                   lastmodelcolor='firebrick',save=True,outfile='spectrum.html')

    Input:
     runpath (string) : the relative path to xmc run folder, which
                        contains the spectrum.* files

     smin/smax (int) : minimum and/or maximum spectrum file include in the
                       averaging of the model spectra. corresponds to the 
                       spectrum* file names, e.g. smax=3 will average the 
                       files spectrum_1.fits, spectrum_2.fits, and 
                       spectrum_3.fits. default is all available spectra.

    nbins (int) : number of bins in the histogram. 
                  defaults to a binsize of 0.015 keV.

    outfile (str) : name of output html file, or 'notebook' if plotting to 
                    an open Jupyter notebook

    save (bool) : if save is False, then will not display the final figure
                  or save it to a file

    Output:
     - plots the spectra for to an interactive plot (if save=True)
     - Returns the figure object

    Usage Notes:
     - must close and save (if desired) the plot manually

    Example:
    """

    #----Import Modules----
    import os
    import astropy.io.fits as fits
    from file_utilities import ls_to_list

    #----Set defaults----
    if smax==None:
        smax = len(ls_to_list(runpath,'spectrum*')) - 1

    # check for MPI file names (if xmc was run with mpi)
    if os.path.isfile(runpath+'/spectrum0_0.fits'):
        specname = 'spectrum0_'
    else:
        specname = 'spectrum_'

    #----Read in data spectrum----
    dataspecfile = runpath+'/'+specname+'0.fits'
    data_table = fits.getdata(dataspecfile,0)
    #field names are 'wave','xdsp',and 'ph'
    data_wave = data_table.field('wave')

    #----Read in model spectra----

    #--read in first model spectrum--
    foundmodel = False
    sm = smin
    #print smin,smax
    while (sm<=smax and foundmodel is False):
        modelspecfile = runpath+'/'+specname+str(sm)+'.fits'
        if os.path.isfile(modelspecfile):
            model_table = fits.getdata(modelspecfile,0)
            model_wave = model_table.field('wave')
            model_wave_avg = model_wave
            foundmodel = True
        else:
            print "Warning: "+modelspecfile+" not found.  Skipping to next spectrum."
        sm = sm+1
            
    if foundmodel is False:
        print "ERROR: no spectrum files found in range."

    #--loop over remaining model spectra--
    for s in range(sm,smax):
        modelspecfile = runpath+'/'+specname+str(s)+'.fits'
        if os.path.isfile(modelspecfile):
            model_table = fits.getdata(modelspecfile,0)
            model_wave = model_table.field('wave')
            model_wave_avg = np.hstack((model_wave_avg,model_wave))
        else:
            print "Warning: "+modelspecfile+" does not exist. Skipping."

    #----Convert to pandas Series----
    data_wave = pd.Series(data_wave,name='Energy (keV)')    
    model_wave_avg = pd.Series(model_wave_avg,name='Energy (keV)')
    model_wave = pd.Series(model_wave,name='Energy (keV)')

    #----Create Histograms----
    if nbins is None:
        nbins = np.ceil((np.max(data_wave.values)-np.min(data_wave.values))
                        /0.015)

    datay,datax = np.histogram(data_wave,bins=nbins,density=True)
    avgmodely,avgmodelx = np.histogram(model_wave_avg,bins=nbins,density=True)
    lastmodely,lastmodelx = np.histogram(model_wave,bins=nbins,density=True)

    #----Set up Plot----
    if save is True:
        if outfile != 'notebook':
            bplt.output_file(outfile)
        else:
            bplt.output_notebook()

    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    fig = bplt.Figure(tools=TOOLS)
    if np.max(datax) < 15.0:
        fig.xaxis.axis_label='Energy (keV)'
    else:
        fig.xaxis.axis_label='Wavelength (Angstroms)'
    fig.yaxis.axis_label=''

    #--format ticks--
    fig.yaxis.formatter=format_ticks(datay)
    fig.xaxis.formatter=format_ticks(datax)

    #----Plot Spectra----
    fig.line(datax,datay,color=datacolor,line_width=2,legend="Data")
    fig.line(avgmodelx,avgmodely,color=modelcolor,legend=
             "Model (Average over Iterations)")
    fig.line(lastmodelx,lastmodely,color=lastmodelcolor,legend=
             "Model (Last Iteration)")

    #----Show the Plot----
    if save is True:
        bplt.show(fig)
        if outfile != 'notebook': bplt.curdoc().clear()

#----Return----
    return fig

#----------------------------------------------------------
def evolution(inframe,iteration_type = 'median',itercol = 'iteration',
              weights=None,save=True,outfile='evolution_plots.html',
              ncols=4,height=300,width=300):
    """
    evolution()
 
   Author: Kari A. Frank
    Date: April 15, 2016
    Purpose: plot parameters vs iteration

   Usage: scatter(inframe,iteration_type='median')

   Input:
 
    inframe (DataFrame):  pandas DataFrame. all columns will be plotted
                          and one column must contain iterations.
             
    iteration_type (string): method used to combine parameter values within
         each iteration
         - 'none' -- plot all blobs separately (very slow!)
         - 'median' -- takes median from the iteration
         - 'average' -- takes average from the iteration
         - 'stdev' -- plots the standard deviation from each iteration
         - 'total' -- sums all blobs within each iteration

    weights (str): column name of weights to apply to each blob 
                   (e.g. emission measure)

    save:        optionally turn off opening and saving the plot as an 
                 html file - returns the figure object only (default=True)

    ncols: optionally specify number of columns in plot grid if plotting 
           more than one parameter. ignored if len(columns)=1

    width/height: specify the size of the individual plots

    outfile: name of output file, or 'notebook' if plotting to a 
             Jupyter notebook

   Output:
   - plots parameter vs iteration to an interactive plot
   - returns the figure object

  Usage Notes:
  - must close and save (if desired) the plot manually
  - may not work properly if trying to plot too many points

  Example:
 
  """

    #----Import Modules----
    from bokeh.models import ColumnDataSource
    from wrangle import weighted_median,weighted_std
    from functools import partial

    inframe = inframe.copy()

    #----Combine blobs in each iteration----

    if weights is None:
        #--group by iteration--
        if iteration_type == 'median':
            iterframe = inframe.groupby(itercol,as_index=False).median()
        elif iteration_type == 'average':
            iterframe = inframe.groupby(itercol,as_index=False).mean()
        elif iteration_type == 'stdev':
            iterframe = inframe.groupby(itercol,as_index=False).std()
        elif iteration_type == 'total':
            iterframe = inframe.groupby(itercol,as_index=False).sum()
        elif iteration_type is None:
            iterframe = inframe
        else:
            print ("Warning: Unrecognized iteration_type. Using"
                   " iteration_type='median'")
            iterframe = inframe.groupby(itercol,as_index=False).median()

            
    else: # weights
        # set up functions to use
        if iteration_type == 'median':
            def func(g):
                w = inframe.ix[g.index][weights]
                return weighted_median(g,weights=w)
        elif iteration_type == 'average':
            def func(g):
                w = inframe.ix[g.index][weights]
                return np.average(g,weights=w)
        elif iteration_type == 'stdev':
            def func(g):
                w = inframe.ix[g.index][weights]
                return weighted_std(g,weights=w)
        elif iteration_type == 'total':
            def func(g):
                w = inframe.ix[g.index][weights]
                return np.sum(g*w)/np.sum(w)
        elif iteration_type is None:
            print ("Warning: iteration_type = None but weights != None. "
                   "If using weights, must aggregrate iterations. Setting"
                   " iteration_type='median'.")
            def func(g):
                w = inframe.ix[g.index][weights]
                return weighted_median(g,weights=w)
        else:
            print ("Warning: Unrecognized iteration_type. Using"
                   " iteration_type='median'")
            def func(g):
                w = inframe.ix[g.index][weights]
                return weighted_median(g,weights=w)
            
        # group and aggregrate
        aggfunc = partial(func)
        inframe_grp = inframe.groupby(itercol,as_index=False)
        iterframe=inframe_grp.agg(aggfunc)#.reset_index()
        # remove extra column index
        # iterframe.columns = iterframe.columns.get_level_values(0)

    #----Set up plot----
    if (save is True):
        if (len(iterframe.columns)>1): 
            itersave=False
            if outfile != 'notebook':
                bplt.output_file(outfile)
            else:
                bplt.output_notebook()
        else:
            itersave=True

    source = ColumnDataSource(iterframe)
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"

    #----Initialize empty figure list----
    figlist=[]

    #----Fill in list of figures----
    for col in iterframe.columns:

        #----Plot parameter vs iteration----
        newfig,newsource = scatter(iterframe,itercol,col,agg=None,source=source,
                         tools=TOOLS,sampling=None,xlog=False,
                         save=itersave,height=height,width=width)
        figlist=figlist+[newfig]    

    #----Plot Grid of Figures----
    p = gridplot(figlist,ncols=ncols,plot_width=width,plot_height=height)

    if save is True:
        bplt.show(p)
        if outfile != 'notebook': bplt.curdoc().clear()

    return p
