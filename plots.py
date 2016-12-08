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
 spectrum
 trace
 plot_lines
"""

#-import common modules-
import pandas as pd
import numpy as np
import bokeh
import bokeh.plotting as bplt
import bokeh.charts as bchart
from bokeh.layouts import gridplot
from wrangle import filterblobs,make_histogram,normalize_histogram

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
def scatter_grid(dframe,sampling=1000.0,agg=None,aggcol=None,
               outfile='scatter_grid.html',save=True,size=None):
# deprecated - columns
    """
    Author: Kari A. Frank
    Date: October 26, 2015
    Purpose: plot interactive matrix of scatter plots from given dataframe

    Usage: scatter_grid(dframe,agg='dscounts,aggcol=None,sampling=1000)

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
            if outfile is None: outfile = 'scatter_grid.html'
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
              density=False,alpha=None,xlog='auto',logbins=None,legend=None,
              norm=False,xmin=None,xmax=None,iterations=None,**kwargs):
    """
    Author: Kari A. Frank
    Date: October 28, 2015
    Purpose: plot a (weighted) histogram from provided series.

    Usage: theplot,bintype = histogram(dataseries,weights=None,bins=30,
                                      width=400,height=300,xlog='auto',
                                      tools=tools,save=True)

    Input:

     dataseries:  a pandas series of data (e.g. column of a pandas 
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

     xlog:        boolean to plot x-axis on log scale. options are
                  - 'auto' (default) -- try to automatically determine
                    which is best based on the data range
                  - True -- force log scale
                  - False -- force linear scale

     logbins:     boolean make x-axis bins uniform in log10 space. default
                  is to be the same as determined by xlog. this argument
                  separated from xlog to allow creation of log-scale bins
                  but plotted on a linear scale, or vice versa.

     xmin,xmax:   explicitly force xaxis min and max

     outfile:     string name of html file to save figure to. set to 
                  outfile='notebook' to plot to jupyter notebook.

     legend:      string to use for legend

     iterations:  optionally provide a series of the iterations that 
                  correspond to each element of dataseries. if provided, 
                  will be used to calculate the errorars for each bin
                  which will be plotted on the histogram.
                
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
        if (infig is not None) or (iterations is not None):
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

    # bin type
    if logbins is None:
        logbins = xlog

    x_axis_type = 'linear'
    if xlog is True:
        x_axis_type = 'log'

    #----Create the weighted histogram and errorbars----
    histy,errors,binedges = make_histogram(dataseries,weights=weights,
                                           density=density,datarange=rng,
                                           bins=bins,logbins=logbins,
                                           iterations=iterations)

#----Normalize and set y-axis range----
    if norm is True:
        histy,errors = normalize_histogram(histy,yerrors=errors)
        # set y axis range
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
                          x_axis_type=x_axis_type,x_range=xaxisrng,
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

#----Plot the histogram----
    h = fig.quad(top=histy,bottom=0,left=binedges[:-1],right=binedges[1:],
                 color=color,alpha=alpha,legend=legend)#,**kwargs)
    
#----Plot Errorbars----
    if iterations is not None:
        xbinsizes = binedges[1:]-binedges[:-1]
        xbins = binedges[:-1]+xbinsizes/2.0
        errorbar(fig, xbins, histy, xerr=None, yerr=errors, color=color, 
                 point_kwargs={}, error_kwargs={})

#----Plot the legend----
    if legend is not None:
        # add some fancier treatment of legend (e.g. move outside the plot)
#        fig.legend.location='top_left'
#        fig.legend.orientation='horizontal'
#        leg = fig.legend
        pass

#----Show the plot----
    if save:
        bplt.show(fig)#,new='window',browser='firefox --no-remote')
        if (outfile != 'notebook'): 
            bplt.curdoc().clear()

#----Return----
    return fig

#----------------------------------------------------------
def histogram_grid(dframes,columns=None,weights=None,bins=100,
                   height=300,width=400,iterations=None,
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

     columns:    list of column names to plot. default=None will
                 plot all columns in the first dframe in list

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
                  
     iterations: optionally provide the iterations associated with each 
                 row. same options and format as weights. if provided,
                 will be used to calculate and plot the errorbars for 
                 each bin in each histogram.

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
     - only the columns specified which exist in at least one of the given 
       dataframes will be plotted as histograms. 
       corresponding columns in other dataframes will be overplotted 
       if they exist. additional columns in any of the dataframes will be 
       ignored.
     - if dframes is a list, and any of the associated arguments which can
       be passed as lists are of longer length than dframes, any
       extra elements in those lists will ignored.

    Example:
     - Plot histograms of the parameters in blobcols list, in grid with 
       1 column, weighted by the blob_weight column, with 300 bins.

         hfigs=xplt.histogram_grid(df,columns=[blobcols],bins=300,
                     ncols=1,height=200,width=400,weights=df.blob_weight)
       
     - Plot histograms of the columns labeled 'blob_kT' and 'blob_mass'
       with the histograms from df2 overplotted on those from df1
       
         hfigs=xplt.histogram_grid([df,df2],columns=['blob_kT','blob_mass'],
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

        if not isinstance(iterations,list):
            if isinstance(iterations,str):
                iterations = [dfr[iterations] for dfr in dframes]
            else:
                iterations = [iterations]*len(dframes)

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
        if isinstance(iterations,str):
            iterations = [dframes[0][iterations]]
        else:
            iterations = [iterations]
        if alphas is None:
            alphas=[1.0]
        bins = [bins]
        legends = [legends]
        if (not isinstance(columns,list)) and (columns is not None):
            columns = [columns]

    if columns is None:
        columns = []
        for df in dframes:
            columns = columns+[c for c in df.columns if c not in columns]

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

    for column in columns:

        # get min and max of xaxis
        xmin = dframes[0][column].min()
        xmax = dframes[0][column].max()
        for dfr in dframes:
            if column in dfr.columns:
                xmin = min(xmin,dfr[column].min())
                xmax = max(xmax,dfr[column].max())

        # get xaxis scale
        xlog = logaxis(xmin,xmax)

        # loop through any remaining dataframes
        newfig = None # reset newfig for each column
        for d in xrange(len(dframes)):
            # proceed only if column exists
            if column in dframes[d].columns:
                if isinstance(weights[d],str):
                    weights[d] = dframes[d][weights[d]]

            # plot histogram
                newfig = histogram(dframes[d][column],
                                     weights=weights[d],bins=bins[d],
                                     save=False,color=colors[d],
                                     alpha=alphas[d],height=height,
                                     width=width,xmin=xmin,xmax=xmax,
                                     norm=norm,#legend=legends[d+1],
                                     infig=newfig,xlog=xlog,
                                   iterations=iterations[d])
        figlist=figlist+[newfig]

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
             modelcolor='steelblue',lastmodelcolor='firebrick',bins=0.03,
             outfile='spectrum.html',ylog=False,xlog=False,logbins=None,
             datarange=None,width=1000,height=500,lines=True,**lineargs):

    """
    Author: Kari A. Frank
    Date: November 18, 2015
    Purpose: Make a plot of the data and (average) model spectra compared

    Usage: spectrum(runpath='./',smin=0,smax=None,nbins=350.0,
                   datacolor='black',modelcolor='steelblue',
                   lastmodelcolor='firebrick',save=True,
                   outfile='spectrum.html')

    Input:
     runpath (string) : the relative path to xmc run folder, which
                        contains the spectrum.* files

     smin/smax (int) : minimum and/or maximum spectrum file include in the
                       averaging of the model spectra. corresponds to the 
                       spectrum* file names, e.g. smax=3 will average the 
                       files spectrum_1.fits, spectrum_2.fits, and 
                       spectrum_3.fits. default is all available spectra.

     bins:        optionally specify the number of bins (int), binsize 
                  (float), or an array containing the bin edge values. 
                  can also accept any of the strings recognized by 
                  np.histogram() for special calculation of bin sizes.
                  - 'auto' (but don't use if using weights)
                  - 'fd' (Freedman Diaconis Estimator)
                  - 'doane'
                  - 'scott'
                  - 'rice' 
                  - 'sturges'
                  - 'sqrt'
                  default is to create equally spaced bins of width 0.03,
                  (30 eV). note that xmc typically uses a binsize of 0.015.

    outfile (str) : name of output html file, or 'notebook' if plotting to 
                    an open Jupyter notebook

    save (bool) : if save is False, then will not display the final figure
                  or save it to a file

    lines (bool) : plot the 10 stongest common emission lines within 
                   the x-axis range. **lineargs will pass any extra 
                   arguments directly to plot_lines()


    logbins, datarange, xlog, ylog : see make_histogram()

    Output:
     - plots the spectra for to an interactive plot (if save=True)
     - Returns the figure object

    Usage Notes:
     - will automatically plot error bars on the average spectrum
       if more than 3 model spectra are used
     - must close and save (if desired) the plot manually

    Example:
    """

    #----Import Modules----
    import os
    import astropy.io.fits as fits
    from bokeh.charts import Step
    from bokeh.models.ranges import Range1d
    from file_utilities import ls_to_list

    #----Set defaults----
    if smax==None:
        smax = len(ls_to_list(runpath,'spectrum*')) - 1

    if logbins is None: logbins=xlog

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
            iters_avg = np.ones_like(model_wave)
            iters_avg.fill(sm)
            foundmodel = True
        else:
            print ("Warning: "+modelspecfile+" not found.  Skipping +"
                   "to next spectrum.")
        sm = sm+1
            
    if foundmodel is False:
        print "ERROR: no spectrum files found in range."

    #--loop over remaining model spectra--
    for s in xrange(sm,smax+1):
        modelspecfile = runpath+'/'+specname+str(s)+'.fits'
        if os.path.isfile(modelspecfile):
            model_table = fits.getdata(modelspecfile,0)
            model_wave = model_table.field('wave')
            model_iters = np.ones_like(model_wave)
            model_iters.fill(s) #assign 'iteration' number
            model_wave_avg = np.hstack((model_wave_avg,model_wave))
            iters_avg = np.hstack((iters_avg,model_iters))
        else:
            print "Warning: "+modelspecfile+" does not exist. Skipping."

    #----Convert to pandas Series----
    data_wave = pd.Series(data_wave,name='Energy (keV)')    
    model_wave_avg = pd.Series(model_wave_avg,name='Energy (keV)')
    model_wave = pd.Series(model_wave,name='Energy (keV)')
    iters_avg = pd.Series(iters_avg,name='iteration')

    #----Create Histograms----
    
    if isinstance(bins,float): # assume binsize, convert to number bins
        bins = np.ceil((np.max(data_wave.values)-
                        np.min(data_wave.values))/bins)    

    datay,dataerrors,dataedges = \
        make_histogram(data_wave,bins=bins,
                       logbins=logbins,
                       datarange=datarange,
                       density=False,iterations=None)
    datay,dataerrors = normalize_histogram(datay,dataerrors)

    avgmodely,avgmodelerrors,avgmodeledges=\
        make_histogram(model_wave_avg,
                       bins=bins,
                       logbins=logbins,datarange=datarange,
                       density=False,iterations=iters_avg)
    avgmodely,avgmodelerrors = normalize_histogram(avgmodely,
                                                   avgmodelerrors)

    lastmodely,lastmodelerrors,lastmodeledges=\
        make_histogram(model_wave,
                       bins=bins,
                       logbins=logbins,datarange=datarange,
                       density=False,iterations=None)
    lastmodely,lastmodelerrors = normalize_histogram(lastmodely,
                                                   lastmodelerrors)

    #----Set up Plot----
    if save is True:
        if outfile != 'notebook':
            bplt.output_file(outfile)
        else:
            bplt.output_notebook()

    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"#,box_select,lasso_select"
    
    x_axis_type = 'linear'
    y_axis_type = 'linear'
    if xlog is True:
        x_axis_type = 'log'
    x_range = (np.min(dataedges),np.max(dataedges))
    if ylog is True:
        y_axis_type = 'log'
        ymin = np.min(datay[np.where(datay != 0.0)])
        y_range = (ymin,1.5*np.max(datay))
    else: 
        y_range = (0.0,1.1*np.max(datay))
    if np.max(dataedges) < 15.0:
        xlabel = 'Energy (keV)'
    else:
        xlabel = 'Wavelength (Angstroms)'

    avglabel = 'Model (average)'
    lastlabel = 'Model (last iteration)'
    specframes = pd.DataFrame({'Data':datay,xlabel:dataedges[:-1],
                               avglabel:avgmodely,lastlabel:lastmodely})

    #----Plot Spectra as Step Chart----
    step = Step(specframes,x=xlabel,y=['Data',avglabel,lastlabel],
                color=[datacolor,modelcolor,lastmodelcolor],legend=True,
                y_mapper_type=y_axis_type,x_mapper_type=x_axis_type,
                dash=['solid','solid','dashed'],
                plot_width=width,plot_height=height)
    step.x_range=Range1d(*x_range)
    step.y_range=Range1d(*y_range)
    step.legend.location='top_right'
    step.ylabel = 'normalized value'

    #----Plot Errorbars----
    xbinsizes = dataedges[1:]-dataedges[:-1]
    xbins = dataedges[:-1]+xbinsizes/2.0
    errorbar(step, xbins, avgmodely, xerr=None, yerr=avgmodelerrors, 
             color=modelcolor, 
             point_kwargs={}, error_kwargs={})

    #----Plot Emission Lines----
    if lines is True:
        step = plot_lines(step,dataedges,**lineargs)

    #possible attributes to Chart are above, background_fill_alpha, background_fill_color, below, border_fill_alpha, border_fill_color, disabled, extra_x_ranges, extra_y_ranges, h_symmetry, height, hidpi, left, lod_factor, lod_interval, lod_threshold, lod_timeout, min_border, min_border_bottom, min_border_left, min_border_right, min_border_top, name, outline_line_alpha, outline_line_cap, outline_line_color, outline_line_dash, outline_line_dash_offset, outline_line_join, outline_line_width, plot_height, plot_width, renderers, right, sizing_mode, tags, title, title_location, tool_events, toolbar, toolbar_location, toolbar_sticky, v_symmetry, webgl, width, x_mapper_type, x_range, xlabel, xscale, y_mapper_type, y_range, ylabel or yscale


    #----Show the Plot----
    if save is True:
        bplt.show(step)
#        if outfile != 'notebook': bplt.curdoc().clear()

#----Return----
    return step

#----------------------------------------------------------
def trace(inframe,iteration_type = 'median',itercol = 'iteration',
          weights=None,itmin=0,itmax=None,
          save=True,outfile='trace_plots.html',
          ncols=4,height=300,width=300):
    """
    Author: Kari A. Frank
    Date: April 15, 2016
    Purpose: plot parameters vs iteration

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

    itmin/itmax: minimum/maximum iteration to plot

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

    #----Copy frame and apply itmin/itmax----
    inframe = inframe.copy()
    if itmax is None:
        itmax = np.max(inframe[itercol])
    inframe = filterblobs(inframe,itercol,minvals=itmin,maxvals=itmax)

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
        newfig,newsource = scatter(iterframe,itercol,col,agg=None,
                                   source=source,
                         tools=TOOLS,sampling=None,xlog=False,
                         save=itersave,height=height,width=width)
        figlist=figlist+[newfig]    

    #----Plot Grid of Figures----
    p = gridplot(figlist,ncols=ncols,plot_width=width,plot_height=height)

    if save is True:
        bplt.show(p)
        if outfile != 'notebook': bplt.curdoc().clear()

    return p

#----------------------------------------------------------
def errorbar(fig, x, y, xerr=None, yerr=None, color='steelblue', 
             point_kwargs={}, error_kwargs={}):
    """Function to plot symmetric errorbars on top of a figure"""

    """
    From: http://stackoverflow.com/questions/29166353/how-do-you-add-error-bars-to-bokeh-plots-in-python

    
    fig can be either a Bokeh Figure object or bokeh Chart object
    """

    fig.circle(x, y, color=None,fill_alpha=0, **point_kwargs)
    
    if xerr is not None:
        x_err_x = []
        x_err_y = []
        for px, py, err in zip(x, y, xerr):
            x_err_x.append((px - err, px + err))
            x_err_y.append((py, py))
        fig.multi_line(x_err_x, x_err_y, color=color, **error_kwargs)

    if yerr is not None:
        y_err_x = []
        y_err_y = []
        for px, py, err in zip(x, y, yerr):
            y_err_x.append((px, px))
            y_err_y.append((py - err, py + err))
        fig.multi_line(y_err_x, y_err_y, color=color, **error_kwargs)

#----------------------------------------------------------
def plot_lines(fig,bins,kT_range=(0.17,10.0),emissivity_range=(1e-17,1.0),
               energy_range=(0.1,10.0),wavelength_range=None,nlines=50
               ,include_lines=None):
    """
    Plot strong emission lines on top of a spectrum    

    Author: Kari A. Frank
    Date: December 1, 2016
    Purpose: Plot strong emission lines on top of a spectrum figure.

    Input:

     fig (bokeh Chart or bokeh Figure) : chart or figure object 
                       containing the spectrum. must already 
                       have spectrum added to the figure.

     bins (numeric array) : 1d array containing the bin edge values of
                       the spectrum, as output by np.histogram() or 
                       make_histogram() in spectrum()

     kT_range (2d tuple) : range of gas temperatures to include emission
                       lines from

     emissivity_range (2d tuple) : range of line emissivities to include

     wavelength_range (2d tuple) : range of line wavelengths (angstroms) 
                         to include emission lines from

     energy_range (2d tuple) : range of line energies (keV) to include
                         emission lines from
     
     nlines (int) : maximum number of lines to include. if more than
                    nlines (combined) lines meet the other criteria, 
                    then will plot only the nlines with highest emissivity

     include_lines (list of strings) : list of element names to include 
                    emission lines from. will first drop all lines from
                    other elements, then apply remaining criteria (e.g.
                    kT_range).

    Output:

     Add vertical lines (at most one per spectral bin) to the provided
     figure. If more than one line per bin, then a single line will be
     plotted with a label that lists all lines in that bin (grouped by
     ion and ionization stage).


    """

    import os
    from bokeh.models.annotations import Span,Label
    from bokeh.models import HoverTool

    xmin = float(fig.x_range.start)
    xmax = float(fig.x_range.end)
    ymin = float(fig.y_range.start)
    ymax = float(fig.y_range.end)

    #----Read in lines----
    linefile = os.path.dirname(__file__)+'/xraylines_atomdb302.txt'
    linedf = pd.read_table(linefile,sep='\s+',comment='#',engine='python')

    #----Remove extra lines----

    #--remove lines not in include_lines--
    if include_lines is not None:
        linedf = linedf[linedf['ion'].isin(include_lines)]

    #--get only lines from gas with kT in range--
    if kT_range is not None:
        linedf = linedf[kT_range[0] <= linedf['kT']]
        linedf = linedf[linedf['kT'] <= kT_range[1]]

    #--get only lines with emissivity in range--
    if emissivity_range is not None:
        linedf = linedf[emissivity_range[0] <= linedf['emissivity']]
        linedf = linedf[linedf['emissivity'] <= emissivity_range[1]]
    
    #--get only lines with energy in range--
    if energy_range is not None:
        linedf = linedf[energy_range[0] <= linedf['energy']]
        linedf = linedf[linedf['energy'] <= energy_range[1]]

    #--get only lines with wavelength in range--
    if wavelength_range is not None:
        linedf = linedf[wavelength_range[0] <= linedf['wavelength']]
        linedf = linedf[linedf['wavelength'] <= wavelength_range[1]]

    #--get only lines that are within the plot range--

    #-set x units-
    if xmax <= 15.0: # epic data (keV scale)
        x = 'energy'
    else: # rgs data (angstrom scale)
        x = 'wavelength'

    plotdf = agg_lines(linedf.groupby(['ion','ionization_stage']),
                       bins,units=x)

    #-cut out-of-range lines-
#    plotdf = plotdf[plotdf[x] >= xmin]
#    plotdf = plotdf[plotdf[x] <= xmax]

    #--truncate to include no more than nlines lines, keeping those with
    # highest emissivity--
    if (nlines is not None) and (len(plotdf.index)>nlines):
        plotdf = plotdf.nlargest(nlines,'emissivity')
    else:
        nlines = len(plotdf.index)

    #----Print line information----
    print plotdf.head(nlines)

    #----Plot lines----
    
    #--calculate label offset--
    xoffset = (xmax - xmin)/300.0
    labely = (ymax - ymin)/100.0 + ymin

    #--plot lines--
    for i,row in plotdf.iterrows():
#        print row[x],row['label'],row['emissivity']
        lspan = Span(location=row[x],dimension='height',
                     line_color='black',
                     line_dash='dashed',line_width=1)
        fig.add_layout(lspan)
        llabel = Label(x=row[x],y=labely,
                       angle=90.0,angle_units='deg',
                       text_font_size='10px',x_offset=xoffset,
                       text=row['label'])
        fig.add_layout(llabel)

    return fig

#----------------------------------------------------------
def agg_lines(groupeddf,bins,units='energy'):
    """
    Merge emission lines from the same ions/ionization stage if they
    are close enough (where 'close enough'=tolerance). Must provide
    the dataframe from xraylines.txt grouped by both ion and ionization 
    stage columns.

    Returns a 2-column dataframe containing the line energies/wavelengths
    and labels for plotting. All filtering of the lines (for temperature,
    axis range, etc.) must be done prior to calling this function.

    Called by spectrum()
    """

    # empty dataframe for output
    outdf=pd.DataFrame(data=np.zeros((0,3)),columns=[units,'label',
                                                     'emissivity'])
    iondf=pd.DataFrame(data=np.zeros((0,3)),columns=[units,'label',
                                                     'emissivity'])

    # merge lines by element and ionization stage
    for ion,grp in groupeddf:
        # histogram with same bins as spectrum
        y,yerrs,edges,x=make_histogram(grp[units],bins=bins,
                   centers=True,weights=grp.emissivity)
        # combine lines in each (non-empty) bin
        for b in xrange(len(y)):
            if y[b] > 0:
                lab = ion[0]+' '+ion[1]
                iondf=iondf.append(pd.DataFrame([[x[b],lab,y[b]]],
                                          columns=[units,'label',
                                                   'emissivity']))
                
    # check for lines from different ions or ionization stages in 
    # the same energy bin and combine
    y,yerrs,edges,x = make_histogram(iondf[units],bins=bins,
                                     centers=True,
                                     weights=iondf['emissivity'])
    for b in xrange(len(y)):
        if y[b] > 0:
            # find labels of all in bin
            subdf = iondf[iondf[units]==x[b]]
            lab = ''
            for i,l in subdf.iterrows():
                if lab != '':
                    lab = lab+', '+l.label
                else:
                    lab = l.label
            outdf=outdf.append(pd.DataFrame([[x[b],lab,y[b]]],
                                            columns=[units,'label',
                                                     'emissivity']))

    outdf.reset_index(inplace=True,drop=True)
    return outdf
