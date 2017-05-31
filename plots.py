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
 blobs_spectrum
"""

#-import common modules-
import pandas as pd
import numpy as np
import bokeh
import bokeh.plotting as bplt
import bokeh.charts as bchart
from bokeh.layouts import gridplot
from wrangle import filterblobs,make_histogram,normalize_histogram
from xmcfiles import fake_deconvolution,merge_output

#----------------------------------------------------------
def logaxis(minval,maxval,limit=2.0):
    """
    Given minimum and maximum values of axis range, decide if log or linear
    scale should be used.
    
    limit specifies how many orders of magnitude must separate the 
    min and max values before switching to log scale.

    Only allows setting log axis if minval is positive.

    Returns either True (log axis) or False (linear axis).
    """

    # only allow log scale if all positive values
    if (minval > 0.0):
        norders = np.log10(maxval) - np.log10(minval)

        if norders > float(limit):
            return True
        else:
            return False
    else:
        return False
        
#----------------------------------------------------------
def chi2(runpath='./',itmin=0,itmax=None,outfile='chi2_vs_iteration.html',
         display=True):
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

     display (bool) : if False, then will not display the figure. 
                      automatically (re)set to True if outfile='notebook'.j
                      ignored if save=False

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
    statframe = merge_output(runpath,filetype='statistic',save=False,
                             itmin=itmin,itmax=itmax)
    if itmax is None:
        itmax = np.max(statframe['iteration'])
    statframe = filterblobs(statframe,'iteration',minvals=itmin,maxvals=itmax)

#----Calculate reduced chi2----
    statframe['redchi2'] = statframe['chi2']/statframe['dof']

#----Set up Plot----
    if outfile != 'notebook':    
        bplt.output_file(outfile)
    else:
        display = True
        bplt.output_notebook()
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    fig = bplt.Figure(tools=TOOLS)
    fig.xaxis.axis_label='iteration'
    fig.yaxis.axis_label='chi2/dof'

#----Plot chi2 vs iteration----
#    statframe.plot(x='iteration',y='redchi2',kind='scatter')
#    plt.show()

    fig.circle(x=statframe['iteration'],y=statframe['redchi2'])
    if display is True:
        bplt.show(fig)#,new='window',browser='firefox --no-remote')
    else: 
        bplt.save(fig)
    if outfile != 'notebook': bplt.curdoc().clear()

#----Return----
    return statframe

#----------------------------------------------------------
def scatter(inframe,x,y,sampling=2000.,agg=None,aggcol=None,save=True,
            display=True,width=600,height=600,source=None,tools=None,size=5,
            xlog='auto',ylog='auto',outfile=None,returnfunc=False,
            span=None,cscale='eq_hist'):
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

    save:   if True, save the plot as an html file (if not using Datashader)
            Returns the figure object only if save=False (default=True)

    display (bool) : if False, then will not display the figure. 
                     automatically (re)set to True if outfile='notebook'
                     ignored if save=False

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
                    from datashader.bokeh_ext import InteractiveImage
                    InteractiveImage(fig,image_callback) 
                 on the python
                 command line. Only enable this if calling 
                 scatter() directly.

      span : [zmin,zmax], to specify the range to use for the colorscale
             only used datashader and aggcol is used.
             
      cscale : string specify type of scale to use for colorscale,
               'log', 'linear', or 'eq_hist' (default). only used if
               datashader is used.


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
            display=True
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
            if display is True:
                bplt.show(fig)#,new='window',browser='firefox --no-remote')
            else:
                bplt.save(fig)
            if outfile != 'notebook': bplt.curdoc().clear()

    else:        
        if (span is None) and (aggcol is not None):
            span = [df[aggcol].min(),df[aggcol].max()]
        else:
            span = None
        def image_callback(x_range,y_range,w,h):
            # re-draw canvas
#            cvs = ds.Canvas(plot_width=w,plot_height=h,x_range=x_range,
#                           y_range=y_range)
            cvs = ds.Canvas(x_range=x_range,y_range=y_range)
            # re-aggregate
            aggr = cvs.points(df,x,y,agg=fdict[agg](aggcol))
            img = tf.interpolate(aggr,cmap=cm,span=span,
                                 how=cscale)
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
                 outfile='scatter_grid.html',save=True,display=True,
                 size=None):
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
            save=False

      display (bool) : if False, then will not display the figure
               automatically (re)set to True if outfile='notebook'
               ignored if save=False

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
            display=True
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
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
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
            if display is True:
                bplt.show(bigfig)
            else:
                bplt.save(bigfig)
            bplt.curdoc().clear()
    else:
        if save is True: 
            if display is True:
                bplt.show(bigfig)
            else:
                bplt.save(bigfig)

#----Return----
#    return figarr
    return bigfig

#----------------------------------------------------------
def histogram(dataseries,weights=None,bins=100,save=True,display=True,
              height=600,
              width=800,tools="pan,wheel_zoom,box_zoom,reset,save",
              infig=None,color='steelblue',outfile='histogram.html',
              density=False,alpha=None,xlog='auto',logbins=None,legend=None,
              norm=False,xmin=None,xmax=None,iterations=None,
              ytitle=None,**kwargs):
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

     save:        optionally turn off saving the plot,
                  returns the figure object only (default=True)

     display (bool) : if False, then will not display the figure
                  automatically (re)set to True if outfile='notebook'
                   ignored if save=False

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

     ytitle:      string specifiying label of y-axis. default is the 
                  column name of the provided data series.

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
            display = True
            bplt.output_notebook()

    if infig is None:
        fig = bplt.Figure(tools=tools,plot_width=width,plot_height=height,
                          x_axis_type=x_axis_type,x_range=xaxisrng,
                          y_range=yaxisrng)
        fig.xaxis.axis_label=dataseries.name
            
        if weights is not None:
            if ytitle is not None:
                fig.yaxis.axis_label=ytitle
            else:
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
                 color=color,alpha=alpha,legend=legend,**kwargs)
    
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
        if display is True:
            bplt.show(fig)#,new='window',browser='firefox --no-remote')
        else:
            bplt.save(fig)
        if (outfile != 'notebook'): 
            bplt.curdoc().clear()

#----Return----
    return fig

#----------------------------------------------------------
def histogram_grid(dframes,columns=None,weights=None,bins=100,
                   height=300,width=400,iterations=None,display=True,
                   ncols=2,outfile='histogram_grid.html',
                   colors=['steelblue','darkolivegreen',
                  'mediumpurple','darkorange','firebrick','gray'],
                   alphas=None,norm=False,legends=None,**histargs):
    """
    Create html grid of (weighted) histograms from a dataframe.

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
                  in order from ['steelblue','darkolivegreen',
                  'mediumpurple','darkorange','firebrick','gray']

     legends:     optional string or list of strings (list if 
                  len(dframes)>1) to use as legend labels.

     ncols:       number of columns in the grid of histograms (default=4)

     outfile:     string name of output html plot file. if plotting to 
                  a jupyter notebook, use outfile='notebook'.

     display (bool) : if False, then will not display the figure
                  automatically (re)set to True if outfile='notebook'
                  ignored if save=False

     **histargs:  pass any number of extra keyword arguments that are 
                  accepted by histogram()

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
        display = True

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
                                   iterations=iterations[d],**histargs)
        figlist=figlist+[newfig]

#----Plot histograms----
#    p = gridplot(figarr)
    p = gridplot(figlist,ncols=ncols,plot_width=width,plot_height=height)
    if display is True:
        bplt.show(p)
    else:
        bplt.save(p)
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
def spectrum(runpath='../',smin=0,smax=None,datacolor='black',
             datalabel='Data',save=True,scale = 1.0,infig=None,
             model=True,display=True,modellabel='Model (average)',
             modelcolor='steelblue',lastmodelcolor='firebrick',
             lastmodellabel='Model (last iteration)',bins=0.03,
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
                        contains the spectrum.* files. alternatively,
                        can pass a dataframe with the x and y values, e.g.
                        as would be read from a text file written by
                        xspec with the wd command.

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

    scale (float or list of floats) : scale the each spectrum 
                    histogram by a constant. if a single float is provided,
                    all spectra will be scaled by the same value. if a list
                    is provided, it must have one element per spectrum to
                    to be plotted (1, 2, or 3), and must be in the 
                    following order - [data,(avg)model,lastmodel]

    outfile (str) : name of output html file, or 'notebook' if plotting to 
                    an open Jupyter notebook

    save (bool) : if save is False, then will not display the final figure
                  or save it to a file

    display (bool) : if False, then will not display the figure
                  automatically (re)set to True if outfile='notebook'
                  ignored if save=False                  

    infig (bokeh Figure) : if provided, then the spectrum will be
                   overplotted on the provided figure.

    model (bool) : if False, then will plot only the data spectrum.

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

    if isinstance(scale,float):
        scale = [scale,scale,scale]
        
    #----Read in data spectrum----
    dataspecfile = runpath+'/'+specname+'0.fits'
    data_table = fits.getdata(dataspecfile,0)
    #field names are 'wave','xdsp',and 'ph'
    data_wave = data_table.field('wave')

    #----Read in model spectra----
    
    if model is True:
    
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
#            print modelspecfile
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
    data_wave = pd.Series(data_wave.byteswap().newbyteorder(),name='Energy (keV)')    
    if model is True:
        model_wave_avg = pd.Series(model_wave_avg+0,name='Energy (keV)')
        model_wave = pd.Series(model_wave+0,name='Energy (keV)')
        iters_avg = pd.Series(iters_avg+0,name='iteration')

    #----Create Histograms----
    # y units are counts per bin (per observation exposure or iteration)
    
    if isinstance(bins,float): # assume binsize, convert to number bins
        bins = np.ceil((np.max(data_wave.values)-
                        np.min(data_wave.values))/bins)    

    datay,dataerrors,dataedges = \
        make_histogram(data_wave,bins=bins,
                       logbins=logbins,
                       datarange=datarange,
                       density=False,iterations=None)
    datay = datay*scale[0]
    if dataerrors is not None:
        dataerrors = dataerrors*scale[0]
#    datay,dataerrors = normalize_histogram(datay,dataerrors)

    if model is True:
        avgmodely,avgmodelerrors,avgmodeledges=\
            make_histogram(model_wave_avg,
                           bins=bins,
                           logbins=logbins,datarange=datarange,
                           density=False,iterations=iters_avg)
        nspecs = len(np.unique(iters_avg))
        print 'nspecs = ',nspecs
        # scale by number of iterations
        avgmodely = avgmodely*scale[1]/float(nspecs)
        avgmodelerrors = avgmodelerrors*scale[1]/float(nspecs)
#        avgmodely,avgmodelerrors = normalize_histogram(avgmodely,
#                                                       avgmodelerrors)


        if smin != smax: # skip if only plotting one iteration
            lastmodely,lastmodelerrors,lastmodeledges=\
            make_histogram(model_wave,
                           bins=bins,
                           logbins=logbins,datarange=datarange,
                           density=False,iterations=None)
            lastmodely=lastmodely*scale[2]
#            lastmodelerrors = lastmodelerrors*scale[2]
#            lastmodely,lastmodelerrors = normalize_histogram(lastmodely,
#                                                       lastmodelerrors)
    print 'max data, model = ',max(datay),max(avgmodely)

    #----Set up Plot----
    if save is True:
        if outfile != 'notebook':
            bplt.output_file(outfile)
        else:
            display = True
            bplt.output_notebook()

    if infig is None:
        TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    
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

    avglabel = modellabel #'Model (average)'
    lastlabel = lastmodellabel #'Model (last iteration)'

    #----Plot Spectra as Step Chart----
    if model is True:
        if smin != smax: # plotting more than one iteration
            specframes = pd.DataFrame({'Data':datay,xlabel:dataedges[:-1],
                               avglabel:avgmodely,lastlabel:lastmodely})

            step = Step(specframes,x=xlabel,y=[datalabel,avglabel,lastlabel],
                color=[datacolor,modelcolor,lastmodelcolor],legend=True,
                y_mapper_type=y_axis_type,x_mapper_type=x_axis_type,
                dash=['solid','solid','dashed'],
                plot_width=width,plot_height=height)
        else:
            specframes = pd.DataFrame({'Data':datay,xlabel:dataedges[:-1],
                               avglabel:avgmodely})

            step = Step(specframes,x=xlabel,y=[datalabel,avglabel],
                color=[datacolor,modelcolor],legend=True,
                y_mapper_type=y_axis_type,x_mapper_type=x_axis_type,
                dash=['solid','solid'],
                plot_width=width,plot_height=height)

            
    else:
        specframes = pd.DataFrame({datalabel:datay,xlabel:dataedges[:-1]})

        step = Step(specframes,x=xlabel,y=[datalabel],
                color=[datacolor],legend=True,
                y_mapper_type=y_axis_type,x_mapper_type=x_axis_type,
                dash=['solid'],
                plot_width=width,plot_height=height)

    step.x_range=Range1d(*x_range)
    step.y_range=Range1d(*y_range)
    step.legend.location='top_right'
    step.ylabel = 'counts'

    #----Plot Errorbars----
    if model is True:
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
        if display is True:
            bplt.show(step)
        else:
            bplt.save(step)
#        if outfile != 'notebook': bplt.curdoc().clear()
    

#----Return----
    return step

#----------------------------------------------------------
def trace(inframe,iteration_type = 'median',itercol = 'iteration',
          weights=None,itmin=0,itmax=None,display=True,
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

    save (bool): optionally turn off saving the plot as an 
                 html file - returns the figure object only (default=True)

    display (bool) : if False, then will not display the figure
                     automatically (re)set to True if outfile='notebook'
                     ignored if save=False

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
                display = True
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
        if display is True:
            bplt.show(p)
        else:
            bplt.save(p)
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
def plot_lines(fig,bins,nlines=50,show=False,**fetchargs):
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

     nlines (int) : maximum number of (grouped) lines to include.
                    if more than nlines (grouped) lines meet the 
                    other criteria, then will plot only the nlines 
                    with highest emissivity

     show (bool) : specify if the line information should be printed to
                   screen as well as plotted

     fetch_lines() input options:

         kT_range (2d tuple) : range of gas temperatures to include 
                               emission lines from

         emissivity_range (2d tuple) : range of line emissivities 
                                       to include

         wavelength_range (2d tuple) : range of line wavelengths 
                             (angstroms) to include emission lines from

         energy_range (2d tuple) : range of line energies (keV) to include
                             emission lines from

         include_lines (list of strings) : list of element names to 
                        include emission lines from. will first drop 
                        all lines from other elements, then apply 
                        remaining criteria (e.g. kT_range).

         redshift (float) : redshift of the source. if non-zero, lines 
                            will be redshifted appropriately before
                            applying any criteria or plotting.


    Output:

     Add vertical lines (at most one per spectral bin) to the provided
     figure. If more than one line per bin, then a single line will be
     plotted with a label that lists all lines in that bin (grouped by
     ion and ionization stage).


    """

    import os
    from bokeh.models.annotations import Span,Label
    from bokeh.models import HoverTool
    from astro_utilities import fetch_lines

    #--fetch lines from atomdb file--
    linedf = fetch_lines(**fetchargs)

    #--get only lines that are within the plot range--

    #-get plot range-
    xmin = float(fig.x_range.start)
    xmax = float(fig.x_range.end)
    ymin = float(fig.y_range.start)
    ymax = float(fig.y_range.end)

    #-set x units-
    if xmax <= 15.0: # epic data (keV scale)
        x = 'energy'
    else: # rgs data (angstrom scale)
        x = 'wavelength'

    #-cut lines-
    linedf = linedf[linedf[x] <= xmax]
    linedf = linedf[linedf[x] >= xmin]

    #-group lines by ion and ionization stage
    plotdf = agg_lines(linedf.groupby(['ion','ionization_stage']),
                       bins,units=x)

    #--truncate to include no more than nlines (grouped) lines, 
    #   keeping those with highest emissivity--
    if (nlines is not None) and (len(plotdf.index)>nlines):
        plotdf = plotdf.nlargest(nlines,'emissivity')
    else:
        nlines = len(plotdf.index)

    #----Print line information----
    if show is True: print plotdf.head(nlines)

    #----Plot lines----
    
    #--calculate label offset--
    xoffset = (xmax - xmin)/300.0
    labely = (ymax - ymin)/100.0 + ymin

    #--plot lines--
    for i,row in plotdf.iterrows():
#        print row[x],row['label'],row['emissivity']
        lspan = Span(location=row[x],dimension='height',
                     line_color='gray',
                     line_dash='dashed',line_width=1)
        fig.add_layout(lspan)
        llabel = Label(x=row[x],y=labely,text_color='gray',
                       angle=90.0,angle_units='deg',
                       text_font_size='10px',x_offset=xoffset,
                       text=row['label'])
        fig.add_layout(llabel)

    return fig

#----------------------------------------------------------
def agg_lines(groupeddf,bins,units='energy'):
    """
    Merge emission lines from the same ions/ionization stage if they
    are close enough (where 'close enough'=in same bin). Must provide
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

#----------------------------------------------------------
def spectrum_from_blobs(df,runpath='../',datacolor='black',save=True,
                        display=True,suffix='99999',datalabel='Data',
                        modelcolor='steelblue',modellabel='Model',
                        bins=0.03,
                        outfile='spectrum_from_blobs.html',ylog=False,
                        xlog=False,logbins=None,datarange=None,width=1000,
                        height=500,lines=True,**lineargs):
    """
    Create and plot spectrum from the blobs in given dataframe

    Author: Kari A. Frank
    Date: May 1, 2017

    Input:

     df (pd.DataFrame) : blob dataframe. must contain, at minimum,
                   all columns in parameters.txt

     suffix (str) : must be string of an integer. will be appended to end
                    of created file names (e.g. deconvolution.suffix)

     other input passed to spectrum()
     

    Output:

     saves and (optionally) displays spectrum as html file.

    Usage Notes:
     - requires input, parameters.txt, start.xmc, events, and exposure
       map files all to be located in runpath directory
     - all new files are created in the working directory, so it is best
       to NOT run this from a run directory, to avoid the possibility 
       of overwriting any files.

   
    """

    #--imports--
    import os

    pwd = os.getcwd()+'/'

    #--create 'fake' deconvolution and other xmc output files--
    # (file names will end in .99999)
    niters=len(np.unique(df.iteration)) # get number iterations
    oldnorm = df.blob_norm.sum() # net norm
    df = fake_deconvolution(df,runpath=runpath,suffix=suffix)

    #--read fake file into xmc and create spectrum--
    
    #-copy and modify input file-
    with open(runpath+'input','r') as file:
        inputscript = file.readlines()

    # find line with '/null' and replace with '/xw'
    indices = [i for i, s in enumerate(inputscript) if '/null' in s]
    index = indices[-1]
    inputscript[index] = inputscript[index].replace('null','xw')
        
    # find line with 'run' (in case blank lines after run)
    indices = [i for i, s in enumerate(inputscript) if 'run' in s]
    index = indices[-1]

    # get start.xmc file name and modify run line
    for p in inputscript[index].split(' '): # in case extra spaces
        if '.xmc' in p:
            startfile = p.rstrip() # remove any newlines
    
    inputscript[index] = 'run spec_start.xmc'
    with open(pwd+'input_spec','w') as file:
        file.writelines(inputscript)

    #-copy and modify the start.xmc file-
    with open(runpath+startfile,'r') as file:
        runscript = file.readlines()

    # find line with 'deconvolvec' (in case blank lines after)
    indices = [i for i, s in enumerate(runscript) if 'deconvolvec' in s]
    index = indices[-1]
        
    # modify deconvolvec line
    lastlineparts = [ p.rstrip() for p in runscript[index].split(' ')]

    if len(lastlineparts)==3: # add iteration
        lastlineparts.append(suffix+'.')
    else: # replace iteration
        lastlineparts[-1] = suffix+'.'

    lastline=' '.join(lastlineparts)
        
    runscript[index] = lastline

    # add path to data/expo calls
    dindices = [i for i, s in enumerate(runscript) if 'data ' in s]
    eindices = [i for i, s in enumerate(runscript) if 'expo ' in s]
    indices = dindices + eindices
    for i in indices:
        parts = runscript[i].split(' ')
        line = parts[0]+' '+runpath+parts[-1]
        runscript[i] = line

    # write modified file
    with open(pwd+'spec_start.xmc','w') as file:
        file.writelines(runscript)

    #-call xmc to create the spectrum-
    # (for now do it manually)
    raw_input("In another terminal: \n\tFrom this directory"
              " ("+pwd+"),\n"
              "\trun xmc with 'xmc < input_spec'\n"
              "\tWhen spectrum_"+suffix+".fits is produced\n "
              "\tand the *."+suffix+" files are overwritten by xmc, "
              "then kill xmc.\n"
    "\tReturn to this terminal and press Enter to continue.")

    #--read in the new deconvolution file to get xspec normalization--
    if not os.path.isfile(pwd+'parameters.txt'):
        os.link(runpath+'/parameters.txt',pwd+'parameters.txt')
    dfnew = merge_output(runpath=pwd,save=False)
    dfnew = dfnew[dfnew.iteration==int(suffix)]
    newnorm = dfnew.blob_norm.sum()
    xspecscale = oldnorm/newnorm
#    print 'xspecscale = ',xspecscale
    
    #--create spectrum--
    fig = spectrum(runpath=pwd,smin=int(suffix),smax=int(suffix),
                   datacolor=datacolor,datalabel=datalabel,
                   modellabel=modellabel,save=save,
                   model=True,display=display,
                   modelcolor=modelcolor,
                   bins=bins,scale=[1.0,xspecscale/float(niters)],
                   outfile=outfile,ylog=ylog,xlog=xlog,logbins=logbins,
                   datarange=datarange,width=width,height=height,
                   lines=lines,**lineargs)
    
    return fig

#----------------------------------------------------------
def spectra(spectra,colors=['black','steelblue','firebrick'],
            labels = None,dashes=None,
            save=True,display=True,scale=1.0,
            outfile='spectra.html',ylog=False,xlog=False,logbins=None,
            width=1000,height=500,lines=True,**lineargs):

    """
    Plot the given spectra.

    Author: Kari A. Frank
    Date: May 4, 2017

    Input:

     spectra (list of tuples) : should contain a list of tuples, where each
                 tuple represents a histogram, of the form 
                 (y,yerrors,yedges), e.g. as output
                 from either xmcfiles.read_spectra() or 
                 wrangle.make_spectrum()

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

    scale (float or list of floats) : scale the each spectrum 
                    histogram by a constant. if a single float is provided,
                    all spectra will be scaled by the same value. if a list
                    is provided, it must have one element per spectrum to
                    to be plotted (same length as spectra list)

    outfile (str) : name of output html file, or 'notebook' if plotting to 
                    an open Jupyter notebook

    save (bool) : if save is False, then will not display the final figure
                  or save it to a file

    display (bool) : if False, then will not display the figure
                  automatically (re)set to True if outfile='notebook'
                  ignored if save=False                  

    dashes (list of strings) : provide the dash type for each spectrum
                   to be plotted

    lines (bool) : plot the 10 stongest common emission lines within 
                   the x-axis range. **lineargs will pass any extra 
                   arguments directly to plot_lines()

    xlog, ylog (bool) : specify if x and y axes should use log scale

    Output:
     - plots the spectra for to an interactive plot (if save=True)
     - Returns the figure object

    Usage Notes:
     - Will automatically plot error bars if provided, e.g. 
       for an average spectrum returned by read_spectra()
     - For now, all given spectra must have the same binsize and
       start at the same x-value. Upper limit on x-axis will be 
       limited to that of the spectrum with the smallest range.

    Example:
    """

    #----Import Modules----
    import os
    import astropy.io.fits as fits
    from bokeh.charts import Step,color
    from bokeh.models.ranges import Range1d
    from file_utilities import ls_to_list

    #----Set defaults----
    if isinstance(scale,float):
        scale = [scale]*len(spectra)

    if labels is None:
        labels = ['Spec'+str(i) for i in xrange(len(spectra))]

    if dashes is None:
        dashes = ['solid']*len(spectra)
        
    #----Set up Plot----
    if save is True:
        if outfile != 'notebook':
            bplt.output_file(outfile)
        else:
            display = True
            bplt.output_notebook()
            
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
    
    x_axis_type = 'linear'
    y_axis_type = 'linear'
    ymin = 1e10 # dummy ymin
    ymax = 0.0 # dummy ymax
    if xlog is True:
        x_axis_type = 'log'
    xmin = min([np.min(t[2]) for t in spectra])
    xmax = max([np.max(t[2]) for t in spectra])
    x_range = (xmin,xmax)
    if ylog is True:
        y_axis_type = 'log'
        for t in spectra: # calculate non-zero min and max
            t0=t[0]
            ymin = min(ymin,np.min(t0[np.where(t0>0.0)]))
            ymax = max(ymax,np.max(t0[np.where(t0>0.0)]))
        ymax = 1.5*ymax
        y_range = (ymin,ymax)
    else: 
        for t in spectra:
            t0=t[0]
            ymin = min(ymin,np.min(t0))
            ymax = max(ymax,np.max(t0))
        ymax = 1.1*ymax
        y_range = (ymin,ymax)
    if xmax < 15.0:
        xlabel = 'Energy (keV)'
    else:
        xlabel = 'Wavelength (Angstroms)'

    #----Plot Spectra as Step Chart----
    
    # assumes all spectra have the same binsize
    edges = np.round(spectra[0][2],decimals=3)[:-1]
    xbinsizes = edges[1:]-edges[:-1]
    xbinsize = max(xbinsizes)
    for si,s in enumerate(spectra):
        # find matching bins (indices, as Bool array)
        mask = np.isclose( np.round(s[2],decimals=2)[:-1],edges,
                           rtol=0.0,atol=xbinsize/2.0) 
#        mask = np.in1d( np.round(s[2],decimals=3)[:-1],edges) 
        
#        print len(mask), len(edges),len(s[0])
#        print s[2][np.where(mask==False)]
#        print edges[np.where(mask==False)]
        # update spectrum
        newx = s[2][np.where(mask)]
        # update y values
        newy = s[0][np.where(mask)]
#        print 'len newx = ',len(newx)
        # update errors
        if s[1] is not None:
            newyerr = s[1][np.where(mask)]
        else:
            newyerr = s[1]
        spectra[si] = (newy,newyerr,newx)

        print 'total counts = ',newy.sum()
        
    #--create dictionary of spectra histograms--

#    spectra0 = spectra[0]
#    edges = spectra0[2]
    specy = [s[0] for s in spectra]
    specframes = dict(zip(labels,specy)) # add y values to dict
    specframes[xlabel] = edges#[:-1] # add x values to dict
    print len(edges)    
    specframes = pd.DataFrame.from_dict(specframes)
    
    # To Do - figure out how to specify a color for each column
    #         (turns out it is extremely difficult)
    step = Step(specframes,x=xlabel,y=labels,
                legend=True,color=colors,
                y_mapper_type=y_axis_type,x_mapper_type=x_axis_type,
                dash=dashes,ylabel='counts',
                plot_width=width,plot_height=height)
    
    step.x_range=Range1d(*x_range)
    step.y_range=Range1d(*y_range)
    step.legend.location='top_right'
#    step.ylabel = 'counts'

    #----Plot Errorbars----
    xbins = edges[:-1]+xbinsizes/2.0
    for si,s in enumerate(spectra):
        if s[1] is not None:
            errorbar(step, xbins, s[0], xerr=None, yerr=s[1], 
                     color='gray', # colors[si], 
                     point_kwargs={}, error_kwargs={})

    #----Plot Emission Lines----
    if lines is True:
        step = plot_lines(step,edges,**lineargs)

    #possible attributes to Chart are above, background_fill_alpha, background_fill_color, below, border_fill_alpha, border_fill_color, disabled, extra_x_ranges, extra_y_ranges, h_symmetry, height, hidpi, left, lod_factor, lod_interval, lod_threshold, lod_timeout, min_border, min_border_bottom, min_border_left, min_border_right, min_border_top, name, outline_line_alpha, outline_line_cap, outline_line_color, outline_line_dash, outline_line_dash_offset, outline_line_join, outline_line_width, plot_height, plot_width, renderers, right, sizing_mode, tags, title, title_location, tool_events, toolbar, toolbar_location, toolbar_sticky, v_symmetry, webgl, width, x_mapper_type, x_range, xlabel, xscale, y_mapper_type, y_range, ylabel or yscale


    #----Show the Plot----
    if save is True:
        if display is True:
            bplt.show(step)
        else:
            bplt.save(step)
#        if outfile != 'notebook': bplt.curdoc().clear()
    

#----Return----
    return step

#----------------------------------------------------------
def standard_spectra(runpath='../',itmin=1,itmax=None,
                     save=True,display=True,
                     outfile='spectra.html',ylog=False,xlog=False,
                     logbins=None,bins=0.03,lastiter=True,
                     legacy=False,
                     width=1000,height=500,lines=True,**lineargs):

    """
    Plot standard data, average model, and most recent model spectra

    Author: Kari A. Frank
    Date: May 5, 2017

    Input:

     spectra (list of tuples) : should contain a list of tuples, where each
                 tuple represents a histogram, of the form 
                 (y,yerrors,yedges), e.g. as output
                 from either xmcfiles.read_spectra() or 
                 wrangle.make_spectrum()

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

    scale (float or list of floats) : scale the each spectrum 
                    histogram by a constant. if a single float is provided,
                    all spectra will be scaled by the same value. if a list
                    is provided, it must have one element per spectrum to
                    to be plotted (same length as spectra list)

    outfile (str) : name of output html file, or 'notebook' if plotting to 
                    an open Jupyter notebook

    save (bool) : if save is False, then will not display the final figure
                  or save it to a file

    display (bool) : if False, then will not display the figure
                  automatically (re)set to True if outfile='notebook'
                  ignored if save=False                  

    lines (bool) : plot the 10 stongest common emission lines within 
                   the x-axis range. **lineargs will pass any extra 
                   arguments directly to plot_lines()

    xlog, ylog (bool) : specify if x and y axes should use log scale

    lastiter (bool) : if True will plot the spectrum from the 
                   last iteration in the iteration range

    legacy (bool) : assume old-style spectrum file names, i.e.
                     spectrum_iter/100.fits


    Output:
     - plots the spectra for to an interactive plot (if save=True)
     - Returns the figure object

    Usage Notes:
     - Will automatically plot error bars if provided, e.g. 
       for an average spectrum returned by read_spectra()
     - For now, all given spectra must have the same bins

    Example:
    """
    # - Imports -
    from xmcfiles import read_spectra
    
    # - Defaults -
    if logbins is None: logbins = xlog

    # - read model spectra -
    if itmin == 0: itmin = 1
    shists = read_spectra(runpath=runpath,itmin=itmin,itmax=itmax,
                          average=True,bins=bins,logbins=logbins,
                          legacy=legacy)
    lasthist = shists[-2]
    avghist = shists[-1]
    
    # - read data spectrum -
    dhist = read_spectra(runpath=runpath,itmin=0,itmax=0,
                             average=False,bins=bins,logbins=logbins)

    datahist = dhist[0]

    # - plot spectra -
    if lastiter is True:
        speclist = [datahist,avghist,lasthist]
        labs = ['Data','Model (average)','Model (latest)']
    else:
        speclist = [datahist,avghist]
        labs = ['Data','Model (average)']
    sfig = spectra(speclist,
                        labels=labs,
                        display=display,
                        outfile=outfile,
                        ylog=ylog,xlog=xlog,
                        lines=lines,**lineargs)

    return sfig
