"""
Module of functions for plotting xmc results.

Contains the following functions:

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
def chi2(runpath='./',itmin=0,itmax=None):
    """
    Author: Kari A. Frank
    Date: October 20, 2015
    Purpose: plot simple scatter plot of chi2 vs iteration, directly from 
             the statistic.* files output by xmc.  This is mainly intended
             as a quick way to check the convergence of an ongoing xmc run.

    Usage: chi2(runpath='./')

    Input:

     runpath:     string containing the relative path to xmc run folder, which
                  contains the statistic.* files.

     itmin/itmax: minimum/maximum iteration to plot

    Output:
     - plots chi2 vs iteration for the specified iteration to an interactive plot

    Usage Notes:
     - must close and save (if desired) the plot manually

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
    bplt.output_file('chi2_vs_iteration.html')
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    fig = bplt.Figure(tools=TOOLS)
    fig.xaxis.axis_label='iteration'
    fig.yaxis.axis_label='chi2/dof'

#----Plot chi2 vs iteration----
#    statframe.plot(x='iteration',y='redchi2',kind='scatter')
#    plt.show()

    fig.circle(x=statframe['iteration'],y=statframe['redchi2'])
    bplt.show(fig)#,new='window',browser='firefox --no-remote')
    bplt.curdoc().clear()

#----Return----
    return statframe

#----------------------------------------------------------
def scatter(inframe,x,y,sampling=2000.,agg=None,aggcol=None,save=True,
            width=600,height=600,source=None,tools=None,size=5,
            xlog='auto',ylog='auto'):
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

   Output:
   - plots x vs y to an interactive plot (if Datashader used or
     save=True)
   - returns the figure object and the ColumnDataSource object

  Usage Notes:
  - must close and save (if desired) the plot manually
  - may not work properly if trying to plot too many points, unless
    using Datashader
  - if using Datashader agg options:
    - MUST BE RUNNING FROM A JUPYTER NOTEBOOK!! bokeh plot output 
      will be set to output_notebook()and an html file will not be
      created or saved.
    - it is necessary to call bplt.output_notebook() in the 
      current notebook session prior to calling this function. otherwise
      the figure will not automatically be displayed inline.
    - a jupyter notebook can be started by simply typing 
      'jupyter notebook' on the command line (assuming you have 
      jupyter notebooks installed)

  Example:
 
  """
#    print len(inframe[x]),len(inframe[y])
    
#----Import Modules and Set Output Location----
    # switch to specify if using Datashader

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
        bplt.output_notebook()
#    else:
    if (save is True) and (agg is None): 
        bplt.output_file(x+'_vs_'+y+'.html')

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
        norders = np.log10(x_rng[1]) - np.log10(x_rng[0])
        if norders > 2.0:
            xlog = True
        else:
            xlog = False
    y_rng = (df[y].min(),df[y].max())
    if ylog == 'auto':
        norders = np.log10(y_rng[1]) - np.log10(y_rng[0])
        if norders > 2.0:
            ylog = True
        else:
            ylog = False

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
    if agg is None:
        fig.circle(x,y,source=source,size=size)
        if save is True: 
            bplt.show(fig)#,new='window',browser='firefox --no-remote')
            bplt.curdoc().clear()

    else:
        def image_callback(x_range,y_range,w,h):
            # re-draw canvas
#            cvs = ds.Canvas(plot_width=w,plot_height=h,x_range=x_range,
#                           y_range=y_range)
            cvs = ds.Canvas(x_range=x_range,y_range=y_range)
            # re-aggregate
            aggr = cvs.points(df,x,y,agg=fdict[agg](aggcol))
            img = tf.interpolate(aggr,cmap=['lightblue','darkblue'],
                                 how='log')
            return tf.dynspread(img,shape='circle',max_px=size,
                                threshold=0.85,how='saturate')
        interimg = InteractiveImage(fig,image_callback)
        if save is True: bplt.show(fig)

#----Return----
    return fig,source


#----------------------------------------------------------
def traceplots(dframe,sampling=1000.0,agg=None,aggcol=None,
               outfile='traceplots.html',save=True,size=50):
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
             - 'sampling' (take random subsample of rows, default)
             - 'contour' (plot density contours instead of points - does not
                use linked brushing) -- not yet implemented

      columns: DEPRECATED list of dataframe column names to include in 
               the plots. default is to use all columns

      save: optionally specify to return the figure list, without 
            saving the plot file or plotting in a browser. ignored 
            if agg != None

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
          
    outfile: name of output html file (ignored if agg!=None)
    
    size: size of data points to plot

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

#----Sample the Data----
    if sampling is None: sampling = len(dframe.index)
    if len(dframe.index) > sampling:
        df = dframe.sample(sampling)
    else:
        df = dframe

#----Set up plot----
    if (agg is None) and (save is True): bplt.output_file(outfile)
    source = ColumnDataSource(df)
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
              density=False,xlog='auto',**kwargs):
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

     bins:        optionally specify the number of bins (default=30)

     height,width: height and width of each histogram

     save:        optionally turn off opening and saving the plot,
                  returns the figure object only (default=True)

     infig:       optionally pass an initialized figure object to plot on 
                  (allows plotting multiple dataseries on the same figure) 
                  (default=None)

     density:     passed to histogram. if True, then histogram is 
                  normalized.

     xlog:        make x-axis bins uniform in log10 space. options are:
                  - 'auto' (default) -- try to automatically determine
                    which is best based on the data range
                  - True -- force log bins
                  - False -- force linear bins

     outfile:     string name of html file to save figure to. set to 
                  outfile='notebook' to plot to jupyter notebook.

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

#----Set up Log Bins----
    rng = (dataseries.min(),dataseries.max())
    if xlog == 'auto':
        norders = np.log10(rng[1]) - np.log10(rng[0])
        if norders > 2.0:
            xlog = True
        else:
            xlog = False
#    print 'xlog = ',xlog

    if xlog == True: # set up log bins
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

#----Set up Plot----
    if save:
        if (outfile != 'notebook'): 
            bplt.output_file(outfile)
    #    TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
        else:
            bplt.output_notebook()

    if infig is None:
        fig = bplt.Figure(tools=tools,plot_width=width,plot_height=height,
                          x_axis_type=bintype,x_range=(rng[0],rng[1]))
        fig.xaxis.axis_label=dataseries.name
    else:
        fig = infig
    if weights is not None:
        fig.yaxis.axis_label=weights.name
        if np.log10(max(histy))>3: 
            fig.yaxis.formatter=PrintfTickFormatter(format = "%1.1e")

#----Format Ticks----
#    fig.yaxis.formatter=PrintfTickFormatter(format="%4.1e")
#    fig.xaxis.formatter=PrintfTickFormatter(format="%4.1e")
#    if infig is None:
#        fig.yaxis.formatter=format_ticks(histy)
#        fig.xaxis.formatter=format_ticks(binedges)

#---Plot the histogram----
    h = fig.quad(top=histy,bottom=0,left=binedges[:-1],right=binedges[1:],
                 color=color,**kwargs)

    if save:
        bplt.show(fig)#,new='window',browser='firefox --no-remote')
        if (outfile != 'notebook'): 
            bplt.curdoc().clear()

#----Return----
    return fig

#----------------------------------------------------------
def histogram_grid(dframe,weights=None,bins=100,height=300,width=400,
                   ncols=2,outfile='histogram_grid.html',**kwargs):
    """
    Author: Kari A. Frank
    Date: October 28, 2015
    Purpose: plot interactive matrix of weighted histograms from 
             given dataframe

    Usage: histogram_grid(dframe,weights=None,bins=30,**kwargs):

    Input:

     dframe:      input dataframe

     weights:     optionally provided a pandas series of weights which 
                  correspond to the values in datacolumn (e.g. emission 
                  measure)

     bins:        optionally specify the number of bins (default=30)

     height,width: height and width of each histogram, passed to histogram()

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

    Example:

    """

#----Import Modules----
    import math

#----Set up plot----
    if outfile != 'notebook':
        bplt.output_file(outfile)
    else:
        bplt.output_notebook()

#----Initialize empty figure list----
    figlist=[]

#----Fill in list of figures----
    for column in dframe:
        newfig = histogram(dframe[column],weights=weights,
                           save=False,height=height,width=width,
                           bins=bins,**kwargs)
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
def spectra(runpath='./',smin=0,smax=None,datacolor='black',
            modelcolor='cornflowerblue',lastmodelcolor='crimson',nbins=None):
    """
    Author: Kari A. Frank
    Date: November 18, 2015
    Purpose: Make a plot of the data and (average) model spectra compared

    Usage: spectra(runpath='./',smin=0,smax=None,nbins=350.0,
                   datacolor='black',modelcolor='steelblue',
                   lastmodelcolor='darkred')

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

    Output:
     - plots the spectra for to an interactive plot
     - Returns the data spectrum numpy array

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
    bplt.output_file('spectrum.html')
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
    bplt.show(fig)
    bplt.curdoc().clear()

#----Return----
    return data_wave
#----------------------------------------------------------
def evolution(inframe,columns=None,iteration_type = 'median',weights=None,save=True,
              ncols=4):
    """
    evolution()
 
   Author: Kari A. Frank
    Date: April 15, 2016
    Purpose: plot parameters vs iteration

   Usage: scatter(inframe,columns=None,iteration_type='median')

   Input:
 
    inframe (DataFrame):  pandas DataFrame containing the data columns to be plotted   

    columns (string or list of strings): names of columns to be plotted
             
    iteration_type (string): method used to combine parameter values within
         each iteration
         - 'none' -- plot all blobs separately (very slow!)
         - 'median' -- takes median from the iteration
         - 'average' -- takes average from the iteration
         - 'stdev' -- plots the standard deviation from each iteration
         - 'total' -- sums all blobs within each iteration

    weights (pd.Series): weights to apply to each blob (e.g. emission measure)

    save:        optionally turn off opening and saving the plot as an 
                  html file - returns the figure object only (default=True)

    ncols: optionally specify number of columns in plot grid if plotting more than
           one parameter. ignored if len(columns)=1

   Output:
   - plots parameter vs iteration to an interactive plot
   - returns the figure object

  Usage Notes:
  - must close and save (if desired) the plot manually
  - may not work properly if trying to plot too many points

  Example:
 
  """
    

    #----Combine blobs in each iteration----

    
    #----Set up plot----
    if (save is True):
        if (len(columns)>1): 
            bplt.output_file('evolution_grid.html')
            itersave=False
        else:
            itersave=True

    source = ColumnDataSource(iterframe)
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    w = 300 # plot width 
    h = 300 # plot height

    #----Initialize empty figure list----
    figlist=[]

    #----Fill in list of figures----
    for col in columns:

        #----Plot parameter vs iteration----
        newfig = scatter(iterframe,'iteration',col,agg='none',save=itersave)
        figlist=figlist+[newfig]    

    #----Reshape list into a 4 column array---

    #--define new shape--
    nfigs = len(figlist)
    nrows = int(math.ceil(float(nfigs)/float(ncols)))

    #--pad list with None to have nrows*ncols elements--
    figlist = figlist+[None]*(nrows*ncols-nfigs)

    #--reshape list--
    figarr = [figlist[ncols*i:ncols*(i+1)] for i in range(nrows)]

    #----Plot histograms----
    p = bplt.gridplot(figarr)
    if save is True: 
        bplt.show(p)
        bplt.curdoc().clear()

    return fig
