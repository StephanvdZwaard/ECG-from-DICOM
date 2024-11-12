'''
Tools to plot ECG signals.

The primary aim is not to make nice ECG visuals but to create input data for
image analysis.

The code heavily borrows from the following brilliant github
`repository <https://github.com/marcodebe/dicom-ecg-plot>`_.
'''

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# imports
import warnings
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from typing import (
    Callable, List, Literal, Type, Union, Tuple, Self, Dict, Optional, Any,
)
from ecgprocess.errors import (
    NotCalledError,
    is_type,
    Error_MSG,
)
from ecgprocess.constants import (
    ProcessDicomNames as PDNames,
)
from ecgprocess.utils.general import (
    assign_empty_default,
    _update_kwargs,
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ECGDrawing(object):
    '''
    Takes a called `ECGDICOMReader` instance and plots the lead-specific
    ECG signals. Includes a method to map the figure to a 3-dimensional numpy
    array.
    
    Attributes
    ----------
    mm_mv : `float`, default 10
        The scaling factor applied to the ECG signals, mapping these from
        1mm/mV to mm_vv * mm/mv.
    paper_w : `float`, default 297 mm
        The figure width.
    paper_h : `float`
        The figure height, default 210 mm
    width : `float`
        The plotting area width, default 250 mm
    height : `float`
        The plotting area height, default 170 mm
    grid_color : `dict` {`minor`: colour_1, `major`: color_2}
        The grid line colours.
    grid_linewidth : `dict`
        The grid line width.
    text_pad_x : `float`, default 40
        Padding of the text x-axis coordinate.
    
    Methods
    -------
    to_numpy(crop)
        maps a matplotlib image to a numpy array.
    
    Notes
    -----
    Calling the class instance will check the ECG unit, and if needed convert
    microvolts (µV) to millivolts (mV).
    
    '''
    
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def __init__(self, update_keys:dict[str,str]|None=None) -> None:
        """
        Initialises a new instance of `ECGDrawing`.
        
        Parameters
        ----------
        update_keys: dict [`str`, `str`], default `NoneType`
            A dictionary to remap lead names: [`old`, `new`]
        """
        # scaling factor for 1/mV to FACTOR/mV
        setattr(self, 'update_keys', update_keys)
        # scaling factor for 1/mV to FACTOR/mV
        setattr(self, PDNames.WAVE_SCALING, 10)
        # #### set sensible plotting defaults matching an actual ECG printout
        self.inch_mm = 24.5
        # standard A4 size in mm - landscape!
        self.paper_w = 297.0
        self.paper_h = 210.0
        # Dimensions in mm of plot area
        self.width = 250.0
        self.height = 170.0
        # The grid constants
        self.grid_color = {'minor': '#ff5333', 'major': '#d43d1a'}
        self.grid_linewidth = {'minor': .1, 'major': .2}
        # the text constants
        self.text_pad_x = 40
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def __str__(self):
        CLASS_NAME = type(self).__name__
        return (f"{CLASS_NAME}"
                )
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def __repr__(self):
        CLASS_NAME = type(self).__name__
        return (f"{CLASS_NAME}()"
                )
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def __call__(self, ecgreader:Callable,
                 wave_type:Literal['rhythm', 'median']='rhythm',
                 start_pos:Literal['first', 'continues']= 'first',
                 add_grid:bool=False, minor_axis:bool=True,
                 image_layout:list[list[str]] | None=None,
                 auto_margins:bool=True,
                 verbose:bool=True,
                 ax:plt.Axes | None=None,
                 ) -> Self:
        '''
        Creates an ECG drawing using either the waveforms `rhythm` or the median
        beats signals.
        
        Parameters
        ----------
        ecgreader : Callable
            An instance of the `ECGDICOMReader`, or a similar, data class that
            has already been called to process the ECG data.
        wave_type : {'rhythm', 'median'}, default `rhythm`
            The type of ECG signal to plot.
        start_pos: {`first`, `continues`}, default `first`
            Whether multiple sequential leads should start at the first
            time point, or whether the next lead should start at the time
            the previous lead stopped.
        add_grid : `bool`, default `False`
            Whether to annotate the ECG image with a canonical grid (with small
            squares indicating 0.04 seconds (x-axis) and 0.1 mV (y-axis) and
            the larger squares a 5 multiple. Depending on the intended
            deployment setting and the applied pre-processing you may want to
            remove the grid.
        minor_axis: `bool`, default `True`
            Whether to add the minor grid axes.
        image_layout : `list` [`list` [`str`]]
            A nested lists organising the image by ECG lead names. When
            wave_type is `median` please ensure the supplied list has the same
            number of columns in each row.
        ax : plt.Axes, default `NoneType`
            An optional `matplotlib.axes.Axes` instance on which the ECG
            signals are plotted. If ommited will simply take an A4 x/y axes
            aspect ratio and mimic a cononical ECG image.
            canonical ECG image (using an A4 size y and x-axis aspect ratio).
        verbose : `bool`, default `False`
            Prints missing files if skip_missing is set to `True`.
        
        Attributes
        ----------
        ecg_signal : `dict` [`str`, `np.ndarray`]
            The ECG signals depicted in the figure. Organised as a dictionary
            with keys matching the lead name and np.ndarray as values recording
            the ECG measurement.
        ecg_unit : `str`
            The measurement unit - should be `millivolt`. Included here for
            checking/debugging.
        image_layout : `list` [`list` [`str`]]
            A nested lists organising the image by ECG lead names.
        wave_type : `str`
            The type of ECG signal.
        fig : `plt.figure`
            The matplotlib figure.
        ax : `plt.axes`
            The matplotlib axes.
        
        Returns
        -------
        self : `ECGDICOMDraw` instance
            Returns the class instance with updated attributes.
        
        Notes
        -----
        To replicate an actual ECG image the signals are plotted on a single
        x-axis. When plotting waveforms (i.e. wave_type=`rhythm`) which are
        periodic the x-axis is simply divided in equally spaced sections and
        the ECG signals are sampled based on the x-axis coordinates which
        match the specific section. Median beats are not periodic and hence
        the class simply plots the entire signal. This does mean that the
        image_layout should have the same number of elements/columns for each
        row.
        
        Example
        -------
        >>> layout = [
        >>>     ['I', 'aVR', 'V1', 'V4'],
        >>>     ['II', 'aVL', 'V2', 'V5'],
        >>>     ['III', 'aVF', 'V3', 'V6'],
        >>> ]
        >>> ecgdicomreader = ECGDICOMReader()
        >>> ecgdicomreader = ecgdicomreader(path=path, skip_empty=skip_empty,
        >>>                                verbose=verbose)
        >>> plt.ion()
        >>> artist = ECGDrawing(add_grid=True, image_layout=layout)
        >>> plt.close()
        '''
        # constants
        WARN1 = ('The Waveform unit attribute is unavailable. Drawing assumes '
                 f'measurements are in {PDNames.MILLIVOLT}.')
        WARN2 = ('The Waveform unit  attribute is not recognised: `{}`. Drawing '
                 f'assumes measurements are in {PDNames.MILLIVOLT}.')
        # check input and assign to self
        is_type(verbose, bool)
        is_type(start_pos, str)
        is_type(wave_type, str)
        is_type(add_grid, bool)
        is_type(minor_axis, bool)
        is_type(ax, (type(None), plt.Axes))
        is_type(image_layout, (type(None), list))
        self.verbose = verbose
        # confirm wave_type
        setattr(self, PDNames.WAVE_TYPE, wave_type)
        WAVE_TYPE = [PDNames.WAVETYPE_RHYTHM,
                    PDNames.WAVETYPE_MEDIAN,
                    ]
        if not getattr(self, PDNames.WAVE_TYPE) in WAVE_TYPE:
            raise ValueError(Error_MSG.CHOICE_PARM.\
                             format('wave_type', ', '.join(WAVE_TYPE)))
        # set layout if needed
        if image_layout is None:
            image_layout = [
                ['I', 'aVR', 'V1', 'V4'],
                ['II', 'aVL', 'V2', 'V5'],
                ['III', 'aVF', 'V3', 'V6'],
            ]
        setattr(self, PDNames.PLOT_LAYOUT, image_layout)
        # confirm start pos
        START_POS = [PDNames.START_POS_CONT, PDNames.START_POS_FIRST]
        if not start_pos in START_POS:
                raise ValueError(Error_MSG.CHOICE_PARM.\
                                 format('start_pos', ', '.join(START_POS)))
        # #### get ecg data
        # confirm class has been called - result_dict whould always be present
        if not hasattr(ecgreader, PDNames.RESULTS_DICT):
            raise NotCalledError(f"`{PDNames.ECG_READER}` __call__ method has "
                                 "not been invoked.")
        setattr(self, PDNames.ECG_READER, ecgreader)
        if getattr(self, PDNames.WAVE_TYPE) == PDNames.WAVETYPE_RHYTHM:
            setattr(self, PDNames.ECG_SIGNAL,
                    getattr(getattr(self, PDNames.ECG_READER),
                            PDNames.LEAD_VOLTAGES)
                    )
            setattr(self, PDNames.PLOT_SAMPLING_NUMBER,
                    getattr(getattr(self, PDNames.ECG_READER),
                               PDNames.RESULTS_DICT)[PDNames.SAMPLING_NUMBER]
                    )
        else:
            setattr(self, PDNames.ECG_SIGNAL,
                    getattr(getattr(self, PDNames.ECG_READER), PDNames.LEAD_VOLTAGES2)
                    )
            # multiplying this by the maximum number of columns in any row
            # of plot_layout
            setattr(self, PDNames.PLOT_SAMPLING_NUMBER,
                    getattr(
                        getattr(self, PDNames.ECG_READER),
                        PDNames.RESULTS_DICT)[PDNames.SAMPLING_NUMBER_M]*\
                    max(len(l) for l in getattr(self, PDNames.PLOT_LAYOUT))
                    )
        # #### convert microvolts (µV) to millivolts (mV)
        try:
            unit = getattr(getattr(self, PDNames.ECG_READER),
                           PDNames.RESULTS_DICT)[PDNames.LEAD_UNITS]
        except KeyError:
            try:
                unit = getattr(getattr(self, PDNames.ECG_READER),
                               PDNames.RESULTS_DICT)[PDNames.LEAD_UNITS2]
            except KeyError:
                warnings.warn(WARN1)
                unit = PDNames.MILLIVOLT
        if unit == PDNames.MICROVOLT:
            setattr(self, PDNames.ECG_SIGNAL,
                    {k:v/1000 for k,v in\
                     getattr(self, PDNames.ECG_SIGNAL).items()}
                    )
            pass
        elif unit == PDNames.MILLIVOLT:
            # already correct units
            pass
        else:
            warnings.warn(WARN2.format(unit))
        # raise warning if unit is not helpful
        if unit is None or pd.isna(unit):
            warnings.warn(WARN2.format(unit))
        # add unit
        setattr(self, PDNames.PLOT_UNIT, unit)
        # #### optionally updated the lead names
        if self.update_keys is not None:
            setattr(self, PDNames.ECG_SIGNAL,
                    {self.update_keys.get(k, k): v for k, v in\
                  getattr(self, PDNames.ECG_SIGNAL).items()}
                    )
        # #### create figure
        self._set_canvas(auto_margins=auto_margins, ax=ax)
        if add_grid == True:
            self._draw_grid(minor_axis=minor_axis)
        # #### draw ecg signal
        self._draw_signal(layout=getattr(self, PDNames.PLOT_LAYOUT),
                          start_pos=start_pos,
                          )
        # #### return self
        return self
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def _set_canvas(self, auto_margins:bool=True,
                    ax:plt.Axes | None=None,
                    kwargs_subplots:dict[Any, Any] | None=None,
                    kwargs_subplots_adjust:dict[Any, Any] | None=None,
                    ) -> None:
        '''
        Creates the matplotlib figure and axes objects and assigns this to the
        class instance.
        
        Parameters
        ----------
        auto_margins : `bool`, default `True`
            Whether the margins should be calculated based on the supplied
            figure and plotting area dimensions.
        ax : plt.Axes, default `NoneType`
            An optional `matplotlib.axes.Axes` instance on which the ECG
            signals are plotted. If ommited will simply take an A4 x/y axes
            aspect ratio and mimic a cononical ECG image.
            canonical ECG image (using an A4 size y and x-axis aspect ratio).
        kwargs_*_dict : `dict` [`any`, `any`], default `NoneType`
            Optional arguments supplied to the various plotting functions:
                kwargs_subplots         --> plt.subplots
                kwargs_subplots_adjust  --> figure.subplots_adjust
        
        Attributes
        ----------
        fig : `plt.figure`
            The matplotlib figure.
        ax : `plt.axes`
            The matplotlib axes.
        '''
        is_type(ax, (type(None), plt.Axes))
        is_type(auto_margins, bool)
        # ### settings empty dict defaults
        kwargs_subplots, kwargs_subplots_adjust = assign_empty_default(
            [kwargs_subplots, kwargs_subplots_adjust], dict)
        # #### set canvas
        if ax is None:
            # using a default A4 landscape page
            f, axes = plt.subplots(**kwargs_subplots)
        else:
            # assign ax and its figure
            axes = ax
            f = ax.figure
            # update height and width with the ax properties
            pos = axes.get_position()
            _, _, width_n, height_n = pos.bounds
            # Get the figure dimensions in inches
            fig_width, fig_height = f.get_size_inches()
            # Convert normalized dimensions to inches
            setattr(self, PDNames.PLOT_WIDTH,
                    np.round(width_n * fig_width * self.inch_mm))
            setattr(self, PDNames.PLOT_HEIGHT,
                    np.round(height_n * fig_height * self.inch_mm))
            setattr(self, PDNames.PAPER_WIDTH,
                    np.round(fig_width * self.inch_mm))
            setattr(self, PDNames.PAPER_HEIGHT,
                    np.round(fig_height * self.inch_mm))
        # calculate marings
        if auto_margins == True:
            margin_left =  .5 * (self.paper_w -
                                 getattr(self, PDNames.PLOT_WIDTH))
            # margin_right = .5 * (paper_w - width)
            margin_bottom = 10.0
            # Normalized in [0, 1]
            self.left = margin_left / self.paper_w
            self.right = self.left +\
                getattr(self, PDNames.PLOT_WIDTH)/self.paper_w
            self.bottom = margin_bottom / self.paper_h
            self.top = self.bottom + self.height / self.paper_h
            # apply margins
            new_kwargs_adjust = _update_kwargs(
                update_dict=kwargs_subplots_adjust,
                left=self.left,
                right=self.right,
                top=self.top,
                bottom=self.bottom,
            )
            f.subplots_adjust(**new_kwargs_adjust)
        # set the plotting axis
        axes.set_ylim([0, getattr(self, PDNames.PLOT_HEIGHT)])
        # axes.set_ylim([0, 5000])
        axes.set_xlim([0, getattr(self, PDNames.PLOT_SAMPLING_NUMBER) - 1])
        # set to self
        setattr(self, PDNames.PLOT_FIG, f)
        setattr(self, PDNames.PLOT_AXES, axes)
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def _draw_grid(self, minor_axis:bool=True) -> None:
        '''
        Draws the minor and major grid lines of an ECG image.
        
        Parameters
        ----------
        minor_axis : `bool`, default `True`
            Whether the minor grid lines should be drawn.
        '''
        # do we want to add minor axes - these should have a stepsize of 1 in
        # ECG images
        if minor_axis == True:
            getattr(self, PDNames.PLOT_AXES).xaxis.set_minor_locator(
                plt.LinearLocator(int(self.width + 1))
            )
            getattr(self, PDNames.PLOT_AXES).yaxis.set_minor_locator(
                plt.LinearLocator(int(getattr(self, PDNames.PLOT_HEIGHT)+ 1))
            )
        # annotating the major axes - these should have a stepsize of 5 in
        # ECG images
        getattr(self, PDNames.PLOT_AXES).xaxis.set_major_locator(
            plt.LinearLocator(int(self.width / 5 + 1))
        )
        getattr(self, PDNames.PLOT_AXES).yaxis.set_major_locator(
            plt.LinearLocator(int(getattr(self, PDNames.PLOT_HEIGHT)/ 5 + 1))
        )
        # set the aesthetics
        for which in 'major', 'minor':
            getattr(self, PDNames.PLOT_AXES).grid(
                which=which,
                linestyle='-',
                linewidth=self.grid_linewidth[which],
                color=self.grid_color[which]
            )
            getattr(self, PDNames.PLOT_AXES).tick_params(
                which=which,
                color=self.grid_color[which],
                bottom=False,
                top=False,
                left=False,
                right=False
            )
    # /////////////////////////////////////////////////////////////////////////
    def _draw_signal(self, layout:list[list[str]],
                     figsize:tuple[float, float] | None=None,
                     wave_type:Literal['rhythm', 'media'] | None=None,
                     start_pos:Literal['first', 'continues']= 'first',
                     kwargs_signal:dict[Any, Any] | None=None,
                     kwargs_text:dict[Any, Any] | None=None,
                     ) -> None:
        '''
        Draws ECG signals on a matplotlib axis.
        
        Parameters
        ----------
        layout : `list` [`list` [`str`]]
            A nested lists organising the image by ECG lead names.
        figsize : `tuple` (`float`, `float`), default `NoneType`
            The width by height figure size in inches. If `NoneType` it will
            use the self.paper_h and self_paper_w attributes which are
            defined in millimetres divided by self.inch_mm.
        start_pos: {`first`, `continues`}, default `first`
            Whether multiple sequential leads should start at the first
            time point, or whether the next lead should start at the time
            the previous lead stopped.
        wave_type : {'rhythm', 'median'}, default `rhythm`
            The type of signal to plot, will update the `wave_type` attribute.
        kwargs_*_dict : `dict` [`any`, `any`], default `NoneType`
            Optional arguments supplied to the various plotting functions:
                kwargs_signal        --> ax.plot
                kwargs_text          --> ax.text
        '''
        # #### check input
        is_type(start_pos, str)
        is_type(wave_type, (type(None),str))
        is_type(layout, (type(None), list))
        is_type(figsize, (type(None), tuple))
        setattr(self, PDNames.PLOT_LAYOUT, layout)
        # set wave type
        if not wave_type is None:
            setattr(self, PDNames.WAVE_TYPE, wave_type)
            WAVE_TYPE = [PDNames.WAVETYPE_RHYTHM,
                        PDNames.WAVETYPE_MEDIAN,
                        ]
            if not getattr(self, PDNames.WAVE_TYPE) in WAVE_TYPE:
                raise ValueError(Error_MSG.CHOICE_PARM.\
                                 format('wave_type', ', '.join(WAVE_TYPE)))
        # confirm start pos
        START_POS = [PDNames.START_POS_CONT, PDNames.START_POS_FIRST]
        if not start_pos in START_POS:
                raise ValueError(Error_MSG.CHOICE_PARM.\
                                 format('start_pos', ', '.join(START_POS)))
        # check an axis is available
        if hasattr(self, PDNames.PLOT_AXES) == False:
            raise NotCalledError('The `_set_canvas` method has not been '
                                 'Called. There are no axes to draw on.')
        # confirm layout matches the signal.keys
        flat_layout = [l for s in getattr(self, PDNames.PLOT_LAYOUT) for l in s]
        msg_elements = [l for l in flat_layout if l not in\
                        getattr(self, PDNames.ECG_SIGNAL)]
        if msg_elements:
            raise KeyError(f"The following `layout` entries are unavailable "
                           f"in the waveform dictionary: {msg_elements}. "
                           f"The dictionary keys are: "
                           f"{list(getattr(self, PDNames.ECG_SIGNAL))}.")
        # set some default
        if figsize is None:
            figsize=(self.paper_w/self.inch_mm, self.paper_h/self.inch_mm)
        setattr(self, PDNames.PLOT_FIGSIZE, figsize)
        # map None to dict
        kwargs_signal, kwargs_text = assign_empty_default(
            [kwargs_signal, kwargs_text], dict)
        # #### sort out the ploting areas (note not using gridspecs currently)
        rows = len(getattr(self, PDNames.PLOT_LAYOUT))
        for numrow, row in enumerate(getattr(self, PDNames.PLOT_LAYOUT)):
            columns = len(row)
            row_height = getattr(self, PDNames.PLOT_HEIGHT) / rows
            # Horizontal shift for lead labels and separators
            h_delta = getattr(self, PDNames.PLOT_SAMPLING_NUMBER)/ columns
            # Vertical shift of the origin
            v_delta = round(
                getattr(self, PDNames.PLOT_HEIGHT) * (1.0 - 1.0 / (rows * 2)) -
                numrow * (getattr(self, PDNames.PLOT_HEIGHT) / rows)
            )
            # Let's shift the origin on a multiple of 5 mm
            v_delta = (v_delta + 2.5) - (v_delta + 2.5) % 5
            # Lenght of a signal chunk
            chunk_size =\
                int(getattr(self, PDNames.PLOT_SAMPLING_NUMBER) / len(row))
            # ### actually plot the signals
            for numcol, k in enumerate(row):
                left = numcol * chunk_size
                right = (1 + numcol) * chunk_size
                sign_start = left
                sign_stop = right
                if start_pos == PDNames.START_POS_FIRST:
                    sign_start = 0
                    sign_stop = chunk_size
                # get the signal
                signal_temp =\
                    getattr(self, PDNames.ECG_SIGNAL)[k][sign_start:sign_stop]
                # scaled by mm/mV factor
                signal = v_delta + getattr(self, PDNames.WAVE_SCALING) *\
                    signal_temp
                # update kwargs
                new_kwargs_signal = _update_kwargs(
                    update_dict=kwargs_signal,
                    clip_on=False,
                    linewidth=0.6,
                    color='black',
                    zorder=2,
                )
                # plot the signal
                getattr(self, PDNames.PLOT_AXES).plot(
                    list(range(left, right)),
                    signal,
                    **new_kwargs_signal,
                    )
                # update kwargs
                new_kwargs_text = _update_kwargs( update_dict=kwargs_text,
                    zorder=3,
                    fontsize=8,
                )
                # plot the lead name
                h = h_delta * numcol
                # v = v_delta + row_height / 2.6
                getattr(self, PDNames.PLOT_AXES).text(
                    x=h + self.text_pad_x,
                    y=v_delta + row_height / 3,
                    s=k,
                    **new_kwargs_text,
                )
            # remove the tick labels
            getattr(self, PDNames.PLOT_AXES).set_xticklabels([])
            getattr(self, PDNames.PLOT_AXES).set_yticklabels([])
            # resize figure
            getattr(self, PDNames.PLOT_FIG).set_size_inches(figsize)
    # /////////////////////////////////////////////////////////////////////////
    def to_numpy(self, crop:bool=False, close:bool=True) -> np.ndarray:
        '''
        Maps a matplotlib image to a numpy array.
        
        Parameters
        ----------
        crop : `bool`, default `False`
            Whether the image should be cropped to focus on the data within
            the axes spines. False simply maps the entire figure object to a
            numpy array.
        close : `bool`, default `True`
            Whether to call `plt.close` after extracting the numpy array data.
        
        Returns
        -------
        array : np.ndarray
            3-dimensional array:
                x: the pixels along the vertical axis.
                y: the pixels along the horizontal axis.
                z: the number of channels (4 for an RGBA image).
        
        Notes
        -----
        Crop does not really crop the image it simply increases the axes to
        cover the entire figure and resizes the figure to the original axes
        size.
        '''
        # #### check input and ensure a figure is present
        is_type(crop, bool)
        is_type(close, bool)
        # check a figure is present (should have a fig and unit attribute)
        if hasattr(self, PDNames.PLOT_FIG) == False or\
                hasattr(self, PDNames.PLOT_UNIT) == False:
            raise NotCalledError()
        # ##### Does the figure need to be cropped
        if crop == True:
            # remove the entire axis
            plt.axis('off')
            # increase the axes size to cover the entire figure
            getattr(self, PDNames.PLOT_AXES).set_position([0, 0, 1, 1])
            # Crop by simply calling the height and width from self
            getattr(self, PDNames.PLOT_FIG).set_size_inches(
                    getattr(self, PDNames.PLOT_WIDTH)/self.inch_mm,
                    getattr(self, PDNames.PLOT_HEIGHT)/self.inch_mm,
                )
        # #### the actual mapping to numpy
        # render the figure buffer
        getattr(self, PDNames.PLOT_FIG).canvas.draw()
        # The buffer is in RGBA format, so we need to use the correct
        # shape. Note np.uint8 is sufficient for RBGA data and
        # stores each pixel channel (without loss) in 1 byte of memory.
        # Using a different (larger) data type would be inefficient and
        # require the data to be scaled to an [0.0, 1.0] interval.
        img_array = np.frombuffer(getattr(self, PDNames.PLOT_FIG).\
                                  canvas.tostring_argb(),
                                  dtype=np.uint8)
        # reshaping it to a 3D array with 4 RGBA channels
        img_array = img_array.reshape(
            getattr(self, PDNames.PLOT_FIG).\
            canvas.get_width_height()[::-1] + (4,)
        )
        # reorder the channels from ARGB to RGBA
        img_array = np.roll(img_array, 3, axis=2)
        # #### close figure
        if close == True:
            plt.close(getattr(self, PDNames.PLOT_FIG))
        # #### return
        return img_array

