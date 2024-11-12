'''
A module to extract median beats and raw waveforms from DICOM ECGs wrapping
pydicom.
'''

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# imports
import os
import re
import sys
import copy
import uuid
import pathlib
import warnings
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from pydicom import dcmread
from pydicom.dataset import FileDataset as DCM_Class
from datetime import datetime
from scipy import signal
from dataclasses import dataclass
from typing import (
    List, Type, Union, Tuple, Self, Dict, Optional, Any, Literal,
)
from ecgprocess.errors import (
    NotCalledError,
    is_type,
    MissingDICOMTagError,
    Error_MSG,
    STDOUT_MSG,
    _check_readable,
    _check_presence,
)
from ecgprocess.constants import (
    ProcessDicomNames as PDNames,
    DICOMTags,
)
from ecgprocess.utils.general import (
    replace_with_tar,
    assign_empty_default,
)
from ecgprocess.plot_ecgs import (
    ECGDrawing,
)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Base class with slots
class BaseECGDICOMReader:
    '''
    An ECGDICOMReader base class implementing the more efficient __slots__ for
    the waveform arrays, while still retaining __dict__ dynamic attribute
    creation.
    '''
    __slots__ = (PDNames.LEAD_VOLTAGES, PDNames.LEAD_VOLTAGES2,
                 PDNames.RESULTS_DICT)
    # /////////////////////////////////////////////////////////////////////////
    def __init__(self, lead_voltages=None, lead_voltages2=None,
                 results_dict=None):
        '''
        Initialises slots entries to `None`.
        '''
        setattr(self, PDNames.LEAD_VOLTAGES, lead_voltages)
        setattr(self, PDNames.LEAD_VOLTAGES2, lead_voltages2)
        setattr(self, PDNames.RESULTS_DICT, results_dict)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
@dataclass
class ECGDICOMReader(BaseECGDICOMReader):
    """
    Takes an ECG DICOM file and extracts metadata, median beats (if available)
    and raw waveforms.
    
    Parameters
    ----------
    augment_leads : bool, default `False`
        Whether the augmented leads are available in the DICOM, if not
        these are calculated.
    resample : bool, default `True`
        Whether to resample the ECG to a frequency of 500 Hertz.
    retain_raw : bool, default `False`
        Whether the raw pydicom instance and raw waveforms should be retained.
        Set to `False` to decrease memory usage. Set to `True` to explore the
        orignal pydicom instance. For example, use this one a few files to
        identify none-standard information to extract.
    
    Attributes
    ----------
    augment_leads : bool
        Whether the augmented leads were calculated if these were unavailable.
    resample : bool
        Whether the ECG was resampled to a 500 Hertz frequency.
    retain_raw : bool
        Whether the raw pydicom data was retained.
    METADATA : dict [`str`, `str`]
        A dictionary describing the metadata one wants to extract from a
        DICOM. The dictionary keys represents the `target` (new) name and the
        dictionary values the `source` (old) names.
    ECG_TRAIT_DICT : dict [`str`, `list[str]`]
        A dictionary with the keys reflecting the desired name of the ECG trait.
        Each key will have a list of strings as a value. These strings will be
        compared to the names in `WaveformAnnotationSequence` attribute.
        Matching is done without case-sensitivity.  If for any key there are
        multiple matching strings the algorithm will check if the extracted
        values are all the same, if not multiple entries will be returned for
        the user to decide what to do next. The extracted ECG traits will be
        included with the extracted METADATA.
    STUDY_DATE_SOURCE : str
        The study date format of the source DICOM.  This will be used in
        `datetime.strptime` as format argument.
    STUDY_DATE_TARGET : str
        The target formatting of the study date.  This will be used in
        `datetime.strptime` as format argument.
    ACQUISITION_DATE_SOURCE : str
        The acquisition date and time format of the source DICOM.  This will be
        used in `datetime.strptime` as format argument.
    ACQUISITION_DATE_TARGET : str
        The target formatting of the acquisition date and time.  This will be
        used in `datetime.strptime` as format argument.
    
    Methods
    -------
    get_metadata(path, skip_empty)
        Extract the dicom metadata.
    make_leadvoltages(waveform_array, lead_info, augment_leads)
        Extracts the voltages from a DICOM file. Will automatically extract the
        limb leads if missing.
    
    Notes
    -----
    The waveforms are stored using `__slots__` to decrease memory usage, which
    is further improved by setting `retain_raw` to `False` - the default
    behaviour.
    """
    # #### parameters, with defaults
    augment_leads: bool=False
    resample_500:bool=True
    retain_raw:bool=False
    # #### check input
    is_type(augment_leads, bool, 'augment_leads')
    is_type(resample_500, bool, 'resample_500')
    # #### default tags - hacking about to make these non-persistance
    # NOTE this is probably related to DICOMTags being a @dataclass
    METADATA = copy.deepcopy(DICOMTags().METADATA)
    WAVE_FORMS = copy.deepcopy(DICOMTags().WAVE_FORMS)
    WAVE_FORMS_DICT = copy.deepcopy(DICOMTags().WAVE_FORMS_DICT)
    MEDIAN_BEATS = copy.deepcopy(DICOMTags().MEDIAN_BEATS)
    MEDIAN_BEATS_DICT = copy.deepcopy(DICOMTags().MEDIAN_BEATS_DICT)
    ECG_TRAIT_DICT = copy.deepcopy(DICOMTags().ECG_INTERPERTATION_DICT)
    # #### default dates to extract and format
    STUDY_DATE_SOURCE = "%Y%m%d"
    STUDY_DATE_TARGET = "%Y-%m-%d"
    ACQUISITION_DATE_SOURCE = "%Y%m%d%H%M%S"
    ACQUISITION_DATE_TARGET = "%Y-%m-%d %H:%M:%S"
    # #### Error MSG
    __MSG1=('Please supply either `path` or `dicom_instance` but not both.')
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def __call__(self, path:str, skip_empty:bool=True, verbose:bool=False,
                 ) -> Self:
        """
        Read a `.dcm` DICOM file and extracts metadata, raw waveforms, and
        median beats.
        
        see `constants.DICOMTags` for the `METADATA`, `WAVE_FORMS`, and
        `MEADIAN_BEATS` tags looked for.
        
        Parameters
        ----------
        path : str
            The path to the .dcm file.
        skip_empty : bool, default `True`
            Whether empty tags should be skipped or throw an error.
        verbose : bool, default `False`
            Prints missing tags if skip_empty is set to `True`.
        
        Attributes
        ----------
        GeneralInfo : list [`str`]
            A list of dcmread extracted attributes.
        Waveforms : dict [`str`, `np.array`]
            The lead specific ECG waveforms.
        MedianWaveforms : dict [`str`, `np.array`]
            The lead specific ECG median beats.
        
        Returns
        -------
        self : `ECGDICOMReader` instance
            Returns the class instance with updated attributes extracted
            from `dcmread`.
        """
        is_type(path, (pathlib.PosixPath, str), 'path')
        is_type(skip_empty, bool, 'skip_empty')
        is_type(verbose, bool, 'verbose')
        # confirm file is readable
        _check_readable(path)
        # #### Read DICOM
        ECG, results_dict, empty_metadata = self.get_metadata(
            path, skip_empty=skip_empty)
        # #### Extract waveforms
        _, wave_dict, empty_wave_forms = self.get_waveforms(
            dicom_instance=ECG, skip_empty=skip_empty)
        results_dict.update(wave_dict)
        # #### Extract Median beats
        _, median_dict, empty_median_beats = self.get_median_beats(
            dicom_instance=ECG, skip_empty=skip_empty)
        results_dict.update(median_dict)
        # #### do we need to resample
        if self.resample_500 == True: 
            # the function will internally adjust
            # LEAD_VOLTAGES, and LEAD_VOLTAGES2
            # Nothing is returned
            self._resampling_500hz(frequency=results_dict[PDNames.SF])
            setattr(self, PDNames.RESAMPLED, True)
            results_dict[PDNames.RESAMPLED] = True
            results_dict[PDNames.SF] = 500
            results_dict[PDNames.SAMPLING_FREQ] = 500
        else:
            setattr(self, PDNames.RESAMPLED, False)
            results_dict[PDNames.RESAMPLED] = False
        # #### remove raw data
        if self.retain_raw == False:
            try:
                delattr(self, PDNames.WAVE_ARRAY)
                delattr(self, PDNames.MEDIAN_ARRAY)
            except AttributeError:
                pass
        ## Remove SF column
        #results_dict.drop[PDNames.SF]
        # #### Extract standard ECG measurements
        #_, ecg_traits, missing_ecg_traits = self._get_waveform_annotation(
        #    dicom_instance=ECG, skip_empty=skip_empty)
        #results_dict.update(ecg_traits)
        #empty_metadata = empty_metadata + missing_ecg_traits
        # #### end of extractions, optionally printing tags which were missing
        if verbose == True:
            if len(empty_metadata) + len(empty_wave_forms) +\
                    len(empty_median_beats)> 0:
                warnings.warn(
                    'The following DICOM tags could not be found: {}.'.format(
                        empty_metadata + empty_wave_forms +\
                    empty_median_beats))
        # #### add duration,
        # NOTE both values are either int/float or NA so no need
        # for a try except
        results_dict[PDNames.DURATION] = (
            results_dict[PDNames.SAMPLING_NUMBER] /
            results_dict[PDNames.SAMPLING_FREQ]
        )
        # #### add time variables
        raw_AD = results_dict[PDNames.ACQUISITION_DATE]
        raw_SD = results_dict[PDNames.STUDY_DATE]
        results_dict[PDNames.ACQUISITION_DATE] = \
            datetime.strptime(raw_AD, self.ACQUISITION_DATE_SOURCE).\
            strftime(self.ACQUISITION_DATE_TARGET)
        results_dict[PDNames.STUDY_DATE] = \
            datetime.strptime(raw_SD, self.STUDY_DATE_SOURCE).\
            strftime(self.STUDY_DATE_TARGET)
        # #### assign results_dict
        setattr(self, PDNames.RESULTS_DICT, results_dict)
        # assign unique identifier
        setattr(self, PDNames.SOP_UID, results_dict[PDNames.SOP_UID])
        # add the original dcmread instance
        if self.retain_raw == True:
            setattr(self, PDNames.ORIG_DCMREAD_INST, ECG)
        # #### return stuff
        return self
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def get_metadata(self, path:str|None=None, dicom_instance: DCM_Class|None=None,
                     skip_empty:bool=True
                     ) -> tuple[DCM_Class, dict[str, Any], list[str]]:
        '''
        Takes a dicom file and extracts its metedata
        
        Parameters
        ----------
        path : str, default `NoneType`.
            The path to the .dcm file.
        dicom_instance : DCM_Class, default `NoneType`.
            A DCM_Class instance.
        
        Returns
        -------
        results : dict,DCM_Class
            - A `DCM_Class` instance.
            - A dictionary with extracted metadata.
            - A list of missing `DCM_Class` attribute names.
        
        Notes
        -----
        Either supply a path to a dicom file or a DCM_Class instance
        '''
        # #### check input
        is_type(path, (type(None), pathlib.PosixPath, str))
        is_type(dicom_instance, (type(None), DCM_Class))
        is_type(skip_empty, bool)
        results_dict = {}
        if (not dicom_instance is None) and (not path is None):
            raise ValueError(self.__MSG1)
        # #### Read DICOM
        # NOTE `with` closes automatically if an error is raised
        if not path is None:
            with open(path, 'rb') as dicom:
                # reads standard dicom content
                ECG=dcmread(dicom)
        else:
            ECG=dicom_instance
        # #### extract metadata
        empty_metadata = []
        for t, s in self.METADATA.items():
            if hasattr(ECG, s):
                # Assign if present
                results_dict[t] = getattr(ECG, s)
            elif skip_empty == False:
                # Should an Error be returned
                raise MissingDICOMTagError(s)
            else:
                # assign NA and append missing metadata
                results_dict[t] = np.nan
                empty_metadata.append(s)
        # return
        return ECG, results_dict, empty_metadata
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def get_waveforms(self, path:str|None=None,
                     dicom_instance: DCM_Class|None=None,
                     skip_empty:bool=True
                     ) -> tuple[DCM_Class, dict[str, Any], list[str]]:
        '''
        Takes a dicom file and extracts the waveforms and waveform metadata.
        
        Parameters
        ----------
        path : str, default `NoneType`.
            The path to the .dcm file.
        dicom_instance : DCM_Class, default `NoneType`.
            A DCM_Class instance.
        
        Returns
        -------
        results : dict,DCM_Class
            - A `DCM_Class` instance.
            - A dictionary with extracted metadata.
            - A list of missing `DCM_Class` attribute names.
        
        Notes
        -----
        Either supply a path to a dicom file or a DCM_Class instance
        '''
        # #### check input
        is_type(path, (type(None), pathlib.PosixPath, str))
        is_type(dicom_instance, (type(None), DCM_Class))
        is_type(skip_empty, bool)
        temp_results_dict = {}
        if (not dicom_instance is None) and (not path is None):
            raise ValueError(self.__MSG1)
        # #### Read DICOM
        # NOTE `with` closes automatically if an error is raised
        if not path is None:
            with open(path, 'rb') as dicom:
                # reads standard dicom content
                ECG=dcmread(dicom)
        else:
            ECG=dicom_instance
        # #### extract waveforms
        try:
            nam, el = self.WAVE_FORMS_DICT[PDNames.WAVE_ARRAY]
            setattr(self, PDNames.WAVE_ARRAY,
                    getattr(ECG, nam)(el).T)
        except:
            raise AttributeError(Error_MSG.MISSING_DICOM.format(nam))
        WAVE = getattr(ECG, PDNames.WAVE_FORM_SEQ)[0]
        SETTINGS = getattr(WAVE, PDNames.CHANNEL_DEF_SEQ)[0]
        # loop over the remaining wave form data
        empty_wave_forms = []
        for t, s in self.WAVE_FORMS.items():
            # if present in WAVE or SETTINGS assign
            if hasattr(WAVE, s):
                temp_results_dict[t] = getattr(WAVE, s)
            elif hasattr(SETTINGS, s):
                temp_results_dict[t] = getattr(SETTINGS, s)
            elif skip_empty == False:
                # Should an Error be returned
                raise MissingDICOMTagError(s)
            else:
                # assign NA and append missing metadata
                temp_results_dict[t] = np.nan
                empty_wave_forms.append(s)
        # Check number of channels in ECG waveform
        if temp_results_dict[PDNames.CHANNEL_NUMBER] >= 8:
            channel_seq = getattr(getattr(ECG, PDNames.WAVE_FORM_SEQ)[0],
                                  PDNames.CHANNEL_DEF_SEQ)
            lead_info_waveform, lead_units=self._get_lead_info(channel_seq)
            setattr(self, PDNames.LEAD_VOLTAGES,
                    self.make_leadvoltages(
                        lead_info=lead_info_waveform,
                        waveform_array=getattr(self, PDNames.WAVE_ARRAY),
                        augment_leads=self.augment_leads,
                        ))
            temp_results_dict[PDNames.LEAD_UNITS] = lead_units
            temp_results_dict[PDNames.SAMPLING_FREQ] =\
                temp_results_dict[PDNames.SF_ORIGINAL]
        else:
            # set to NA
            setattr(self, PDNames.LEAD_VOLTAGES, np.nan)
            empty_wave_forms.append(PDNames.LEAD_VOLTAGES)
            temp_results_dict[PDNames.SAMPLING_FREQ] =\
                temp_results_dict[PDNames.SF_ORIGINAL]
            temp_results_dict[PDNames.LEAD_UNITS] = np.nan
        # return
        return ECG, temp_results_dict, empty_wave_forms
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def get_median_beats(self, path:str|None=None,
                     dicom_instance: DCM_Class|None=None,
                     skip_empty:bool=True
                     ) -> tuple[DCM_Class, dict[str, Any], list[str]]:
        '''
        Takes a dicom file and extracts the median beats and its metadata.
        
        Parameters
        ----------
        path : str, default `NoneType`.
            The path to the .dcm file.
        dicom_instance : DCM_Class, default `NoneType`.
            A DCM_Class instance.
        
        Returns
        -------
        results : dict,DCM_Class
            - A `DCM_Class` instance.
            - A dictionary with extracted metadata.
            - A list of missing `DCM_Class` attribute names.
        
        Notes
        -----
        Either supply a path to a dicom file or a DCM_Class instance
        '''
        # #### check input
        is_type(path, (type(None), pathlib.PosixPath, str))
        is_type(dicom_instance, (type(None), DCM_Class))
        is_type(skip_empty, bool)
        temp_results_dict = {}
        if (not dicom_instance is None) and (not path is None):
            raise ValueError(self.__MSG1)
        # #### Read DICOM
        # NOTE `with` closes automatically if an error is raised
        if not path is None:
            with open(path, 'rb') as dicom:
                # reads standard dicom content
                ECG=dcmread(dicom)
        else:
            ECG=dicom_instance
        # ##### extract the median beats data
        # check if the wave form is present and if the index is correct
        # the median beats should be index 1 (starting at 0)
        temp_results_dict = {k: np.nan for k in self.MEDIAN_BEATS}
        empty_median_beats = list(self.MEDIAN_BEATS.keys())
        sccss = True
        try:
            nam, el = self.MEDIAN_BEATS_DICT[PDNames.MEDIAN_ARRAY]
            setattr(self, PDNames.MEDIAN_ARRAY,
                    getattr(ECG, nam)(el).T)
        except:
            if skip_empty == True:
                sccss = False
                setattr(self, PDNames.MEDIAN_ARRAY, np.nan)
                pass
            else:
                raise AttributeError(Error_MSG.MISSING_DICOM.format(nam))
        # getting WaveformSequence and ChannelDefinitionSequence attributes
        if sccss == True:
            WAVE_M = getattr(ECG, PDNames.WAVE_FORM_SEQ)[1]
            SETTINGS_M = getattr(WAVE_M, PDNames.CHANNEL_DEF_SEQ)[1]
            for t, s in self.MEDIAN_BEATS.items():
                # if present in WAVE or SETTINGS assign
                if hasattr(WAVE_M, s):
                    temp_results_dict[t] = getattr(WAVE_M, s)
                    empty_median_beats.remove(t)
                elif hasattr(SETTINGS_M, s):
                    temp_results_dict[t] = getattr(SETTINGS_M, s)
                    empty_median_beats.remove(t)
                elif skip_empty == False:
                    # Should an Error be returned
                    raise MissingDICOMTagError(s)
                else:
                    # the dict and list have already been pre-populated
                    pass
            # Update the leads if needed
            channel_seq_median = getattr(
                getattr(ECG, PDNames.WAVE_FORM_SEQ)[1],
                PDNames.CHANNEL_DEF_SEQ)
            lead_info_median, lead_units2=self._get_lead_info(channel_seq_median)
            # assign the median beats
            setattr(self, PDNames.LEAD_VOLTAGES2,
                    self.make_leadvoltages(
                        lead_info=lead_info_median,
                        waveform_array=\
                        getattr(self, PDNames.MEDIAN_ARRAY),
                        augment_leads=self.augment_leads,
                    ))
            setattr(self, PDNames.MEDIAN_PRESENT, True)
            temp_results_dict[PDNames.MEDIAN_PRESENT] = True
            temp_results_dict[PDNames.LEAD_UNITS2] = lead_units2
        else:
            # set to missing
            setattr(self, PDNames.LEAD_VOLTAGES2, np.nan)
            setattr(self, PDNames.MEDIAN_PRESENT, False)
            temp_results_dict[PDNames.MEDIAN_PRESENT] = False
            temp_results_dict[PDNames.LEAD_UNITS2] = np.nan
        # return
        return ECG, temp_results_dict, empty_median_beats
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def _get_waveform_annotation(
        self, path:str|None=None, dicom_instance: DCM_Class|None=None,
        skip_empty:bool=True,) -> tuple[DCM_Class, dict[str, Any], list[str]]:
        '''
        Extract information from the `WaveformAnnotationSequence` attribute
        of an dicom file.
        
        Parameters
        ----------
        path : str, default `NoneType`.
            The path to the .dcm file.
        dicom_instance : DCM_Class, default `NoneType`.
            A DCM_Class instance.
        
        Returns
        -------
        results : dict,DCM_Class
            - A `DCM_Class` instance.
            - A dictionary with the extracted data.
            - A list of missing `DCM_Class` attribute names.
        
        Notes
        -----
        Either supply a path to a dicom file or a DCM_Class instance
        '''
        # #### check input
        is_type(path, (type(None), pathlib.PosixPath, str))
        is_type(dicom_instance, (type(None), DCM_Class))
        is_type(skip_empty, bool)
        default_results_dict = {}
        temp_results_dict = {}
        temp_unit_dict = {}
        if (not dicom_instance is None) and (not path is None):
            raise ValueError(self.__MSG1)
        # #### Read DICOM
        # NOTE `with` closes automatically if an error is raised
        if not path is None:
            with open(path, 'rb') as dicom:
                # reads standard dicom content
                ECG=dcmread(dicom)
        else:
            ECG=dicom_instance
        # first set everything to NA
        for e in self.ECG_TRAIT_DICT:
           default_results_dict[e] = np.nan
           default_results_dict[e+PDNames.ECG_UNIT_STRING] = np.nan
           free_text = np.nan
           default_results_dict[PDNames.PACEMAKER_SPIKE]=np.nan
        # next see if we can extract some interpretations.
        if hasattr(ECG, PDNames.ECG_INTERPERTATION):
            free_text = ''
            # get the values lower-case strings (which are lists)
            ECG_CMPR_LWR = [it.lower() for sl\
                            in self.ECG_TRAIT_DICT.values() for it in sl]
            for w in getattr(ECG, PDNames.ECG_INTERPERTATION):
                if hasattr(w, PDNames.ECG_CONCEPTNAME):
                    ecg_int = getattr(w, PDNames.ECG_CONCEPTNAME)[0]
                    try:
                        ecg_trait = getattr(ecg_int, PDNames.CODE_MEANING)
                        # find matching element - using lower case again
                        if ecg_trait.lower() in ECG_CMPR_LWR:
                            # get the measurement and set to float
                            temp_results_dict[ecg_trait.lower()] = float(
                                        getattr(w, PDNames.ECG_TRAIT_VALUE)
                                    )
                            # get the unit
                            try:
                                temp_unit_dict[
                                    ecg_trait.lower()+\
                                    PDNames.ECG_UNIT_STRING]=getattr(
                                        getattr(w, PDNames.ECG_UNIT)[0],
                                        PDNames.CODE_MEANING)
                            except (AttributeError, IndexError):
                                pass
                    except AttributeError:
                        pass
                    # #### ancillary info
                    # see if there is a PaceMakerSpike
                    if ecg_trait.lower() == PDNames.PACEMAKER_SPIKE.lower():
                        try:
                            default_results_dict[PDNames.PACEMAKER_SPIKE]=\
                                getattr(w, PDNames.REFERENCED_POS)
                        except AttributeError:
                            pass
                    # get free text
                    if ecg_trait.lower() == PDNames.FREE_TEXT.lower():
                        try:
                            free_text = "%s\n%s" % (
                                free_text,
                                getattr(w, PDNames.FREE_TEXT)
                            )
                        except AttributeError:
                            pass
        # assing the final free_text object
        default_results_dict[PDNames.FREE_TEXT]=free_text
        # now assign temp_results_dict todefault_results_dict dealing with
        # keys with more than one matching string.
        if len(temp_results_dict) > 0:
            for k, idx in self.ECG_TRAIT_DICT.items():
                # Check if idx has more than one entry
                # extract all entries and make sure they are the
                # same  - NOTE using set comprehension to get the
                # unique elements
                unique_set = list({temp_results_dict[el.lower()] for el in\
                                   idx if el.lower() in temp_results_dict})
                if len(unique_set) == 1:
                    # if only one unique entry simply assign this to k
                    for e in idx:
                        try:
                           default_results_dict[k] =\
                                temp_results_dict[e.lower()]
                           default_results_dict[k+PDNames.ECG_UNIT_STRING] =\
                                temp_unit_dict[
                                    e.lower()+PDNames.ECG_UNIT_STRING]
                        except KeyError:
                            # NOTE the KeyError is expected behaviour,
                            # some of the `e`'s will not be in the
                            # temp_results_dict
                            pass
                elif len(unique_set) > 1:
                    # given that the results are not unique
                    # we will return all using the individual `idx` elements
                    # instead of `k` - NOTE the index without call to lower
                    # is intended.
                    for e in idx:
                        try:
                           default_results_dict[e] =\
                                temp_results_dict[e.lower()]
                           default_results_dict[e+PDNames.ECG_UNIT_STRING] =\
                                temp_unit_dict[
                                    e.lower()+PDNames.ECG_UNIT_STRING]
                        except KeyError:
                            pass
        # which ECG traits are still nan
        missing_ecg_traits =\
            [k for k,v in default_results_dict.items() if pd.isna(v)]
        if skip_empty == False and len(missing_ecg_traits) > 0:
            # Should an Error be returned
            raise ValueError('The following ECG measurments are '
                             'unavailable: {}'.format(missing_ecg_traits))
        # return
        return ECG, default_results_dict, missing_ecg_traits
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def _get_lead_info(self, channel_seq: DCM_Class
                       ) -> tuple[dict[int, str], str]:
        """
        Extracts the lead names from a `dcmread` instance.
        
        Parameters
        ----------
        channel_seq : DCM_Class
            A `pydicom.sequence.Sequence` instance.
        
        Returns
        -------
        tuple of a dict and a string:
                
                - leadnames : dict [`int`, `str`]
                    A dictionary with numerical (interval) keys and the lead names as
                    values.
                - leadunit : str
                    The measurment unit of the ECG leads.
        """
        # #### extracting the lead names
        leadnames = {}
        leadunits = {}
        for k, channel in enumerate(channel_seq):
            # extracting lead names
            source = getattr(getattr(channel,PDNames.CHANNEL_SOURCE_SEQ)[0],
                             PDNames.CHANNEL_CODE_MEANING)
            # extracting units
            if hasattr(channel, PDNames.CHANNEL_SENS_UNIT):
                unit = getattr(
                    getattr(channel, PDNames.CHANNEL_SENS_UNIT)[0],
                    PDNames.CODE_MEANING)
            else:
                unit = np.nan
            # assign lead names to numericals
            leadnames[k] = source.replace(PDNames.LEAD, '').strip()
            leadunits[leadnames[k]] = unit
        # confirm the units are all the same
        unique_unit = list(set(leadunits.values()))
        if len(unique_unit) != 1:
            raise ValueError('The ECG leads were measured using different '
                             'units: `{}`.'.format(unique_unit))
        # return stuff
        return leadnames, unique_unit[0]
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def make_leadvoltages(self, waveform_array: np.ndarray,
                          lead_info:Dict[int, str],
                          augment_leads:bool,
                          ) -> Dict[str, np.ndarray]:
        """
        Extracts the voltages from a DICOM file. Will automatically extract the
        limb leads if missing, please see:
        `url <https://ecgwaves.com/topic/ekg-ecg-leads-electrodes-systems-limb-chest-precordial/>`_.
        
        Parameters
        ----------
        waveform_array : np.ndarray
            An array with the waveforms.
        lead_info: dict [`int`, `str`]
            A dictionary where the keys match the `waveform_array` rows and
            the dictionary values contain the lead name strings.
        augment_leads : bool
            Whether to calculate the additional augmented leads. Will only work
            if there are exactly 8 leads available.
        
        Returns
        -------
        leads: dict [`str`, `np.ndarray`]
            A dictionary with the lead name string as keys and leads as
            np.ndarray values.
        """
        # #### check input and set constants
        is_type(waveform_array, np.ndarray, 'waveform_array')
        is_type(lead_info, dict, 'lead_info')
        leads = {}
        # #### the main function
        for k, lead in enumerate(waveform_array):
            leads[lead_info[k]] = lead
        if (k+1) == 8 and augment_leads == True:
            # Calculate limb leads
            leads[PDNames.LEAD_III] = np.subtract(leads[PDNames.LEAD_II],
                                                  leads[PDNames.LEAD_I])
            leads[PDNames.LEAD_aVR] = np.add(leads[PDNames.LEAD_I],
                                  leads[PDNames.LEAD_II]) * (-0.5)
            leads[PDNames.LEAD_aVL] = np.subtract(leads[PDNames.LEAD_I],
                                       0.5 * leads[PDNames.LEAD_II])
            leads[PDNames.LEAD_aVF] = np.subtract(leads[PDNames.LEAD_II],
                                       0.5 * leads[PDNames.LEAD_I])
        # return
        return leads
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def _resampling_500hz(self, frequency:int|float):
        """
        Re-sample the frequency to 500 hz.
        
        Parameters
        ----------
        frequency : int or float
            The original sampling frequency. Will simply check that this is
            not 500 and if not proceed with the resampling.
        
        """
        if int(frequency) != 500:
            # #### wave forms
            lead_volt_temp = getattr(self, PDNames.LEAD_VOLTAGES)
            for i in lead_volt_temp:
                lead_volt_temp[f"{i}"] =\
                    signal.resample(lead_volt_temp[f"{i}"], 5000)
            setattr(self, PDNames.LEAD_VOLTAGES, lead_volt_temp)
            # #### median beats
            if getattr(self, PDNames.MEDIAN_PRESENT) == True:
                lead_volt_temp2 = getattr(self, PDNames.LEAD_VOLTAGES2)
                for i in lead_volt_temp2:
                    lead_volt_temp2[f"{i}"] = \
                        signal.resample(lead_volt_temp2[f"{i}"], 600)
                setattr(self, PDNames.LEAD_VOLTAGES2, lead_volt_temp2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ECGDICOMTable(object):
    '''
    Takes an `ECGDICOMReader` instance and loops over a list of dicom paths
    and maps these to which can be used in analyses or saved to disk.
    
    Attributes
    ----------
    RawPathList : list [`str`]
        A list of dcm file paths.
    info_type : str
        The type of information one wants to extract from the dicom files.
    
    Methods
    -------
    get_table(update_keys,**kwargs)
        extracts data from multiple dicom files and maps these to
        pandas.DataFrames.
    write_ecg(target_tar, target_path, sep, mode, compression, update_keys,
    **kwargs)
        writes each dicom file to a single sets of files (metadata,
        waveforms, median beats). This is done by appending the extracted
        data from each file the target file set, minimising the memory
        footprint.
    write_pdf(ecgdrawing, target_path, write_failed, kwargs_reader,
    kwargs_drawing, kwargs_savefig)
        writes dicom files to pdfs using the dicom unique id as file name.
    '''
    
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def __init__(self, ecgdicomreader:ECGDICOMReader, path_list:List[str],
                 info_type:Literal['all', 'rhythm', 'median', 'meta'] = 'all'
                 ) -> None:
        """
        Initialises a new instance of `ECGDICOMTable`.
        
        Parameters
        ----------
        ecgdicomreader : ECGDICOMReader
            An instance of the ECGDICOMReader data class.
        path_list : list [`str`]
            A list of paths to one or more .dcm files.
        info_type : {`all`, `rhythm`, `median`, `meta`}
            Which information should be extracted.
        """
        EXP_INFO=[PDNames.INFO_TYPE_ALL, PDNames.INFO_TYPE_RTM,
                  PDNames.INFO_TYPE_MED, PDNames.INFO_TYPE_MET,
                  ]
        self.INFO_MET = [PDNames.INFO_TYPE_ALL, PDNames.INFO_TYPE_MET]
        self.INFO_MED = [PDNames.INFO_TYPE_ALL, PDNames.INFO_TYPE_MED]
        self.INFO_RTM = [PDNames.INFO_TYPE_ALL, PDNames.INFO_TYPE_RTM]
        # #### check input
        is_type(ecgdicomreader, ECGDICOMReader, 'ecgdicomreader')
        is_type(path_list, list, 'path_list')
        if not info_type in EXP_INFO:
            raise ValueError(f'`info_type` is restricted to `{EXP_INFO}`.')
        self.ecgdicomreader = ecgdicomreader
        setattr(self, PDNames.INFO_TYPE, info_type)
        setattr(self, PDNames.RPATH_L, path_list)
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def __str__(self):
        CLASS_NAME = type(self).__name__
        return (f"{CLASS_NAME} instance with "
                f"ecgdicomreader={self.ecgdicomreader}, "
                f"path_list={getattr(self, PDNames.RPATH_L)}."
                )
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def __repr__(self):
        CLASS_NAME = type(self).__name__
        return (f"{CLASS_NAME}(ecgdicomreader={self.ecgdicomreader}, "
                f"augment_leads={getattr(self, PDNames.RPATH_L)})"
                )
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    def __call__(
        self, skip_missing:Literal['Permissions', 'Data', 'None']='Permissions',
        verbose:bool=False,
    ) -> Self:
        """
        Will take a ECGDICOMReader and loops over a list of dcm file paths and
        and confirm the files exist and have appropriate read permission
        
        Parameters
        ----------
        skip_missing : {'Permissions', 'Data', 'None'}, default `Permissions`
            Whether files with a `PermissionError` or `AttributeError` should
            be skipped, or whether all errors should be raised.
            A PermissionError implies either the file is absent or it cannot
            be read. An AttributeError indicates the waveform_array
            attribute is absent from the pydicom instance.
            `Permissions` will only skip PermissionErrors, while `Data` will
            skip both PermssionErrors and AttributeErrors (due to a missing
            waveform_array), `None` will not skip errors.
        verbose : bool, default `False`
            Prints missing files if skip_missing is set to `True`.
        
        Attributes
        ----------
        FailedPathList : list [`str`]
            File paths which were either absent or without read permission.
        CuratedPathList : list [`str`]
            File paths which are readable.
        
        Returns
        -------
        self : `ECGDICOMTable` instance
            Returns the class instance with updated attributes.
        """
        is_type(skip_missing, str)
        is_type(verbose, bool)
        self.skip_missing = skip_missing
        self.verbose = verbose
        SKIP_MISSING = [PDNames.SKIP_PERMISSIONS,
                        PDNames.SKIP_DATA,
                        PDNames.SKIP_NONE,
                        ]
        if not skip_missing in SKIP_MISSING:
            raise ValueError(Error_MSG.CHOICE_PARM.\
                             format('skip_missing', ', '.join(SKIP_MISSING)))
        # #### loop over path
        empty_list = []
        curated_list = []
        # #### loop over individual dcm files and assign to self
        for p in getattr(self, PDNames.RPATH_L):
            try:
                # add p and remove if skipp_missing == True
                curated_list.append(p)
                _check_readable(p)
            except PermissionError as PE:
                empty_list.append(p)
                # try next
                if self.skip_missing != PDNames.SKIP_NONE:
                    curated_list.pop()
                    pass
                else:
                    raise PE
        setattr(self, PDNames.FPATH_L, empty_list)
        setattr(self, PDNames.CPATH_L, curated_list)
        # #### do we want to print empty_list
        if self.verbose == True:
            if len(getattr(self, PDNames.FPATH_L))>0:
                warnings.warn(
                    'The following files could not be accessed or found: {}.'.\
                    format(getattr(self, PDNames.FPATH_L))
                )
        # #### Return
        return self
    # /////////////////////////////////////////////////////////////////////////
    def get_table(self, update_keys:Optional[Dict[str, str]]=None,
                  **kwargs:Optional[Any],
                  ) -> Self:
        """
        Will extracted dicom data to tables.
        
        Parameters
        ----------
        update_keys: dict [`str`, `str`], default `NoneType`
            A dictionary to remap lead names: [`old`, `new`]
        **kwargs: optional
            Keyword arguments used in the call method of a `ECGDICOMReader`
            instance
        
        Attributes
        ----------
        IndexList: list [`str`]
            A list of unique identifiers.
        GeneralInfoTable: pandas.DataFrame
            A table of the dicom metadata.
        WaveFormsTable: pandas.DataFrame
            A  long-format table with waveforms.
        MedianWaveTable: pandas.DataFrame
            A long-format table with the median beat waveforms.
        NoDataList : `list` [`str`]
            A list of dicom files without a waveform_array.
        
        Returns
        -------
        self : `ECGDICOMTable` instance
            Returns the class instance with updated attributes.
        """
        self.kwargs = kwargs
        # #### check if __call__ has been run
        if not hasattr(self, PDNames.CPATH_L):
            raise NotCalledError()
        # #### extract dicom data
        no_data_list, key_list, info_list, wave_list, median_list =\
            [[] for _ in range(5)]
        # loop over individual dcm files
        for p in getattr(self, PDNames.CPATH_L):
            if self.verbose == True:
                print(STDOUT_MSG.PROCESSING_PATH.format(p), file=sys.stdout)
            # get instance
            try:
                ecg_inst = self.ecgdicomreader(p, **self.kwargs)
            except AttributeError as AE:
                if self.skip_missing == PDNames.SKIP_DATA:
                    no_data_list.append(p)
                    # moving to the next path
                    continue
                else:
                    raise AE
            # extract unique identifier and check if it has been used before
            if hasattr(ecg_inst, PDNames.SOP_UID) == False:
                raise AttributeError(Error_MSG.MISSING_ATTR.format(
                    PDNames.SOP_UID, 'ecg_inst'))
            key = str(getattr(ecg_inst, PDNames.SOP_UID))
            if key in key_list:
                raise IndexError('{0}:{1} was already extracted before. Please '
                                 'ensure the supplied files are unique.'.\
                                 format(PDNames.SOP_UID, key))
            else:
                key_list.append(key)
            # extract the remaining
            if getattr(self, PDNames.INFO_TYPE) in self.INFO_MET:
                info_list.append(getattr(ecg_inst, PDNames.RESULTS_DICT))
            if getattr(self, PDNames.INFO_TYPE) in self.INFO_RTM:
                wave_list.append(getattr(ecg_inst, PDNames.LEAD_VOLTAGES))
            if getattr(self, PDNames.INFO_TYPE) in self.INFO_MED:
                median_list.append(getattr(ecg_inst, PDNames.LEAD_VOLTAGES2))
        # #### make tables
        setattr(self, PDNames.FAILED_DATA_L, no_data_list)
        setattr(self, PDNames.KEY_L, key_list)
        # general info
        if getattr(self, PDNames.INFO_TYPE) in self.INFO_MET:
            setattr(self, PDNames.INFO_T, pd.DataFrame(
                info_list, index=getattr(self, PDNames.KEY_L)
            ))
        # wave forms, and median beats
        if getattr(self, PDNames.INFO_TYPE) in self.INFO_RTM:
            setattr(self, PDNames.WAVE_T, self._get_long_table(
                wave_list, wave_type=PDNames.WAVETYPE_RHYTHM,
                update_keys=update_keys,
                purge_header=True,
            ))
        if getattr(self, PDNames.INFO_TYPE) in self.INFO_MED:
            setattr(self, PDNames.MEDIAN_T, self._get_long_table(
                median_list, wave_type=PDNames.WAVETYPE_MEDIAN,
                update_keys=update_keys,
                purge_header=True,
            ))
        # #### Return
        setattr(self, PDNames.TABLE_CALLED, 'True')
        return self
    # /////////////////////////////////////////////////////////////////////////
    def _get_long_table(self, lead_list:List[Dict[str, np.ndarray]],
                        wave_type:str,
                        update_keys:Optional[Dict[str,str]]=None,
                        purge_header:bool=True,
                        **kwargs,
                        ) -> pd.DataFrame:
        '''
        Mapping lists of dictionaries containing lead specific numpy arrays
        to a long-formatted pandas table.
        
        Parameters
        ----------
        lead_list : list [`dict`]
            A list of dictionaries with the lead mapped to the keys and the
            voltage mapped to the values.
        wave_type : str,
            Adds a column `Waveform type` containing this string as a constant.
        update_keys : dict [`str`, `str`], default `NoneType`
            A dictionary to remap lead names: [`old`, `new`]
        purge_header : bool, default True
            Set to `False` to make sure the file header persists between calls.
            This is used to ensure the headers are the same between files.
        **kwargs : optional
            keyword arguments passed to pd.DataFrame.
        
        Attributes
        ----------
        kwargs_pdd : optional
            The keyword arguments supplied to pd.DataFrame.
        
        Returns
        -------
        self : pd.DataFrame
            A long-formatted table with lead, voltage, sampling indicator
            columns grouped by file (including an unique file indicator column)
        '''
        # #### check input and set constants
        is_type(lead_list, list, 'lead_list')
        is_type(wave_type, str, 'wave_type')
        self.kwargs_pdd = kwargs
        PREV = 'previous'
        # do we need to purge the header information
        # otherwise this persists throughout calls
        if purge_header == True:
            setattr(self, PREV, None)
        else:
            if not hasattr(self, PREV):
                setattr(self, PREV, None)
        # #### initiate table and map lists
        table = pd.DataFrame()
        for w, k in zip(lead_list, getattr(self, PDNames.KEY_L), strict=True):
            # do we need to remap key names
            if update_keys is not None:
                w = {update_keys.get(k, k): v for k, v in w.items()}
            # confirm the dictionary keys are identical
            current_keys = list(w.keys())
            if getattr(self, PREV) is not None:
                if current_keys != getattr(self, PREV):
                    raise KeyError('The dictionaries contain distinct keys. '
                                   'The last set of valid keys was {}, '
                                   'compared to {}.'.\
                                   format(getattr(self, PREV), current_keys))
            # if the same map to dataframe
            temp_df = pd.DataFrame(w, **self.kwargs_pdd)
            temp_df[PDNames.SOP_UID] = k
            table = pd.concat([table, temp_df], ignore_index=True)
            # clean
            setattr(self, PREV, current_keys)
            del current_keys
        # #### formatting to long
        table[PDNames.SAMPLING_SEQ] = table.groupby(
            PDNames.SOP_UID).cumcount()
        long_table = pd.melt(table, id_vars=\
                             [PDNames.SOP_UID, PDNames.SAMPLING_SEQ],
                             value_vars=getattr(self, PREV),
                             var_name=PDNames.COL_LEAD,
                             value_name=PDNames.COL_VOLTAGE
                             )
        long_table[PDNames.COL_WAVETYPE] = wave_type
        # #### return
        return long_table
    # /////////////////////////////////////////////////////////////////////////
    def write_ecg(self, target_tar:Union[None,str]=None, target_path:str='.',
                  table_prefix:str='',
                  sep:str='\t', mode:str='w:gz', compression:str='gzip',
                  update_keys:Optional[Dict[str,str]]=None,
                  write_failed:bool=True,
                  **kwargs:Optional[Any],
                  ) -> Self:
        '''
        Extracts dicom files, and write these one by one to a set of target
        files. By appending each dicom individually to the target files the
        memory footprint of this function will be minimal.
        
        Parameters
        ----------
        target_tar : str, default `NoneType`
            The `name` of an optional tarfile where the individual files will
            be written to. The target_tar will be concatenated to
            `target_path` and depending on `mode` this directory will be
            tar.gz compressed. Set `target_tar` to `NoneType` to simply add
            the files directly to `target_path`. Note this will overwrite
            any potential directory or files with the identical names.
        target_path : str, default '.'
            The full path where the files should be written to. if provided
            `target_tar` will be created underneath this path, otherwise the
            files will be directly written to the `target_path` terminal
            directory (assuming this is writable).
        table_prefix : str, default ''
            Prefix for filenames of all tables, e.g. to add context information 
        sep : str, default '\t'`
            The file separator, which will be passed to
            pandas.DataFrame.to_csv.
        mode : str, default `w:gz`
            The tarfile.open mode.
        compression : str, default `gzip`
            The file compression passed to pandas.DataFrame.write_csv.
        update_keys : dict [`str`, `str`], default `NoneType`
            A dictionary to remap lead names: [`old`, `new`]
        write_failed : bool, default `True`
            Whether to write a text file to disk containing the failed file
            names.
        **kwargs : Optional[Any],
            Keyword arguments used in the call method of a `ECGDICOMReader`
            instance.
        
        Attributes
        ----------
        target_path : str
            The directory or tar file path were the files are written to.
        NoDataList : `list` [`str`]
            A list of dicom files without a waveform_array.
        
        Returns
        -------
        self : `ECGDICOMTable` instance
            Returns the class instance with updated attributes.
        
        Notes
        -----
        This method writes the following files to disk:
        - `GeneralInfoTable.tsv`
        - `WaveFormsTable.tsv`
        - `MedianWaveTable.tsv`
        - `FailedFiles.txt`
        
        Raises
        ------
        NotADirectoryError or PermissionError
            If the target directory does not exist or is not writable.
        '''
        self.kwargs = kwargs
        # #### check input and set constants
        is_type(target_tar, (type(None), str), 'target_tar')
        is_type(target_path, (pathlib.PosixPath, str), 'target_path')
        is_type(sep, str, 'sep')
        is_type(table_prefix, str, 'table_prefix')
        is_type(mode, str, 'mode')
        is_type(compression, (type(None), str), 'compression')
        # check readability
        _check_presence(target_path)
        _check_readable(target_path)
        # #### create target path
        # get the current wd if requested
        if target_path == '.':
            target_path = os.getcwd()
        # create a temp directory under target_path or otherwise simply use
        # target_path
        if target_tar is not None:
            # adding a temp dir
            target = os.path.join(target_path, str(uuid.uuid4()))
            # make the new child dir
            # NOTE if target_tar has extensions such as `tar.gz` this will be
            # included but the directly is not really compressed (will be done
            # below).
            os.makedirs(target, exist_ok=True)
        else:
            target = target_path
            setattr(self, PDNames.WRITE_ECG_PATH, target)
        # #### check if __call__ has been run
        if not hasattr(self, PDNames.CPATH_L):
            raise NotCalledError()
        # #### extract dicom data
        first=True
        key_list, no_data_list = [[] for _ in range(2)]
        # loop over individual dcm files
        for p in getattr(self, PDNames.CPATH_L):
            no_data_list, key_list, ecg_inst = self._write_internal(
                path=p, no_data_list=no_data_list, key_list=key_list,
                **self.kwargs,
                )
            if isinstance(ecg_inst, str):
                continue
            # #### extract data from ecg_inst
            if getattr(self, PDNames.INFO_TYPE) in self.INFO_MET:
                info=getattr(ecg_inst, PDNames.RESULTS_DICT)
            if getattr(self, PDNames.INFO_TYPE) in self.INFO_RTM:
                wave=getattr(ecg_inst, PDNames.LEAD_VOLTAGES)
            if getattr(self, PDNames.INFO_TYPE) in self.INFO_MED:
                median=getattr(ecg_inst, PDNames.LEAD_VOLTAGES2)
            # assign key to self for use in `_get_long_table`
            setattr(self, PDNames.KEY_L, [key_list[-1]])
            # #### write to disk
            if first==True:
                first = False
                # metadata
                print([key_list[-1]])
                if getattr(self, PDNames.INFO_TYPE) in self.INFO_MET:
                    pd.DataFrame([info], index=[key_list[-1]]).to_csv(
                        os.path.join(target, table_prefix + PDNames.INFO_FILE), sep=sep,
                        header=True, index=False, compression=compression)
                # waveforms
                if getattr(self, PDNames.INFO_TYPE) in self.INFO_RTM:
                    self._get_long_table(
                        [wave], wave_type=PDNames.WAVETYPE_RHYTHM,
                        update_keys=update_keys,
                        purge_header=True,
                    ).to_csv(
                            os.path.join(target, table_prefix + PDNames.WAVE_FILE), sep=sep,
                            header=True, index=False, compression=compression)
                # median beats
                if getattr(self, PDNames.INFO_TYPE) in self.INFO_MED:
                    self._get_long_table(
                        [median], wave_type=PDNames.WAVETYPE_MEDIAN,
                        update_keys=update_keys,
                        purge_header=True,
                    ).to_csv(
                            os.path.join(target, table_prefix + PDNames.MEDIAN_FILE), sep=sep,
                            header=True, index=False, compression=compression)
            else:
                # appending using mode = 'a'
                # metadata
                if getattr(self, PDNames.INFO_TYPE) in self.INFO_MET:
                    pd.DataFrame([info], index=[key_list[-1]]).to_csv(
                        os.path.join(target, table_prefix + PDNames.INFO_FILE), sep=sep,
                        header=False, index=False, mode='a', compression=compression)
                # waveforms
                if getattr(self, PDNames.INFO_TYPE) in self.INFO_RTM:
                    self._get_long_table(
                        [wave], wave_type=PDNames.WAVETYPE_RHYTHM,
                        update_keys=update_keys,
                        purge_header=False,
                    ).to_csv(
                            os.path.join(target, table_prefix + PDNames.WAVE_FILE), sep=sep,
                            header=False, index=False, mode='a', compression=compression)
                # median beats
                if getattr(self, PDNames.INFO_TYPE) in self.INFO_MED:
                    self._get_long_table(
                        [median], wave_type=PDNames.WAVETYPE_MEDIAN,
                        update_keys=update_keys,
                        purge_header=False,
                    ).to_csv(
                            os.path.join(target, table_prefix + PDNames.MEDIAN_FILE), sep=sep,
                            header=False, index=False, mode='a', compression=compression)
            # delete key
            delattr(self, PDNames.KEY_L)
        # #### write failed files, note not compressing these
        DELIM = '\t'
        if write_failed == True:
            # adding the reason for failing
            total_failures = [
                (p, PDNames.SKIP_PERMISSIONS) for p in\
                getattr(self, PDNames.FPATH_L) ] + [
                (p, PDNames.SKIP_DATA) for p in\
                    getattr(self, PDNames.FAILED_DATA_L)]
            # writing to text file
            with open(os.path.join(target, table_prefix + PDNames.FAILED_FILE), 'w') as file:
                for p, cause in total_failures:
                    file.write(p + DELIM + cause + "\n")
        # #### if needed replace directory by tar.gz version
        if target_tar is not None:
            # create the final target path
            target_final = os.path.join(target_path, target_tar)
            replace_with_tar(target, target_final, mode=mode)
            setattr(self, PDNames.WRITE_ECG_PATH, target_final)
        # #### return
        return self
    # /////////////////////////////////////////////////////////////////////////
    def write_pdf(self, ecgdrawing:ECGDrawing,
                  target_path:str='.', write_failed:bool=True,
                  kwargs_reader:Dict[Any,Any] | None=None,
                  kwargs_drawing:Dict[Any,Any] | None=None,
                  kwargs_savefig:Dict[Any,Any] | None=None,
                  ) -> Self:
        '''
        Extracts dicom files, and write these one by one to a pdf files
        using a supplied `ECGDrawing` instance.
        
        Parameters
        ----------
        ecgdrawing : ECGDrawing
            An instance of the ECGDrawing data class.
        target_path : `str`, default '.'
            The full path where the pdfs should be written to.
        write_failed : `bool`, default `True`
            Whether to write a text file to disk containing the failed file
            names.
        kwargs_*: dict [`any`, `any`], default `NoneType`
            dictionaries with keyword arguments for the `plt.savefig`,
            ECGDrawing, or ECGDICOMReader instances.
        
        Attributes
        ----------
        target_path : `str`
            The directory or tar file path were the files are written to.
        NoDataList : `list` [`str`]
            A list of dicom files without a waveform_array.
        
        Returns
        -------
        self : `ECGDICOMTable` instance
            Returns the class instance with updated attributes.
        
        Notes
        -----
        The dicom UID instance will be used as file name for the pdfs.
        
        Raises
        ------
        NotADirectoryError or PermissionError
            If the target directory does not exist or is not writable.
        '''
        # #### check input and set constants
        is_type(kwargs_reader, (type(None), dict))
        is_type(kwargs_drawing, (type(None), dict))
        is_type(kwargs_savefig, (type(None), dict))
        is_type(ecgdrawing, ECGDrawing)
        is_type(target_path, (pathlib.PosixPath, str))
        is_type(write_failed, bool)
        # check readability
        _check_presence(target_path)
        _check_readable(target_path)
        # map None to dict
        kwargs_savefig, kwargs_reader, kwargs_drawing= assign_empty_default(
                [kwargs_savefig, kwargs_reader, kwargs_drawing], dict)
        # map info type to type of plot
        if getattr(self, PDNames.INFO_TYPE) == PDNames.INFO_TYPE_MET:
            raise AttributeError('`info` should be {0} or {1}.'.format(
                PDNames.INFO_TYPE_MED, PDNames.INFO_TYPE_RTM))
        elif getattr(self, PDNames.INFO_TYPE) == PDNames.INFO_TYPE_ALL:
            signal_type = PDNames.WAVETYPE_RHYTHM
        else:
            signal_type = getattr(self, PDNames.INFO_TYPE)
        # #### create target path
        if target_path == '.':
            target_path = os.getcwd()
        # assign to self
        target = target_path
        setattr(self, PDNames.WRITE_ECG_PATH, target)
        # #### check if __call__ has been run
        if not hasattr(self, PDNames.CPATH_L):
            raise NotCalledError()
        # #### extract dicom data
        key_list, no_data_list = [[] for _ in range(2)]
        # loop over individual dcm files
        for p in getattr(self, PDNames.CPATH_L):
            no_data_list, key_list, ecg_inst = self._write_internal(
                path=p, no_data_list=no_data_list, key_list=key_list,
                **kwargs_reader,
                )
            if isinstance(ecg_inst, str):
                continue
            # #### updated leads
            # if update_keys is not None:
            #     w = {update_keys.get(k, k): v for k, v in w.items()}
            # ##### draw and write figure
            artist = ecgdrawing(ecgreader=ecg_inst, wave_type=signal_type,
                                **kwargs_drawing)
            filename_pdf = re.sub(r"[ ,\-\(\)\{\}]", '_', key_list[-1]) + '.pdf'
            plt.savefig(fname=os.path.join(target, filename_pdf), **kwargs_savefig)
            plt.close(artist.fig)
        # #### write failed files, note not compressing these
        DELIM = '\t'
        if write_failed == True:
            # adding the reason for failing
            total_failures = [
                (p, PDNames.SKIP_PERMISSIONS) for p in\
                getattr(self, PDNames.FPATH_L) ] + [
                (p, PDNames.SKIP_DATA) for p in\
                    getattr(self, PDNames.FAILED_DATA_L)]
            # writing to text file
            with open(os.path.join(target, PDNames.FAILED_FILE), 'w') as file:
                for p, cause in total_failures:
                    file.write(p + DELIM + cause + "\n")
        # #### return
        return self
    # /////////////////////////////////////////////////////////////////////////
    def _write_internal(self, path:str, no_data_list:list[str],
                        key_list:list[str],
                        **kwargs,
                        ) -> tuple[list, list, ECGDICOMReader|str]:
        '''
        An internal function to read ECG dicom data, check whether there these
        data may have been extracted before (compared to key_list) and record
        whether some files did not contain waveform_array attributes.
        
        Parameters
        ----------
        no_data_list : list [`str`]
            A list of file names without an waveform_array attribute.
        key_list : list [`str`]
            A list of dicom UIDs which were processed before.
        
        Returns
        -------
        `tuple`
            A tuple with `no_data_list`, `key_list`, and an ECGDICOMReader
            instance.
        
        Raises
        ------
        AttributeError
            raised if waveform_array or SOPinstanceUID attributes are absent.
        IndexError
            raised if a dicom with the same SOPinstanceUID is processed
        '''
        
        if self.verbose == True:
            print(STDOUT_MSG.PROCESSING_PATH.format(path), file=sys.stdout)
        # get instance
        try:
            ecg_inst = self.ecgdicomreader(path, **kwargs)
        except AttributeError as AE:
            if self.skip_missing == PDNames.SKIP_DATA:
                no_data_list.append(path)
                # moving to the next path - cannot use continue here
                return no_data_list, key_list, 'continue'
            else:
                raise AE
        # extract unique identifier and check if it has been used before
        if hasattr(ecg_inst, PDNames.SOP_UID) == False:
            raise AttributeError(Error_MSG.MISSING_ATTR.format(
                PDNames.SOP_UID, 'ecg_inst'))
        key = str(getattr(ecg_inst, PDNames.SOP_UID))
        if key in key_list:
            raise IndexError('{0}:{1} was already extracted before. Please '
                             'ensure the supplied files are unique.'.\
                             format(PDNames.SOP_UID, key))
        else:
            key_list.append(key)
        # extract the remaining
        setattr(self, PDNames.FAILED_DATA_L, no_data_list)
        # return
        return no_data_list, key_list, ecg_inst



