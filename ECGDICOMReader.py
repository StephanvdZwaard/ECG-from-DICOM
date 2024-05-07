# Define class for reading ECGs from DICOM files.

import numpy as np;
import pandas as pd;
from pydicom import dcmread
from pydicom.waveforms import multiplex_array
from datetime import datetime
from scipy import signal

class ECGDICOMReader:
    """ Extract voltage data from a ECG DICOM file
        Authors: Stephan van der Zwaard, s.vanderzwaard@amsterdamumc.nl and Philip Croon, p.croon@amsterdamumc.nl
        for questions feel free to email.    
        
        -- First version based on code from Philip Croon, p.croon@amsterdamumc.nl
    """

    def __init__(self, augmentLeads=False, resample_500=True):
        """ 
        Initialize class. If resample_500 is True ECGs with sampling frequency that are not 500 will be resampled to 500. 
        If AugmentLeads = True and the augmented leads are not available, they are calculated. 
        """
        self.augmentLeads = augmentLeads
        self.resample_500 = resample_500
        print("Initialization succesfull")

    def __call__(self, path, verbose=False):
        
        try:
            with open(path, 'rb') as DICOM:
                
                # Open DICOM file
                self.ECG                   = dcmread(DICOM)
                
                # Add DICOM information by relevant tags
                # General information
                
                # Check if key DICOM tags are present
                if set(['SOPInstanceUID','SeriesInstanceUID','StudyInstanceUID','PatientID','AccessionNumber','StudyDate','StudyTime']).issubset(self.ECG.dir()) == True:
                    self.SOPinstanceUID        = self.ECG.SOPInstanceUID
                    self.SERIESinstanceUID     = self.ECG.SeriesInstanceUID
                    self.STUDYinstanceUID      = self.ECG.StudyInstanceUID
                    self.PatientID             = self.ECG.PatientID
                    self.StudyDate             = self.ECG.StudyDate
                    self.StudyTime             = self.ECG.StudyTime
                    self.AccessionNumber       = self.ECG.AccessionNumber
                else:
                    if (verbose==True): 
                        print('Essential DICOM tags are missing: inspect DICOM file')
                    return(pd.DataFrame({'filetype': [self.ECG.data_element('SOPClassUID').repval], 'error': 'Essential DICOM tags missing'}) )
                    
                # Fill general info with other DICOM tags (if present)                
                self.PatientBirthDate      = self.ECG.PatientBirthDate      if set(['PatientBirthDate']).issubset(self.ECG.dir()) else ''
                self.PatientName           = self.ECG.PatientName           if set(['PatientName']).issubset(self.ECG.dir()) else ''
                self.PatientSex            = self.ECG.PatientSex            if set(['PatientSex']).issubset(self.ECG.dir()) else ''
                self.StudyDescription      = self.ECG.StudyDescription      if set(['StudyDescription']).issubset(self.ECG.dir()) else ''
                self.AcquisitionDateTime   = self.ECG.AcquisitionDateTime   if set(['AcquisitionDateTime']).issubset(self.ECG.dir()) else ''
                self.AcquisitionTimeZone   = self.ECG.TimezoneOffsetFromUTC if set(['TimezoneOffsetFromUTC']).issubset(self.ECG.dir()) else ''
                self.Manufacturer          = self.ECG.Manufacturer          if set(['Manufacturer']).issubset(self.ECG.dir()) else ''
                self.ManufacturerModelName = self.ECG.ManufacturerModelName if set(['ManufacturerModelName']).issubset(self.ECG.dir()) else ''
                self.SoftwareVersions      = self.ECG.SoftwareVersions      if set(['SoftwareVersions']).issubset(self.ECG.dir()) else ''
                self.DataExportedBy        = self.ECG.IssuerOfPatientID     if set(['IssuerOfPatientID']).issubset(self.ECG.dir()) else ''

                
                # Waveform information
                # Check existance of raw Waveform
                try:
                    self.Waveforms             = self.ECG.waveform_array(0).T
                except:
                    if (verbose==True): 
                        print('No ECG waveform present: inspect DICOM file')
                    return(pd.DataFrame({'filetype': [self.ECG.data_element('SOPClassUID').repval], 'error': 'No ECG waveform present'}) )
                
                # Add channel settings and lead information 
                wave                       = self.ECG.WaveformSequence[0]
                settings                   = self.ECG.WaveformSequence[0].ChannelDefinitionSequence[0]
                self.TimeOffset            = wave.MultiplexGroupTimeOffset if set(['MultiplexGroupTimeOffset']).issubset(wave.dir()) else '' 
                self.ChannelNumber         = wave.NumberOfWaveformChannels if set(['NumberOfWaveformChannels']).issubset(wave.dir()) else '' 
                self.ChannelSensitivity    = settings.ChannelSensitivity   if set(['ChannelSensitivity']).issubset(settings.dir()) else '' 
                self.ChannelBaseline       = settings.ChannelBaseline      if set(['ChannelBaseline']).issubset(settings.dir()) else '' 
                self.ChannelSampleSkew     = settings.ChannelSampleSkew    if set(['ChannelSampleSkew']).issubset(settings.dir()) else '' 
                self.FilterLowFrequency    = settings.FilterLowFrequency   if set(['FilterLowFrequency']).issubset(settings.dir()) else '' 
                self.FilterHighFrequency   = settings.FilterHighFrequency  if set(['FilterHighFrequency']).issubset(settings.dir()) else '' 
                self.NotchFilterFrequency  = settings.NotchFilterFrequency if set(['NotchFilterFrequency']).issubset(settings.dir()) else ''
                self.sf                    = wave.SamplingFrequency
                self.sf_original           = wave.SamplingFrequency
                # Check number of channels in ECG waveform
                if self.ChannelNumber >= 8: 
                    self.lead_info_final       = self.lead_info(0)
                    self.LeadVoltages          = self.make_leadvoltages(0)
                else: 
                    self.LeadVoltages          = np.zeros((0,0))
                    if (verbose==True): 
                        print('Limited number of ECG waveform channels present: inspect DICOM file')
                    #return(pd.DataFrame({'filetype': [self.ECG.data_element('SOPClassUID').repval], 'error': 'Less than 8 or 12-leads ECG present:'+str(self.ECG.WaveformSequence[0][0x003a0005].value)}) )

                # Check existance of MedianWaveform
                try: 
                    self.MedianWaveforms       = self.ECG.waveform_array(1).T
                    self.lead_info_final       = self.lead_info(1)
                    self.LeadVoltages2         = self.make_leadvoltages(1)
                    self.MedianWaveformPresent = 'Yes'
                except: 
                    self.MedianWaveformPresent = 'No'
                    self.MedianWaveforms       = ''
                    self.LeadVoltages2         = np.zeros((0,0))
                    if (verbose==True): 
                        print('No Median Waveform present')

                self.samplingfrequency     = self.resampling_500hz()

                # Create dictionary from the above 
                self.read_dict_final       = self.readable_dict()
            
            return self.read_dict_final

        except Exception as e:
            print(str(e))
            pass

    def readable_dict(self):
        """Make a readable dict"""
        read_dict                             = {}
        read_dict["SOPinstanceUID"]           = self.SOPinstanceUID
        read_dict["SERIESinstanceUID"]        = self.SERIESinstanceUID
        read_dict["STUDYinstanceUID"]         = self.STUDYinstanceUID
        read_dict["PatientID"]                = self.PatientID
        read_dict["PatientBirthDate"]         = self.PatientBirthDate
        read_dict["PatientName"]              = str(self.PatientName)
        read_dict["PatientSex"]               = self.PatientSex
        read_dict["StudyDate"]                = datetime.strptime(self.StudyDate, "%Y%m%d").strftime('%Y-%m-%d') #if your date is different format adapt
        read_dict["StudyTime"]                = self.StudyTime[:6]
        read_dict["StudyDescription"]         = self.StudyDescription
        read_dict["AcquisitionDateTime"]      = datetime.strptime(self.AcquisitionDateTime[:14], "%Y%m%d%H%M%S").strftime('%Y-%m-%d %H:%M:%S') 
        read_dict["AcquisitionTimeZone"]      = self.AcquisitionTimeZone
        #read_dict["DatapointsWaveform"]       = len(list(self.LeadVoltages.values())[0])
        #read_dict["DatapointsMedianWaveform"] = len(list(self.LeadVoltages2.values())[0])
        read_dict["AccessionNumber"]          = self.AccessionNumber
        read_dict["SamplingFrequency"]        = self.samplingfrequency if self.ChannelNumber >= 8 else self.sf_original
        read_dict["OriginalSamplingFrequency"]= self.sf_original
        read_dict["ChannelNumber"]            = self.ChannelNumber
        read_dict["ChannelSensitivity"]       = self.ChannelSensitivity
        read_dict["ChannelBaseline"]          = self.ChannelBaseline
        read_dict["ChannelSampleSkew"]        = self.ChannelSampleSkew
        read_dict["FilterLowFrequency"]       = self.FilterLowFrequency
        read_dict["FilterHighFrequency"]      = self.FilterHighFrequency
        read_dict["NotchFilterFrequency"]     = self.NotchFilterFrequency
        read_dict["Manufacturer"]             = self.Manufacturer
        read_dict["ManufacturerModelName"]    = self.ManufacturerModelName
        read_dict["SoftwareVersions"]         = self.SoftwareVersions
        read_dict["DataExportedBy"]           = self.DataExportedBy
        read_dict["Waveforms"]                = self.LeadVoltages
        read_dict["MedianWaveforms"]          = self.LeadVoltages2
        return read_dict

    def make_leadvoltages(self,nr):
        """Extracts the voltages out of the DICOM"""
        num_leads = 0
        leads = {}

        for i, lead in enumerate(self.ECG.waveform_array(nr).T):
            num_leads += 1
            leads[self.lead_info_final[i]] = lead
        if num_leads == 8 and self.augmentLeads:
            leads['III'] = np.subtract(leads['II'], leads['I'])
            leads['aVR'] = np.add(leads['I'], leads['II']) * (-0.5)
            leads['aVL'] = np.subtract(leads['I'], 0.5 * leads['II'])
            leads['aVF'] = np.subtract(leads['II'], 0.5 * leads['I'])
        return leads

    def lead_info(self,nr):
        """returns the names of the channels from the DICOM"""
        leadnames = {}
        for ii, channel in enumerate(self.ECG.WaveformSequence[nr].ChannelDefinitionSequence):
            source = channel.ChannelSourceSequence[0].CodeMeaning
            units = "unitless"
            if "ChannelSensitivity" in channel:
                units = channel.ChannelSourceSequence[0].CodeMeaning
            leadnames[ii] = source.replace('Lead','').strip()
        return leadnames

    def resampling_500hz(self):
        """In case sf is 250, make 500"""
        if self.resample_500 is False:
            return int(self.sf)
        else:
            if int(self.sf) != 500:
                for i in self.LeadVoltages:
                    self.LeadVoltages[f"{i}"] = signal.resample(self.LeadVoltages[f"{i}"], 5000)
                    self.oversampled = "Yes"
                if self.MedianWaveformPresent == 'Yes':
                    for i in self.LeadVoltages:
                        self.LeadVoltages2[f"{i}"] = signal.resample(self.LeadVoltages2[f"{i}"], 600)
                self.sf = 500
                return self.sf
            else:
                return self.sf
        return self