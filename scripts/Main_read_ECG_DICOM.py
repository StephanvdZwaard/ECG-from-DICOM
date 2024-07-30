# -------------------------------------------------------------------------------------
#                   Code for converting ECG DICOM's into tabular data
# -------------------------------------------------------------------------------------
#
# Description:      This script converts ECG data stored in DICOM into tabular data.
# Authors:          Stephan van der Zwaard
# Date:             18-06-2024
# Python.version:   3.9.13

# -------------------------------------------------------------------------------------
#                                  Settings & dependencies 
# -------------------------------------------------------------------------------------

# Install required packages
# pip install -r requirements.txt

# Import required libraries
import os;
import shutil;
import sys;
import numpy as np;
import pandas as pd;
#import matplotlib.pyplot as plt;
import progressbar; 
from datetime import datetime
from scipy import signal
from pydicom import dcmread
from pydicom.waveforms import multiplex_array
from functions import recursive_copy
from functions import get_batch_prefix

#Set options
np.set_printoptions(threshold = 500)

# Import and initialize class to parse DICOM file with ECG 
from ECGDICOMReader import ECGDICOMReader
ecgreader = ECGDICOMReader()
#verbose = True

# Define required paths
base_path      = "E:/AnacondaData/SvanderZwaard/Python/ecg-pipeline/"
path_to_dicom  = "E:/DataExchange/Hartcentrum/ECG_input/" #base_path+"Python/ECG/dicom/"
path_to_archive= path_to_dicom+"processed/"
path_to_export = "E:/DataExchange/Hartcentrum/ECG_output/" #"D:/AnacondaData/Stephan/"
path_to_logs   = "E:/DataExchange/Hartcentrum/ECG_archive/logs/"

# Define current date for filename of output
today = datetime.today().strftime('%Y%m%d')

# -------------------------------------------------------------------------------------
#                       Data pipeline: conversion from DICOM to CSV 
# -------------------------------------------------------------------------------------

# Retrieve all DICOM files to be converted
ECG_files = list() 
for path, subdirs, files in os.walk(path_to_dicom):
    subdirs[:] = [d for d in subdirs if d not in path_to_archive]
    for filename in files:
        if path.endswith('/'):
            ECG_files.append(path+filename)
        else:
            ECG_files.append(path+'/'+filename)
#file = recursive_copy(path_to_dicom+'ECG-DICOM-DT4H-AMC/')
#print(ECG_files)
  
# Loop over all DICOM files to convert using ECGDICOMReader()

# Preallocate dataframes()
general_info   = pd.DataFrame()
summary        = pd.DataFrame()
median_waves   = pd.DataFrame()
original_waves = pd.DataFrame()
error_dicom    = pd.DataFrame()

# Set-up progressbar
i_start = 0
i_end   = len(ECG_files) #account for Python indexing
offset  = 0 #number of scans already processed today.
print('Number of DICOMs: '+str(len(range(i_start,i_end))))

pbar = progressbar.ProgressBar(widgets = [progressbar.Percentage(), " ", progressbar.GranularBar(), " ", progressbar.ETA()], 
                               maxval = len(ECG_files[i_start:i_end])-1, 
                               redirect_stdout=True);

for i in range(i_start,i_end) : 
    
    try:
        dicom = ecgreader(ECG_files[i], verbose=False)

        # Generate Table 1: median waveform
        mw = pd.DataFrame(dicom['MedianWaveforms'])
        mw = mw.add_prefix('lead_')
        mw["id"] = mw.index
        mw = pd.wide_to_long(mw, stubnames ='lead_', i="id", j="lead",suffix = r'\w+').sort_index(level=0)
        mw["SOPinstanceUID"] = dicom["SOPinstanceUID"]
        mw["waveform"] = "median_beat"
        mw = mw.reset_index()
        mw = mw[["SOPinstanceUID", "waveform", "lead","id","lead_"]]
        mw = mw.rename(columns = {'lead_':'voltage', 'id':'sample_id', 'SOPinstanceUID':'record_id_ecg'})

        # Generate Table 2: waveform rhythm
        w = pd.DataFrame(dicom['Waveforms'])
        w = w.add_prefix('lead_')
        w["id"] = w.index
        w = pd.wide_to_long(w, stubnames ='lead_', i="id", j="lead",suffix = r'\w+').sort_index(level=0)
        w["SOPinstanceUID"] = dicom["SOPinstanceUID"]
        w["waveform"] = "rhythm"
        w = w.reset_index()
        w = w[["SOPinstanceUID", "waveform","lead", "id","lead_"]]
        w = w.rename(columns = {'lead_':'voltage', 'id':'sample_id', 'SOPinstanceUID':'record_id_ecg'})

        # Generate Table 3: summary
        #summs = dicom['Summary']
        #result = [summs[key] for key in summs if key.startswith('Summary')]
        #s_text = ' \n'.join([str(item) for item in result])
        #s_values = dict(filter(lambda item: not item[0].startswith('Summary'),
        #          summs.items()))
        #s = dict({'Summary':s_text}|s_values)
        #s = pd.DataFrame.from_dict(s, orient = 'index').transpose()
        #s["RECORD_ID_ECG"] = dicom["SOPinstanceUID"]
        #print(s.keys())

        # Generate Table 4: general info
        wave = dicom.pop('Waveforms')
        mbeat= dicom.pop('MedianWaveforms')
        #summ = dicom.pop('Summary')
        info = pd.DataFrame.from_dict(dicom, orient = 'index').transpose()
        info = info.rename(columns = {'SOPinstanceUID':'RECORD_ID_ECG'})

        # Combine data with previous records
        general_info    = pd.concat([general_info,info], axis=0)
        #summary         = pd.concat([summary,s], axis=0)
        median_waves    = pd.concat([median_waves,mw], axis=0)
        original_waves  = pd.concat([original_waves,w], axis=0)
    
    except: # Retrieve relevant information when DICOM could not be read including the error
        dicom['file_no']  = i+1 #account for python indexing
        dicom['filename'] = ECG_files[i].replace(path_to_dicom,'')
        error_dicom       = pd.concat([error_dicom,dicom], axis=0)

    # Move processed ECG DICOMs to archive to distinguish processed from unread files.
    if not os.path.exists(path_to_archive):
           os.makedirs(path_to_archive)
    os.rename(ECG_files[i], ECG_files[i].replace(path_to_dicom,path_to_archive))

    # Save to CSV-files (at end of query or for every X records defined by batch size)
    batch_size = 1000
    if ((i+1) == i_end or (i+1)%batch_size == 0): #account for python indexing

        if ((i+1) == i_end):
            
            batch_pre  = i_start if (i_end)<=batch_size else int(((i//batch_size))*batch_size)
            batch_post = i_end

        elif ((i+1)%batch_size == 0):
            
            batch_pre  = i_start if (i+1)==batch_size else int((((i+1)/batch_size)-1)*batch_size)
            batch_post = int(((i+1)/batch_size)*batch_size)

        # Summarise batch range
        batch_prefix  = get_batch_prefix(batch_pre+offset)
        batch_postfix = get_batch_prefix(batch_post+offset)
        batch         = batch_prefix+str(batch_pre+offset)+'_'+batch_postfix+str(batch_post+offset) 

        # Save separate CSV-files for each batch
        error_dicom.to_csv(path_to_logs+today+'_'+'DICOM_error_'+batch+'.csv', index=False)
        general_info.to_csv(path_to_logs+today+'_'+'DICOM_ECG_GENERALINFO_'+batch+'.csv', index=False)
        general_info.to_csv(path_to_export+today+'_'+'DICOM_ECG_GENERALINFO_'+batch+'.csv', index=False)
        #summary.to_csv(path_to_export+today+'_'+'DICOM_ECG_SUMMARY_'+batch+'.csv', index=False)
        median_waves.to_csv(path_to_export+today+'_'+'DICOM_ECG_WAVEFORM_MEDIANBEAT_'+batch+'.csv', index=False)
        original_waves.to_csv(path_to_export+today+'_'+'DICOM_ECG_WAVEFORM_RHYTHM_'+batch+'.csv', index=False)

        # Preallocate dataframe after saving
        error_dicom    = pd.DataFrame()
        general_info   = pd.DataFrame()
        #summary        = pd.DataFrame()
        median_waves   = pd.DataFrame()
        original_waves = pd.DataFrame()

    # Update progressbar
    pbar.update(i-i_start)
    

# If finished
pbar.finished()
print('\nConversion finished! --- ')