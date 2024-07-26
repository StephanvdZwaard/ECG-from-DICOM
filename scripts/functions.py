## Contains helper functions for ECG conversion from DICOM to CSV

# Import required libraries
import os
import shutil

# Helper function to obtain all files within folder (incl subfolders)
def recursive_copy(path):

    for f in sorted(os.listdir(os.path.join(os.getcwd(), path))):

        file = os.path.join(path, f)

        if os.path.isfile(file):

            temp = os.path.split(path)
            f_name = '_'.join(temp)
            file_name = f_name + '_' + f
            shutil.move(file, file_name)

        else:

            recursive_copy(file)

# Helper function for filename when DICOMs are processed in batches (so that of equal length using prefilling with zero's).
def get_batch_prefix(batch):
            
    if batch < 0 or type(batch) is not int:
        raise ValueError('Incorrect batch index: inspect selected start/end index')
    elif batch == 0:
        prefix = '00000'
    elif batch < 100:
        prefix = '0000' # if batch size is 100, prefill with four zeros, so that 0000100
    elif batch < 1000:
        prefix = '000'
    elif batch < 10000:
        prefix = '00'
    elif batch < 100000:
        prefix = '0'
    else:
        prefix = ''

    return(prefix)