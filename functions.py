import os
import shutil

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

def get_batch_prefix(batch):
            
    if batch < 0 or type(batch) is not int:
        raise ValueError('Incorrect batch index: inspect selected start/end index')
    elif batch < 100:
        prefix = '0000'
    elif batch <1000:
        prefix = '000'
    elif batch <10000:
        prefix = '00'
    elif batch <100000:
        prefix = '0'
    else:
        prefix = ''

    return(prefix)