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