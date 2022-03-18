from tqdm import tqdm
import pandas as pd
from urllib.request import urlretrieve
import sys, os
import progressbar
import shutil

"""
Small script to download ENCODE files
5 cols = SRS sequencing
12 cols = LRS sequencing
"""


pbar = None


def show_progress(block_num, block_size, total_size):
    global pbar
    if pbar is None:
        pbar = progressbar.ProgressBar(maxval=total_size)
        pbar.start()

    downloaded = block_num * block_size
    if downloaded < total_size:
        pbar.update(downloaded)
    else:
        pbar.finish()
        pbar = None


# def dl_file(file, output_dir, output_dir_cp):
def dl_file(file, output_dir):
	"""
	Download fct for pandas apply => for each row of a pandas df, fct will download file associated to id
	"""
	
    encode_prefix = "https://www.encodeproject.org"
    file_lite = file.split("/")[-1]
    if os.path.isfile(output_dir + file_lite) is False:
        print(file_lite)
        urlretrieve(encode_prefix + file, output_dir + file_lite, show_progress)
    else:
        print("File already exists : {}, copying ...".format(file_lite))

        # CP in a 2nd version in an another dir
        # shutil.copyfile(output_dir + file_lite, output_dir_cp + file_lite)


encode_list_files = sys.argv[1]
output_dir = sys.argv[2]
# output_dir_cp = sys.argv[3]
output_dir = output_dir + "/" if output_dir.endswith("/") is False else output_dir

df = pd.read_csv(encode_list_files, sep="\t", skiprows=1)

print(df)
# df["Download URL"].apply(lambda r: dl_file(r, output_dir, output_dir_cp))

# Core part to iterate on dataframe
df["Download URL"].apply(lambda r: dl_file(r, output_dir))


