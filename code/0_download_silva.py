# Campen_Brine_2018
# 0_download_silva.py

"""
This file downloads the mothur compatible SILVA v128 reference files, and
trims it to the V4 region.

By using the SILVA files you should read and agree to their dual-use license
available at https://www.arb-silva.de/silva-license-information.

"""

import os
import shutil
import tarfile
import wget

from mothur_py import Mothur

# configure directories
base_dir = os.path.join('..')
base_dir = os.path.abspath(base_dir)
data_dir = os.path.join(base_dir, 'data')

silva_dir = os.path.join(data_dir, 'silva')
if not os.path.isdir(silva_dir):
	os.mkdir(silva_dir)

# configure path to mothur executeable
# some example paths provided
MOTHUR_PATH = 'mothur'  # mothur is installed in PATH environment variable
# MOTHUR_PATH = './mothur'  # mothur is in current dir (linux/mac)
# MOTHUR_PATH = 'mothur.exe'  # mothur is in current dir (windows)

# setup source url and target file
url = 'https://mothur.org/w/images/b/b4/Silva.nr_v128.tgz'
compressed_file_path = os.path.join(silva_dir, 'silva_v128.tgz')
file_path = os.path.join(silva_dir, 'silva_v128')

# conditionally download Silva reference
print('-Downloading Silva v128 reference files...')
if not os.path.isfile(compressed_file_path):
	wget.download(url=url, out=compressed_file_path)
	print(' ...done!')
else:
	print(' ...Silva v128 already downloaded, skipping')

# conditionally unpack Silva reference
print('-Unpacking Silva v128 reference files...')
if not os.path.isdir(file_path):
	with tarfile.open(compressed_file_path, mode='r:gz') as tar:
		tar.extractall(path=file_path)
	print(' ...done!')
else:
	print(' ...Silva v128 already unpacked, skipping')

# trim Silva reference alignment to V4 region
v4_file = os.path.join(file_path, 'silva.nr_v128.V4.align')
print('-Trimming Silva v128 reference alignment to V4 region...')
if not os.path.isfile(v4_file):

	# configure mothur to do the trimming
	m = Mothur()
	m.mothur_path = MOTHUR_PATH
	m.suppress_logfile = True
	m.verbostiy = 1
	m.current_dirs['input'] = m.current_dirs['output'] = file_path

	# trim full alignment to specific region
	m.pcr.seqs(
		fasta='silva.nr_v128.align', 
		start=13862, 
		end=23445, 
		keepdots=False, 
		processors=2
	)

	# rename file
	shutil.move(
		src=os.path.join(file_path, 'silva.nr_v128.pcr.align'),
		dst=v4_file
	)

	print(' ...done!')
else:
	print(' ...Silva v128 trimmed to V4 region already exists, skipping')

print('\n\nDownloading and Trimming of Silva v128 to V4 region complete!!')
