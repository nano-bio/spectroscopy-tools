#!/usr/bin/python

# labbook filelist
lb_fl = 'input_864nm.txt'
# download folder
dl_folder = '864nm'
# download url
dl_url = 'http://138.232.72.25/clustof/export/'

import urllib
import os

# create folder
try:
    os.stat(dl_folder)
except:
    os.mkdir(dl_folder)

# file handle for labbook filelist
fh_lb_fl = open(lb_fl, 'r')

# go through file list
for entry in fh_lb_fl:
    entrylist = entry.split(' ')
    complete_filename = dl_folder + '/' + entrylist[1] + ' ' + entrylist[2].replace('\n','')
    complete_filename.replace('  ', ' ')
    print('Downloading measurement ID: ' + entrylist[0])
    urllib.urlretrieve(dl_url + entrylist[0], complete_filename)
