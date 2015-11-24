#!/usr/bin/python

import re
import os

# filenames to look at
wavemeter_dir = '21.22.november'

# output file to write to
output_file = '21and22nov_wavelengths.txt'

# timespan over which should be averaged
averaging_time = 120

# this matches lines like 613969	942,751045
wavelength_pattern = r'([0-9]{2,7})\t([0-9]{3,4}[,|.][0-9]{6})'

wl_regex = re.compile(wavelength_pattern)

# loop through all tof files and find the corresponding wavemeter file
# once the corresponding wavemeter file is found, we average every 
# $timespan and write the wavelength to a new file.
# this new file contains the tof filename and the corresponding
# wavelength for each $timespan

output_fh = open(output_file, 'w+')

all_wl_files = sorted(os.listdir(wavemeter_dir))

for wl_file in all_wl_files:
    fh = open(wavemeter_dir + '/' + wl_file, 'r')
    print wl_file

    linecounter = 0
    averagecounter = 0
    wl_sum = 0
    printcounter = 0

    for line in fh:
        values = wl_regex.match(line)
        # only use lines that actually match the above pattern
        if values:
            # get the first column (time) and the wavelength
            time = values.group(1)
            wl = values.group(2)
            # have we just crossed the 60 seconds? if so, reset and print avg
            if int(time)-averagecounter*averaging_time*1000 > averaging_time*1000:
                averagecounter = averagecounter + 1
                output_fh.write(str(wl_sum/linecounter) + '\r\n')
                printcounter = printcounter + 1
                linecounter = 0
                wl_sum = 0

            # count the line and sum up the wavelength
            wl_sum = wl_sum + float(wl.replace(',', '.'))
            linecounter = linecounter + 1
    fh.close()

output_fh.close()
