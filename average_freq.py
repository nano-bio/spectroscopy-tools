#!/usr/bin/python

import re
import os

# filename to look at
wavemeter_dir = 'wl-files'
tof_file_list = 'filelist_bereich_3_2_fine_960nm.txt'

# output file to write to
output_file = 'output_bereich_3_2_fine_960nm.txt'

# timespan over which should be averaged
averaging_time = 60

# this matches lines like 613969	942,751045
wavelength_pattern = r'^([0-9]{2,7})\t([0-9]{3,4},[0-9]{6})$'

# this matches lines like:
#935,818731 DataFile_2015.10.10-22h47m39s_AS.h5 
tof_file_pattern = r'^([0-9]{3},[0-9]{4,6})\s(DataFile_20[0-9]{2}.[0-9]{2}.[0-9]{2}-[0-9]{2}h[0-9]{2}m[0-9]{2}s_AS.h5)$'

wl_regex = re.compile(wavelength_pattern)
tof_regex = re.compile(tof_file_pattern)

# loop through all tof files and find the corresponding wavemeter file
# once the corresponding wavemeter file is found, we average every 
# $timespan and write the wavelength to a new file.
# this new file contains the tof filename and the corresponding
# wavelength for each $timespan

tof_fh = open(tof_file_list)
output_fh = open(output_file, 'w+')

all_wl_files = os.listdir(wavemeter_dir)

for tof_file in tof_fh:
    # see if the line is actually what we expect
    tof_file_details = tof_regex.match(tof_file)
    if tof_file_details:
        # it is!
        rough_wavelength = tof_file_details.group(1)
        tof_file_name = tof_file_details.group(2)
        # now lets find a wavemeter file, where the filename contains the wavelength
        found_wl_file = False
        for wl_file in all_wl_files:
            if wl_file.find(rough_wavelength) is not -1:
                found_wl_file = True
                fh = open(wavemeter_dir + '/' + wl_file, 'r')

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
                while printcounter < 10:
                    output_fh.write(str(wl_sum/linecounter) + '\r\n')
                    print('repeated last value, because otherwise there would only be 9 values: ' + wl_file)
                    printcounter = printcounter + 1

                fh.close()
        if found_wl_file is False:
            print('No WL-File found for: ' + tof_file_name)
            print('Corresponding rough wavelength: ' + rough_wavelength)

output_fh.close()
