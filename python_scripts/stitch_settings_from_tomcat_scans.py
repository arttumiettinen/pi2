#!/usr/bin/env python

import argparse
import glob
import re
import os
from logFileParser import logFileParser

from default_stitch_settings import *

def natural_sort_key(s, nsre=re.compile('([0-9]+)')):
    """
    Key to use in various sorting commands in order to sort alphanumeric strings to natural order,
    i.e. str1 before str10 and str11.
    """
    return [int(text) if text.isdigit() else text.lower() for text in nsre.split(s)]

def main():


    try:
        argsparser = argparse.ArgumentParser(description='Creates input file for nr_stitcher.py program from TOMCAT scan log files. Saves the log file to the current directory and assumes that the input sub-images are also saved to the current folder.')
        argsparser.add_argument('scan_folders', nargs='*', type=str, help='Specify path that corresponds to base folder of each sub-scan. The base folder is the folder that contains sub-folders tif, rec_16bit_Paganin_0 etc. Use glob syntax to select multiple directories or pass multiple entries separated by space. Specify folders that contain the original scan log file in subdirectory tif.')
        argsparser.add_argument('-b', '--binning', default=1, type=int, help='Binning that is to be applied before stitching.')
        argsparser.add_argument('-n', '--name', default='stitched', type=str, help='Sample name.')
        argsparser.add_argument('-r', '--recdir', default='', type=str, help='Try to find reconstructions from this directory. Use e.g. if the reconstructions were done in Ra and the original log files do not contain reconstruction folder. String %%s will be replaced by sub-scan name, e.g. /das/work/p1234/Data10/disk1/%%s/ may expand to /das/work/p1234/Data10/disk1/01_BigSample_B7/')

        args = argsparser.parse_args()

        if len(args.scan_folders) <= 0:
            args.scan_folders = ['./*']

        if args.recdir:
            print(f"Reconstruction directory override: {args.recdir}")
            if '%s' not in args.recdir:
                print(f"Reconstruction directory override does not contain %s. Use %s to mark the location where the scan name is to be inserted.")
                return 1

        logs = []
        for folder in args.scan_folders:
            if not os.path.isfile(folder):
                logglob = f"{folder}/tif/*.log"
                files = glob.glob(logglob)
                if len(files) > 0:
                    logs.extend(files)
                else:
                    print(f"Warning: No log files found from {folder}")
            else:
                # folder is actually file name. We assume that it is name of log file to be stitched.
                # Skip non-supported log files (like .h5 files that often get into the folder list)
                if folder.endswith('xml') or folder.endswith('log'):
                    logs.append(folder)

        # Normalize paths
        logs = [os.path.normpath(p) for p in logs]

        # Remove duplicates
        logs = list(set(logs))

        # Sort to natural order
        logs.sort(key=natural_sort_key)
        

        if len(logs) <= 0:
            raise SystemExit(f"No log files found.")

        first_x = 0
        first_y = 0
        first_z = 0

        good_count = 0
        positions = ""
        for logfile in logs:
            # Read log file
            print(f"Reading log from {logfile}...")
            parser = logFileParser()
            parser.parser_input(logfile)

            scan_name = parser.logDict['Scan Settings File Prefix']
            pixel_size = parser.logDict['Detector Settings Actual pixel size [um]']
            xx = parser.logDict['Sample coordinates XX-coordinate']
            zz = parser.logDict['Sample coordinates ZZ-coordinate']
            y = parser.logDict['Sample coordinates Y-coordinate']
            gigafrost = parser.logDict['Detector Settings Camera'] == 'GigaFRoST'
            
            rec_dir = ""
            if args.recdir:
                # Use reconstruction directory override
                rec_dir = args.recdir % scan_name
            elif 'Reconstruction Parameters Reconstruction Dir' in parser.logDict:
                rec_dir = parser.logDict['Reconstruction Parameters Reconstruction Dir']
            else:
                #raise Exception(f"No reconstruction directory found from log file {logfile}, and no reconstruction directory override specified. Consider using --recdir command line argument.")
                print(f"No reconstruction directory found from log file {logfile}, and no reconstruction directory override specified. Consider using --recdir command line argument.")
            
            
            if rec_dir != "":
                #xroi = parser.logDict['Detector Settings X-ROI']
                #yroi = parser.logDict['Detector Settings Y-ROI']

                # Get image size from camera ROI
                #numbers = re.findall('\d+', xroi)
                #wh = int(numbers[1]) - int(numbers[0]) + 1
                #numbers = re.findall('\d+', yroi)
                #d = int(numbers[1]) - int(numbers[0]) + 1

                # Convert from beamline coordinates to image coordinates
                X = zz / pixel_size
                Y = -xx / pixel_size
                Z = y / pixel_size

                # GigaFRoST images use different coordinate system -> fix that here
                if gigafrost:
                    X *= -1
                    Y *= -1

                # Subtract coordinates of the first image so that it is always at (0, 0, 0).
                # This is not needed for stitching but in case of problems it is easier to
                # check the locations of the images visually if the first image is the reference point.
                #if logfile == logs[0]:
                if good_count == 0:
                    first_x = X
                    first_y = Y
                    first_z = Z

                X -= first_x
                Y -= first_y
                Z -= first_z

                # Normalize rec_dir path and make sure it exists
                rec_dir = os.path.normpath(rec_dir)

                rec_dir2 = os.path.dirname(logfile) + os.path.sep + ".." + os.path.sep + os.path.basename(rec_dir)
                rec_dir2 = os.path.normpath(rec_dir2)

                if not os.path.isdir(rec_dir) and os.path.isdir(rec_dir2):
                    # The original rec directory does not exist, try the secondary one
                    print(f"Warning: Original reconstruction directory {rec_dir} does not exist, using directory {rec_dir2}")
                    rec_dir = rec_dir2

                if os.path.isdir(rec_dir):

                    # Create an entry to the stitching settings file
                    line = f"{rec_dir} = {X}, {Y}, {Z}"
                    positions += line + "\n"
                    good_count += 1

                else:
                    print(f"Warning: reconstruction directory is not found from {rec_dir} or from {rec_dir2}")

        if good_count > 0:
            print("Locations of reconstructions that were found:")
            print(positions)

        if good_count > 1:
            print("Writing output to stitch_settings.txt...")

            write_stitch_settings(args.name, args.binning, positions)
            
            print('All done. Consider running nr_stitcher.py now.')
        else:
            print('Less than two images were found. There is no need for stitching. No output written.')

        return 0

    except SystemExit as se:
        if len(se.args) > 0:
            print(se.args[0])
        return se.code


if __name__ == '__main__':
    main()
