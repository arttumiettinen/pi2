#!/usr/bin/env python3

import argparse
import glob
import json
import re
import os
from logFileParser import logFileParser

from default_stitch_settings import write_stitch_settings

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
        argsparser.add_argument('-j', '--json', action='store_true', help='Set to true to use the logs from json files. Default is False.')
        argsparser.add_argument('-n', '--name', default='stitched', type=str, help='Sample name.')
        argsparser.add_argument('-r', '--recdir', default='', type=str, help='Try to find reconstructions from this directory. Use e.g. if the reconstructions were done in Ra and the original log files do not contain reconstruction folder. String %%s will be replaced by sub-scan name, e.g. /das/work/p1234/Data10/disk1/%%s/ may expand to /das/work/p1234/Data10/disk1/01_BigSample_B7/. Can be used multiple times, e.g. /das/work/p1234/Data10/disk1/%%s/%%s/ expands to /das/work/p1234/Data10/disk1/01_BigSample_B7/01_BigSample_B7/')
        argsparser.add_argument('-m', '--mask', action='store_true', help='Set to true to mask the input images to the maximum inscribed circle in each cross-section. Use this setting to erase possible bad data outside of the well-reconstructed region.')
        argsparser.add_argument('-u', '--units', default='um', type=str, help='Units used for sample coordinates in the log files. Can be mm or um.', choices=['mm', 'um'])
        argsparser.add_argument('-xx', '--xx_sign', default='positive', type=str, help='Sign of the XX coordinate. Can be "positive" or "negative"', choices=['positive', 'negative'])
        argsparser.add_argument('-zz', '--zz_sign', default='positive', type=str, help='Sign of the ZZ coordinate. Can be "positive" or "negative"', choices=['positive', 'negative'])
        argsparser.add_argument('--valid_scans', nargs='*', default='', type=str, help='Specify valid scans to be stitched.')


        args = argsparser.parse_args()

        coord_factor = 1
        if args.units == 'mm':
            coord_factor = 1000

        xx_fac = 1
        # If XX axis needs to be inverted
        if args.xx_sign == 'positive':
            xx_fac = 1
        elif args.xx_sign == 'negative':
            xx_fac = -1

        zz_fac = 1
        # If ZZ axis needs to be inverted
        if args.zz_sign == 'positive':
            zz_fac = 1
        elif args.zz_sign == 'negative':
            zz_fac = -1

        if len(args.scan_folders) <= 0:
            args.scan_folders = ['./*']

        if args.recdir:
            if '%s' not in args.recdir:
                args.recdir = args.recdir + '%s'
                #print(f"Reconstruction directory override does not contain %s. Use %s to mark the location where the scan name is to be inserted.")
                #return 1
            print(f"Reconstruction directory override: {args.recdir}")

        logs = []
        for folder in args.scan_folders:
            if not os.path.isfile(folder):
                logglob = f"{folder}/tif/*.json" if args.json else f"{folder}/tif/*.log"
                files = glob.glob(logglob)
                if len(files) > 0:
                    logs.extend(files)
                else:
                    logglob = f"{folder}/*.json" if args.json else f"{folder}/*.log"
                    files = glob.glob(logglob)
                    if len(files) > 0:
                        logs.extend(files)
                    else:
                        print(f"Warning: No log files found from {folder}")
            else:
                # folder is actually file name. We assume that it is name of log file to be stitched.
                # Skip non-supported log files (like .h5 files that often get into the folder list)
                if folder.endswith('xml') or folder.endswith('json' if args.json else 'log'):
                    logs.append(folder)

        # Remove '_dacatXXX.log' entries
        logs = [p for p in logs if '_dacat' not in p]

        if args.json:
            # Remove '_config.json entries
            logs = [p for p in logs if '_config' not in p]

        # Normalize paths
        logs = [os.path.normpath(p) for p in logs]

        # Remove duplicates
        logs = list(set(logs))

        # Sort to natural order
        logs.sort(key=natural_sort_key)
        

        if len(logs) <= 0:
            raise SystemExit("No log files found.")

        first_x = 0
        first_y = 0
        first_z = 0

        good_count = 0
        positions = ""

        def is_valid_scan(log_file_path):
            for scan in args.valid_scans:
                if scan in log_file_path:
                    return True
            return False
        
        for logfile in logs:
            # Check if the scan is valid
            if args.valid_scans and not is_valid_scan(logfile):
                continue
            # Read log file
            print(f"Reading log from {logfile}...")
            if args.json:
                with open(logfile, 'r') as f:
                    logDict = json.load(f)
                scan_name = logDict['scientificMetadata']['scanParameters']['File Prefix']
                pixel_size = logDict['scientificMetadata']['detectorParameters']['Actual pixel size']['v']
                xx = logDict['scientificMetadata']['scanParameters']['Sample holder XX-position']['v']
                zz = logDict['scientificMetadata']['scanParameters']['Sample holder ZZ-position']['v']
                y = logDict['scientificMetadata']['scanParameters']['Sample holder Y-position']['v']
                gigafrost = logDict['scientificMetadata']['detectorParameters']['Camera'] == 'GigaFRoST'
            else:
                parser = logFileParser()
                parser.parser_input(logfile)

                logDict = parser.logDict

                scan_name = logDict['Scan Settings File Prefix']
                pixel_size = logDict['Detector Settings Actual pixel size [um]']
                xx = logDict['Sample user coordinates XX-coordinate']
                zz = logDict['Sample user coordinates ZZ-coordinate']
                y = logDict['Sample user coordinates Y-coordinate']
                gigafrost = logDict['Detector Settings Camera'] == 'GigaFRoST'
            
            rec_dir = ""
            if args.recdir:
                # Use reconstruction directory override
                rec_dir = args.recdir % tuple(scan_name for _ in range(str(args.recdir).count('%s'))) # Allow the user to pass multiple %s in the string, happens sometimes
            elif 'Reconstruction Parameters Reconstruction Dir' in logDict:
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
                X = zz * coord_factor * zz_fac / pixel_size
                Y = -xx * coord_factor * xx_fac / pixel_size
                Z = y * coord_factor / pixel_size

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
                rec_dir = os.path.normpath(rec_dir)#+ f"/rec_16bit_Paganin"

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

            mask_d = -1 # Rectangle mask
            if args.mask:
                mask_d = 0 # Automatic circle mask
            write_stitch_settings(args.name, args.binning, positions, cluster_name='Slurm', max_circle_diameter=mask_d)
            
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
