#!/usr/bin/env python

import argparse
import glob
import os

from default_stitch_settings import *


def get_value(lines, key):
    """
    Finds line from the lines list that starts with given key, and returns what is at the line after the key + one character.
    """

    for line in lines:
        if line.startswith(key):
            line = line[len(key)+1:]
            return line.strip()
    raise SystemExit(f"No key {key} found.");

def main():


    try:
        argsparser = argparse.ArgumentParser(description='Creates input file for nr_stitcher.py program from Phoenix scan log files. Saves the log file to the current directory and assumes that the input sub-images are also saved to the current folder.')
        argsparser.add_argument('pcm_name', nargs='?', default='', type=str, help='Specify name of .pcm file that contains information about the tiled scan. If no name is specified, the first .pcm file found from the current directory is used.')
        argsparser.add_argument('-b', '--binning', default=1, type=int, help='Binning that is to be applied before stitching.')
        argsparser.add_argument('-n', '--name', default='', type=str, help='Filename prefix for the stitched file and temporary files. If omitted or empty, the pcm file name is used.')

        args = argsparser.parse_args()

        # Check that we have .pcm file or try to find it
        if len(args.pcm_name) <= 0:
            lst = glob.glob('*.pcm')
            if len(lst) > 0:
                args.pcm_name = lst[0]
            else:
                raise SystemExit(f"No .pcm file specified and no .pcm file found from the current directory.")

        # Get info from the .pcm file
        with open(args.pcm_name, 'r') as file:
            lines = file.readlines() 

        image_prefix = os.path.basename(args.pcm_name)
        image_prefix = os.path.splitext(image_prefix)[0]

        if len(args.name) <= 0:
            args.name = image_prefix

        count = int(get_value(lines, 'NrScans'))
        axis = int(get_value(lines, 'MScanAxis'))

        print(f"The scan contains {count} tiles.")

        # axis == 1 corresponds to shift in z-direction in image coordinates
        if axis != 1:
            raise SystemExit(f"This script has been tested only for the case where MScanAxis=1, but in the present .pcm file MScanAxis={axis}.")

        positions = ''
        curr_pos = 0
        for n in range(0, count):

            X = 0
            Y = 0
            Z = curr_pos
            line = f"{image_prefix}s{n+1}.pcr = {X}, {Y}, {Z}"
            positions += line + "\n"

            line = get_value(lines, f"{n}")
            parts = line.split()
            if len(parts) == 3:
                pos = float(parts[0])
                distance = float(parts[1])
                voxshift = float(parts[2])
                curr_pos += voxshift

        print("Locations of tiles in pixels:")
        print(positions)

        print("Writing output to stitch_settings.txt...")
        write_stitch_settings(args.name, args.binning, positions, point_spacing=300, coarse_block_radius=200, coarse_binning=4, use_cluster=False, normalize_while_stitching=True)

        print('All done. Consider running nr_stitcher.py now.')



    except SystemExit as se:
        if len(se.args) > 0:
            print(se.args[0])
        return se.code


if __name__ == '__main__':
    main()