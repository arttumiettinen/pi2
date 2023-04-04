#!/usr/bin/env python3


# Non-rigid stitching script.
# In order to use this script, please make a stitch_settings.txt file according to
# the provided template. Place it into the directory where the output and temporary files
# can be placed and run this script (without command line arguments).


from base import *

import time
import configparser
import argparse


class NonGridStitchSettings:
    """
    Stores values of all the stitching settings.
    Values are read from input file or set to default values.
    """

    def __init__(self):
        self.sample_name = ''
        self.binning = 1
        self.dimensions = np.array([0, 0, 0])
        self.point_spacing = 0
        self.coarse_block_radius = ['', '', '']
        self.coarse_binning = [0, 0, 0]
        self.fine_block_radius = ['', '', '']
        self.fine_binning = [0, 0, 0]
        self.normalize_in_blockmatch = 0
        self.normalize_while_stitching = 0
        self.filter_threshold = '0'
        self.global_optimization = True
        self.allow_rotation = True
        self.allow_local_shifts = True
        self.allow_local_deformations = True




def get_position(config, key):
    """
    Gets a position value from configuration dictionary if the key is set. If not, returns default.
    """

    if 'stitch' in config:
        if key in config['stitch']:
            return config['stitch'][key]

    return default



def overlaps_range(x1, x2, y1, y2):
    """
    Tests if [x1, x2] overlaps with [y1, y2].
    """

    return max(x1, y1) <= min(x2, y2)


def overlaps(scan1, scan2):
    """
    Tests if two scans overlap.
    """

    return overlaps_range(scan1.position[0], scan1.position[0] + scan1.dimensions[0], scan2.position[0], scan2.position[0] + scan2.dimensions[0]) and overlaps_range(scan1.position[1], scan1.position[1] + scan1.dimensions[1], scan2.position[1], scan2.position[1] + scan2.dimensions[1]) and overlaps_range(scan1.position[2], scan1.position[2] + scan1.dimensions[2], scan2.position[2], scan2.position[2] + scan2.dimensions[2])

def to_bool(val):
    """
    Convert val to boolean.
    """

    return val in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']


def main():
    """
    Main entry point for elastic stitching.
    Reads configuration from stitch_settings.txt and runs the stitching process.
    """

    argsparser = argparse.ArgumentParser(description='Perform non-rigid stitching process. Requires stitch settings file and tile images (as specified in the stitch settings file) as input, and produces the stitched image. Optionally runs individual computation jobs using a cluster. After running this script, consider using rm-stitch-temp shell script to remove all the temporary files generated.')
    argsparser.add_argument('settings_file', nargs='?', default="stitch_settings.txt", type=str, help='The stitch settings file.')
    argsparser.add_argument('--sample_name', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--binning', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--point_spacing', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--coarse_block_radius', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--coarse_binning', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--fine_block_radius', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--fine_binning', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--normalize_in_blockmatch', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--normalize_while_stitching', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--mask_to_max_circle', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--global_optimization', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--allow_rotation', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--allow_local_deformations', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--displacement_filter_threshold', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--create_goodness', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--max_block_size', type=str, help='Overrides corresponding setting read from the configuration file.')
    argsparser.add_argument('--cluster', type=str, help='Overrides corresponding setting read from the configuration file.')
    
    args = argsparser.parse_args()


    settings_file = args.settings_file
    #if len(sys.argv) > 1:
    #    settings_file = sys.argv[1]
    

    print(f"Reading stitch settings from {settings_file}")

    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(settings_file)

    read_global_settings(config, args)


    # Stitching settings
    settings = NonGridStitchSettings()
    settings.sample_name = get(config, 'sample_name', 'sample')
    if args.sample_name:
        settings.sample_name = args.sample_name

    # Output binning
    sval = get(config, 'binning', '1')
    if args.binning:
        sval = args.binning
    settings.binning = int(sval)

    # Space between reference points in the deformation field in pixels
    settings.point_spacing = 1
    sval = get(config, 'point_spacing', '20')
    if args.point_spacing:
        sval = args.point_spacing
    settings.point_spacing = int(sval)
    
    # Size of calculation block in coarse matching
    # Indicates also maximum local shift between images
    sval = get(config, 'coarse_block_radius', '[25, 25, 25]')
    if args.coarse_block_radius:
        sval = args.coarse_block_radius
    settings.coarse_block_radius = from_string(sval)

    # Amount of binning for coarse matching
    settings.coarse_binning = 0
    sval = get(config, 'coarse_binning', '4')
    if args.coarse_binning:
        sval = args.coarse_binning
    settings.coarse_binning = int(sval)

    # Size of calculation block in fine matching
    # Indicates maximum shift on top of shift determined in coarse matching.
    sval = get(config, 'fine_block_radius', '[25, 25, 25]')
    if args.fine_block_radius:
        sval = args.fine_block_radius
    settings.fine_block_radius = from_string(sval)
    
    # Amount of binning for fine matching
    settings.fine_binning = 0
    sval = get(config, 'fine_binning', '1')
    if args.fine_binning:
        sval = args.fine_binning
    settings.fine_binning = int(sval)

    # Normalization flags
    sval = get(config, 'normalize_in_blockmatch', 'True')
    if args.normalize_in_blockmatch:
        sval = args.normalize_in_blockmatch
    settings.normalize_in_blockmatch = to_bool(sval)

    sval = get(config, 'normalize_while_stitching', 'True')
    if args.normalize_while_stitching:
        sval = args.normalize_while_stitching
    settings.normalize_while_stitching = to_bool(sval)

    # Mask each input image to maximum inscribed circle?
    sval = get(config, "mask_to_max_circle", 'False')
    if args.mask_to_max_circle:
        sval = args.mask_to_max_circle
    settings.max_circle = to_bool(sval)

    # Allow global optimization of transformations?
    sval = get(config, 'global_optimization', 'True')
    if args.global_optimization:
        sval = args.global_optimization
    settings.global_optimization = to_bool(sval)

    # Allow rotational transformation
    sval = get(config, 'allow_rotation', 'True')
    if args.allow_rotation:
        sval = args.allow_rotation
    settings.allow_rotation = to_bool(sval)

    # Allow local tranformations
    sval = get(config, 'allow_local_deformations', 'True')
    if args.allow_local_deformations:
        sval = args.allow_local_deformations
    settings.allow_local_deformations = to_bool(sval)

    # Threshold for displacement filtering
    sval = get(config, 'displacement_filter_threshold', '3')
    if args.displacement_filter_threshold:
        sval = args.displacement_filter_threshold
    settings.filter_threshold = float(sval)

    # Create stitch goodness output?
    sval = get(config, 'create_goodness', 'False')
    if args.create_goodness:
        sval = args.create_goodness
    settings.create_goodness_file = to_bool(sval)

    # Read image names and locations
    if not ('positions' in config):
        raise RuntimeError("No 'positions' section found in the stitching settings file.")
    
    settings_contents = ""
    data = config.items('positions')
    index = 0
    relations = nx.DiGraph()
    for line in data:
        file = line[0]
        position = np.fromstring(line[1], sep=',')

        sc = Scan()
        sc.rec_file = file
        sc.dimensions = get_image_size(file)
        sc.index = index
        sc.position = position

        relations.add_node(sc, scan=sc)

        out_line = f"{file} at position {position}, size = {sc.dimensions}"
        settings_contents += out_line + "\n"
        print(out_line)

        index = index + 1


    # Find if the positions section has changed
    settingsfile = f"{settings.sample_name}_position_settings.txt"
    current_contents = get_contents(settingsfile)
    redo_all = False
    if settings_contents != current_contents:
        # Tile positions or sizes have changed. Redo everything.
        print("Tile positions or sizes have changed.")
        redo_all = True

    write_contents(settingsfile, settings_contents)

    # Perform binning
    redo_all_binning = redo_all
    while True:
        if auto_binning(relations, settings.binning, redo_all_binning):
            break;

        redo_all_binning = False;
        wait_for_cluster_jobs()
        
    if settings.binning != 1:
        for node in relations.nodes():
            node.rec_file = node.binned_file
            node.dimensions = node.dimensions / settings.binning
            node.position = node.position / settings.binning

    settings.point_spacing = int(np.round(settings.point_spacing / settings.binning))
    settings.coarse_block_radius = np.round(settings.coarse_block_radius / settings.binning).astype(int)
    settings.coarse_binning = np.maximum(1, np.round(settings.coarse_binning / settings.binning).astype(int))
    settings.fine_block_radius = np.round(settings.fine_block_radius / settings.binning).astype(int)
    settings.fine_binning = np.maximum(1, np.round(settings.fine_binning / settings.binning).astype(int))

    settings.filter_threshold /= settings.binning

    if settings.binning != 1:
        settings.sample_name = f"bin{settings.binning}_{settings.sample_name}"



    # Find overlapping images
    for node1 in relations.nodes():
        for node2 in relations.nodes():
            if node1.index < node2.index and overlaps(node1, node2):
                relations.add_edge(node1, node2)

    # Show the relations
    #import matplotlib
    #matplotlib.use('agg')
    #import matplotlib.pyplot as plt
    #positions = {scan: [scan.position[1], scan.position[2]] for scan in relations.nodes()}
    #nx.draw(relations, positions)
    #plt.show()
    #plt.savefig('network.png')

    # Calculate displacement fields
    redo_all_displacements = redo_all
    while True:
        if calculate_displacement_fields(settings.sample_name, relations, settings.point_spacing, settings.coarse_block_radius, settings.coarse_binning, settings.fine_block_radius, settings.fine_binning, settings.normalize_in_blockmatch, settings.filter_threshold, redo_all_displacements):
            break

        redo_all_displacements = False
        wait_for_cluster_jobs()
        
    read_displacement_fields(settings.sample_name, relations, settings.allow_rotation)

    run_stitching_for_all_connected_components(relations, settings.sample_name, settings.normalize_while_stitching, settings.max_circle, settings.global_optimization, settings.allow_rotation, settings.allow_local_deformations, settings.create_goodness_file)

    wait_for_cluster_jobs()


if __name__ == "__main__":
    start = time.time()
    main()
    elapsed = time.time() - start
    print(f"Operation took {elapsed} seconds.")
