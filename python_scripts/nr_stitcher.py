#!/usr/bin/env python3


# Non-rigid stitching script that does not assume that the sub-images are arranged in a grid.
# In order to use this script, please make a stitch_settings.txt file according to
# the provided template. Place it into the directory where the output and temporary files
# can be placed and run this script (without command line arguments).


from base import *

import sys
import configparser


class NonGridStitchSettings:
    """
    Stores values of all the stitching settings.
    Values are read from input file or set to default values.
    """

    def __init__(self):
        sample_name = ''
        binning = 1
        dimensions = np.array([0, 0, 0])
        point_spacing = 0
        coarse_block_radius = ['', '', '']
        coarse_binning = [0, 0, 0]
        fine_block_radius = ['', '', '']
        fine_binning = [0, 0, 0]
        normalize_in_blockmatch = 0
        normalize_while_stitching = 0
        filter_threshold = '0'
        global_optimization = True
        allow_rotation = True
        allow_local_shifts = True
        allow_local_deformations = True
        force_redo = False




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



def main():
    """
    Main entry point for elastic stitching.
    Reads configuration from stitch_settings.txt and runs the stitching process.
    """

    settings_file = 'stitch_settings.txt'
    if len(sys.argv) >= 1:
        settings_files = sys.argv[0]

    print(f"Reading stitch settings from {settings_file}")

    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(settings_file)

    read_global_settings(config)


    # Stitching settings
    settings = NonGridStitchSettings()
    settings.sample_name = get(config, 'sample_name', 'sample')

    # Output binning
    settings.binning = int(get(config, 'binning', '1'))

    # Space between reference points in the deformation field in pixels
    settings.point_spacing = 1
    settings.point_spacing = int(get(config, 'point_spacing', 20))

    # Size of calculation block in coarse matching
    # Indicates also maximum local shift between images
    sval = get(config, 'coarse_block_radius', '[25, 25, 25]')
    settings.coarse_block_radius = from_string(sval)

    # Amount of binning for coarse matching
    settings.coarse_binning = 0
    settings.coarse_binning = int(get(config, 'coarse_binning', 4))

    # Size of calculation block in fine matching
    # Indicates maximum shift on top of shift determined in coarse matching.
    sval = get(config, 'fine_block_radius', '[25, 25, 25]')
    settings.fine_block_radius = from_string(sval)

    # Amount of binning for fine matching
    settings.fine_binning = 0
    settings.fine_binning = int(get(config, 'fine_binning', 1))

    # Normalization flags
    settings.normalize_in_blockmatch = get(config, 'normalize_in_blockmatch', True) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']
    settings.normalize_while_stitching = get(config, 'normalize_while_stitching', True) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Allow global optimization of transformations?
    settings.global_optimization = get(config, 'global_optimization', True) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Allow rotational transformation
    settings.allow_rotation = get(config, 'allow_rotation', True) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Allow local tranformations
    settings.allow_local_deformations = get(config, 'allow_local_deformations', True) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Threshold for displacement filtering
    settings.filter_threshold = float(get(config, 'displacement_filter_threshold', '3'))

    # Redo displacement fields?
    settings.force_redo = get(config, 'redo', False) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Create stitch goodness output?
    settings.create_goodness_file = get(config, 'create_goodness', False) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Read image names and locations
    if not ('positions' in config):
        raise RuntimeError("No 'positions' section found in the stitching settings file.")

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

        print(f"{file} at position {position}, size = {sc.dimensions}")

        index = index + 1


    # Perform binning
    while True:
        if auto_binning(relations, settings.binning):
            break;

        if is_use_cluster():
            if is_wait_for_jobs():
                wait_for_cluster_jobs()
            else:
                print("Please run this program again when all the cluster jobs have finished.")
                exit()

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
    while True:
        if calculate_displacement_fields(settings.sample_name, relations, settings.point_spacing, settings.coarse_block_radius, settings.coarse_binning, settings.fine_block_radius, settings.fine_binning, settings.normalize_in_blockmatch, settings.filter_threshold, settings.force_redo):
            break

        if is_use_cluster():
            if is_wait_for_jobs():
                wait_for_cluster_jobs()
            else:
                print("Please run this program again when all the cluster jobs have finished.")
                exit()

    read_displacement_fields(settings.sample_name, relations, settings.allow_rotation)

    run_stitching_for_all_connected_components(relations, settings.sample_name, settings.normalize_while_stitching, settings.global_optimization, settings.allow_rotation, settings.allow_local_deformations, settings.create_goodness_file)


    if is_use_cluster() and is_wait_for_jobs():
        wait_for_cluster_jobs()


if __name__ == "__main__":
    start = time.time()
    main()
    elapsed = time.time() - start
    print(f"Operation took {elapsed} seconds.")
