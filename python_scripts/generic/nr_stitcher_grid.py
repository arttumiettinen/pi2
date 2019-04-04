#!/usr/bin/env python3

# Non-rigid stitching script that does assumes that the sub-images are arranged in a grid.
# In order to use this script, please make a stitch_settings.txt file according to
# the provided template. Place it into the directory where the output and temporary files
# can be placed and run this script (without command line arguments).

import configparser
import time
import copy
import numpy as np
from base import *


class StitchSettings:
    """
    Stores values of all the stitching settings.
    Values are read from input file or set to default values.
    """
    binning = 1
    sample_name = ''
    input_file_template = ''
    grid_width = 0
    grid_height = 0
    grid_depth = 0
    overlap_x = 0
    overlap_y = 0
    overlap_z = 0
    point_spacing = 0
    coarse_block_radius = ['', '', '']
    coarse_binning = [0, 0, 0]
    fine_block_radius = ['', '', '']
    fine_binning = [0, 0, 0]
    normalize_in_blockmatch = 0
    normalize_while_stitching = 0
    filter_threshold = '0'
    xdir = 1
    ydir = 1
    zdir = 1
    global_optimization = True
    allow_rotation = True
    allow_local_shifts = True
    allow_local_deformation = True
    force_redo = [False, False, False]


class GridScan(Scan):
    """
    Holds information about one sub-scan that is placed in a grid.
    """

    # Position of the scan in the grid
    grid_pos = np.array([-1, -1, -1])







def add_edges(relations, direction):
    """
    Adds relevant edges to the relations tree given the directions that have been stitched (xdone, ydone, zdone) and
    the direction that should be stitched (0 = x, 1 = y, 2 = z)
    """

    for node1 in relations.nodes():
        x = node1.grid_pos[0]
        y = node1.grid_pos[1]
        z = node1.grid_pos[2]

        for node2 in relations.nodes():

            if direction == 0:
                xok = node2.grid_pos[0] == x + 1
            else:
                xok = node2.grid_pos[0] == x

            if direction == 1:
                yok = node2.grid_pos[1] == y + 1
            else:
                yok = node2.grid_pos[1] == y

            if direction == 2:
                zok = node2.grid_pos[2] == z + 1
            else:
                zok = node2.grid_pos[2] == z

            if xok and yok and zok:
                relations.add_edge(node1, node2)


def xyz_to_index(x, y, z, settings):
    """
    Converts image coordinates (image indices in coordinate directions) to image index (used in image file name).
    """

    w = settings.grid_width
    h = settings.grid_height
    d = settings.grid_depth

    n = 1
    if settings.xdir > 0:
        n = n + x
    else:
        n = n + (w - 1 - x)

    if settings.ydir > 0:
        n = n + y * w
    else:
        n = n + (h - 1 - y) * w

    if settings.zdir > 0:
        n = n + z * w * h
    else:
        n = n + (d - 1 - y) * w * h

    return n



def init_grid(settings, show_warnings):
    """
    Initializes stiching grid and returns the relations graph (containing no edges) and settings adjusted for binning.
    """

    # Create grid of all scans
    scans = [[[GridScan() for j in range(settings.grid_depth)] for i in range(settings.grid_height)] for j in range(settings.grid_width)]
    relations = nx.DiGraph()
    for z in range(settings.grid_depth):
        for y in range(settings.grid_height):
            for x in range(settings.grid_width):
                sc = GridScan()
                sc.index = xyz_to_index(x, y, z, settings)
                sc.grid_pos = np.array([x, y, z])
                sc.rec_file = settings.input_file_template.format(index=sc.index, x=x, y=y, z=z)
                sc.dimensions = get_image_size(sc.rec_file)
                sc.position = np.array([x * (sc.dimensions[0] - settings.overlap_x), y * (sc.dimensions[1] - settings.overlap_y), z * (sc.dimensions[2] - settings.overlap_z)])
                scans[x][y][z] = sc

                if sc.is_rec_ok():
                    relations.add_node(sc)

                    if show_warnings:
                        print(f"{sc.rec_file} at position {sc.position}, size = {sc.dimensions}")
                else:
                    if show_warnings:
                        print(f"Warning: Reconstruction {sc.rec_file} is missing.")

    # Adjust for binning
    settings2 = copy.deepcopy(settings)

    auto_binning(relations, settings2.binning)

    if is_use_cluster():
        if is_wait_for_jobs():
            wait_for_cluster_jobs()
        else:
            # TODO: It is not good to exit() here. This should be moved to the top-level function.
            print("Please run this program again when all the cluster jobs have finished.")
            exit()


    settings2.overlap_x /= settings2.binning
    settings2.overlap_y /= settings2.binning
    settings2.overlap_z /= settings2.binning

    settings2.point_spacing[0] = int(np.round(settings2.point_spacing[0] / settings2.binning))
    settings2.point_spacing[1] = int(np.round(settings2.point_spacing[1] / settings2.binning))
    settings2.point_spacing[2] = int(np.round(settings2.point_spacing[2] / settings2.binning))

    settings2.coarse_block_radius[0] = np.round(settings2.coarse_block_radius[0] / settings2.binning).astype(int)
    settings2.coarse_block_radius[1] = np.round(settings2.coarse_block_radius[1] / settings2.binning).astype(int)
    settings2.coarse_block_radius[2] = np.round(settings2.coarse_block_radius[2] / settings2.binning).astype(int)

    settings2.coarse_binning[0] = np.maximum(1, np.round(settings2.coarse_binning[0] / settings2.binning).astype(int))
    settings2.coarse_binning[1] = np.maximum(1, np.round(settings2.coarse_binning[1] / settings2.binning).astype(int))
    settings2.coarse_binning[2] = np.maximum(1, np.round(settings2.coarse_binning[2] / settings2.binning).astype(int))

    settings2.fine_block_radius[0] = np.round(settings2.fine_block_radius[0] / settings2.binning).astype(int)
    settings2.fine_block_radius[1] = np.round(settings2.fine_block_radius[1] / settings2.binning).astype(int)
    settings2.fine_block_radius[2] = np.round(settings2.fine_block_radius[2] / settings2.binning).astype(int)

    settings2.fine_binning[0] = np.maximum(1, np.round(settings2.fine_binning[0] / settings2.binning).astype(int))
    settings2.fine_binning[1] = np.maximum(1, np.round(settings2.fine_binning[1] / settings2.binning).astype(int))
    settings2.fine_binning[2] = np.maximum(1, np.round(settings2.fine_binning[2] / settings2.binning).astype(int))

    settings2.filter_threshold /= settings2.binning

    if settings2.binning != 1:
        settings2.sample_name = f"bin{settings2.binning}_{settings2.sample_name}"



    return relations, settings2







def stitch_all(settings):
    """
    Main entry point for stiching.
    Main method reads parameter file and passes control to this function.
    """

    must_wait = False
    for direction in [0, 1, 2]:
        print(f"Calculating displacement fields in dimension {direction}...")
        relations, adj_settings = init_grid(settings, direction == 0)
        add_edges(relations, direction)
        if not calculate_displacement_fields(adj_settings.sample_name, relations, adj_settings.point_spacing[direction], adj_settings.coarse_block_radius[direction], adj_settings.coarse_binning[direction], adj_settings.fine_block_radius[direction], adj_settings.fine_binning[direction], adj_settings.normalize_in_blockmatch, adj_settings.filter_threshold, adj_settings.force_redo[direction]):
            must_wait = True

    if must_wait:
        return False # We have to wait until jobs are finished etc.

    relations, adj_settings = init_grid(settings, False)
    add_edges(relations, 0)
    add_edges(relations, 1)
    add_edges(relations, 2)
    read_displacement_fields(adj_settings.sample_name, relations, adj_settings.allow_rotation)

    # Process all subgraphs separately
    run_stitching_for_all_connected_components(relations, adj_settings.sample_name, adj_settings.normalize_while_stitching, adj_settings.global_optimization, adj_settings.allow_rotation, adj_settings.allow_local_deformations)

    # TODO: Waiting here is a small hack... We should return flag and wait in main method, but this way we avoid reading displacement fields again.
    wait_for_cluster_jobs()

    return True







def main():
    """
    Main entry point for elastic stitching.
    Reads configuration from stitch_settings.txt and runs the stitching process.
    """

    config = configparser.ConfigParser()
    config.read('stitch_settings.txt')

    read_global_settings(config)


    # Stitching settings
    settings = StitchSettings()
    settings.sample_name = get(config, 'sample_name', 'sample')
    settings.input_file_template = get(config, 'input_file_template', 'tile{index}')

    # Output binning
    settings.binning = int(get(config, 'binning', '1'))


    # Counts of sub-scans in each coordinate direction.
    settings.grid_width = int(get(config, 'grid_width', 1))
    settings.grid_height = int(get(config, 'grid_height', 1))
    settings.grid_depth = int(get(config, 'grid_depth', 1))

    # Direction of increasing image number
    settings.xdir = int(get(config, 'xdir', 1))
    settings.ydir = int(get(config, 'ydir', 1))
    settings.zdir = int(get(config, 'zdir', 1))

    # Approximate overlap between scans in each coordinate direction.
    settings.overlap_x = int(get(config, 'overlap_x', 100))
    settings.overlap_y = int(get(config, 'overlap_y', 100))
    settings.overlap_z = int(get(config, 'overlap_z', 100))

    # Space between reference points in the deformation field in pixels
    # Different value can be used when stitching in different coordinate directions
    settings.point_spacing = [1, 1, 1]
    settings.point_spacing[0] = int(get(config, 'point_spacing_x_stitch', 20))
    settings.point_spacing[1] = int(get(config, 'point_spacing_y_stitch', 20))
    settings.point_spacing[2] = int(get(config, 'point_spacing_z_stitch', 20))

    # Size of calculation block in coarse matching
    # Different value can be used when stitching in different coordinate directions
    # Indicates also maximum local shift between images.
    sval = get(config, 'coarse_block_radius_x_stitch', '50')
    settings.coarse_block_radius[0] = from_string(sval)
    sval = get(config, 'coarse_block_radius_y_stitch', '50')
    settings.coarse_block_radius[1] = from_string(sval)
    sval = get(config, 'coarse_block_radius_z_stitch', '50')
    settings.coarse_block_radius[2] = from_string(sval)

    # Amount of binning for coarse matching
    settings.coarse_binning = [0, 0, 0]
    settings.coarse_binning[0] = int(get(config, 'coarse_binning_x_stitch', 4))
    settings.coarse_binning[1] = int(get(config, 'coarse_binning_y_stitch', 4))
    settings.coarse_binning[2] = int(get(config, 'coarse_binning_z_stitch', 4))

    # Size of calculation block in fine matching
    # Different value can be used when stitching in different coordinate directions
    # Indicates maximum shift on top shift determined in coarse matching.
    sval = get(config, 'fine_block_radius_x_stitch', '20')
    settings.fine_block_radius[0] = from_string(sval)
    sval = get(config, 'fine_block_radius_y_stitch', '20')
    settings.fine_block_radius[1] = from_string(sval)
    sval = get(config, 'fine_block_radius_z_stitch', '20')
    settings.fine_block_radius[2] = from_string(sval)

    # Amount of binning for fine matching
    settings.fine_binning = [0, 0, 0]
    settings.fine_binning[0] = int(get(config, 'fine_binning_x_stitch', 1))
    settings.fine_binning[1] = int(get(config, 'fine_binning_y_stitch', 1))
    settings.fine_binning[2] = int(get(config, 'fine_binning_z_stitch', 1))

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

    # Redo some direction?
    settings.force_redo[0] = get(config, 'redo_x', False) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']
    settings.force_redo[1] = get(config, 'redo_y', False) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']
    settings.force_redo[2] = get(config, 'redo_z', False) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Run
    while True:
        if stitch_all(settings):
            break
        # Set force_redo to false so that we don't calculate the same things again and again!
        settings.force_redo = [False, False, False]

        if not is_use_cluster():
            break

        if is_wait_for_jobs():
            wait_for_cluster_jobs()
        else:
            print("Run this script again when all the results of the jobs are available.")
            exit()



if __name__ == "__main__":
    start = time.time()
    main()
    elapsed = time.time() - start
    print(f"Operation took {elapsed} seconds.")
