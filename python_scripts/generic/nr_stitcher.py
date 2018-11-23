#!/usr/bin/env python3

import os.path
import subprocess
import struct
import configparser
import time
import networkx as nx
import numpy as np
from pyquaternion import Quaternion


# Path to pi2 program.
pi_path = "."

# Indicates if cluster should be used for compute-intensive tasks
use_cluster = False

# Partition to use if calculations are performed on cluster
cluster_partition = ""

# Init commands that are run before pi2 program when a cluster job starts.
cluster_job_init_commands = ""

# Maximum image dimension to use while stitching. 2500 corresponds to ~120 GB memory requirement.
max_block_size = 2500

# Should the script wait until all cluster jobs are finished?
wait_for_jobs = True

# Extra parameters for sbatch
sbatch_extra_params = ""

# List of submitted but non-finished jobs
submitted_jobs = []


class StitchSettings:
    """
    Stores values of all the stitching settings.
    Values are read from input file or set to default values.
    """
    sample_name = ''
    width = 0
    height = 0
    depth = 0
    data_type = ''
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
    allow_rotation = True
    force_redo = [False, False, False]


class Scan:
    """
    Holds information about one sub-scan.
    """

    # File that contains the image data
    rec_file = ''

    # Index of the scan (first scan = 1, second scan = 2, ...)
    index = -1

    # Position of the scan in the grid
    grid_x = -1
    grid_y = -1
    grid_z = -1

    # Size of the image
    width = 0
    height = 0
    depth = 0

    # Similarity transformation (rotation matrix, scaling factor, translation)
    R = np.eye(3)
    a = 1
    c = np.zeros((3, 1))

    # Indicates if the similarity transformation has been set
    transformation_set = False

    # Name of file where the total transformation of this scan has been stored
    transformation_file = ''

    # Filename prefix used while storing world to local transformation grid.
    world_to_local_prefix = ''

    # Approximate position of the first pixel of the image in world coordinates.
    world_x = 0
    world_y = 0
    world_z = 0

    # Normalization factor for gray values
    norm_fact = 0

    def is_rec_ok(self):
        """
        Tests if reconstructed image file exists.
        """
        return os.path.isfile(self.rec_file)

    def is_at(self, pos):
        """
        Tests if this scan it at given grid position.
        """
        return (pos[0] == self.grid_x) and (pos[1] == self.grid_y) and (pos[2] == self.grid_z)




def run_pi2(pi_script, output_prefix):
    """
    Runs pi2 job either locally or on cluster.
    """
    global submitted_jobs

    if use_cluster:

        job_cmdline = ""
        if len(cluster_job_init_commands) > 0:
            job_cmdline = cluster_job_init_commands + ";"
        job_cmdline = job_cmdline + f"{pi_path}/pi2 '{pi_script}'"

        sbatch_params = f"--job-name=stitch --partition={cluster_partition} --output={output_prefix}-out.txt --error={output_prefix}-err.txt {sbatch_extra_params} --wrap=\"{job_cmdline}\""

        # For testing
        #cmd = "echo"
        # Real command line
        cmd = "sbatch"

        res = subprocess.check_output(cmd + " " + sbatch_params, shell=True)
        job_id = int(res.split()[-1])
        submitted_jobs.append(job_id)

    else:
        print(f"Executing pi2...")
        subprocess.call([pi_path + "/pi2", pi_script])





def displacement_file_prefix(sample_name, scan1, scan2):
    """
    Gets name prefix for displacement field files for given two scans.
    """
    return f"{sample_name}_{scan1.index}-{scan2.index}"






def is_displacement_ok(sample_name, scan1, scan2):
    """
    Checks if displacement field has been calculated between the two given scans.

    Returns true if displacement field has been calculated, and false otherwise.
    """
    filename = "%s_refpoints.txt" % displacement_file_prefix(sample_name, scan1, scan2)
    return os.path.isfile(filename)




def is_filtered_displacement_ok(sample_name, scan1, scan2):
    """
    Checks if filtered displacement field has been calculated between the two given scans.

    Returns true if filtered displacement field has been calculated, and false otherwise.
    """
    filename = "%s_filtered_refpoints.txt" % displacement_file_prefix(sample_name, scan1, scan2)
    return os.path.isfile(filename)




def calculate_displacement_field(sample_name, scan1, scan2, image_data_type, point_spacing, coarse_block_radius, coarse_binning, fine_block_radius, fine_binning, normalize, filter_threshold):
    """
    Calculates displacement field between scan1 and scan2 locally or submits displacement field calculation job to cluster.
    Calculates also filtered displacement field.
    """

    # Calculate initial minimum and maximum for the reference grid.
    xmin = max([0, scan2.world_x - scan1.world_x])
    ymin = max([0, scan2.world_y - scan1.world_y])
    zmin = max([0, scan2.world_z - scan1.world_z])
    xmax = min([scan1.width, scan2.world_x - scan1.world_x + scan2.width])
    ymax = min([scan1.height, scan2.world_y - scan1.world_y + scan2.height])
    zmax = min([scan1.depth, scan2.world_z - scan1.world_z + scan2.depth])
    shiftx = -(scan2.world_x - scan1.world_x)
    shifty = -(scan2.world_y - scan1.world_y)
    shiftz = -(scan2.world_z - scan1.world_z)
    #shiftx = -((scan2.world_x + scan2.width / 2) - (scan1.world_x + scan1.width / 2));

    file_prefix = displacement_file_prefix(sample_name, scan1, scan2)

    # This script loads or maps the .raw files and runs normal block match. Currently processing mapped files is very slow.
    #params = (f"echo;"
    #          f"mapraw(ref, {image_data_type}, {scan1.rec_file});"
    #          f"mapraw(def, {image_data_type}, {scan2.rec_file});"
    #          f"blockmatch(ref, def, {xmin}, {xmax}, {point_spacing}, {ymin}, {ymax}, {point_spacing}, {zmin}, {zmax}, {point_spacing}, {shiftx}, {shifty}, {shiftz}, {file_prefix}, {block_radius});"
    #         )
    # This script uses ad-hoc implementation of blockmatch to save memory.
    params = (f"echo;"
              f"blockmatchmemsave({scan1.rec_file}, {scan2.rec_file}, {xmin}, {xmax}, {point_spacing}, {ymin}, {ymax}, {point_spacing}, {zmin}, {zmax}, {point_spacing}, [{shiftx}, {shifty}, {shiftz}], {file_prefix}, {coarse_block_radius}, {coarse_binning}, {fine_block_radius}, {fine_binning}, {normalize});"
              f"filterdisplacements({file_prefix}, {filter_threshold});"
             )

    run_pi2(params, file_prefix)



def filter_displacement_field(sample_name, scan1, scan2, filter_threshold):
    """
    Calculates filtered displacement field from non-filtered one.
    This is mainly used if the user deleted filtered displacement fields (but not unfiltered ones) to tune filtering parameters.
    """

    file_prefix = displacement_file_prefix(sample_name, scan1, scan2)

    params = (f"echo;"
              f"filterdisplacements({file_prefix}, {filter_threshold});"
             )

    run_pi2(params, f"{file_prefix}_filter")




def read_displacement_field(sample_name, scan1, scan2):
    """
    Reads displacement field between scan1 and scan 2 from files saved by pi2 program.

    Returns (x0, y0, z0, x1, y1, z1, gof, norm_fact), where
    x0, y0, z0 contain points in reference image, and
    x1, y1, z1 contain corresponding points in deformed image.
    norm_fact is value that must be added to deformed image to make mean of deformed image and reference image the same.
    """

    # Read reference point definition file
    refpoints_file = displacement_file_prefix(sample_name, scan1, scan2) + "_filtered_refpoints.txt"
    with open(refpoints_file, 'r') as f:
        lines = f.readlines()
        parts = lines[0].split(",")
        xmin = int(parts[0])
        xmax = int(parts[1])
        xstep = int(parts[2])
        parts = lines[1].split(",")
        ymin = int(parts[0])
        ymax = int(parts[1])
        ystep = int(parts[2])
        parts = lines[2].split(",")
        zmin = int(parts[0])
        zmax = int(parts[1])
        zstep = int(parts[2])
        norm_fact = float(lines[3])


    xcount = int(np.floor((xmax - xmin) / xstep)) + 1
    ycount = int(np.floor((ymax - ymin) / ystep)) + 1
    zcount = int(np.floor((zmax - zmin) / zstep)) + 1

    total_count = xcount * ycount * zcount

    # Read goodness of fit
    gof_file = displacement_file_prefix(sample_name, scan1, scan2) + f"_filtered_gof_{xcount}x{ycount}x{zcount}.raw"
    gof_data = read_floats(gof_file, total_count)
    defpoints_file = displacement_file_prefix(sample_name, scan1, scan2) + f"_filtered_defpoints_{xcount}x{ycount}x{zcount}.raw"
    defpoints_data = read_doubles(defpoints_file, total_count * 3)

    # Create reference points
    x = np.zeros([zcount, ycount, xcount])
    y = np.zeros([zcount, ycount, xcount])
    z = np.zeros([zcount, ycount, xcount])
    x1 = np.zeros([zcount, ycount, xcount])
    y1 = np.zeros([zcount, ycount, xcount])
    z1 = np.zeros([zcount, ycount, xcount])
    gof = np.zeros([zcount, ycount, xcount])
    for zi in range(zcount):
        for yi in range(ycount):
            for xi in range(xcount):
                x[zi, yi, xi] = xmin + xi * xstep
                y[zi, yi, xi] = ymin + yi * ystep
                z[zi, yi, xi] = zmin + zi * zstep

                ind = zi * xcount * ycount + yi * xcount + xi
                gof[zi, yi, xi] = gof_data[ind]
                x1[zi, yi, xi] = defpoints_data[3 * ind]
                y1[zi, yi, xi] = defpoints_data[3 * ind + 1]
                z1[zi, yi, xi] = defpoints_data[3 * ind + 2]



    return x, y, z, x1, y1, z1, gof, norm_fact



def read_floats(filename, count):
    """
    Reads count 32-bit float values from given file. The values must be in native byte order.

    Returns values that were read.
    """

    data = []
    with open(filename, 'rb') as f:
        for n in range(count):
            val = struct.unpack('f', f.read(4))
            data.append(val[0])

    return data



def read_doubles(filename, count):
    """
    Reads count 64-bit float values from given file. The values must be in native byte order.

    Returns values that were read.
    """

    data = []
    with open(filename, 'rb') as f:
        for n in range(count):
            val = struct.unpack('d', f.read(8))
            data.append(val[0])

    return data




def rigid_body_transformation_from_displacement_field(x_orig, y_orig, z_orig, u_orig, v_orig, w_orig, gof_orig, allow_rotation):
    """
    Calculates similarity transformation that matches the given displacement field as well as possible.
    The transformation is
    X1 = a * R * (X2 - c) + c + delta


    Returns (delta, c, a, R)


    Algorithm is from
    Sorkine-Hornung - Least-Squares Rigid Motion Using SVD
    but it is probably unpublished.
    Other sources are
    https://stackoverflow.com/questions/13432805/finding-translation-and-scale-on-two-sets-of-points-to-get-least-square-error-in
    https://gist.github.com/nh2/bc4e2981b0e213fefd4aaa33edfb3893#file-rigid-transform-with-scale-py-L20-L39
    https://ch.mathworks.com/matlabcentral/fileexchange/25746-kabsch-algorithm?focused=3790280&tab=function
    Umeyama - Least-Squares Estimation of Transformation Parameters Between Two Point Patterns

    TODO: I'm not completely sure if the algorithm actually works as in the comment above, or if it returns the inverse of the transformation!
    """

    # Remove bad values
    good_ind = ~(np.isnan(u_orig) | (gof_orig <= 0))

    x = x_orig[good_ind]
    y = y_orig[good_ind]
    z = z_orig[good_ind]
    u = u_orig[good_ind]
    v = v_orig[good_ind]
    w = w_orig[good_ind]
    gof = gof_orig[good_ind]

    # Build data matrices
    U = np.array([u.flatten(), v.flatten(), w.flatten()]).T
    P = np.array([x.flatten(), y.flatten(), z.flatten()]).T
    Q = P + U
    W = gof.flatten()
    W = W / sum(W)

    n, dim = P.shape

    # Rotation center and translation
    c = np.reshape(np.average(P, axis=0, weights=W), (1, 3))
    delta = np.reshape(np.average(Q, axis=0, weights=W), (1, 3)) - c

    # Centered data
    centeredP = P - c
    centeredQ = Q - (delta + c)

    # Weighted covariance matrix
    tmp = np.zeros_like(centeredP)
    for i in range(len(W)):
        tmp[i, :] = W[i] * centeredP[i, :]
    C = np.dot(np.transpose(centeredP), centeredQ) / n

    V, S, Ww = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(Ww)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    R = np.dot(V, Ww)
    R = R.T

    # Scaling factor
    # TODO: This scaling factor calculation using weighted variance does not seem to work.
    #varP = np.average(centeredP**2, axis=0, weights=W)
    #a = 1/varP.sum() * np.sum(S)
    #a = np.sum(S) / np.average(centeredP**2, axis=0).sum()
    # Disable scaling altogether
    a = 1

    # Reshape the vectors for output
    delta = delta.T
    c = c.T

    # Disable rotation if requested to do so
    if not allow_rotation:
        R = np.eye(3, 3)

    return delta, c, a, R





def calculate_displacement_fields(sample_name, relations, data_type, point_spacing, coarse_block_radius, coarse_binning, fine_block_radius, fine_binning, normalize, filter_threshold, force_redo):
    """
    Calculate displacement fields and load results to relations network.
    Returns true if all displacement fields have been calculated and read, and false otherwise.
    """

    # Determine which pairwise displacement fields should be calculated (i.e. don't calculate for non-existing images or if the calculation has been made already)
    # and start calculations or submit jobs
    jobs_started = 0
    for edge in relations.edges(data=True):
        scan1 = edge[2]["scan1"]
        scan2 = edge[2]["scan2"]

        if force_redo or (not is_displacement_ok(sample_name, scan1, scan2)):
            print("Calculate displacement field %i -> %i" % (scan1.index, scan2.index))
            calculate_displacement_field(sample_name, scan1, scan2, data_type, point_spacing, coarse_block_radius, coarse_binning, fine_block_radius, fine_binning, normalize, filter_threshold)
            jobs_started = jobs_started + 1


    # If using a cluster, tell the user to wait until the results of the calculations are available.
    if (jobs_started > 0) and use_cluster:
        return False

    # Check that filtered displacement fields are also ok.
    jobs_started = 0
    for edge in relations.edges(data=True):
        scan1 = edge[2]["scan1"]
        scan2 = edge[2]["scan2"]

        if not is_filtered_displacement_ok(sample_name, scan1, scan2):
            print("Filter displacement field %i -> %i" % (scan1.index, scan2.index))
            filter_displacement_field(sample_name, scan1, scan2, filter_threshold)
            jobs_started = jobs_started + 1

    if (jobs_started > 0) and use_cluster:
        return False

    return True





def read_displacement_fields(sample_name, relations, allow_rotation):
    """
    Reads all displacement fields and determines similarity transformation for each field.
    """

    print("Reading displacement fields to calculate average similarity transformations...")

    # Read displacement fields (including goodness-of-fit information and assig weights for the edges of the relations graph)
    counter = 0
    for edge in relations.edges(data=True):
        scan1 = edge[2]["scan1"]
        scan2 = edge[2]["scan2"]

        x, y, z, x1, y1, z1, gof, norm_fact = read_displacement_field(sample_name, scan1, scan2)
        u = x1 - x
        v = y1 - y
        w = z1 - z

        delta, c, a, R = rigid_body_transformation_from_displacement_field(x, y, z, u, v, w, gof, allow_rotation)

        # TODO: It looks like rigid_body_transformation_from_displacement_field gives inverse of what we want. I have not dug into details yet.
        a = 1 / a
        R = R.T

        cn = delta + c
        #deltan + cn = c ==>
        #deltan = c - cn = c - delta - c = -delta
        deltan = delta

        delta = deltan
        c = cn

        edge[2]["weight"] = 1.0 - np.mean(gof)
        edge[2]["delta"] = delta
        edge[2]["c"] = c
        edge[2]["a"] = a
        edge[2]["R"] = R
        edge[2]["norm_fact"] = norm_fact

        counter = counter + 1
        print(f"{counter} / {relations.number_of_edges()}\r", end="")




def convert_to_world(node_pos, child_pos, tree):
    """
    Determines similarity transformation of child node in world coordinates, based on world-coordinate similarity transformation of parent node and parent to child transformation.
    """

    parent_scan = tree[node_pos][child_pos]["scan1"]
    child_scan = tree[node_pos][child_pos]["scan2"]
    delta = tree[node_pos][child_pos]["delta"]
    c = tree[node_pos][child_pos]["c"]
    a = tree[node_pos][child_pos]["a"]
    R = tree[node_pos][child_pos]["R"]
    norm_fact = tree[node_pos][child_pos]["norm_fact"]

    delta = -delta
    deltadot = c + delta

    if not child_scan.is_at(child_pos):
        # We are moving from scan2 to scan1
        raise RuntimeError("Invalid scan1 position.")

    if not parent_scan.transformation_set:
        raise RuntimeError(f"Parent scan at ({node_pos[0]}, {node_pos[1]}, {node_pos[2]}) has no transformation set yet, and its transformation is required for child scan at ({child_pos[0]}, {child_pos[1]}, {child_pos[2]}).")

    # Convert transformation from T(n-1, n) to T(0, n)
    wa = parent_scan.a * a
    wR = np.matmul(parent_scan.R, R)
    wc = 1 / a * np.matmul(np.linalg.inv(R), parent_scan.c - deltadot) + c
    wnorm_fact = parent_scan.norm_fact + norm_fact

    return wa, wR, wc, wnorm_fact




def determine_average_transformations(relations):
    """
    Sets node's similarity transformation to average of transformations between scans connected to it,
    and converted to world coordinates.
    """


    ordered = list(nx.topological_sort(relations))

    for node_pos in ordered:

        scan = nx.get_node_attributes(relations, 'scan')[node_pos]

        # Process all neighbours that have been scanned earlier
        count = 0
        atot = 0
        ctot = 0
        qtot = Quaternion(axis=[1, 0, 0], angle=0)
        normtot = 0

        parent_scans = relations.predecessors(node_pos)
        for parent_pos in parent_scans:
            #parent_scan = nx.get_node_attributes(relations, 'scan')[parent_pos]
            # print("\tProcessing neighbour {0}, {1}, {2}".format(nb_pos[0], nb_pos[1], nb_pos[2]))

            a, R, c, nf = convert_to_world(parent_pos, node_pos, relations)
            atot = atot + a
            ctot = ctot + c
            # TODO: This quaternion averaging is just approximation that might work if all the rotations are almost the same...
            qtot = qtot + Quaternion(matrix=R)
            normtot = normtot + nf
            count = count + 1


        if count > 0:
            scan.a = atot / count
            scan.c = ctot / count
            qtot = qtot / count
            qtot = qtot.normalised
            scan.R = qtot.rotation_matrix
            scan.norm_fact = normtot / count


        scan.transformation_set = True






def determine_bounds(comp):
    """
    Determines approximate bounds for the stitched image.
    Returns xmin, ymin, zmin, xmax, ymax, zmax
    """


    xmin = float('inf')
    ymin = float('inf')
    zmin = float('inf')
    xmax = -float('inf')
    ymax = -float('inf')
    zmax = -float('inf')
    for node_pos in comp:
        scan = nx.get_node_attributes(comp, "scan")[node_pos]
        a = scan.a
        R = scan.R
        c = scan.c
        width = scan.width
        height = scan.height
        depth = scan.depth

        # Corners of the image 1 in image n coordinates
        points = np.array([[0, 0, 0], [width, 0, 0], [width, height, 0], [0, height, 0],
                           [0, 0, depth], [width, 0, depth], [width, height, depth], [0, height, depth]])

        # Convert from coordinates of image n to coordinates of image 1 using the similarity transformation
        for n in range(0, points.shape[0]):
            p = a * np.matmul(R, points[n, :] - c.flatten())

            # Determine bounds
            xmin = min(xmin, p[0])
            ymin = min(ymin, p[1])
            zmin = min(zmin, p[2])

            xmax = max(xmax, p[0])
            ymax = max(ymax, p[1])
            zmax = max(zmax, p[2])


    return xmin, ymin, zmin, xmax, ymax, zmax




def find_first_node(comp):
    """
    Finds the root scan in the given tree.
    """

    ordered = list(nx.topological_sort(comp))
    return ordered[0]



def save_transformation(sample_name, node_pos, relations):
    """
    Saves world to local similarity tranformation of node at given position.
    """

    scan = nx.get_node_attributes(relations, 'scan')[node_pos]
    a = scan.a
    R = scan.R
    c = scan.c
    norm_fact = scan.norm_fact

    scan_name = scan.rec_file[0:scan.rec_file.rfind('_')]
    file = f"{scan_name}_transformation.txt"
    with open(file, 'wb') as f:
        # Write world to local similarity transformation of this node
        np.savetxt(f, np.array([a]))
        np.savetxt(f, c)
        np.savetxt(f, R)
        np.savetxt(f, np.array([norm_fact]))

        # Save information about parent node indices
        parents = list(relations.predecessors(node_pos))

        # Write count of parent nodes
        np.savetxt(f, np.array([len(parents)]))
        for parent in parents:
            parent_scan = nx.get_node_attributes(relations, 'scan')[parent]
            parent_scan_name = parent_scan.rec_file[0:parent_scan.rec_file.rfind('_')]

            prefix = displacement_file_prefix(sample_name, parent_scan, scan)

            np.savetxt(f, np.array([f"{parent_scan_name}_world_to_local"]), fmt="%s")
            np.savetxt(f, np.array([f"{prefix}_filtered"]), fmt="%s")


    scan.transformation_file = file




def find_roots(tree):
    """
    Creates list of nodes with no predecessors.
    """

    roots = []

    for node_pos in tree:
        if len(list(tree.predecessors(node_pos))) <= 0:
            roots.append(node_pos)

    return roots


def calculate_world_to_local(tree):
    """
    Calculates world to local transformations for all nodes in the given tree, starting from the root nodes that have no incoming connections.
    """

    curr_level = find_roots(tree)
    next_level = []

    while len(curr_level) > 0:
        while len(curr_level) > 0:
            curr_node = curr_level.pop()

            scan = nx.get_node_attributes(tree, 'scan')[curr_node]
            scan_name = scan.rec_file[0:scan.rec_file.rfind('_')]

            # Calculate world to local grid transform using pi
            scan.world_to_local_prefix = f"{scan_name}_world_to_local"
            script = (f"echo;"
                      f"determine_world_to_local({scan.transformation_file}, [{scan.width}, {scan.height}, {scan.depth}], {scan.world_to_local_prefix});"
                     )
            run_pi2(script, scan.world_to_local_prefix)

            # Insert all successors to the next level
            children = tree.successors(curr_node)
            for child in children:
                if next_level.count(child) <= 0:
                    next_level.append(child)

        wait_for_cluster_jobs()
        curr_level, next_level = next_level, curr_level




def run_stitching(comp, sample_name, data_type, normalize):
    """
    Prepares and runs pi2 stitching process for connected component 'comp' of scan relations tree 'tree'.
    - determines final world to image transformations
    - saves them
    - determines bounds of the stitched image
    - makes index file for pi2, runs pi2 stitching in blocks
    """

    if not nx.is_directed_acyclic_graph(comp):
        raise RuntimeError("Connected component is not directed acyclic graph.")

    # Determine similarity transformations
    determine_average_transformations(comp)

    # Save the transformations
    for node_pos in comp.nodes:
        save_transformation(sample_name, node_pos, comp)

    # Calculate world to local grid transformations
    print("Calculating world to local transformation for each image...")
    calculate_world_to_local(comp)


    print("Stitching...")

    # Determine approximate bounds for the stitched image
    minx, miny, minz, maxx, maxy, maxz = determine_bounds(comp)

    out_width = int(round(maxx - minx))
    out_height = int(round(maxy - miny))
    out_depth = int(round(maxz - minz))
    minx = int(round(minx))
    miny = int(round(miny))
    minz = int(round(minz))

    start_node_pos = find_first_node(comp)
    out_template = f"{sample_name}_{start_node_pos[0]}_{start_node_pos[1]}_{start_node_pos[2]}"
    out_file = f"{out_template}_{out_width}x{out_height}x{out_depth}.raw"

    # Make index file
    index_file = out_template + "_index.txt"
    with open(index_file, 'w') as f:
        for node_pos in comp:
            scan = nx.get_node_attributes(comp, 'scan')[node_pos]
            f.write(f"{scan.rec_file}\n")
            f.write(f"{scan.world_to_local_prefix}\n")


    # Cut big sample to small regions and process them separately to save memory
    block_size = int(max_block_size)

    jobs_started = 0
    for zstart in range(minz, out_depth + minz, block_size):
        for ystart in range(miny, out_height + miny, block_size):
            for xstart in range(minx, out_width + minx, block_size):

                curr_width = min(block_size, out_width - xstart)
                curr_height = min(block_size, out_height - ystart)
                curr_depth = min(block_size, out_depth - zstart)

                pi_script = (f"echo;"
                             f"newimage(outimg, {data_type});"
                             f"stitch_ver2(outimg, {index_file}, {xstart}, {ystart}, {zstart}, {curr_width}, {curr_height}, {curr_depth}, {normalize});"
                             f"writerawblock(outimg, {out_file}, {xstart - minx}, {ystart - miny}, {zstart - minz}, {out_width}, {out_height}, {out_depth});"
                            )

                run_pi2(pi_script, f"{out_template}_{jobs_started}")
                jobs_started = jobs_started + 1

    # This version does the stitching without cutting the sample into blocks.
    # NOTE: Uses old version of stitching function.
    #pi_script = (f"echo;"
    #             f"newimage(outimg, {data_type});"
    #             f"stitch(outimg, {index_file}, {minx}, {miny}, {minz}, {out_width}, {out_height}, {out_depth});"
    #             f"writeraw(outimg, {out_template});"
    #            )
    #subprocess.call([pi_path + "/pi2", pi_script])

    return jobs_started


def run_stitching_for_all_connected_components(relations, sample_name, data_type, normalize):
    """
    Calls run_stitching for each connected component in relations network.
    """

    jobs_started = 0
    comps = (relations.subgraph(c) for c in nx.weakly_connected_components(relations))
    for comp in comps:
        jobs_started = jobs_started + run_stitching(comp, sample_name, data_type, normalize)

    if (jobs_started > 0) and use_cluster:
        return False

    return True






def add_edges(relations, direction):
    """
    Adds relevant edges to the relations tree given the directions that have been stitched (xdone, ydone, zdone) and
    the direction that should be stitched (0 = x, 1 = y, 2 = z)
    """

    for node1 in relations.nodes(True):
        node1_pos = node1[0]
        node1_scan = node1[1]["scan"]
        x = node1_scan.grid_x
        y = node1_scan.grid_y
        z = node1_scan.grid_z

        for node2 in relations.nodes(True):
            node2_pos = node2[0]
            node2_scan = node2[1]['scan']

            if direction == 0:
                xok = node2_pos[0] == x + 1
            else:
                xok = node2_pos[0] == x

            if direction == 1:
                yok = node2_pos[1] == y + 1
            else:
                yok = node2_pos[1] == y

            if direction == 2:
                zok = node2_pos[2] == z + 1
            else:
                zok = node2_pos[2] == z

            if xok and yok and zok:
                relations.add_edge(node1_pos, node2_pos, scan1=node1_scan, scan2=node2_scan)


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
    Initializes stiching grid and returns the relations graph (containing no edges)
    """
    filename_template = "%s%d_%dx%dx%d.raw"

    # Create grid of all scans
    scans = [[[Scan() for j in range(settings.grid_depth)] for i in range(settings.grid_height)] for j in range(settings.grid_width)]
    relations = nx.DiGraph()
    for z in range(settings.grid_depth):
        for y in range(settings.grid_height):
            for x in range(settings.grid_width):
                sc = Scan()
                sc.index = xyz_to_index(x, y, z, settings)
                sc.grid_x = x
                sc.grid_y = y
                sc.grid_z = z
                sc.world_x = sc.grid_x * (settings.width - settings.overlap_x)
                sc.world_y = sc.grid_y * (settings.height - settings.overlap_y)
                sc.world_z = sc.grid_z * (settings.depth - settings.overlap_z)
                sc.width = settings.width
                sc.height = settings.height
                sc.depth = settings.depth
                sc.rec_file = filename_template % (settings.sample_name, sc.index, settings.width, settings.height, settings.depth)
                scans[x][y][z] = sc

                if sc.is_rec_ok():
                    relations.add_node((x, y, z), scan=sc)
                else:
                    if show_warnings:
                        print(f"Warning: Reconstruction {sc.rec_file} is missing.")

    return relations



def get(config, key, default):
    """
    Gets a value from configuration dictionary if the key is set. If not, returns default.
    """

    if 'stitch' in config:
        if key in config['stitch']:
            return config['stitch'][key]

    return default



def stitch_all(settings):
    """
    Main entry point for stiching.
    Main method reads parameter file and passes control to this function.
    """

    must_wait = False
    for direction in [0, 1, 2]:
        print(f"Calculating displacement fields in dimension {direction}...")
        relations = init_grid(settings, direction == 0)
        add_edges(relations, direction)
        if not calculate_displacement_fields(settings.sample_name, relations, settings.data_type, settings.point_spacing[direction], settings.coarse_block_radius[direction], settings.coarse_binning[direction], settings.fine_block_radius[direction], settings.fine_binning[direction], settings.normalize_in_blockmatch, settings.filter_threshold, settings.force_redo[direction]):
            must_wait = True

    if must_wait:
        return False # We have to wait until jobs are finished etc.

    relations = init_grid(settings, False)
    add_edges(relations, 0)
    add_edges(relations, 1)
    add_edges(relations, 2)
    read_displacement_fields(settings.sample_name, relations, settings.allow_rotation)

    # Process all subgraphs separately
    run_stitching_for_all_connected_components(relations, settings.sample_name, settings.data_type, settings.normalize_while_stitching)

    # TODO: Waiting here is a small hack... We should return flag and wait in main method, but this way we avoid reading displacement fields again.
    wait_for_cluster_jobs()

    return True




def is_job_running(job_id):
    """
    Tests if a cluster job with given job id is not finished (running or waiting to be run).
    """

    cmd = "squeue"
    params = f"--noheader --jobs={job_id}"
    res = subprocess.run(cmd + " " + params, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
    res = res.decode("utf-8")
    res = res.strip()
    if res.startswith(str(job_id)):
        return True

    return False



def wait_for_cluster_jobs():
    """
    Waits until all jobs in the submitted jobs array are finished.
    """

    global submitted_jobs

    if not use_cluster:
        return

    if len(submitted_jobs) <= 0:
        return

    print(f"Waiting for {len(submitted_jobs)} cluster jobs to finish...")

    time.sleep(2)

    while True:
        if len(submitted_jobs) <= 0:
            return

        job_id = submitted_jobs.pop(0)

        while is_job_running(job_id):
            time.sleep(2)






def main():
    """
    Main entry point for elastic stitching.
    Reads configuration from stitch_settings.txt and runs stitching process.
    """

    config = configparser.ConfigParser()
    config.read('stitch_settings.txt')

    # Global settings
    global pi_path
    global use_cluster
    global cluster_partition
    global cluster_job_init_commands
    global max_block_size
    global wait_for_jobs
    global sbatch_extra_params

    # Path to pi2 program. By default ./
    pi_path = get(config, 'pipath', pi_path)

    # Should we wait for the cluster jobs to complete?
    wait_for_jobs = get(config, 'wait_for_jobs', wait_for_jobs) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Indicates if the calculations should be performed on a cluster
    use_cluster = get(config, 'use_cluster', use_cluster) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Cluster partition to use (if cluster calculation is enabled). Empty value means default partition.
    cluster_partition = get(config, 'cluster_partition', cluster_partition)

    # Job init commands
    cluster_job_init_commands = get(config, 'cluster_job_init_commands', cluster_job_init_commands)

    # Maximum stitching block size
    max_block_size = get(config, 'max_block_size', max_block_size)

	# Extra parameters for sbatch command
    sbatch_extra_params = get(config, 'cluster_extra_params', sbatch_extra_params)


    # Stitching settings
    settings = StitchSettings()
    settings.sample_name = get(config, 'sample_name', 'sample')
    settings.width = int(get(config, 'width', 100))
    settings.height = int(get(config, 'height', 100))
    settings.depth = int(get(config, 'depth', 100))

    # Pixel data type. Can be uint8, uint16, or float32
    settings.data_type = get(config, 'data_type', 'uint16')

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
    # Indicates also maximum local shift between images
    settings.coarse_block_radius = ['', '', '']
    settings.coarse_block_radius[0] = get(config, 'coarse_block_radius_x_stitch', '50')
    settings.coarse_block_radius[1] = get(config, 'coarse_block_radius_y_stitch', '50')
    settings.coarse_block_radius[2] = get(config, 'coarse_block_radius_z_stitch', '50')

    # Amount of binning for coarse matching
    settings.coarse_binning = [0, 0, 0]
    settings.coarse_binning[0] = int(get(config, 'coarse_binning_x_stitch', 4))
    settings.coarse_binning[1] = int(get(config, 'coarse_binning_y_stitch', 4))
    settings.coarse_binning[2] = int(get(config, 'coarse_binning_z_stitch', 4))

    # Size of calculation block in fine matching
    # Different value can be used when stitching in different coordinate directions
    # Indicates maximum shift on top shift determined in coarse matching.
    settings.fine_block_radius = ['', '', '']
    settings.fine_block_radius[0] = get(config, 'fine_block_radius_x_stitch', '20')
    settings.fine_block_radius[1] = get(config, 'fine_block_radius_y_stitch', '20')
    settings.fine_block_radius[2] = get(config, 'fine_block_radius_z_stitch', '20')

    # Amount of binning for fine matching
    settings.fine_binning = [0, 0, 0]
    settings.fine_binning[0] = int(get(config, 'fine_binning_x_stitch', 1))
    settings.fine_binning[1] = int(get(config, 'fine_binning_y_stitch', 1))
    settings.fine_binning[2] = int(get(config, 'fine_binning_z_stitch', 1))

    # Normalization flags
    settings.normalize_in_blockmatch = get(config, 'normalize_in_blockmatch', True) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']
    settings.normalize_while_stitching = get(config, 'normalize_while_stitching', True) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Allow rotational transformation
    settings.allow_rotation = get(config, 'allow_rotation', True) in [True, 'true', 'True', 'TRUE', '1', 't', 'y', 'yes']

    # Threshold for displacement filtering
    settings.filter_threshold = get(config, 'displacement_filter_threshold', '3')

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

        if not use_cluster:
            break

        if wait_for_jobs:
            wait_for_cluster_jobs()
        else:
            print("Run this script again when all the results of the jobs are available.")
            exit



if __name__ == "__main__":
    start = time.time()
    main()
    elapsed = time.time() - start
    print(f"Operation took {elapsed} seconds.")
