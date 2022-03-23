

import os.path
import subprocess
import struct
import networkx as nx
from pyquaternion import Quaternion
import time
import math
import glob
import numpy as np
import scipy
import scipy.optimize

from pi2py2 import *
pi = Pi2()


# Path to pi2 program.
pi_path = "."

# Indicates which cluster system should be used for compute-intensive tasks, if any
cluster = ""

# Maximum image dimension to use while stitching. 2500 corresponds to ~120 GB memory requirement.
max_block_size = 2500



def is_use_cluster():
    return cluster != ""





def get(config, key, default):
    """
    Gets a value from configuration dictionary if the key is set. If not, returns default.
    """

    if 'stitch' in config:
        if key in config['stitch']:
            return config['stitch'][key]

    return default



def read_global_settings(config):
    """
    Reads cluster and pi2 related settings from a ConfigParser.
    """

    # Global settings
    global pi_path
    global cluster
    global max_block_size

    # Path to pi2 program. By default directory of running script.
    pi_path = os.path.dirname(os.path.realpath(__file__))
    pi_path = get(config, 'pipath', pi_path)

    # Indicates if the calculations should be performed on a cluster
    cluster = get(config, 'cluster', cluster)
    if cluster == "None" or cluster == "none":
        cluster = ""

    if cluster != "":
        pi.distribute(cluster)

    # Maximum stitching block size
    max_block_size = get(config, 'max_block_size', max_block_size)




def from_string(str):

    if len(str) > 0 and str[0] == '[':
        str = str[1:]

    if len(str) > 0 and str[-1] == ']':
        str = str[:-1]

    return np.fromstring(str, dtype=int, sep=',')



def run_pi2_locally(pi_script):
    """
    Runs pi2 on local computer. Returns output as string.
    """

    return subprocess.check_output([pi_path + "/pi2", pi_script])




def run_pi2(pi_script, output_prefix):
    """
    Runs pi2 job either locally or on cluster.
    output_prefix is a prefix for error and output log files. (currently not in use)
    """

    pi.submitjob(pi_script, "normal")



def wait_for_cluster_jobs():
    """
    Waits until all jobs in the submitted jobs array are finished.
    """

    if not is_use_cluster():
        return

    print("Waiting for cluster jobs to finish...")
    pi.waitforjobs()



class Scan:
    """
    Holds information about one sub-scan.
    """

    def __init__(self):
        # File that contains the image data
        rec_file = ''

        # Index of the scan (first scan = 1, second scan = 2, ...)
        index = -1

        # Size of the image
        dimensions = np.array([0, 0, 0])

        # Similarity transformation (rotation matrix, scaling factor, negative of position in stitched image)
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
        # This member is used as initial guess of the position to determine approximately overlapping regions
        # of the sub-images.
        position = np.array([0, 0, 0])

        # Normalization factors for gray values
        norm_fact = 0
        norm_fact_std = 1
        mean_def = 0

    def is_rec_ok(self):
        """
        Tests if reconstructed image file exists.
        """
        
        # TODO: Replace this hack with pi2py.
        pi_script = f"fileinfo({self.rec_file});"
        s = run_pi2_locally(pi_script)
        s = s.decode('ASCII')
        lines = s.splitlines()
        return len(lines) == 3 # Three lines means the file was found and can be read.







def get_image_size(filename):
    """
    Finds out size of given image and returns it as numpy array.
    """

    # TODO: Replace this hack with pi2py.
    pi_script = f"fileinfo({filename});"
    s = run_pi2_locally(pi_script)
    s = s.decode('ASCII')
    lines = s.splitlines()
    if len(lines) == 3:
        return from_string(lines[1])

    raise RuntimeError("Unable to read dimensions from image file " + filename)




def raw_exists(prefix):
    """
    Tests if a .raw file with given prefix exists.
    """

    return len(glob.glob(f"{prefix}*.raw")) > 0




def auto_binning(relations, binning):
    """
    Makes binned versions of the original input files if they do not exist yet.
    Returns true if all binned files are already done.
    """

    if binning < 1:
        raise RuntimeError("Value of binning must be greater than or equal to 1.")

    result = True
    if binning != 1:

        # Find common part in all the input file names so that we can remove that
        filenames = [node.rec_file for node in relations.nodes()]
        common_prefix = os.path.commonpath(filenames)

        # Create binned input files and adjust scan settings
        for node in relations.nodes():

            orig_file = node.rec_file

            # The orig_file may be in different folder!
            #binned_file = f"bin{binning}_{orig_file}"
            binned_file = orig_file[len(common_prefix):]
            if binned_file[0] == '\\' or binned_file[0] == '/':
                binned_file = binned_file[1:]
            binned_file = binned_file.replace('\\', '-').replace('/', '-')
            binned_file = f"bin{binning}_{binned_file}"

            node.binned_file = binned_file

            # TODO: If dimensions are used in non-computational context (file names etc.) they should be rounded/read from file again.
            #node.dimensions = node.dimensions / binning
            #node.position = node.position / binning

            if not raw_exists(binned_file):

                print(f"Binning {orig_file} to {binned_file}")

                params = (f"echo;"
                          f"read(img, {orig_file});"
                          f"maskedbin(img, binned, {binning}, 0, 0);"
                          f"writeraw(binned, {binned_file});"
                         )

                run_pi2(params, f"binning_{binned_file}")

                result = False

    return result


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



def overlap_region(scan1, scan2):
    """
    Calculates coordinates of edges of the overlap region between scan1 and scan2 in scan1 coordinates.
    Returns (xmin, ymin, zmin, xmax, ymax, zmax, shiftx, shifty, shiftz),
    where shift[xyz] are shift vector from scan2 to scan1.
    """

    xmin = int(round(max([0, scan2.position[0] - scan1.position[0]])))
    ymin = int(round(max([0, scan2.position[1] - scan1.position[1]])))
    zmin = int(round(max([0, scan2.position[2] - scan1.position[2]])))
    xmax = int(round(min([scan1.dimensions[0], scan2.position[0] - scan1.position[0] + scan2.dimensions[0]])))
    ymax = int(round(min([scan1.dimensions[1], scan2.position[1] - scan1.position[1] + scan2.dimensions[1]])))
    zmax = int(round(min([scan1.dimensions[2], scan2.position[2] - scan1.position[2] + scan2.dimensions[2]])))
    shiftx = int(round(-(scan2.position[0] - scan1.position[0])))
    shifty = int(round(-(scan2.position[1] - scan1.position[1])))
    shiftz = int(round(-(scan2.position[2] - scan1.position[2])))

    return xmin, ymin, zmin, xmax, ymax, zmax, shiftx, shifty, shiftz



def calculate_displacement_field(sample_name, scan1, scan2, point_spacing, coarse_block_radius, coarse_binning, fine_block_radius, fine_binning, normalize, filter_threshold):
    """
    Calculates displacement field between scan1 and scan2 locally or submits displacement field calculation job to cluster.
    Calculates also filtered displacement field.
    """

    # Calculate initial minimum and maximum for the reference grid.
    xmin, ymin, zmin, xmax, ymax, zmax, shiftx, shifty, shiftz = overlap_region(scan1, scan2)

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

    Returns (x0, y0, z0, x1, y1, z1, gof, norm_fact, norm_fact_std, mean_def), where
    x0, y0, z0 contain points in reference image, and
    x1, y1, z1 contain corresponding points in deformed image.
    norm_fact is value that must be added to deformed image to make mean of deformed image and reference image the same.
    norm_fact_std is multiplicative normalization factor. Values of scan2 must be multiplied by this value after subtraction of mean_def.
    mean_def is mean of scan2 in the overlapping region.
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
        norm_fact_std = float(lines[4])
        mean_def = float(lines[5])


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



    return x, y, z, x1, y1, z1, gof, norm_fact, norm_fact_std, mean_def



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

    TODO: I'm not completely sure if the algorithm actually works as in the comment above, or if it returns the inverse of the transformation.
    """

    # Remove bad values
    good_ind = ~(np.isnan(u_orig) | (gof_orig <= 0))

    if not good_ind.any():
        # No good values. This has been taken care by c++ code, so proceed with "all good".
        print("Warning: no usable points in displacement field. Consider using smaller grid step, larger filtering threshold, and larger overlap between images.")
        good_ind = gof_orig == gof_orig
        allow_rotation = False

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
    if sum(W) > 0:
        W = W / sum(W)
    else:
        W = np.ones(gof.shape)

    n, dim = P.shape

    # Rotation center and translation
    c = np.reshape(np.average(P, axis=0, weights=W), (1, 3))
    delta = np.reshape(np.average(Q, axis=0, weights=W), (1, 3)) - c

    if allow_rotation:
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
    else:
        R = np.eye(3, 3)

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
    #if not allow_rotation:
    #    R = np.eye(3, 3)

    return delta, c, a, R





def calculate_displacement_fields(sample_name, relations, point_spacing, coarse_block_radius, coarse_binning, fine_block_radius, fine_binning, normalize, filter_threshold, force_redo):
    """
    Calculate displacement fields and load results to relations network.
    Returns true if all displacement fields have been calculated and read, and false otherwise.
    """

    # Determine which pairwise displacement fields should be calculated (i.e. don't calculate for non-existing images or if the calculation has been made already)
    # and start calculations or submit jobs
    jobs_started = 0
    for edge in relations.edges(data=True):
        scan1 = edge[0]
        scan2 = edge[1]

        if force_redo or (not is_displacement_ok(sample_name, scan1, scan2)):
            print("Calculate displacement field %i -> %i" % (scan1.index, scan2.index))
            calculate_displacement_field(sample_name, scan1, scan2, point_spacing, coarse_block_radius, coarse_binning, fine_block_radius, fine_binning, normalize, filter_threshold)
            jobs_started = jobs_started + 1


    # If using a cluster, tell the user to wait until the results of the calculations are available.
    if (jobs_started > 0) and is_use_cluster():
        return False


    for edge in relations.edges(data=True):
        scan1 = edge[0]
        scan2 = edge[1]

        if not is_displacement_ok(sample_name, scan1, scan2):
            raise RuntimeError("Displacement field calculation has failed, unable to continue. See error messages above.")


    # Check that filtered displacement fields are also ok.
    jobs_started = 0
    for edge in relations.edges(data=True):
        scan1 = edge[0]
        scan2 = edge[1]

        if not is_filtered_displacement_ok(sample_name, scan1, scan2):
            print("Filter displacement field %i -> %i" % (scan1.index, scan2.index))
            filter_displacement_field(sample_name, scan1, scan2, filter_threshold)
            jobs_started = jobs_started + 1

    if (jobs_started > 0) and is_use_cluster():
        return False

    return True





def read_displacement_fields(sample_name, relations, allow_rotation):
    """
    Reads all displacement fields and determines similarity transformation for each field.
    """

    print("Reading displacement fields...")

    # Read displacement fields and store all data in relations graph.
    counter = 0
    for edge in relations.edges(data=True):
        scan1 = edge[0]
        scan2 = edge[1]

        print(f"{scan1.index} -> {scan2.index} ({scan1.rec_file} -> {scan2.rec_file})")

        x, y, z, x1, y1, z1, gof, norm_fact, norm_fact_std, mean_def = read_displacement_field(sample_name, scan1, scan2)
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
        deltan = -delta

        delta = deltan
        c = cn

        edge[2]["weight"] = 1.0 - np.mean(gof)
        edge[2]["delta"] = delta
        edge[2]["c"] = c
        edge[2]["a"] = a
        edge[2]["R"] = R
        edge[2]["norm_fact"] = norm_fact
        edge[2]["norm_fact_std"] = norm_fact_std
        edge[2]["mean_def"] = mean_def

        # Store point correspondences for global optimization
        # r0 = (x0, y0, z0) are points in scan1 coordinates,
        # r1 = (x1, y1, z1) are points in scan2 coordinates,
        edge[2]["r0"] = np.array([x.flatten(), y.flatten(), z.flatten()])
        edge[2]["r1"] = np.array([x1.flatten(), y1.flatten(), z1.flatten()])

        counter = counter + 1
        print(f"{counter} / {relations.number_of_edges()}\r", end="")




def convert_to_world(parent_scan, child_scan, tree):
    """
    Determines similarity transformation of child node in world coordinates, based on world-coordinate similarity transformation of parent node and parent to child transformation.
    """

    delta = tree[parent_scan][child_scan]["delta"]
    c = tree[parent_scan][child_scan]["c"]
    a = tree[parent_scan][child_scan]["a"]
    R = tree[parent_scan][child_scan]["R"]
    norm_fact = tree[parent_scan][child_scan]["norm_fact"]
    norm_fact_std = tree[parent_scan][child_scan]["norm_fact_std"]
    mean_def = tree[parent_scan][child_scan]["mean_def"]

    deltadot = c + delta

    # Sanity check
    if not parent_scan.transformation_set:
        raise RuntimeError(f"Parent scan at ({parent_scan.position[0]}, {parent_scan.position[1]}, {parent_scan.position[2]}) has no transformation set yet, and its transformation is required for child scan at ({child_scan.position[0]}, {child_scan.position[1]}, {child_scan.position[2]}).")

    # Convert transformation from T(n-1, n) to T(0, n)
    wa = parent_scan.a * a
    wR = np.matmul(parent_scan.R, R)
    wc = 1 / a * np.matmul(np.linalg.inv(R), parent_scan.c - deltadot) + c
    wnorm_fact = parent_scan.norm_fact + norm_fact
    wnorm_fact_std = parent_scan.norm_fact_std * norm_fact_std
    wmean_def = mean_def

    return wa, wR, wc, wnorm_fact, wnorm_fact_std, wmean_def




def determine_average_transformations(relations):
    """
    Sets node's transformation to average of transformations between scans connected to it,
    and converted to world coordinates.
    Does not change transformation of the root node (the node that has no predecessors).
    """


    ordered = list(nx.topological_sort(relations))

    for scan in ordered:

        # Process all neighbours that have been scanned earlier
        count = 0
        atot = 0
        ctot = 0
        qtot = Quaternion(0, 0, 0, 0)
        normtot = 0
        normtot_std = 0
        mean_tot = 0

        parent_scans = relations.predecessors(scan)
        for parent in parent_scans:

            a, R, c, nf, nf_std, mean_def = convert_to_world(parent, scan, relations)
            atot = atot + a
            ctot = ctot + c
            # TODO: This quaternion averaging is just approximation that might work if all the rotations are almost the same...
            qtot = qtot + Quaternion(matrix=R)
            normtot = normtot + nf
            normtot_std = normtot_std + nf_std
            mean_tot = mean_tot + mean_def
            count = count + 1


        if count > 0:
            scan.a = atot / count
            scan.c = ctot / count
            qtot = qtot / count
            qtot = qtot.normalised
            scan.R = qtot.rotation_matrix
            scan.norm_fact = normtot / count
            scan.norm_fact_std = normtot_std / count
            scan.mean_def = mean_tot / count

        scan.transformation_set = True




def check_result(result):
    """
    Checks if minimization results indicates success and prints warning message if not.
    """

    if not result.success:
        print("Warning: Unable to find globally optimal poses for sub-images. Minimization process returned error message")
        print(result.message)
        print("It is possible that the result is still good. If not, consider switching to non-globally optimal solution.")


#def add_corner_points_3d(relations):
#    """
#    For each edge in relations, adds list of corner points of overlapping region between the scans in older scan coordinates to edge['corner_points'].
#    """

#    # Determine corner points of overlapping region between each pair of scans
#    for edge in relations.edges():
#        parent = edge[0]
#        scan = edge[1]

#        # Find corner points of initial overlap region between the two scans in parent scan coordinates
#        xmin, ymin, zmin, xmax, ymax, zmax, shiftx, shifty, shiftz = overlap_region(parent, scan)

#        points = np.array([np.array([[xmin, ymin, zmin]]).T,
#                np.array([[xmax, ymin, zmin]]).T,
#                np.array([[xmin, ymax, zmin]]).T,
#                np.array([[xmax, ymax, zmin]]).T,
#                np.array([[xmin, ymin, zmax]]).T,
#                np.array([[xmax, ymin, zmax]]).T,
#                np.array([[xmin, ymax, zmax]]).T,
#                np.array([[xmax, ymax, zmax]]).T])

#        relations[parent][scan]["corner_points"] = points



#def add_corner_points_2d(relations):
#    """
#    For each edge in relations, adds list of corner points of overlapping region between the scans in older scan coordinates to edge['corner_points'].
#    Discards z dimension
#    """

#    # Determine corner points of overlapping region between each pair of scans
#    for edge in relations.edges():
#        parent = edge[0]
#        scan = edge[1]

#        # Find corner points of initial overlap region between the two scans in parent scan coordinates
#        xmin, ymin, zmin, xmax, ymax, zmax, shiftx, shifty, shiftz = overlap_region(parent, scan)

#        points = np.array([np.array([[xmin, ymin, 0]]).T,
#                           np.array([[xmax, ymin, 0]]).T,
#                           np.array([[xmin, ymax, 0]]).T,
#                           np.array([[xmax, ymax, 0]]).T])

#        relations[parent][scan]["corner_points"] = points


#def weight(cs, rotmats, edgelist):
#    """
#    Calculates minimization cost function for global optimization of sub-image locations and orientations.
#    This version uses only corner points of the overlapping region in the weight calculation.
#    """

#    value = 0

#    # Process all overlapping scans
#    for edge in edgelist:
#        parent = edge[0]
#        scan = edge[1]

#        # Get pairwise transformation from parent to child (determined from point correspondences)
#        delta = edge[2]
#        c = edge[3]
#        a = edge[4]
#        R = edge[5]

#        # Get corner points of overlapping region in parent scan coordinates
#        points = edge[6]



#        # Get current optimized world to local transformations of parent and child
#        parent_c = cs[parent]
#        parent_R = rotmats[parent]
#        parent_a = 1

#        scan_c = cs[scan]
#        scan_R = rotmats[scan]
#        scan_a = 1



#        # Transform the points to child scan coordinates using data from pose array
#        # Notice:
#        # world to scan conversion
#        # scanP = 1 / a * (Rinv * worldP) + c
#        # scan to world conversion
#        # worldP = a * R * (scanP - c)
#        # where
#        # c = -position

#        # Transform first to world coordinates using parent scan transformation...
#        transf1 = parent_a * np.matmul(parent_R, points - parent_c)

#        # ...and then to child scan coordinates using child scan transformation
#        transf1 = 1 / scan_a * np.matmul(scan_R.T, transf1) + scan_c


#        # Transform points to child scan coordinates using data from relations[parent][child] (determined from point correspondences)
#        transf2 = 1 / a * np.matmul(R.T, points - (c + delta)) + c


#        # Calculate sum of magnitude of difference between the points
#        #d = np.linalg.norm(transf1 - transf2, axis=1)
#        #value = value + np.sum(d**2)

#        D = transf1 - transf2
#        d2 = D[:, 0]**2 + D[:, 1]**2 + D[:, 2]**2
#        value = value + np.sum(d2)
#        #value = value + sum(d2)


#    return value




def weight(cs, rotmats, edgelist):
    """
    Calculates minimization cost function for global optimization of sub-image locations and orientations.
    This version uses all the point correspondences determined with phase correlation.
    """

    value = 0

    # Process all overlapping scans
    for edge in edgelist:
        parent = edge[0]
        scan = edge[1]

        # Get corresponding points in the two scans
        parent_points = edge[7]
        scan_points = edge[8]

        # Get current optimized world to local transformations of parent and child
        parent_c = cs[parent]
        parent_R = rotmats[parent]
        parent_a = 1

        scan_c = cs[scan]
        scan_R = rotmats[scan]
        scan_a = 1


        # Transform the points to world coordintes
        # Notice:
        # world to scan conversion
        # scanP = 1 / a * (Rinv * worldP) + c
        # scan to world conversion
        # worldP = a * R * (scanP - c)
        # where
        # c = -position

        parent_points = parent_a * np.matmul(parent_R, parent_points - parent_c)

        scan_points = scan_a * np.matmul(scan_R, scan_points - scan_c)


        # Calculate sum of magnitude of difference between the points
        D = parent_points - scan_points
        d2 = D[0, :]**2 + D[1, :]**2 + D[2, :]**2
        value = value + np.sum(d2)


    return value



def angle_to_rot_mat(angle):

    return np.array([[math.cos(angle), -math.sin(angle), 0],
                     [math.sin(angle), math.cos(angle), 0],
                     [0, 0, 1]])


def build_edge_data_list(edges):
    """
    Converts relations.edges() type object to simple edge list.
    """

    # Version to be used with corner point-based optimization
    #return [[edge[0], edge[1], edge[2]["delta"], edge[2]["c"], edge[2]["a"], edge[2]["R"], edge[2]["corner_points"], edge[2]["r0"], edge[2]["r1"]] for edge in edges]
    return [[edge[0], edge[1], edge[2]["delta"], edge[2]["c"], edge[2]["a"], edge[2]["R"], 0, edge[2]["r0"], edge[2]["r1"]] for edge in edges]


def get_edge_list(relations):
    """
    Creates edge list for weight functions.
    """

    return build_edge_data_list(relations.edges(data=True))



def weight_and_derivative_generic(pose, variables_per_scan, pose_to_cs, pose_to_rotmat, relations, edges, scan_to_index):
    """
    Calculates value of cost function and its (numerical) derivative.
    """


    # Extract -position and rotation matrix for each sub-image from the pose array
    cs = {scan: pose_to_cs(pose, scan, scan_to_index[scan]) if scan_to_index[scan] >= 0 else scan.c for scan in relations}
    rotmats = {scan: pose_to_rotmat(pose, scan, scan_to_index[scan]) if scan_to_index[scan] >= 0 else scan.R for scan in relations}

    # Calculate value of the weight function
    total_weight = weight(cs, rotmats, edges)


    # Calculate gradient vector, \Delta_i = dF / dx_i
    eps = 0.000001
    der = np.zeros(pose.shape)
    for scan in relations:
        # Only process scans that are not fixed
        if scan_to_index[scan] >= 0:
            relevant_edges = list(relations.in_edges(scan, data=True))
            relevant_edges.extend(relations.out_edges(scan, data=True))

            relevant_edge_list = build_edge_data_list(relevant_edges)

            # For each variable n in the scan
            for n in range(0, variables_per_scan):
                i = variables_per_scan * scan_to_index[scan] + n

                # Change pose
                old_pose_elem = pose[i]
                pose[i] = pose[i] + eps

                # Update cs and rotmat of the scan
                cs[scan] = pose_to_cs(pose, scan, scan_to_index[scan])
                rotmats[scan] = pose_to_rotmat(pose, scan, scan_to_index[scan])

                # Calculate weight with current pose, with only those edges that connect to or from the current scan
                w1 = weight(cs, rotmats, relevant_edge_list)


                # Change pose
                pose[i] = old_pose_elem - eps

                # Update cs and rotmat of the scan
                cs[scan] = pose_to_cs(pose, scan, scan_to_index[scan])
                rotmats[scan] = pose_to_rotmat(pose, scan, scan_to_index[scan])

                # Calculate weight with current pose, with only those edges that connect to or from the current scan
                w0 = weight(cs, rotmats, relevant_edge_list)



                # Restore original pose
                pose[i] = old_pose_elem


                # Calculate derivative
                der[i] = (w1 - w0) / (2 * eps)




    print(f"Weight = {total_weight}", end='\r')
    return total_weight, der





def weight_3d_trans_rot(pose, relations, edges, scan_to_index):

    def pose_to_cs(pose, scan, index):
        return np.array([[pose[7 * index + 0] * scan.dimensions[0], pose[7 * index + 1] * scan.dimensions[1], pose[7 * index + 2] * scan.dimensions[2]]]).T

    def pose_to_rotmat(pose, scan, index):
        return Quaternion(pose[7 * index + 3],
                          pose[7 * index + 4],
                          pose[7 * index + 5],
                          pose[7 * index + 6]).normalised.rotation_matrix

    return weight_and_derivative_generic(pose, 7,
                                         pose_to_cs,
                                         pose_to_rotmat,
                                         relations, edges, scan_to_index)



def weight_3d_trans(pose, relations, edges, scan_to_index):

    def pose_to_cs(pose, scan, index):
        return np.array([[pose[3 * index + 0] * scan.dimensions[0], pose[3 * index + 1] * scan.dimensions[1], pose[3 * index + 2] * scan.dimensions[2]]]).T

    def pose_to_rotmat(pose, scan, index):
        return scan.R

    return weight_and_derivative_generic(pose, 3,
                                         pose_to_cs,
                                         pose_to_rotmat,
                                         relations, edges, scan_to_index)




def weight_2d_trans_rot(pose, relations, edges, scan_to_index):

    def pose_to_cs(pose, scan, index):
        return np.array([[pose[3 * index + 0] * scan.dimensions[0], pose[3 * index + 1] * scan.dimensions[1], 0]]).T

    def pose_to_rotmat(pose, scan, index):
        return angle_to_rot_mat(pose[3 * index + 2])

    return weight_and_derivative_generic(pose, 3,
                                         pose_to_cs,
                                         pose_to_rotmat,
                                         relations, edges, scan_to_index)




def weight_2d_trans(pose, relations, edges, scan_to_index):

    def pose_to_cs(pose, scan, index):
        return np.array([[pose[2 * index + 0] * scan.dimensions[0], pose[2 * index + 1] * scan.dimensions[1], 0]]).T

    def pose_to_rotmat(pose, scan, index):
        return scan.R

    return weight_and_derivative_generic(pose, 2,
                                         pose_to_cs,
                                         pose_to_rotmat,
                                         relations, edges, scan_to_index)





def optimize_transformations_3d_trans_rot(relations):
    """
    Optimizes transformations of each sub-image so that they are globally optimal in some sense.
    Uses current transformations as initial guess.
    Does not change transformation of the root node that does not have any predecessors, i.e. the node
    that has no predecessors is fixed.
    """

    root = find_first_node(relations)

    #add_corner_points_3d(relations)

    # Initial guess for pose = (3 components of negation of translation vector + 4 components of rotation quaternion) * number of sub-images
    # Pose of the root image is not estimated.
    initial_pose = np.zeros(((3 + 4) * (relations.number_of_nodes() - 1)))
    scan_to_index = {}
    scan_to_index[root] = -1
    scan_index = 0
    for scan in relations:
        if scan != root:
            scan_to_index[scan] = scan_index

            q = Quaternion(matrix=scan.R)

            initial_pose[7 * scan_index + 0] = scan.c[0] / scan.dimensions[0]
            initial_pose[7 * scan_index + 1] = scan.c[1] / scan.dimensions[1]
            initial_pose[7 * scan_index + 2] = scan.c[2] / scan.dimensions[2]
            initial_pose[7 * scan_index + 3] = q[0]
            initial_pose[7 * scan_index + 4] = q[1]
            initial_pose[7 * scan_index + 5] = q[2]
            initial_pose[7 * scan_index + 6] = q[3]

            scan_index = scan_index + 1

    result = scipy.optimize.minimize(weight_3d_trans_rot, initial_pose, (relations, get_edge_list(relations), scan_to_index), tol=500, jac=True)
    check_result(result)
    optimal_pose = result.x

    # Create transformations from the optimization result
    for scan, index in scan_to_index.items():

        if scan != root:
            c = np.array([optimal_pose[7 * index + 0] * scan.dimensions[0],
                          optimal_pose[7 * index + 1] * scan.dimensions[1],
                          optimal_pose[7 * index + 2] * scan.dimensions[2]])
            R = Quaternion(optimal_pose[7 * index + 3],
                           optimal_pose[7 * index + 4],
                           optimal_pose[7 * index + 5],
                           optimal_pose[7 * index + 6]).normalised.rotation_matrix

            scan.c = c
            scan.R = R




def optimize_transformations_3d_trans(relations):
    """
    Optimizes translation of each sub-image so that they are globally optimal in some sense.
    Uses current transformations as initial guess.
    Does not change transformation of the root node that does not have any predecessors, i.e. the node
    that has no predecessors is fixed.
    """

    root = find_first_node(relations)

    #add_corner_points_3d(relations)


    # Initial guess for pose = (3 components of negation of translation vector) * number of sub-images
    # Pose of the root image is not estimated.
    initial_pose = np.zeros((3 * (relations.number_of_nodes() - 1)))
    scan_to_index = {}
    scan_to_index[root] = -1
    scan_index = 0
    for scan in relations:
        if scan != root:
            scan_to_index[scan] = scan_index

            initial_pose[3 * scan_index + 0] = scan.c[0] / scan.dimensions[0]
            initial_pose[3 * scan_index + 1] = scan.c[1] / scan.dimensions[1]
            initial_pose[3 * scan_index + 2] = scan.c[2] / scan.dimensions[2]

            scan_index = scan_index + 1

    result = scipy.optimize.minimize(weight_3d_trans, initial_pose, (relations, get_edge_list(relations), scan_to_index), tol=500, jac=True)
    check_result(result)
    optimal_pose = result.x


    # Create transformations from the optimization result
    for scan, index in scan_to_index.items():

        if scan != root:
            c = np.array([optimal_pose[3 * index + 0] * scan.dimensions[0],
                          optimal_pose[3 * index + 1] * scan.dimensions[1],
                          optimal_pose[3 * index + 2] * scan.dimensions[2]])

            scan.c = c





def optimize_transformations_2d_trans_rot(relations):
    """
    Optimizes transformations of each sub-image so that they are globally optimal in some sense.
    Uses current transformations as initial guess.
    Does not change transformation of the root node that does not have any predecessors, i.e. the node
    that has no predecessors is fixed.
    """

    root = find_first_node(relations)

    #add_corner_points_2d(relations)

    # Initial guess for pose = (2 components of negation of translation vector + 1 rotation angle) * number of sub-images
    # Pose of the root image is not estimated.
    initial_pose = np.zeros((3 * (relations.number_of_nodes() - 1)))
    scan_to_index = {}
    scan_to_index[root] = -1
    scan_index = 0
    for scan in relations:
        if scan != root:
            scan_to_index[scan] = scan_index

            angle = math.atan2(scan.R[1, 0], scan.R[0, 0])

            initial_pose[3 * scan_index + 0] = scan.c[0] / scan.dimensions[0]
            initial_pose[3 * scan_index + 1] = scan.c[1] / scan.dimensions[1]
            initial_pose[3 * scan_index + 2] = angle

            scan_index = scan_index + 1


    result = scipy.optimize.minimize(weight_2d_trans_rot, initial_pose, (relations, get_edge_list(relations), scan_to_index), tol=500, jac=True)

    check_result(result)
    optimal_pose = result.x

    # Create transformations from the optimization result
    for scan, index in scan_to_index.items():

        if scan != root:
            c = np.array([optimal_pose[3 * index + 0] * scan.dimensions[0],
                          optimal_pose[3 * index + 1] * scan.dimensions[1],
                          0])
            angle = optimal_pose[3 * index + 2]
            R = angle_to_rot_mat(angle)

            scan.c = c
            scan.R = R




def optimize_transformations_2d_trans(relations):
    """
    Optimizes translation of each sub-image so that they are globally optimal in some sense.
    Uses current transformations as initial guess.
    Does not change transformation of the root node that does not have any predecessors, i.e. the node
    that has no predecessors is fixed.
    """

    root = find_first_node(relations)

    #add_corner_points_2d(relations)


    # Initial guess for pose = (2 components of negation of translation vector) * number of sub-images
    # Pose of the root image is not estimated.
    initial_pose = np.zeros((3 * (relations.number_of_nodes() - 1)))
    scan_to_index = {}
    scan_to_index[root] = -1
    scan_index = 0
    for scan in relations:
        if scan != root:
            scan_to_index[scan] = scan_index

            initial_pose[2 * scan_index + 0] = scan.c[0] / scan.dimensions[0]
            initial_pose[2 * scan_index + 1] = scan.c[1] / scan.dimensions[1]

            scan_index = scan_index + 1

    result = scipy.optimize.minimize(weight_2d_trans, initial_pose, (relations, get_edge_list(relations), scan_to_index), tol=500, jac=True)
    check_result(result)
    optimal_pose = result.x


    # Create transformations from the optimization result
    for scan, index in scan_to_index.items():

        if scan != root:
            c = np.array([optimal_pose[2 * index + 0] * scan.dimensions[0],
                          optimal_pose[2 * index + 1] * scan.dimensions[1],
                          0])

            scan.c = c





def optimize_transformations(relations, allow_rotation):

    if relations.number_of_nodes() <= 1:
        return

    node = find_first_node(relations)
    if node.dimensions[2] <= 1:
        if allow_rotation:
            optimize_transformations_2d_trans_rot(relations)
        else:
            optimize_transformations_2d_trans(relations)
    else:
        if allow_rotation:
            optimize_transformations_3d_trans_rot(relations)
        else:
            optimize_transformations_3d_trans(relations)


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
    for scan in comp:
        a = scan.a
        R = scan.R
        c = scan.c
        width = scan.dimensions[0]
        height = scan.dimensions[1]
        depth = scan.dimensions[2]

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



def fix_directories(path):
    """
    Removes directory separator characters from the path and returns the result.
    """

    while len(path) > 0 and (path[0] == '/' or path[0] == '\\'):
        path = path[1:]

    path = path.replace('.', '')
    path = path.replace('/', '-')
    path = path.replace('\\', '-')

    return path


def save_transformation(sample_name, scan, relations):
    """
    Saves world to local similarity tranformation of node at given position.
    """

    a = scan.a
    R = scan.R
    c = scan.c
    norm_fact = scan.norm_fact
    norm_fact_std = scan.norm_fact_std
    mean_def = scan.mean_def

    scan_name = fix_directories(scan.rec_file)
    file = f"{scan_name}_transformation.txt"
    with open(file, 'wb') as f:
        # Write world to local similarity transformation of this node
        np.savetxt(f, np.array([a]))
        np.savetxt(f, c)
        np.savetxt(f, R)
        np.savetxt(f, np.array([norm_fact]))
        np.savetxt(f, np.array([norm_fact_std]))
        np.savetxt(f, np.array([mean_def]))

        # Save information about parent node indices
        parents = list(relations.predecessors(scan))

        # Write count of parent nodes
        np.savetxt(f, np.array([len(parents)]))
        for parent_scan in parents:
            parent_scan_name = fix_directories(parent_scan.rec_file)

            prefix = displacement_file_prefix(sample_name, parent_scan, scan)

            np.savetxt(f, np.array([f"{parent_scan_name}_world_to_local"]), fmt="%s")
            np.savetxt(f, np.array([f"{prefix}_filtered"]), fmt="%s")


    scan.transformation_file = file




def find_roots(tree):
    """
    Creates list of nodes with no predecessors.
    """

    roots = []

    for node in tree:
        if len(list(tree.predecessors(node))) <= 0:
            roots.append(node)

    return roots


def is_world_to_local_ok(prefix):
    """
    Tests if world to local transformation has been calculated.
    """
    
    filename = f"{prefix}_refpoints.txt"
    return os.path.isfile(filename)
        

def calculate_world_to_local(tree, allow_local_deformations, force_redo):
    """
    Calculates world to local transformations for all nodes in the given tree, starting from the root nodes that have no incoming connections.
    """

    # We don't want to use topological_sort here as we can process multiple images at once as they
    # depend only on images that have been processed already.

    done = {scan: False for scan in tree}

    while not all(done[scan] for scan in tree):
        done_prev = done.copy()
        for scan in tree:
            if not done[scan]:
                pred = tree.predecessors(scan)
                all_done = all([done_prev[scan] for scan in pred])

                if all_done:
                    #print(scan.rec_file)
                    done[scan] = True

                    scan_name = fix_directories(scan.rec_file)

                    # Calculate world to local grid transform using pi
                    scan.world_to_local_prefix = f"{scan_name}_world_to_local"
                    script = (f"echo;"
                              f"determine_world_to_local({scan.transformation_file}, [{scan.dimensions[0]}, {scan.dimensions[1]}, {scan.dimensions[2]}], {scan.world_to_local_prefix}, {allow_local_deformations});"
                             )
                             
                    if force_redo or (not is_world_to_local_ok(scan.world_to_local_prefix)):     
                        run_pi2(script, scan.world_to_local_prefix)

        wait_for_cluster_jobs()
        
    


def run_stitching(comp, sample_name, normalize, max_circle, global_optimization, allow_rotation, allow_local_deformations, create_goodness_file, force_redo):
    """
    Prepares and runs pi2 stitching process for connected component 'comp' of scan relations tree 'tree'.
    - determines final world to image transformations
    - saves them
    - determines bounds of the stitched image
    - makes index file for pi2, runs pi2 stitching in blocks
    """

    if not nx.is_directed_acyclic_graph(comp):
        raise RuntimeError("Connected component is not directed acyclic graph.")

    # Position root node
    root = find_first_node(comp)
    root.R = np.eye(3, 3)
    root.a = 1
    root.c = -root.position.reshape(-1, 1)

    # Determine transformations based on transformations of older connected nodes
    determine_average_transformations(comp)

    # This line can be used to rotate the first image to test the global optimization process
    #root.R = Quaternion(axis=[0, 0, 1], angle=20 / 180.0 * 3.14).rotation_matrix

    # Alternatively one can use estimated node position as initial guess for node location
    #for node in comp.nodes:
    #    node.c = -node.position.reshape(-1, 1)

    if global_optimization:
        print("Finding globally optimal locations and orientations for the sub-images...")
        optimize_transformations(comp, allow_rotation)
        # TODO: Here we could do similar global optimization process for intensities, too.

    # Save the transformations
    for node in comp.nodes:
        save_transformation(sample_name, node, comp)

    # Calculate world to local grid transformations
    print("Calculating world to local transformation for each image...")
    calculate_world_to_local(comp, allow_local_deformations, force_redo)


    print("Stitching...")

    # Determine approximate bounds for the stitched image
    minx, miny, minz, maxx, maxy, maxz = determine_bounds(comp)

    out_width = int(round(maxx - minx))
    out_height = int(round(maxy - miny))
    out_depth = int(round(maxz - minz))
    minx = int(round(minx))
    miny = int(round(miny))
    minz = int(round(minz))

    out_width = max(1, out_width)
    out_height = max(1, out_height)
    out_depth = max(1, out_depth)

    start_node = find_first_node(comp)
    out_template = f"{sample_name}_{start_node.position[0]}_{start_node.position[1]}_{start_node.position[2]}"
    out_file = f"{out_template}_{out_width}x{out_height}x{out_depth}.raw"
    out_goodness_file = f"{out_template}_goodness_{out_width}x{out_height}x{out_depth}.raw"

    # Make index file
    index_file = out_template + "_index.txt"
    with open(index_file, 'w') as f:
        for scan in comp:
            f.write(f"{scan.rec_file}\n")
            f.write(f"{scan.world_to_local_prefix}\n")

    # Find name of one of the input files so that we can copy data type info from that.
    for scan in comp:
        first_file_name = scan.rec_file
        break

    # Cut big sample to small regions and process them separately to save memory
    block_size = int(max_block_size)

    jobs_started = 0
    for zstart in range(minz, out_depth + minz, block_size):
        for ystart in range(miny, out_height + miny, block_size):
            for xstart in range(minx, out_width + minx, block_size):

                curr_width = min(block_size, out_width - (xstart - minx))
                curr_height = min(block_size, out_height - (ystart - miny))
                curr_depth = min(block_size, out_depth - (zstart - minz))

                if not create_goodness_file:
                    pi_script = (f"echo;"
                                 f"newlikefile(outimg, {first_file_name}, Unknown, 1, 1, 1);"
                                 f"stitch_ver2(outimg, {index_file}, {xstart}, {ystart}, {zstart}, {curr_width}, {curr_height}, {curr_depth}, {normalize}, {max_circle});"
                                 f"writerawblock(outimg, {out_file}, [{xstart - minx}, {ystart - miny}, {zstart - minz}], [{out_width}, {out_height}, {out_depth}]);"
                                 f"newimage(marker, uint8, 1, 1, 1);"
                                 f"writetif(marker, {out_template}_{jobs_started}_done);"
                                )
                else:
                    pi_script = (f"echo;"
                             f"newlikefile(outimg, {first_file_name}, Unknown, 1, 1, 1);"
                             f"newlikefile(goodnessimg, {first_file_name}, Unknown, 1, 1, 1);"
                             f"stitch_ver3(outimg, goodnessimg, {index_file}, {xstart}, {ystart}, {zstart}, {curr_width}, {curr_height}, {curr_depth}, {normalize}, {max_circle});"
                             f"writerawblock(outimg, {out_file}, [{xstart - minx}, {ystart - miny}, {zstart - minz}], [{out_width}, {out_height}, {out_depth}]);"
                             f"writerawblock(goodnessimg, {out_goodness_file}, [{xstart - minx}, {ystart - miny}, {zstart - minz}], [{out_width}, {out_height}, {out_depth}]);"
                             f"newimage(marker, uint8, 1, 1, 1);"
                             f"writetif(marker, {out_template}_{jobs_started}_done);"
                            )

                if force_redo or (not os.path.isfile(f"{out_template}_{jobs_started}_done.tif")):
                     run_pi2(pi_script, "")
                else:
                     print(f"Stitch job {jobs_started} has been done already. Skipping it.")
                    
                # jobs_started must be increased even if no job started due to result of above test.
                jobs_started = jobs_started + 1

    return jobs_started


def run_stitching_for_all_connected_components(relations, sample_name, normalize, max_circle, global_optimization, allow_rotation, allow_local_deformations, create_goodness_file, force_redo):
    """
    Calls run_stitching for each connected component in relations network.
    """

    jobs_started = 0
    comps = (relations.subgraph(c) for c in nx.weakly_connected_components(relations))
    for comp in comps:
        jobs_started = jobs_started + run_stitching(comp, sample_name, normalize, max_circle, global_optimization, allow_rotation, allow_local_deformations, create_goodness_file, force_redo)

    if (jobs_started > 0) and is_use_cluster():
        return False

    return True
