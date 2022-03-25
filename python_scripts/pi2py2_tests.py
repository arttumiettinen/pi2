"""
This file contains some of the tests of the pi2py2 library.
"""



import sys
from pi2py2 import *
import numpy as np
pi2 = Pi2()


# Define convenience functions that return input and output file names.
# This is just to avoid copying the paths to all the examples in case they change.
def input_file(filename='t1-head_256x256x129.raw'):
    return '../../testing/input_data/' + filename

def input_file_bin():
    return '../../testing/input_data/t1-head_bin_256x256x129.raw'

def output_file(name):
    return '../../testing/pi2py2/' + name

total_tests = 0
failed_tests = 0
fail_fast = True

def check_result(is_ok, msg, print_if_ok = True):
    """
    If is_ok is false, shows msg on screen and increments failed test count.
    Increments total test count on every call.
    """

    global total_tests
    total_tests = total_tests + 1

    global failed_tests

    
    if not is_ok:
        print('-------')
        print('Result:')
        print(msg)
        print('-------')

        failed_tests = failed_tests + 1

        if fail_fast:
            print('Stopping further testing as requested by the fail_fast flag.')
            sys.exit(1)
    else:
        if print_if_ok:
            print('-------')
            print('Result:')
            print('OK')
            print('-------')
    



def calc_difference(img1, img2):
    """
    Determines maximum absolute difference between img1 and img2.
    """

    imgf = pi2.newimage()
    img2f = pi2.newimage()

    pi2.convert(img1, imgf, ImageDataType.FLOAT32)
    pi2.convert(img2, img2f, ImageDataType.FLOAT32)

    pi2.subtract(imgf, img2f)
    pi2.abs(imgf)
    M = pi2.newimage(ImageDataType.FLOAT32)
    pi2.maxval(imgf, M)
    M = M.get_value()

    return M



def check_distribution_test_result(file1, file2, operation_name, compname, tolerance=0.00001, data_type=ImageDataType.UNKNOWN):
    """
    Helper function used to check results of distributed processing tests.
    """

    # Load both results and calculate maximum absolute difference
    img1 = pi2.read(file1, data_type)
    img2 = pi2.read(file2, data_type)

    M = calc_difference(img1, img2)

    # Check that the difference is zero
    check_result(M <= tolerance, f"ERROR: Found difference in results of {compname} processing while testing {operation_name}. Max absolute difference = {M}.")
    


def test_difference_delaying(testname, script, resultname='result', tolerance=0.00001, maxmem=30, chunk_size=[30, 32, 33]):
    """
    Calculates given script using distributed processing with delaying enabled and disabled.
    Compares the results and reports errors.
    """

    outfile_normal = output_file(f"{testname}_distributed_normal")
    outfile_delayed = output_file(f"{testname}_distributed_delayed")

    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(maxmem)
    pi2.chunksize(chunk_size)

    pi2.delaying(False)
    pi2.run_script(script)
    pi2.writeraw(resultname, outfile_normal)

    pi2.delaying(True)
    pi2.run_script(script)
    pi2.writeraw(resultname, outfile_delayed)

    pi2.distribute(Distributor.NONE)

    check_distribution_test_result(outfile_normal, outfile_delayed, testname, 'normal distributed and delayed distributed', tolerance)




def test_difference_normal_distributed(opname, args, resultname='result', infile=input_file(), tolerance=0.00001, convert_to_type=ImageDataType.UNKNOWN, maxmem=15, chunk_size=[100, 101, 102]):
    """
    Calculates operation normally and using distributed processing.
    Calculates difference between the results of the two versions and prints a message if the results do not match.
    opname is name of the function to call, e.g. 'gaussfilter'.
    args is array of parameters to the function. The source image is loaded to variables 'img' and 'result',
    and name of variable containing the result is stored in 'resultname' argument.
    If 'resultname' is empty, variable 'result' is assumed to store the result.
    """

    print(f"Testing {opname}...")
    print('-' * (len(opname) + 8 + 3))

    argstr = '_'.join(str(e) for e in args)
    argstr = argstr.replace('*', '_') # Remove non-filename characters
    outfile_normal = output_file(f"head_{opname}_{argstr}_normal")
    outfile_distributed = output_file(f"head_{opname}_{argstr}_distributed")

    


    # Run the operation using distributed processing
    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(maxmem)
    pi2.chunksize(chunk_size)
    img = pi2.read(infile)
    if convert_to_type != ImageDataType.UNKNOWN:
        pi2.convert(img, convert_to_type)

    result = pi2.newlike(img)
    pi2.set(result, img)

    temp1 = pi2.newimage(ImageDataType.FLOAT32)
    temp2 = pi2.newimage(ImageDataType.FLOAT32)

    args_tmp = args.copy()
    for n, i in enumerate(args_tmp):
      if i == 'img':
         args_tmp[n] = img
      elif i == resultname:
          args_tmp[n] = result
      elif i == 'temp1':
          args_tmp[n] = temp1
      elif i == 'temp2':
          args_tmp[n] = temp2

    pi2.run_command(opname, args_tmp)

    for n, i in enumerate(args):
      if i == resultname:
          pi2.writeraw(args_tmp[n], outfile_distributed)


    pi2.distribute(Distributor.NONE)



    # Run the operation without distributed processing
    img = pi2.read(infile)
    if convert_to_type != ImageDataType.UNKNOWN:
        pi2.convert(img, convert_to_type)

    result = pi2.newlike(img)
    pi2.set(result, img)

    temp1 = pi2.newimage(ImageDataType.FLOAT32)
    temp2 = pi2.newimage(ImageDataType.FLOAT32)

    args_tmp = args.copy()
    for n, i in enumerate(args_tmp):
      if i == 'img':
         args_tmp[n] = img
      elif i == resultname:
          args_tmp[n] = result
      elif i == 'temp1':
          args_tmp[n] = temp1
      elif i == 'temp2':
          args_tmp[n] = temp2

    pi2.run_command(opname, args_tmp)

    for n, i in enumerate(args):
      if i == resultname:
          pi2.writeraw(args_tmp[n], outfile_normal)



    # Compare results
    check_distribution_test_result(outfile_normal, outfile_distributed, opname, 'distributed and local', tolerance, convert_to_type)




def create_particle_labels_test():
    """
    Creates initial data for particle label growing test.
    """

    img = pi2.read(input_file('complicated_particles_1_38x36x21.raw'))

    pi2.set(img, [7, 7, 10], 2)
    pi2.set(img, [20, 18, 10], 3)
    pi2.set(img, [25, 25, 10], 4)

    pi2.writeraw(img, output_file('complicated_particles_point_labels'))




def trace_skeleton_test():
    """
    Tests skeleton tracing in normal and distributed mode.
    """
    

    # Normal mode
    skele = pi2.read(input_file("real_skele_200x200x200.raw"))
    vertices = pi2.newimage()
    edges = pi2.newimage()
    measurements = pi2.newimage()
    points = pi2.newimage()
    pi2.tracelineskeleton(skele, vertices, edges, measurements, points)

    pi2.writeraw(vertices, output_file('vertices'))
    pi2.writeraw(edges, output_file('edges'))
    pi2.writeraw(measurements, output_file('measurements'))
    pi2.writeraw(points, output_file('points'))

    skele = pi2.read(input_file("real_skele_200x200x200.raw"))
    pi2.fillskeleton(skele, vertices, edges, measurements, points, 1)
    pi2.writeraw(skele, output_file('filled_skeleton'))
    pi2.clear()




    # Distributed mode
    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(5)
    pi2.chunksize([100, 100, 100])

    skele = pi2.read(input_file("real_skele_200x200x200.raw"))
    vertices = pi2.newimage()
    edges = pi2.newimage()
    measurements = pi2.newimage()
    points = pi2.newimage()
    pi2.tracelineskeleton(skele, vertices, edges, measurements, points)

    pi2.writeraw(vertices, output_file('vertices_distributed'))
    pi2.writeraw(edges, output_file('edges_distributed'))
    pi2.writeraw(measurements, output_file('measurements_distributed'))
    pi2.writeraw(points, output_file('points_distributed'))

    skele = pi2.read(input_file("real_skele_200x200x200.raw"))
    pi2.fillskeleton(skele, vertices, edges, measurements, points, 1)
    pi2.writeraw(skele, output_file('filled_skeleton_distributed'))
    pi2.clear()

    pi2.distribute(Distributor.NONE)

    check_distribution_test_result(output_file('vertices'), output_file('vertices_distributed'), "tracelineskeleton", "vertices")
    # NOTE: Only vertices can be compared; the edge order might be different and probably the edges might be a bit different, too,
    # although bot are valid representations of the skeleton.
    #check_distribution_test_result(output_file('edges'), output_file('edges_distributed'), "tracelineskeleton", "edges")
    #check_distribution_test_result(output_file('measurements'), output_file('measurements_distributed'), "tracelineskeleton", "measurements")
    #check_distribution_test_result(output_file('points'), output_file('points_distributed'), "tracelineskeleton", "points", data_type='int32')
    # Skeleton filling should give the same result
    check_distribution_test_result(output_file('filled_skeleton'), output_file('filled_skeleton_distributed'), "tracelineskeleton", "filled skeleton")

def fill_skeleton_test():
    """
    Tests distributed and normal skeleton filling.
    """

    outfile_normal = output_file('filled_skeleton_length_normal')
    outfile_distributed = output_file('filled_skeleton_length_distributed')

    # First trace the skeleton.
    skele = pi2.read(input_file("real_skele_200x200x200.raw"))
    vertices = pi2.newimage()
    edges = pi2.newimage()
    measurements = pi2.newimage()
    points = pi2.newimage()
    pi2.tracelineskeleton(skele, vertices, edges, measurements, points)
    
    pi2.writeraw(vertices, output_file('vertices'))
    pi2.writeraw(edges, output_file('edges'))
    pi2.writeraw(measurements, output_file('measurements'))
    pi2.writeraw(points, output_file('points'))

    # Fill in normal mode
    skele = pi2.read(input_file("real_skele_200x200x200.raw"))
    pi2.fillskeleton(skele, vertices, edges, measurements, points, 1)
    pi2.writeraw(skele, outfile_normal)
    pi2.clear()

    # Fill in distributed mode
    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(5)

    vertices = pi2.read(output_file('vertices'))
    edges = pi2.read(output_file('edges'))
    measurements = pi2.read(output_file('measurements'))
    points = pi2.read(output_file('points'), 'int32') # NOTE: We need to specify the data type here as otherwise it defaults to uint32. This is required only for .raw files.
    skele = pi2.read(input_file("real_skele_200x200x200.raw"))
    pi2.fillskeleton(skele, vertices, edges, measurements, points, 1)
    pi2.writeraw(skele, outfile_distributed)

    pi2.distribute(Distributor.NONE)

    check_distribution_test_result(outfile_normal, outfile_distributed, 'fillskeleton', 'distributed and local', 0)



def morphorec_test():
    """
    Tests normal and distributed morphological reconstruction.
    """

    geom = pi2.read(input_file('complicated_particles_1_38x36x21.raw'))
    
    img = pi2.newlike(geom)
    pi2.set(img, [7, 7, 10], 2)
    pi2.set(img, [20, 18, 10], 3)
    pi2.set(img, [25, 25, 10], 4)
    pi2.writeraw(img, output_file('morphorec_seeds'))

    pi2.morphorec(img, geom)
    pi2.writeraw(img, output_file('morphorec_normal'))

    pi2.clear()
    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(0.025)
    geom = pi2.read(input_file('complicated_particles_1_38x36x21.raw'))
    img = pi2.read(output_file('morphorec_seeds'))
    pi2.morphorec(img, geom)
    pi2.writeraw(img, output_file('morphorec_distributed'))

    pi2.distribute(Distributor.NONE)

    check_distribution_test_result(output_file('morphorec_normal'), output_file('morphorec_distributed'), 'morphorec', 'distributed and local', 0)




def multimax_test(direction):
    """
    Tests normal and distributed max projection using second image.
    """


    img = pi2.read(input_file())
    imgVal = pi2.read(input_file_bin())

    maxproj = pi2.newimage(img.get_data_type())
    maxval = pi2.newimage(imgVal.get_data_type())
    pi2.maxproject(img, maxproj, direction, imgVal, maxval)

    pi2.writeraw(maxproj, output_file(f"max_projection_{direction}_normal"))
    pi2.writeraw(maxval, output_file(f"max_projection_value_{direction}_normal"))

    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(15)

    img = pi2.read(input_file())
    imgVal = pi2.read(input_file_bin())

    maxproj = pi2.newimage(img.get_data_type())
    maxval = pi2.newimage(imgVal.get_data_type())
    pi2.maxproject(img, maxproj, direction, imgVal, maxval)

    pi2.writeraw(maxproj, output_file(f"max_projection_{direction}_distributed"))
    pi2.writeraw(maxval, output_file(f"max_projection_value_{direction}_distributed"))

    pi2.distribute(Distributor.NONE)

    check_distribution_test_result(output_file(f"max_projection_{direction}_normal"), output_file(f"max_projection_{direction}_distributed"), "multimax projection", "distributed and local", 0)
    check_distribution_test_result(output_file(f"max_projection_value_{direction}_normal"), output_file(f"max_projection_value_{direction}_distributed"), "multimax value", "distributed and local", 0)




def generate_particles():
    """
    analyze_particles function uses this to generate an image containing some random particles.
    """

    import random

    pi2.echo(False, False)

    img = pi2.newimage(ImageDataType.UINT16, 500, 500, 500)

    # Generate some particles

    # Spheres
    #for i in range(0, 1000):
    #    pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
    #    r = random.randint(1, 20)
    #    col = random.randint(1, 255)
    #    pi2.sphere(img, pos, r, col)

    # Generic ellipsoids
    for i in range(0, 1000):
        pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        r = [random.randint(1, 20), random.randint(1, 20), random.randint(1, 20)]
        col = random.randint(1, 255)

        u1 = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        tmp = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        u2 = np.cross(u1, tmp)

        pi2.ellipsoid(img, pos, r, col, u1, u2)

    # Axis-aligned boxes
    #for i in range(0, 1000):
    #    pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
    #    size = [random.randint(1, 20), random.randint(1, 20), random.randint(1, 20)]
    #    col = random.randint(1, 255)
    #    pi2.box(img, pos, size, col)

    # Generic boxes
    for i in range(0, 1000):
        pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        r = [random.randint(1, 20), random.randint(1, 20), random.randint(1, 20)]
        col = random.randint(1, 255)

        u1 = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        tmp = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        u2 = np.cross(u1, tmp)

        pi2.box(img, pos, r, col, u1, u2)


    pi2.writeraw(img, output_file('particles'))

    pi2.clear()


def analyze_particles_local(analyzers):
    """
    Helper for analyze_particles demo.
    """

    # Read generated data file
    img = pi2.read(output_file('particles'))
    pi2.threshold(img, 0)

    # Analyze particles
    result = pi2.newimage(ImageDataType.FLOAT32)
    pi2.analyzeparticles(img, result, analyzers)

    # Draw ellipsoids
    img = pi2.read(output_file('particles'))
    pi2.threshold(img, 0)
    pi2.drawellipsoids(img, analyzers, result, 2)
    pi2.writeraw(img, output_file('ellipsoids_local'))

    pa = result.get_data()
    

    return pa

def analyze_particles_distributed(analyzers):
    """
    Helper for analyze_particles demo.
    """

    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(50)

    # Read generated data file
    img = pi2.read(output_file('particles'))
    pi2.threshold(img, 0)

    # Analyze particles
    result = pi2.newimage(ImageDataType.FLOAT32)
    pi2.analyzeparticles(img, result, analyzers)

     # Draw ellipsoids
    img = pi2.read(output_file('particles'))
    pi2.threshold(img, 0)
    pi2.drawellipsoids(img, analyzers, result, 2)
    pi2.writeraw(img, output_file('ellipsoids_distr'))

    pa = result.get_data()

    pi2.distribute(Distributor.NONE)

    return pa


def check_column(a, b, n, msg, tol):

    max_diff = 0
    for i in range(0, a.shape[0]):
        if np.isnan(a[i, n]) and np.isnan(b[i, n]):
            # OK
            pass
        elif np.isposinf(a[i, n]) and np.isposinf(b[i, n]):
            # OK
            pass
        elif np.isneginf(a[i, n]) and np.isneginf(b[i, n]):
            # OK
            pass
        else:
            diff = np.abs(a[i, n] - b[i, n])
            if diff > max_diff:
                max_diff = diff


    check_result(max_diff <= tol, f"Local and distributed mode particle analysis differs in column {msg}, max absolute difference = {max_diff}")


def analyze_particles():
    """
    Demonstrates particle (i.e, region or blob) analysis, and some drawing commands.
    """

#    generate_particles()

    # Show analyzer names
    pi2.listanalyzers()

    # Analyze particles locally and using distributed computing mode
    analyzers = 'volume coordinates bounds boundingsphere isonedge pca'

    # Show titles of data columns
    print('Titles of columns in data table:')
    pi2.headers(analyzers)

    pa_local = analyze_particles_local(analyzers)
    

    # Flood fill all particles so that only background should be left
    img = pi2.read(output_file('particles'))
    pi2.threshold(img, 0)
    V = pa_local[:, 0]
    x = pa_local[:, 1]
    y = pa_local[:, 2]
    z = pa_local[:, 3]

    for i in range(0, len(x)):
        pi2.floodfill(img, [x[i], y[i], z[i]], 0)

    M = img.get_data().max()

    pi2.writeraw(img, output_file('black'))

    check_result(M == 0, "Not all particles are in the results table.")



    print(f"Local processing returned results array of size {pa_local.shape}")

    pa_dist = analyze_particles_distributed(analyzers)
    print(f"Distributed processing returned results array of size {pa_dist.shape}")

    # Compare analysis results
    
    if not np.isclose(pa_local.shape, pa_dist.shape).all():
        check_result(False, 'Different number of results in local and distributed particle analysis.')
    else:
        # Particle order may be different in local and distributed processing, so sort both result arrays
        # before comparing.
        pa_local.sort(0)
        pa_dist.sort(0)

        # We need different tolerances for different columns!
        check_column(pa_local, pa_dist, 0, 'volume', 0)
        # Don't check X, Y, Z as those might differ as different points of the particle may be encountered first
        #check_column(pa_local, pa_dist, 1, 'X')
        #check_column(pa_local, pa_dist, 2, 'Y')
        #check_column(pa_local, pa_dist, 3, 'Z')
        check_column(pa_local, pa_dist, 4, 'minx', 0)
        check_column(pa_local, pa_dist, 5, 'maxx', 0)
        check_column(pa_local, pa_dist, 6, 'miny', 0)
        check_column(pa_local, pa_dist, 7, 'maxy', 0)
        check_column(pa_local, pa_dist, 8, 'minz', 0)
        check_column(pa_local, pa_dist, 9, 'maxz', 0)
        check_column(pa_local, pa_dist, 10, 'bounding sphere X', 0.01)
        check_column(pa_local, pa_dist, 11, 'bounding sphere Y', 0.01)
        check_column(pa_local, pa_dist, 12, 'bounding sphere Z', 0.01)
        check_column(pa_local, pa_dist, 13, 'bounding sphere radius', 0.01)
        check_column(pa_local, pa_dist, 14, 'is on edge', 0)
        check_column(pa_local, pa_dist, 15, 'CX', 0.01)
        check_column(pa_local, pa_dist, 16, 'CY', 0.01)
        check_column(pa_local, pa_dist, 17, 'CZ', 0.01)
        check_column(pa_local, pa_dist, 18, 'e', 0.01)
        check_column(pa_local, pa_dist, 19, 'l1', 0.01)
        check_column(pa_local, pa_dist, 20, 'l2', 0.01)
        check_column(pa_local, pa_dist, 21, 'l3', 0.01)
        # These angles are suprisingly inaccurate although the ellipsoids
        # are very nearly same (see the filling test)
        # 0.2 rad = 11.5 deg
        check_column(pa_local, pa_dist, 22, 'phi1', 0.2)
        check_column(pa_local, pa_dist, 23, 'theta1', 0.2)
        check_column(pa_local, pa_dist, 24, 'phi2', 0.2)
        check_column(pa_local, pa_dist, 25, 'theta2', 0.2)
        check_column(pa_local, pa_dist, 26, 'phi3', 0.2)
        check_column(pa_local, pa_dist, 27, 'theta3', 0.2)
        check_column(pa_local, pa_dist, 28, 'rmax', 0.01)
        check_column(pa_local, pa_dist, 29, 'd1', 1)
        check_column(pa_local, pa_dist, 30, 'd2', 1)
        check_column(pa_local, pa_dist, 31, 'd3', 1)
        check_column(pa_local, pa_dist, 32, 'bounding scale', 0.05)

    # Numerical precision causes small differences in the shapes of the ellipsoids, allow for that
    # This kind of expressions might be easier to write in NumPy as follows,
    # but then distributed processing possibility is of course lost:
    vis_local = pi2.read(output_file('ellipsoids_local')).get_data()
    vis_distr = pi2.read(output_file('ellipsoids_distr')).get_data()
    diff = np.abs(vis_local - vis_distr)
    diff_count = np.sum(diff > 0)
    check_result(diff_count < 100, "More than 100 pixels are different between ellipsoid visualizations made using local and distributed processing.")



def analyze_labels():

    generate_particles()

    img = pi2.read(output_file('particles'))
    analyzers = 'volume coordinates'
    pa_local_img = pi2.newimage(ImageDataType.FLOAT32)
    pi2.analyzelabels(img, pa_local_img, analyzers)
    pa_local = pa_local_img.get_data()

    # Fill all particles so that only background should be left
    V = pa_local[:, 0]
    x = pa_local[:, 1]
    y = pa_local[:, 2]
    z = pa_local[:, 3]

    data = img.get_data()

    for i in range(0, len(x)):
        pix = data[int(y[i]), int(x[i]), int(z[i])] # NOTE: the order of coordinates is y, x, z in the Numpy array!
        pi2.replace(img, pix, 0)

    M = img.get_data().max()

    pi2.writeraw(img, output_file('black'))

    check_result(M == 0, "Not all particles are in the results table.")



def dimension_broadcast():

    img = pi2.read(input_file())
    img2 = pi2.newlike(img, width=img.get_width(), height=img.get_height(), depth=1)
    pi2.set(img2, 2)
    pi2.divide(img, img2, True)

    pi2.writeraw(img, output_file('broadcasted_divide'))

    img_gt = pi2.read(input_file())
    pi2.divide(img_gt, 2)
    pi2.writeraw(img, output_file('normal_divide'))
    
    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(1)

    img = pi2.read(input_file())
    img2 = pi2.newlike(img, width=img.get_width(), height=img.get_height(), depth=1)
    pi2.set(img2, 2)
    pi2.divide(img, img2, True)

    pi2.writeraw(img, output_file('broadcasted_distributed_divide'))

    pi2.distribute(Distributor.NONE)

    check_result(output_file('broadcasted_divide'), output_file('normal_divide'), "Invalid broadcasted division result.")
    check_result(output_file('broadcasted_distributed_divide'), output_file('normal_divide'), "Invalid broadcasted distributed division result.")


def rotate():

    img = pi2.read(input_file())
    out = pi2.newimage()
    pi2.rot90cw(img, out)
    pi2.writeraw(out, output_file('rot90cw'))



def twoimage_distribution():

    img = pi2.read(input_file())
    out = pi2.newimage()

    pi2.finalizetmap(img, out)

    pi2.writeraw(out, output_file('finalizetmap_normal'))

    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(15)

    img = pi2.read(input_file())
    out = pi2.newimage()

    pi2.finalizetmap(img, out)

    pi2.writeraw(out, output_file('finalizetmap_distributed'))
    pi2.distribute(Distributor.NONE)

    check_distribution_test_result(output_file('finalizetmap_normal'), output_file('finalizetmap_distributed'), 'two-image input-output commands distribution', 'distributed and local')



def histogram():

    img = pi2.newimage(ImageDataType.FLOAT32, 200, 200, 200)
    pi2.ramp(img, 0)
    pi2.noise(img, 0, 10000000)
    pi2.writeraw(img, output_file('hist_input'))


    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(15)
    
    img = pi2.read(output_file('hist_input'))

    hst1 = pi2.newimage(ImageDataType.FLOAT32)
    bins1 = pi2.newimage(ImageDataType.FLOAT32)
    pi2.hist(img, hst1, bins1, 0.0, 1000.0, 500)

    hst2 = pi2.newimage(ImageDataType.FLOAT32)
    bins2 = pi2.newimage(ImageDataType.FLOAT32)
    pi2.hist(img, hst2, bins2, 0.0, 1000.0, 250)

    hdata1 = hst1.get_data()
    hdata2 = hst2.get_data()

    print(hdata1.sum())
    print(hdata2.sum())

    check_result(hdata1.sum() == hdata2.sum(), "Histogram normalization does not match.")


    pi2.distribute(Distributor.NONE)




def autothreshold():

    def single_test(method):

        img = pi2.read(input_file())
        pi2.linmap(img, 0, 600, 0, 255)
        pi2.convert(img, ImageDataType.UINT8)
        pi2.writetif(img, output_file('threshold_orig'))

        th = pi2.newlike(img)
        pi2.set(th, img)
        pi2.autothreshold(th, method)
        pi2.multiply(th, 255)
        pi2.writetif(th, output_file(f"threshold_{method}"))

    single_test(AutoThresholdMethod.OTSU)
    single_test(AutoThresholdMethod.HUANG)
    single_test(AutoThresholdMethod.INTERMODES)
    single_test(AutoThresholdMethod.ISODATA)
    single_test(AutoThresholdMethod.LI)
    single_test(AutoThresholdMethod.MAXENTROPY)
    single_test(AutoThresholdMethod.MEAN)
    single_test(AutoThresholdMethod.MINERROR)
    single_test(AutoThresholdMethod.MINIMUM)
    single_test(AutoThresholdMethod.MOMENTS)
    single_test(AutoThresholdMethod.PERCENTILE)
    single_test(AutoThresholdMethod.RENYI)
    single_test(AutoThresholdMethod.SHANBHAG)
    single_test(AutoThresholdMethod.TRIANGLE)
    single_test(AutoThresholdMethod.YEN)
    single_test(AutoThresholdMethod.MEDIAN)
    single_test(AutoThresholdMethod.MIDGREY)
    single_test(AutoThresholdMethod.NIBLACK)
    single_test(AutoThresholdMethod.PHANSALKAR)
    single_test(AutoThresholdMethod.SAUVOLA)
    single_test(AutoThresholdMethod.BERNSEN)



def tif_and_tiff():
    """
    Ensures that both .tif and .tiff sequences can be read.
    """

    img = pi2.read(input_file())

    pi2.writesequence(img, output_file("./sequence_tif/img_@.tif"))
    pi2.writesequence(img, output_file("./sequence_tiff/img_@.tiff"))

    img2 = pi2.read(output_file("./sequence_tif/img_@.tif"))
    img3 = pi2.read(output_file("./sequence_tiff/img_@.tiff"))

    M1 = calc_difference(img, img2)
    M2 = calc_difference(img, img3)

    # Check that the difference is zero
    check_result(M1 <= 0, f"ERROR: Difference in tif and tiff sequence reading.")
    check_result(M2 <= 0, f"ERROR: Difference in tif and tiff sequence reading.")


def get_pixels():
    """
    Checks that various ways of getting pixels from the (normal and distributed) images give the same result.
    """

    pi2.distribute(Distributor.NONE)
    # Generate image
    w = 100
    h = 80
    d = 120
    img = pi2.newimage(ImageDataType.UINT16, w, h, d)
    pi2.ramp(img, 0)
    pi2.noise(img, 0, 25)
    pi2.writeraw(img, output_file("ramp"))

    # Generate some positions
    N = 300
    #positions = pi2.newimage(ImageDataType.FLOAT32, 3, N)
    #pos = positions.get_data()
    pos = np.zeros([N, 3])
    for i in range(0, N):
        pos[i, 0] = random.randint(0, w-2) + 0.5
        pos[i, 1] = random.randint(0, h-1)
        pos[i, 2] = random.randint(0, d-1)
    #pos = np.array([[34.5, 42, 13]])

    # Get pixels at positions (non-distributed)
    img = pi2.read(output_file("ramp"))
    out = pi2.newimage(img.get_data_type())
    pi2.get(img, out, pos)
    data_normal = out.get_data()


    # Get pixels at positions (distributed)
    pi2.distribute(Distributor.LOCAL)
    pi2.maxmemory(0.25)

    img = pi2.read(output_file("ramp"))
    out = pi2.newimage(img.get_data_type())
    pi2.get(img, out, pos)
    data_distributed = out.get_data()

    pi2.distribute(Distributor.NONE)

    # Get pixels at positions (through NumPy)
    img = pi2.read(output_file("ramp"))
    pyimg = img.get_data()

    N = pos.shape[0]
    data_numpy = np.zeros(N)
    for i in range(0, N):
        data_numpy[i] = pyimg[int(pos[i][1] + 0.5), int(pos[i][0] + 0.5), int(pos[i][2] + 0.5)]

    check_result(np.isclose(data_normal, data_distributed).all(), "get pixel normal != distributed")
    check_result(np.isclose(data_normal, data_numpy).all(), "get pixel normal != numpy")

    print("point: normal = numpy = distributed")
    for i in range(0, N):
        print(f"{pos[i]}: {data_normal[i]} = {data_numpy[i]} = {data_distributed[i]}")



def set_pixels():

    img = pi2.newimage(ImageDataType.UINT8, 100, 200, 300)
    p = np.array([[0, 0, 0], [10, 20, 30], [70, 70, 70]], dtype=float)
    v = np.array([[10, 20, 30]], dtype=np.uint8)

    pos = pi2.newimage()
    pos.set_data(p)

    gt = pi2.newimage()
    gt.set_data(v)

    pi2.set(img, pos, gt)

    pi2.writetif(img, output_file("set_pixels"))


    v2 = pi2.newimage()
    pi2.get(img, v2, pos)

    

    M = calc_difference(gt, v2)
    
    
    check_result(M <= 0, "Pixels are not written and read correctly.")



def distributed_numpy():
    """
    Writing numpy arrays to disk caused exception in distributed mode.
    """

    pi2.distribute(Distributor.LOCAL)

    nparr = np.zeros([100, 100, 100])

    pi2.writeraw(nparr, output_file("np_distributed"))

    pi2.distribute(Distributor.NONE)


def memory():
    """
    Checks that image memory is freed when variables are cleared.
    """


    import os
    import psutil
    process = psutil.Process(os.getpid())

    mem = process.memory_info().rss
    check_result(mem < 250 * 1024 * 1024, "consuming too much memory before allocation of an image")

    img = pi2.newimage(ImageDataType.UINT8, 1000, 1000, 1000)

    mem = process.memory_info().rss
    check_result(mem > 500 * 1024 * 1024, "consuming too little memory after allocation of an image")

    pi2.clear(img.name)

    mem = process.memory_info().rss
    check_result(mem < 250 * 1024 * 1024, "consuming too much memory after de-allocation of an image")



def metadata():
    """
    Tests storing and retrieving image meta-data.
    """

    img = pi2.newimage(ImageDataType.UINT8, 1, 1, 1)
    str = pi2.newstring()

    pi2.setmeta(img, "key1", "value1")
    pi2.setmeta(img, "key2", "value2")
    pi2.getmeta(img, "key1", str)

    check_result(str.as_string() == "value1", "value")
    
    pi2.listmeta(img, str)
    check_result(str.as_string() == "key1, key2", "listmeta before")

    pi2.writemeta(img, output_file("metadata"))
    pi2.clearmeta(img)

    pi2.listmeta(img, str)
    check_result(str.as_string() == "", "meta key list not empty")

    pi2.readmeta(img, output_file("metadata"))
    pi2.listmeta(img, str)
    check_result(str.as_string() == "key1, key2", "meta key list contains incorrect values")

    pi2.getmeta(img, "key1", str)
    check_result(str.as_string() == "value1", "key1 after read")

    pi2.getmeta(img, "key2", str)
    check_result(str.as_string() == "value2", "key2 after read")

def named_variables():
    """
    Tests named non-image variables.
    """

    img = pi2.newimage(ImageDataType.UINT8, 1, 1, 1)
    str_key = pi2.newstring()
    str_value = pi2.newstring("VALUE")

    check_result(str_key.as_string() == "", "key before doing anything")
    check_result(str_value.as_string() == "VALUE", "value before doing anything")

    pi2.set(str_key, "KEY")

    check_result(str_key.as_string() == "KEY", "key after setting")

    pi2.setmeta(img, str_key, str_value)

    check_result(str_key.as_string() == "KEY", "key after setmeta")
    check_result(str_value.as_string() == "VALUE", "value after setmeta")

    pi2.set(str_value, "---")
    check_result(str_value.as_string() == "---", "value after set")

    pi2.getmeta(img, str_key, str_value)
    
    check_result(str_key.as_string() == "KEY", "key after getmeta")
    check_result(str_value.as_string() == "VALUE", "value after getmeta")


def set_overloads():

    img1 = pi2.newimage(ImageDataType.UINT16)
    img2 = pi2.newlike(img1)

    pi2.set(img2, img1)


def big_tiff_write():
    img = pi2.newimage(ImageDataType.UINT8, 5*1024, 1024, 1024)
    pi2.ramp(img, 0)
    pi2.writetif(img, output_file("big_tiff.tif"))
    pi2.writeraw(img, output_file("big_raw"))

def big_tiff():

    big_tiff_write()

    # Whoops: I don't have enough memory for this in my computer...
    img1 = pi2.read(output_file("big_tiff.tif"))
    img2 = pi2.read(output_file("big_raw"))
    M = calc_difference(img1, img2)
    check_result(M <= 0, f"ERROR: Images saved as .tif and .raw are not equal.")
    

def lz4_files():

    img1 = pi2.newimage(ImageDataType.FLOAT32, 100, 250, 200)
    pi2.ramp3(img1)
    pi2.writelz4(img1, output_file("lz4test"))
    pi2.writetif(img1, output_file("lz4test.tif"))

    img2 = pi2.read(output_file("lz4test.lz4raw"))
    check_result(calc_difference(img1, img2) == 0, "Image saved and read from .lz4raw dataset changed in the I/O process.")
    

def nn5_files():

    img1 = pi2.newimage(ImageDataType.FLOAT32, 100, 250, 200)
    pi2.ramp3(img1)
    pi2.writenn5(img1, output_file("nn5test"))
    pi2.writenn5(img1, output_file("nn5test_small_chunks"), [20, 40, 60])
    pi2.writetif(img1, output_file("nn5test.tif"))

    img2 = pi2.read(output_file("nn5test"))
    check_result(calc_difference(img1, img2) == 0, "Image saved and read from NN5 dataset changed in the I/O process.")

    img3 = pi2.read(output_file("nn5test_small_chunks"))
    check_result(calc_difference(img1, img3) == 0, "Image saved and read from NN5 small chunks dataset changed in the I/O process.")


def dead_pixels():

    img = pi2.newimage(ImageDataType.FLOAT32, 100, 100)
    x = [20, 20, 0]
    pi2.set(img, x, float('nan'))
    pi2.deadpixelremoval(img)
    val = img.get_data()[x[0], x[1]]
    check_result(val == 0, "dead pixel not removed")



def lsf_cluster():

    pi2.distribute(Distributor.LSF)
    pi2.maxmemory(1)
    img = pi2.newimage(ImageDataType.UINT16, 100, 100, 100)
    pi2.add(img, 10)
    pi2.writeraw(img, output_file("lsf/result"))

    pi2.distribute(Distributor.NONE)



# Enable or disable echoing of commands and timing info on screen
pi2.echo(True, False)



#morphorec_test()
#fill_skeleton_test()
#autothreshold()
#tif_and_tiff()
#get_pixels()
#set_pixels()
#distributed_numpy()
#named_variables()
#metadata()
#set_overloads()
#lz4_files()
#nn5_files()
#dead_pixels()
#trace_skeleton_test()
#multimax_test(0)
#multimax_test(1)
#multimax_test(2)
#generate_particles()
#analyze_particles()
#analyze_labels()
#dimension_broadcast()
#rotate()
#twoimage_distribution()
#histogram()
#for r in range(1, 10):
#    test_difference_normal_distributed('maxfilter', ['img', 'result', r, True, 'ellipsoidal', 'zero'])
#    test_difference_normal_distributed('maxfilter', ['img', 'result', r, True, 'ellipsoidal', 'nearest'])
#    test_difference_normal_distributed('minfilter', ['img', 'result', r, True, 'ellipsoidal', 'zero'])
#    test_difference_normal_distributed('minfilter', ['img', 'result', r, True, 'ellipsoidal', 'nearest'])
#for r in range(1, 10):
#    test_difference_normal_distributed('openingfilter', ['result', r, True, 'ellipsoidal', 'zero'])
#    test_difference_normal_distributed('openingfilter', ['result', r, True, 'ellipsoidal', 'nearest'])
#    test_difference_normal_distributed('closingfilter', ['result', r, True, 'ellipsoidal', 'zero'])
#    test_difference_normal_distributed('closingfilter', ['result', r, True, 'ellipsoidal', 'nearest'])
test_difference_normal_distributed('floodfill', ['img', [0, 0, 0], 100], 'img', input_file_bin(), maxmem=5)
test_difference_normal_distributed('gaussfilter', ['img', 'result', 2])
test_difference_normal_distributed('bin', ['img', 'result', 9])
test_difference_normal_distributed('scale', ['img', 'result', 0], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 0.2, False, Interpolation.NEAREST], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 0.2, True, Interpolation.NEAREST], tolerance=1, maxmem=35)
test_difference_normal_distributed('scale', ['img', 'result', 0.2, False, Interpolation.LINEAR], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 0.2, True, Interpolation.LINEAR], tolerance=1, maxmem=35)
test_difference_normal_distributed('scale', ['img', 'result', 0.2, False, Interpolation.CUBIC], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 0.2, True, Interpolation.CUBIC], tolerance=1, maxmem=35)
test_difference_normal_distributed('scale', ['img', 'result', 0.7, False, Interpolation.NEAREST], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 0.7, True, Interpolation.NEAREST], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 0.7, False, Interpolation.LINEAR], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 0.7, True, Interpolation.LINEAR], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 0.7, False, Interpolation.CUBIC], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 0.7, True, Interpolation.CUBIC], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 2.1, False, Interpolation.NEAREST], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 2.1, True, Interpolation.NEAREST], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 2.1, False, Interpolation.LINEAR], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 2.1, True, Interpolation.LINEAR], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 2.1, False, Interpolation.CUBIC], tolerance=1)
test_difference_normal_distributed('scale', ['img', 'result', 2.1, True, Interpolation.CUBIC], tolerance=1)
test_difference_normal_distributed('scalelabels', ['img', 'result', 2])
test_difference_normal_distributed('scalelabels', ['img', 'result', 3], maxmem=50)
test_difference_normal_distributed('crop', ['img', 'result', '[10, 20, 30]', '[110, 120, 130]'])
test_difference_normal_distributed('sumproject', ['img', 'result', 2])
test_difference_normal_distributed('meanproject', ['img', 'result', 1])
test_difference_normal_distributed('maxproject', ['img', 'result', 0])
test_difference_normal_distributed('sum', ['img', 'result', True], tolerance=0)
test_difference_normal_distributed('sum', ['img', 'result', True], tolerance=0.000001e8, convert_to_type=ImageDataType.FLOAT32)
test_difference_normal_distributed('mean', ['img', 'result'], tolerance=0.001)
test_difference_normal_distributed('maxval', ['img', 'result', True])
test_difference_normal_distributed('satofilter', ['img', 'result', 1], maxmem=100)
test_difference_normal_distributed('dualthreshold', ['img', 300, 500], 'img')
test_difference_normal_distributed('regionremoval', ['img', 500, 'all'], 'img', input_file_bin())
test_difference_normal_distributed('dmap2', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT16)
test_difference_normal_distributed('dmap', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT16)
test_difference_normal_distributed('dmap2', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT32)
test_difference_normal_distributed('dmap', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT32)
test_difference_normal_distributed('dmap2', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.FLOAT32)
test_difference_normal_distributed('dmap', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.FLOAT32)
test_difference_normal_distributed('dmap', ['img', 'img'], 'img', input_file_bin(), convert_to_type=ImageDataType.FLOAT32)
test_difference_normal_distributed('derivative', ['img', 'result', 1,  0], 'result', maxmem=20)
test_difference_normal_distributed('derivative', ['img', 'result', 1,  1], 'result', maxmem=20)
test_difference_normal_distributed('derivative', ['img', 'result', 1,  2], 'result', maxmem=20)
test_difference_normal_distributed('canny', ['img', 1,  10, 100], 'img', maxmem=20)
test_difference_normal_distributed('convert', ['img', 'out', ImageDataType.FLOAT32], 'out')
test_difference_normal_distributed('convert', ['img', 'img', ImageDataType.FLOAT32], 'img')
test_difference_normal_distributed('convert', ['img', ImageDataType.FLOAT32], 'img')
test_difference_normal_distributed('hist', ['img', 'result', 'temp1', 0, 1000, 1000], 'result')
test_difference_normal_distributed('hist', ['img', 'result', 'temp1', 0, 500, 500], 'result')
test_difference_normal_distributed('whist', ['img', 'img', 'result', 'temp1', 0, 500, 500], 'result', convert_to_type=ImageDataType.FLOAT32)
test_difference_normal_distributed('hist2', ['img', 0, 1000, 1000, 'img', 0, 100, 100, 'result', 'temp1', 'temp2'], 'result')
test_difference_normal_distributed('whist2', ['img', 0, 1000, 1000, 'img', 0, 100, 100, 'img', 'result', 'temp1', 'temp2'], 'result', convert_to_type=ImageDataType.FLOAT32)
test_difference_normal_distributed('danielsson2', ['img', 'result'], 'result', input_file('t1-head_bin_dmap_256x256x129.raw'))
test_difference_normal_distributed('maskedmean', ['img', 'result', 0], 'result', input_file_bin(), convert_to_type=ImageDataType.FLOAT32)
test_difference_normal_distributed('tmap', ['img', 'result', 0, False, False], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT16)
test_difference_normal_distributed('tmap', ['img', 'result', 0, False, False, '[50, 50, 50]'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT16)
create_particle_labels_test()
test_difference_normal_distributed('growlabels', ['img', 1, 0], 'img', output_file('complicated_particles_point_labels'), maxmem=0.03)
# We do this test in two parts (we still check only theta or phi and not both, but if one is ok, the other should be ok, too!)
test_difference_normal_distributed('cylinderorientation', ['img', 'result', 'result', 1, 1], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, maxmem=100)
test_difference_normal_distributed('cylinderorientation', ['img', 'result', 'result', 1, 1], 'img', input_file(), convert_to_type=ImageDataType.FLOAT32, maxmem=100)
test_difference_normal_distributed('rot90cw', ['img', 'result'], 'result', input_file(), maxmem=10)
test_difference_normal_distributed('rot90ccw', ['img', 'result'], 'result', input_file(), maxmem=10)
test_difference_normal_distributed('reslice', ['img', 'result', 'top'], 'result', input_file(), maxmem=3)
test_difference_normal_distributed('reslice', ['img', 'result', 'bottom'], 'result', input_file(), maxmem=3)
test_difference_normal_distributed('reslice', ['img', 'result', 'left'], 'result', input_file(), maxmem=3)
test_difference_normal_distributed('reslice', ['img', 'result', 'right'], 'result', input_file(), maxmem=3)
test_difference_normal_distributed('flip', ['img', 0], 'img', input_file(), maxmem=5)
test_difference_normal_distributed('flip', ['img', 1], 'img', input_file(), maxmem=5)
test_difference_normal_distributed('flip', ['img', 2], 'img', input_file(), maxmem=5)
# General rotations involve interpolation so we expect images to match only to a small tolerance.
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=10)
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [0, 0, 1], [128, 128, 64], [128, 128, 64]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=10)
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 0, 0], [128, 128, 64], [128, 128, 64]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=20)
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 1, 1], [128, 128, 64], [128, 128, 64]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=40)
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [0, 0, 1], [10, 10, 0], [10, 10, 0]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=10)
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 0, 0], [10, 10, 0], [10, 10, 0]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=40)
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 1, 1], [10, 10, 0], [10, 10, 0]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=40)
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [0, 0, 1], [10, 10, 0], [30, 20, 40]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=10)
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 0, 0], [10, 10, 0], [30, 20, 40]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=20)
test_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 1, 1], [10, 10, 0], [30, 20, 40]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=40)
test_difference_normal_distributed('meancurvature', ['img', 'result', 1], 'result', maxmem=50)
# Here we use tolerance as curvature seems to be sensitive to the origin of the calculation blocks in the distributed mode.
test_difference_normal_distributed('curvature', ['img', 10, 'result', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.FLOAT32, maxmem=30, tolerance=0.15)
test_difference_normal_distributed('set', ['img', [100, 100, 100], 100], 'img', maxmem=5)
test_difference_normal_distributed('sphere', ['img', [100, 100, 100], 100, 200], 'img', maxmem=5)
test_difference_normal_distributed('box', ['img', [100, 100, 100], [100, 80, 60], 200], 'img', maxmem=5)
test_difference_normal_distributed('ellipsoid', ['img', [100, 100, 100], [100, 80, 60], 200, [1, 1, 1], [1, 1, -1]], 'img', maxmem=5)
test_difference_normal_distributed('box', ['img', [100, 100, 100], [100, 80, 60], 200, [1, 1, 1], [1, 1, -1]], 'img', maxmem=5)
test_difference_normal_distributed('ramp', ['img', 2], 'img', maxmem=5)
test_difference_normal_distributed('ramp', ['img', 1], 'img', maxmem=5)
test_difference_normal_distributed('ramp', ['img', 0], 'img', maxmem=5)
test_difference_normal_distributed('line', ['img', [10, 20, 30], [100, 200, 100], 200], 'img', maxmem=5)
test_difference_normal_distributed('capsule', ['img', [10, 20, 30], [100, 200, 100], 10, 200], 'img', maxmem=5)
test_difference_normal_distributed('stddev', ['img', 'result'], 'result', maxmem=5, tolerance=1e-4)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.OTSU], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.HUANG], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.ISODATA], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.LI], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MAXENTROPY], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MEAN], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MINERROR], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MOMENTS], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.PERCENTILE], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.RENYI], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.SHANBHAG], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.TRIANGLE], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.YEN], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MEDIAN], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MIDGREY], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.NIBLACK], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.PHANSALKAR], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.SAUVOLA], 'img', maxmem=5)
test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.BERNSEN], 'img', maxmem=5)
test_difference_normal_distributed('localthreshold', ['img', 'result', [2, 2, 2], AutoThresholdMethod.OTSU], 'result', maxmem=5)
test_difference_normal_distributed('add', ['result', 'img'], 'result', maxmem=5)
test_difference_normal_distributed('eval', ['77', 'img'], 'img', maxmem=5)
test_difference_normal_distributed('eval', ['x0+x1', 'result', 'img'], 'result', maxmem=5)
test_difference_normal_distributed('eval', ['x0+x1*x2', 'result', 'img', 'img'], 'result', maxmem=5)
test_difference_normal_distributed('ramp3', ['img'], 'img', maxmem=5)
infile = input_file()
test_difference_delaying('delaying_1', f"read(img, {infile}); gaussfilter(img, out, 1); add(img, 100); clear(img); convert(out, conv, uint8);", 'conv');
test_difference_delaying('delaying_2', f"read(img, {infile}); add(img, 100); crop(img, out, [0,0,0], [100, 100, 120]); subtract(out, 100);", 'out');
test_difference_delaying('delaying_3', f"read(img, {infile}); convert(img, img32, float32); clear(img); threshold(img32, 5e-4); convert(img32, cyl, uint8); clear(img32);", 'cyl');
test_difference_delaying('delaying_4', f"read(img, {infile}); convert(img, img32, float32); clear(img); cylindricality(img32, 0.5, 0.5); threshold(img32, 5e-4); convert(img32, cyl, uint8); clear(img32);", 'cyl', maxmem=100);





## These commands are not guaranteed to give the same result in distributed and local processing
## mode so they are not tested here. Despite that, they should give a valid result in both cases, but their
## result is not unique.
##test_difference_normal_distributed('surfacethin', ['result'])
##test_difference_normal_distributed('surfaceskeleton', ['result'])
##test_difference_normal_distributed('linethin', ['result'])
##test_difference_normal_distributed('lineskeleton', ['result'])
## This test involves pretty large dataset
##test_difference_normal_distributed('scalelabels', ['img', 'result', 4])
## This test will never succeed as the approximate bilateral filtering algorithm uses random numbers.
##test_difference_normal_distributed('bilateralfilterapprox', ['img', 'result', 5, 200], 'result') 
## These two tests do not succeed as the default settings assume the data is spread on the whole 16-bit value range
##test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.INTERMODES], 'img', maxmem=5)
###test_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MINIMUM], 'img', maxmem=5)
## This test requires a lot of memory, currently at least 2 * 20 GB
##big_tiff()
## This is not 100 % reliable test for some reason. Sometimes it fails for no reason.
##memory()
## This test requires LSF cluster.
##lsf_cluster()


pi2.timing()

print(f"{total_tests} checks run.")
print(f"{failed_tests} checks failed.")

#input("Press Enter to continue...")
