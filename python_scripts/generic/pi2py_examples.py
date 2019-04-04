
import pi2py

pi2 = pi2py.Pi2()

# Define convenience functions that return input and output file names.
# This is just to avoid copying the paths to all the examples in case they change.
def input_file(filename='t1-head_256x256x129.raw'):
    return '../testing/' + filename

def input_file_bin():
    return '../testing/t1-head_bin_256x256x129.raw'

def output_file(name):
    return '../testing/pi2py/' + name



def information():
    """
    Demonstrates use of info and list commands.
    """

    pi2.info()

    # Create some images and list them on screen
    pi2.newimage('B-img8', 'uint8', 512, 512, 512)
    pi2.newimage('A-img16', 'uint16', 512, 512, 512)
    pi2.newimage('C-img32', 'float32', 512, 512, 512)
    pi2.newimage('E-imgc32', 'complex32', 512, 512, 512)
    pi2.list()


def read_image():

    # .tif image
    pi2.read('img', input_file('t1-head.tif'))

    # .png image
    pi2.read('img', input_file('uint8.png'))

    # .raw file
    pi2.read('img', input_file('t1-head_bin_'))

    # Save sequence sequence
    pi2.writesequence('img', output_file('binary_head/head_@(3)_bin.tif'))

    # Read sequence back. Notice that reading does not support (3) directive in the sequence template!
    pi2.readsequence('img2', output_file('binary_head/head_@_bin.tif'))




def add_and_subtract():
    """
    Demonstrates use of add and subtract commands by adding and subtracting image from itself.
    """

    pi2.readraw('img', input_file())
    pi2.add('img', 'img')
    pi2.writeraw('img', output_file('head_added_to_itself'))

    pi2.readraw('img', input_file())
    pi2.subtract('img', 'img')
    pi2.writeraw('img', output_file('head_subtracted_from_itself'))


def convert_format():

    pi2.readvol('img', input_file('simple_structures.vol'))
    pi2.writeraw('img', output_file('vol_to_raw'))


def gauss():
    """
    Calculates Gaussian blur of image normally and using distributed processing.
    Calculates difference of the two versions.
    """

    # Gaussian filtering (distributed processing)
    pi2.distribute('local')
    pi2.readraw('img', input_file())
    pi2.gaussfilter('img', 'filtered', 5)
    pi2.writeraw('filtered', output_file('head_gauss_distributed'))
    pi2.distribute('')

    # Gaussian filtering (local processing)
    pi2.readraw('img', input_file())
    pi2.gaussfilter('img', 'filtered', 5)
    pi2.writeraw('filtered', output_file('head_gauss_normal'))


    # Compare results of normal and distributed processing
    pi2.readraw('img', output_file('head_gauss_normal'))
    pi2.readraw('img2', output_file('head_gauss_distributed'))

    pi2.convert('img', 'imgf', 'float32')
    pi2.convert('img2', 'img2f', 'float32')

    pi2.subtract('imgf', 'img2f')
    pi2.writeraw('imgf', output_file('gauss_difference'))


def sequence():
    """
    Demonstrates use of image sequence in distributed processing.
    """

    pi2.readraw('img', input_file())
    pi2.writesequence('img', output_file('head_sequence'))
    pi2.clear('img')

    pi2.distribute('local')
    pi2.readsequence('img', output_file('head_sequence'))
    pi2.add('img', 500)
    pi2.writesequence('img', output_file('head_sequence_added'))
    pi2.writeraw('img', output_file('head_sequence_raw_added')) # TODO: At the moment writing to file with same prefix and two formats might not work.
    pi2.distribute('')


def skeletonize():
    """
    Demonstrates skeletonization
    """

    pi2.distribute('local')
    pi2.readraw('img', input_file())
    pi2.hybridskeleton('img')
    pi2.writeraw('img', output_file('head_skeleton'))
    pi2.distribute('')


def calc_difference(img1, img2):
    """
    Determines maximum absolute difference between img1 and img2 in pi system.
    """

    pi2.convert(img1, 'imgf', 'float32')
    pi2.convert(img2, 'img2f', 'float32')

    pi2.subtract('imgf', 'img2f')
    pi2.abs('imgf')
    pi2.maxval('imgf', 'M')
    M = pi2.getvalue('M')

    pi2.clear('imgf')
    pi2.clear('img2f')
    pi2.clear('M')

    return M



def get_set_data():
    """
    Demonstrates how to get access to image data stored in the pi system.
    """

    # Read image
    pi2.readraw('orig', input_file())

    pi2.crop('orig', 'img', [0, 0, 0], [200, 256, 128])
    print(f"Shape of image in pi = [200, 256, 128]")

    # Do processing in pi...

    # Get image data as NumPy array.
    # This does not work for distributed images, though.
    data = pi2.getdata('img')

    print(f"Shape of returned NumPy array = {data.shape}")

    # Now you can do whatever you want with the image.
    # Changes to the data are reflected in pi2.
    # The NumPy array is valid as long as the image is not resized or data type changed pi2,
    # i.e. at least as long as no pi2 commands are run.
    from PIL import Image
    import numpy as np
    slice = data[:, :, 64]
    im = Image.fromarray(slice)
    im.save(output_file("test.tif"))


    # Value of image in pi can be set like this:
    pi2.setdata('slice', slice)
    pi2.writeraw('slice', output_file('slice'))

    # We can also set 3D array
    pi2.setdata('substack', data[50:200, 50:180, 20:100])
    pi2.writeraw('substack', output_file('substack'))

    # Check that getting and setting data without modifications does not change the image in pi:
    data = pi2.getdata('img')
    data = np.copy(data)
    #data[0] = 250 # The error below should be triggered if you enable this line (modification of pixel values).
    pi2.setdata('img2', data)

    diff = calc_difference('img', 'img2')

    if abs(diff) > 0.0001:
        raise RuntimeError('Difference between unmodified images.')



def numpy_arrays_as_images():
    """
    Demonstrates usage of Numpy arrays as images in the pi system.
    """

    import numpy as np

    # Create NumPy array, remember to set data type to desired value (that is compatible with pi2)
    im = np.zeros((100, 200, 10), dtype=np.float32)

    # Call pi2 functions using numpy array as image.
    pi2.ramp(im, 1)
    pi2.writeraw(im, output_file('ramp_to_numpy_array'))


    # 2D arrays can also be used.
    im = np.zeros((100, 200, 1), dtype=np.float32)
    pi2.noise(im, 0, 50)
    pi2.writeraw(im, output_file('noise_to_numpy_array_2d'))


def network_analysis():
    """
    Calculates skeleton, traces it, and shows the result as a Mayavi 3D graph.
    """
    if True:
        # Read input file and calculate skeleton
        pi2.readraw('skeleton', input_file_bin())
        pi2.hybridskeleton('skeleton')

        # Optionally save the skeleton/load it
        pi2.writeraw('skeleton', output_file('skele'))

        # Enable this to use distributed processing
        #pi2.distribute('local')
        pi2.readraw('skeleton', output_file('skele'))

        # Trace skeleton (with branch geometry measurements)
        #pi2.readraw('img', input_file_bin())
        #pi2.tracelineskeleton('skeleton', 'img', 'vertices', 'edges', 'measurements')

        # Trace skeleton (without branch geometry measurements)
        # This is enabled by default because we don't need branch geometry measurements in hte following code.
        pi2.tracelineskeleton('skeleton', 'vertices', 'edges', 'measurements')

        # Save skeleton to disk
        pi2.writeraw('vertices', output_file('verts'))
        pi2.writeraw('edges', output_file('edges'))
        pi2.writeraw('measurements', output_file('measurements'))

        # Switch to local execution
        pi2.distribute('')

    # Load edges and vertices back
    pi2.readraw('vertices', output_file('verts'), 'float32')
    pi2.readraw('edges', output_file('edges'), 'uint64')

    # Extract vertices and edges.
    verts = pi2.getdata('vertices')
    edges = pi2.getdata('edges')

    import plot_network

    # Vertex coordinates are divided by 100 to make the figure smaller without zooming.
    plot_network.draw_graph_simple(verts / 100, edges)


def linefilters():
    """
    Calculates Frangi line-enhancing filter.
    """

    pi2.readraw('img', input_file())

    # Frangi filter, 16-bit result
    pi2.frangifilter('img', 'frangi_16bit', 3.0)
    pi2.writeraw('frangi_16bit', output_file('frangi_16bit'))


    # Frangi filter, 32-bit result. Notice different scaling in the output, compared to the 16-bit version.
    pi2.convert('img', 'img32', 'float32')
    pi2.frangifilter('img32', 'frangi_32bit', 3.0)
    pi2.writeraw('frangi_32bit', output_file('frangi_32bit'))



def binning_scaling():
    """
    Demonstrated distributed binning. Makes an image scaled to one half the size of the original in each dimension.
    """

    pi2.readraw('img', input_file())
    pi2.bin('img', 'img2', 2)
    pi2.writeraw('img2', output_file('binning_2'))



def canny_edge_detection():

    #pi2.distribute('local')

    pi2.readraw('img', input_file())
    pi2.canny('img', 1, 50, 100)
    pi2.writeraw('img', output_file('canny'))

    #pi2.clear()
    #pi2.distribute('')



def generate_particles():
    """
    analyze_particles function uses this to generate an image containing some random particles.
    """

    import random

    pi2.echo(False, False)

    pi2.newimage('img', 'uint8', 500, 500, 500)

    # Generate some particles
    for i in range(0, 1000):
        pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        r = random.randint(1, 20)
        pi2.sphere('img', pos, r, 255)

    for i in range(0, 1000):
        pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        size = [random.randint(1, 20), random.randint(1, 20), random.randint(1, 20)]
        pi2.box('img', pos, size, 255)

    pi2.writeraw('img', output_file('particles'))

    pi2.clear()


def analyze_particles_local(analyzers):
    """
    Helper for analyze_particles demo.
    """

    # Read generated data file
    pi2.readraw('img', output_file('particles'))

    # Analyze particles
    pi2.analyzeparticles('img', 'result', analyzers)

    # Show titles of data columns
    print('Titles of columns in data table:')
    pi2.headers(analyzers)

    import numpy
    pa = pi2.getdata('result')

    # Copy the result to new array as it is a reference to pi2 system and if we run
    # pi2.clear, the reference becomes invalid.
    pa = numpy.copy(pa)

    pi2.clear()

    return pa

def analyze_particles_distributed(analyzers):
    """
    Helper for analyze_particles demo.
    """

    pi2.distribute('local')

    # Read generated data file
    pi2.readraw('img', output_file('particles'))

    # Analyze particles
    pi2.analyzeparticles('img', 'result', analyzers)

    import numpy
    pa = pi2.getdata('result')
    pa = numpy.copy(pa)

    pi2.clear()
    pi2.distribute('off')

    return pa


def analyze_particles():
    """
    Demonstrates particle (i.e, region or blob) analysis, and some drawing commands.
    """

    generate_particles()

    # Show analyzer names
    pi2.listanalyzers()

    # Analyze particles locally and using distributed computing mode
    analyzers = 'volume coordinates bounds boundingsphere isonedge'
    pa_local = analyze_particles_local(analyzers)
    pa_dist = analyze_particles_distributed(analyzers)

    print(f"Local processing returned results array of size {pa_local.shape}")
    print(f"Distributed processing returned results array of size {pa_dist.shape}")

    # Compare analysis results, raise error if they are not the same.
    import numpy

    if not numpy.isclose(pa_local.shape, pa_dist.shape).all():
        raise ValueError('Different number of results in local and distributed particle analysis.')

    # Particle order may be different in local and distributed processing, so sort both result arrays
    # before comparing.
    pa_local.sort(0)
    pa_dist.sort(0)

    if not numpy.isclose(pa_local, pa_dist).all():
        raise ValueError('Difference in results of local and distributed particle analysis.')
    else:
        print(f"Contents of the two results arrays are the same, ok.")

    # Plot volume histogram.
    # Volume is in the first column (see output of pi2.headers command)
    volume = pa_local[:, 0]

    import matplotlib.pyplot as plt
    plt.hist(volume, normed=True, bins=30)
    plt.xlabel('Volume [pixel]')
    plt.ylabel('Probability')
    plt.show(block=False)

    pi2.clear()


def fill_particles():
    """
    Demonstrates analysis and filling of particles.
    """

    analyzers = 'volume coordinates bounds boundingsphere isonedge'

    #generate_particles()


    # Read generated data file
    pi2.readraw('img', output_file('particles'))

    # Analyze particles

    pi2.analyzeparticles('img', 'results', analyzers)

    # Result can be saved as an image
    pi2.writeraw('results', output_file('unfiltered_results'))

    # ...or read to Python
    import numpy
    pa = pi2.getdata('results')
    pa = numpy.copy(pa)

    print(f"Before filtering the shape of the results matrix is {pa.shape}")

    # Filter the measurement results using any Python methods you want
    # Command pi2.headers(analyzers) can be used to check the column names.
    # Particle volume will be in the first column, i.e.
    # pa[:, 0] gives volumes of all particles.
    limit = 20000
    pa = pa[pa[:, 0] > limit, :]

    print(f"After filtering the shape of the results matrix is {pa.shape}")

    # Push the filtered results back to pi2
    pi2.setdata('results', pa)

    # Save again just for fun (and to compare the output files if that is wanted)
    pi2.writeraw('results', output_file('filtered_results'))

    # Fill the particles that we left into the pa array (i.e. the big ones)
    pi2.fillparticles('img', analyzers, 'results', 2)
    pi2.writeraw('img', output_file('large_colored_local'))

    pi2.clear()


    # Now try similar process with distributed processing enabled
    pi2.distribute('local')

    pi2.readraw('img', output_file('particles'))
    pi2.setdata('results', pa)
    pi2.fillparticles('img', analyzers, 'results', 2)
    pi2.writeraw('img', output_file('large_colored_distributed'))
    pi2.clear()
    pi2.distribute('none')

    # Check difference
    check_distribution_test_result(output_file('large_colored_local'), output_file('large_colored_distributed'), "fillparticles")



def test_difference(opname, args, resultname='result', infile=input_file(), tolerance=0.00001):
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
    argstr = argstr.replace(',', '-') # TODO: Currently there is a bug that makes it impossible to have commas in pi2 arguments.
    outfile_normal = output_file(f"head_{opname}_{argstr}_normal")
    outfile_distributed = output_file(f"head_{opname}_{argstr}_distributed")

    # Calculate the operation locally
    pi2.readraw('img', infile)
    pi2.set('result', 'img')
    pi2.run_command(opname, args)
    pi2.writeraw(resultname, outfile_normal)
    pi2.clear()



    # Calculate the operation using distributed processing
    pi2.distribute('local')
    pi2.readraw('img', infile)
    pi2.set('result', 'img')
    pi2.run_command(opname, args)
    pi2.writeraw(resultname, outfile_distributed)
    pi2.clear()
    pi2.distribute('')


    check_distribution_test_result(outfile_normal, outfile_distributed, opname, tolerance)

    pi2.clear()


total_tests = 0
failed_tests = 0

def check_distribution_test_result(file1, file2, operation_name, tolerance=0.00001):

    # Load both results and calculate maximum absolute difference
    pi2.readraw('img', file1)
    pi2.readraw('img2', file2)

    M = calc_difference('img', 'img2')

    global total_tests
    total_tests = total_tests + 1

    global failed_tests
    # Check that the difference is zero
    print('-------')
    print('Result:')
    if M > tolerance:
        print(f"ERROR: Found difference in results of distributed and local processing while testing {operation_name}. Max absolute difference = {M}.")
        failed_tests = failed_tests + 1
        #pi2.writeraw('imgf', output_file(f"head_{operation_name}_difference"))
    else:
        print('OK')
    print('-------')







# Enable or disable echoing of commands and timing info on screen
pi2.echo(True, False)


# Run some examples
information()
#read_image()
#add_and_subtract()
#gauss()
#sequence()
#get_set_data()
#numpy_arrays_as_images()
#skeletonize()
#network_analysis()
#linefilters()
#binning_scaling()
#canny_edge_detection()
#analyze_particles()
#fill_particles()
#convert_format()

# Test that local and distributed processing give the same results
#test_difference('gaussfilter', ['img', 'result', 2])
#test_difference('hybridthin', ['result'])
#test_difference('hybridskeleton', ['result'])
#test_difference('bin', ['img', 'result', 3])
#test_difference('scale', ['img', 'result', 0.7], tolerance = 1)
#test_difference('scale', ['img', 'result', 2.1], tolerance = 1)
#test_difference('scale', ['img', 'result', 2.1, pi2.NEAREST_INTERPOLATION], tolerance = 1)
#test_difference('sumproject', ['img', 'result', 2])
#test_difference('meanproject', ['img', 'result', 1])
#test_difference('maxproject', ['img', 'result', 0])
#test_difference('sum', ['img', 'result'], tolerance=100) # Large tolerance as the value of the sum is approx. 3.5e8
#test_difference('mean', ['img', 'result'], tolerance=0.001)
#test_difference('maxval', ['img', 'result', True])
#test_difference('crop', ['img', 'result', '[10, 20, 30]', '[110, 120, 130]'])
#test_difference('satofilter', ['img', 'result', 1])
#test_difference('canny', ['img', 1,  10, 100], 'img')
#test_difference('dualthreshold', ['img', 300, 500], 'img')
#test_difference('regionremoval', ['img', 500], 'img', input_file_bin())

print(f"{total_tests} tests run.")
print(f"{failed_tests} tests failed.")

#input("Press Enter to continue...")
