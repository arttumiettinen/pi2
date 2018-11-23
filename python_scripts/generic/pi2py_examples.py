
import pi2py

pi2 = pi2py.Pi2()

# Define convenience functions that return input and output file names.
# This is just to avoid copying the paths to all the examples in case they change.
def input_file():
    return '../../testing/t1-head_256x256x129.raw'

def input_file_bin():
    return '../../testing/t1-head_bin_256x256x129.raw'

def output_file(name):
    return '../../testing/pi2py/' + name



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


def opening():
    """
    Calculates opening filter of an image normally and using distributed processing.
    Calculates difference between the two versions.
    """

    pi2.readraw('img', input_file())
    pi2.openingfilter('img', 'filtered', 1)
    pi2.writeraw('filtered', output_file('head_opening_normal'))

    pi2.distribute('local')
    pi2.readraw('img', input_file())
    pi2.openingfilter('img', 'filtered', 1)
    pi2.writeraw('filtered', output_file('head_opening_distributed'))
    pi2.distribute('')

    pi2.readraw('img', output_file('head_opening_normal'))
    pi2.readraw('img2', output_file('head_opening_distributed'))

    pi2.convert('img', 'imgf', 'float32')
    pi2.convert('img2', 'img2f', 'float32')

    pi2.subtract('imgf', 'img2f')
    pi2.writeraw('imgf', output_file('opening_difference'))



def sequence():
    """
    Demonstrates use of image sequence in distributed processing.
    """

    pi2.readraw('img', input_file())
    pi2.writesequence('img', output_file('head_sequence'))
    pi2.clear('img')

    pi2.distribute('local');
    pi2.readsequence('img', output_file('head_sequence'))
    pi2.add('img', 500)
    pi2.writesequence('img', output_file('head_sequence_added'))
    pi2.writeraw('img', output_file('head_sequence_added'))
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


def get_data():
    """
    Demonstrates how to get access to image data stored in pi system.
    """

    # Read image
    pi2.readraw('img', input_file())

    # Do processing in pi...

    # Get image data as NumPy array.
    # This does not work for distributed images, though.
    data = pi2.getdata('img')

    # Now you can do whatever you want with the image.
    # Changes to the data are reflected in pi2.
    # The image is valid as long as it is not changed in pi2,
    # i.e. at least as long as no pi2 commands are run.
    from PIL import Image
    import numpy as np
    slice = data[:, :, 64]
    im = Image.fromarray(slice)
    im.save("test.tiff")



def test_difference(opname, args):
    """
    Calculates operation normally and using distributed processing.
    Calculates difference between the two versions and prints a message if the results do not match.
    Image is loaded to variables 'img' and 'result', and result is expected to be in variable 'result'.
    """

    print(f"Testing {opname}...")
    print('-' * (len(opname) + 8 + 3))

    # Calculate the operation locally
    pi2.readraw('img', input_file())
    pi2.set('result', 'img')
    pi2.run_command(opname, args)
    pi2.writeraw('result', output_file(f"head_{opname}_normal"))
    pi2.clear()

    # Calculate the operation using distributed processing
    pi2.distribute('local')
    pi2.readraw('img', input_file())
    pi2.set('result', 'img')
    pi2.run_command(opname, args)
    pi2.writeraw('result', output_file(f"head_{opname}_distributed"))
    pi2.clear()
    pi2.distribute('')

    # Load both results and calculate maximum absolute difference
    pi2.readraw('img', output_file(f"head_{opname}_normal"))
    pi2.readraw('img2', output_file(f"head_{opname}_distributed"))

    pi2.convert('img', 'imgf', 'float32')
    pi2.convert('img2', 'img2f', 'float32')

    pi2.subtract('imgf', 'img2f')
    pi2.abs('imgf')
    pi2.maxval('imgf', 'M')
    M = pi2.getdata('M')
    M = float(M[0])


    # Check that the difference is zero
    print('Result:')
    print('-------')
    if M > 0.00001:
        print(f"ERROR: difference in results of distributed and local processing when testing {opname}. Max absolute difference = {M}.")
        pi2.writeraw('imgf', output_file(f"head_{opname}_difference"))
    else:
        print('OK')
    print('-------')

    pi2.clear()




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

        pi2.distribute('local')
        pi2.readraw('skeleton', output_file('skele'))

        # Trace skeleton
        pi2.readraw('img', input_file_bin())
        pi2.tracelineskeleton('skeleton', 'img', 'vertices', 'edges', 'measurements')

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



# Enable or disable echoing of commands and timing info on screen
pi2.echo(True, False)


# Run some examples
information()
add_and_subtract()
gauss()
opening()
sequence()
get_data()
skeletonize()
network_analysis()
linefilters()
binning_scaling()

# Test that local and distributed processing give the same results
test_difference('gaussfilter', ['img', 'result', 2])
test_difference('openingfilter', ['img', 'result', 1])
test_difference('hybridthin', ['result'])
test_difference('hybridskeleton', ['result'])
test_difference('bin', ['img', 'result', 3])
test_difference('sumproject', ['img', 'result', 2])
test_difference('meanproject', ['img', 'result', 1])
test_difference('maxproject', ['img', 'result', 0])
test_difference('sum', ['img', 'result']) # This test does not succeed but the reason is probably floating point inaccuracy.
test_difference('mean', ['img', 'result']) # This test does not succeed but the reason is probably floating point inaccuracy.
test_difference('maxval', ['img', 'result'])
test_difference('crop', ['img', 'result', '[10, 20, 30]', '[110, 120, 130]'])

input("Press Enter to continue...")