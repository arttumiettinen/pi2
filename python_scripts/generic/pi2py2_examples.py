"""
This file contains examples of how to use pi2py2 library.
"""








import numpy as np
from pi2py2 import *
pi = Pi2()



# Define convenience functions that return input and output file names.
# This is just to avoid copying the paths to all the examples in case they change.
def input_file(filename='t1-head_256x256x129.raw'):
    return '../../testing/input_data/' + filename

def input_file_bin():
    return '../../testing/input_data/t1-head_bin_256x256x129.raw'

def output_file(name):
    return '../../testing/pi2py2/' + name






def create_and_access_images():
    """
    Demonstrates creation of image and accessing its data.
    """

    # newimage command is used to create new images
    img1 = pi.newimage(ImageDataType.UINT8, 10, 20, 30)
    print(f"Width = {img1.get_width()}")
    print(f"Height = {img1.get_height()}")
    print(f"Depth = {img1.get_depth()}")
    print(f"Data type = {img1.get_data_type()}")
    print(f"All in one: {img1}")

    # newlike command can be used to create new image that is similar to existing image,
    # This line creates second image that has the same dimensions than img1 and pixel
    # data type UINT16.
    img2 = pi.newlike(img1, ImageDataType.UINT16);
    print(f"img2 is {img2}")

    # This line creates second image that has the same pixel data type than img1,
    # but its dimensions are 50x50x50.
    img3 = pi.newlike(img1, ImageDataType.UNKNOWN, 50, 50, 50);
    print(f"img3 is {img3}")

    # The data in the image can be retrieved as a NumPy array
    data = img1.get_data()
    print(f"When converted to a NumPy array, the shape of the image is {data.shape} and data type is {data.dtype}.")

    # Image data can also be set from a NumPy array
    data = np.eye(100, 100)
    img1.set_data(data)

    # This writes img1 to disk as a .raw file.
    # The dimensions of the image and the .raw suffix are automatically appended to
    # the image name given as second argument.
    pi.writeraw(img1, output_file("img1"))

    # NumPy arrays can be used directly as input in commands.
    # Changes made by Pi2 are NOT reflected in the NumPy arrays as
    # in some cases that would require re-shaping of the arrays, and that
    # does not seem to be wise...
    pi.writetif(data, output_file("numpy_array_tif"))




def read_and_write_image():
    """
    Demonstrates reading and writing images.
    """

    # .tif image
    img1 = pi.read(input_file('t1-head.tif'))

    # .png image
    img2 = pi.read(input_file('uint8.png'))

    # .raw file
    img3 = pi.read(input_file('t1-head_bin_'))

    # Save sequence. The @(3) declares that the file name should contain 3 numeric digits in the place of @, e.g. head_001_bin.tif.
    pi.writesequence(img3, output_file('binary_head/head_@(3)_bin.tif'))

    # Read sequence back. Notice that reading does not support (3) directive in the sequence template!
    img4 = pi.read(output_file('binary_head/head_@_bin.tif'))

    # Images can also be read into existing images.
    pi.read(output_file('binary_head/head_@_bin.tif'), target_image=img1)

    # When saving .raw files, the dimensions are appended to the file name
    pi.writeraw(img2, output_file('uint8_as_raw'))

    # When reading .raw files, it is not necessary to give the dimensions and .raw suffix
    # if the given part of the file name identifies an unique .raw file
    img2 = pi.read(output_file('uint8_as_raw'))



def help():
    """
    Demonstrates Pi2 help commands.
    IPython docstrings contain the same information.
    """

    # Show generic info
    pi.info()

    # Show license(s) related to the pi2
    pi.license()

    # List all available commands
    pi.help()

    # Show help for convert command
    pi.help('convert')


    



def math():
    """
    Demonstrates use of simple image math.
    """

    # Add images
    img = pi.read(input_file())
    pi.add(img, img)
    pi.writeraw(img, output_file('head_added_to_itself'))


    # Subtract images
    img = pi.read(input_file())
    pi.subtract(img, img)
    pi.writeraw(img, output_file('head_subtracted_from_itself'))


    # Add constant
    # The math operations are saturating, i.e. if the result of an operation is out of
    # bounds that can be represented with pixel data type, the value is clipped
    # to the bounds.
    # For example,
    # 200 + 200 = 255 for uint8 image,
    # 200 + 200 = 400 for uint16 image,
    # 2*30000 = 65535 for uint16 image,
    # etc.
    img = pi.read(input_file())
    pi.add(img, 65400) # Add large value to partially saturate 16-bit range
    pi.writeraw(img, output_file('head_saturated'))


    # Do you have a 2D mask that you would like to apply to
    # all slices of a 3D stack?
    # No problem, just specify True for 'allow dimension broadcast' parameter:
    img = pi.read(input_file())

    # Create mask whose size is the same than the size of the original but it contains
    # only one slice. Then draw a circle into it, with color 1.
    mask = pi.newimage(img.get_data_type(), img.get_width(), img.get_height())
    pi.sphere(mask, [img.get_width() / 2, img.get_height() / 2, 0], img.get_width() / 4, 1)
    pi.writetif(mask, output_file('mask'))
    
    pi.multiply(img, mask, True)
    pi.writetif(img, output_file('head_masked'))
    


def convert_format():
    """
    Demonstrates file format conversion.
    """

    # Read image
    img = pi.read(input_file('simple_structures.vol'))

    # Write in any supported format
    pi.writeraw(img, output_file('vol_to_raw'))
    pi.writetif(img, output_file('vol_to_tif'))


def filtering():
    """
    Calculates Gaussian blur of image normally and using distributed processing.
    Calculates difference of the two versions.
    The example generalizes to any available filtering procedure like
    vawefilter, bilateralfilter, maxfilter, openingfilter, etc.
    """

    # Gaussian filtering (local sequential 'distributed' processing)
    # --------------------------------------------------------------

    # Enable distributed mode
    pi.distribute(Distributor.LOCAL)

    # For demonstration, set memory per one job to low value.
    # 25 megabytes results in 2 jobs for the default input image in this example.
    # Typically you would set this value in local_config.txt file.
    pi.maxmemory(25)

    # Read image
    img = pi.read(input_file())

    # Create output image
    filtered = pi.newlike(img)

    # Filter
    pi.gaussfilter(img, filtered, 5)

    # Write output to disk.
    # The distributed mode saves internal temporary images as .raw files or .png
    # sequences. Writeraw command in distributed mode therefore often converts
    # into a simple file rename.
    pi.writeraw(filtered, output_file('head_gauss_distributed'))

    # Disable distributed mode
    pi.maxmemory(0) # Sets max memory to automatically determined value
    pi.distribute(Distributor.NONE)


    # Gaussian filtering (local processing)
    # -------------------------------------

    # This code is the same than in distributed case above, but without
    # pi.distribute-commands.
    img = pi.read(input_file())
    filtered = pi.newlike(img)
    pi.gaussfilter(img, filtered, 5)
    pi.writeraw(filtered, output_file('head_gauss_normal'))


    # Calculate difference of results of normal and distributed processing
    # --------------------------------------------------------------------

    # Read both images
    img = pi.read(output_file('head_gauss_normal'))
    img2 = pi.read(output_file('head_gauss_distributed'))

    # Convert them to float32 so that negative values can be represented, too.
    pi.convert(img, ImageDataType.FLOAT32)
    pi.convert(img2, ImageDataType.FLOAT32)

    # Subtract img2 from img
    pi.subtract(img, img2)
    pi.writeraw(img, output_file('gauss_difference'))

    # Calculate absolute value of each pixel
    pi.abs(img)

    # Calculate maximal value and place it to image M.
    # M will be a 1x1x1 image.
    M = pi.newimage(ImageDataType.FLOAT32)
    pi.maxval(img, M)

    # Get the value of the first pixel of image M.
    # In this case M is a 1x1x1 image so we have only one pixel anyway.
    M = M.get_value()
    print(f"Maximal difference = {M}")



def thickmap():
    """
    Demonstrates use of thickness map functions.
    """

    # Read input image. Input must be binary (or it must be made binary,
    # see e.g. 'threshold' command)
    geom = pi.read(input_file_bin())

    # Create output image
    tmap = pi.newimage(ImageDataType.FLOAT32)

    # Convert the input to large enough data type so that it can hold squared distance values.
    # If the data type is too small, an error is raised.
    pi.convert(geom, ImageDataType.UINT32)
    pi.tmap(geom, tmap)
    pi.writeraw(tmap, output_file('head_tmap'))

    # If the default version consumes too much memory, there is also a slower memory-saver version.
    # It is activated by setting last parameter of tmap to True:
    geom2 = pi.read(input_file_bin())
    tmap2 = pi.newimage(ImageDataType.FLOAT32)
    pi.convert(geom2, ImageDataType.UINT32)
    pi.tmap(geom2, tmap2, 0, True)
    pi.writeraw(tmap, output_file('head_tmap_memsave'))



def watershed():
    """
    Demonstrates Meyer's flooding algorithm for calculation of watershed.
    """

    # Read image
    weights = pi.read(input_file('t1-head.tif'))

    # Create new image, taking unspecified properties from the old image
    labels = pi.newlike(weights, 'uint8')

    # Set some seeds
    pi.set(labels, [110, 90, 63], 100) # Brain
    pi.set(labels, [182, 165, 63], 200) # Skull

    # Save seeds so that they can be viewed later
    pi.writetif(labels, output_file('meyer_grow_seeds'))

    # Grow the seeds. (Normally you would use gradient of input image etc. as weight)
    pi.grow(labels, weights)

    # Save result
    pi.writetif(labels, output_file('meyer_grow'))




def skeleton_vtk():
    """
    Creates skeleton, traces it to a graph/network structure, and saves
    the skeleton branches in a .vtk file.
    """

    # Generate a simple image for demonstration
    geometry = pi.newimage(ImageDataType.UINT8, 100, 100)

    p0 = [1, 1, 0]
    p1 = [25, 50, 0]
    p2 = [75, 50, 0]
    p3 = [5, 95, 0]
    p4 = [95, 95, 0]
    
    pi.line(geometry, p0, p1, 255)
    pi.line(geometry, p1, p2, 255)
    pi.line(geometry, p0, p2, 255)
    pi.line(geometry, p1, p3, 255)
    pi.line(geometry, p2, p4, 255)

    # Save the geometry
    pi.writeraw(geometry, output_file("geometry"));

    # Convert geometry to line skeleton
    pi.lineskeleton(geometry);

    # Write the skeleton so that it can be visualized
    pi.writeraw(geometry, output_file("skeleton"));

    # Create images for tracing the skeleton into a graph
    vertices = pi.newimage(ImageDataType.FLOAT32)
    edges = pi.newimage(ImageDataType.UINT64)
    measurements = pi.newimage(ImageDataType.FLOAT32)
    points = pi.newimage(ImageDataType.INT32)

    # Trace the skeleton
    # The last 1 gives count of threads. For images with small number of branches
	# (like here) single-threaded processing is faster.
    pi.tracelineskeleton(geometry, vertices, edges, measurements, points, True, 1)

    # Convert the traced skeleton into points and lines format, as the .vtk files
    # need that format.
    vtkpoints = pi.newimage()
    vtklines = pi.newimage()
    pi.getpointsandlines(vertices, edges, measurements, points, vtkpoints, vtklines);

    pi.writevtk(vtkpoints, vtklines, output_file("vtk_test_file"))



def big_endian_and_little_endian():
    """
    Demonstrates how to change the endianness of pixel values.
    """

    # The images are in the native byte order
    geometry = pi.newimage(ImageDataType.UINT16, 100, 100)
    pi.ramp(geometry, 1)

    # By default, images are written in the native byte order.
    # This usually means little endian at the time this was written, and at the
    # processors this code has ever run.
    pi.writeraw(geometry, output_file("little_endian_ramp"))

    # Change to big endian format can be made by running swapbyteorder command
    pi.swapbyteorder(geometry)
    pi.writeraw(geometry, output_file("big_endian_ramp"))

    # NOTE: Here the geometry image is stored in big endian format, and you cannot do
    # any calculations on it until you run the swapbyteorder command on it again!
    # For example, writing something else than .raw files for image in non-native byte
    # order does not make sense! The .tif file saved below will
    # contain pixel values in wrong byte order.
    pi.writetif(geometry, output_file("big_endian_ramp.tif"))

    



def greedy_coloring():
    """
    Shows how to change colors of regions.
    """

    # Create image
    img = pi.newimage(ImageDataType.UINT8, 200, 200, 1)

    # Plot some overlapping spheres.
    # Each of the spheres has different color
    pi.sphere(img, [100, 100, 0], 50, 7)
    pi.sphere(img, [150, 100, 0], 50, 10)
    pi.sphere(img, [50, 100, 0], 50, 15)
    pi.sphere(img, [100, 50, 0], 50, 20)
    pi.sphere(img, [150, 150, 0], 50, 25)

    # Save the initial image
    pi.writeraw(img, output_file("before_coloring"))

    # Minimize number of colors in the image, still making
    # sure that all neighbouring spheres have different color.
    pi.greedycoloring(img)

    # Save the result
    pi.writeraw(img, output_file("after_coloring"))



def linefilters():
    """
    Calculates Frangi line-enhancing filter.
    """

    # Read image
    img = pi.read(input_file())

    # Frangi filter, 16-bit result
    frangi_16bit = pi.newlike(img, ImageDataType.UINT16)
    pi.frangifilter(img, frangi_16bit, 3.0)
    pi.writeraw(frangi_16bit, output_file('frangi_16bit'))


    # Frangi filter, 32-bit result. Notice different scaling in the output, compared to the 16-bit version.
    frangi_32bit = pi.newlike(img, ImageDataType.FLOAT32)
    pi.convert(img, ImageDataType.FLOAT32)
    pi.frangifilter(img, frangi_32bit, 3.0)
    pi.writeraw(frangi_32bit, output_file('frangi_32bit'))




def binning_scaling():
    """
    Demonstrates binning.
    """

    # Read image
    img = pi.read(input_file())

    # Normal binning
    binned_mean = pi.newlike(img)
    pi.bin(img, binned_mean, 4)
    pi.writeraw(binned_mean, output_file('binning_mean'))

    # Min binning
    binned_min = pi.newlike(img)
    pi.bin(img, binned_min, 4, 'min')
    pi.writeraw(binned_min, output_file('binning_min'))

    # Max binning
    binned_max = pi.newlike(img)
    pi.bin(img, binned_max, 4, 'max')
    pi.writeraw(binned_max, output_file('binning_max'))



def rotations():
    """
    Demonstrates rotations and re-slicing.
    """

    # Read image
    img = pi.read(input_file())

    # Rotate 90 degrees clockwise (around z-axis)
    rot90 = pi.newimage()
    pi.rot90cw(img, rot90)
    pi.writetif(rot90, output_file('rotate_90_clockwise'))

    # Reslice (rotate 90 degrees around x- or y-axis)
    top = pi.newimage()
    pi.reslice(img, top, ResliceDirection.TOP)
    pi.writetif(top, output_file('reslice_top'))

    # General rotation.
    # NOTE:
    # - The size of the output must be set to desired value before the call
    #   to rotate function.
    # - The angle is given in radians. Here it is -60 degrees.
    # - The axis is given as a vector and it does not need to be a unit vector.
    grot = pi.newimage(img.get_data_type(), img.get_dimensions())
    pi.rotate(img, grot, -60/180*3.14, [1, 1, 0])
    pi.writetif(grot, output_file('general_rotation'))






def generate_particles(sphere_count=1000, box_count=1000):
    """
    analyze_particles demo uses this to generate an image containing some random regions.
    """

    import random

    img = pi.newimage(ImageDataType.UINT8, 500, 500, 500)

    # Generate some spheres
    for i in range(0, sphere_count):
        pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        r = random.randint(1, 20)
        pi.sphere(img, pos, r, 255)

    # Generate some boxes
    for i in range(0, box_count):
        pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        size = [random.randint(1, 20), random.randint(1, 20), random.randint(1, 20)]
        pi.box(img, pos, size, 255)

    pi.writeraw(img, output_file('particles'))




def analyze_particles():
    """
    Demonstrates particle (i.e, region or blob) analysis, and some drawing commands.
    """

    # Generate particle image
    generate_particles()

    # Show analyzer names.
    # This does not do anything else than shows which analyzers are available.
    pi.listanalyzers()

    # Make a list of analyzers that we want to use
    analyzers = 'volume coordinates bounds boundingsphere isonedge'


    # Read generated data file
    img = pi.read(output_file('particles'))

    # Analyze particles
    result = pi.newimage(ImageDataType.FLOAT32)
    pi.analyzeparticles(img, result, analyzers)

    # Show titles of data columns
    print('Titles of columns in data table:')
    pi.headers(analyzers)

    # Get result data from the pi2 system
    pa = result.get_data()

    # Volume is in the first column (see output of 'headers' command)
    volume = pa[:, 0]

    # Plot volume histogram
    import matplotlib.pyplot as plt
    plt.hist(volume, density=True, bins=30)
    plt.xlabel('Volume [pixel]')
    plt.ylabel('Probability')
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig(output_file('volume_distribution.png'))





def fill_particles():
    """
    Demonstrates analysis and filling of particles.
    """

    # Generate particle image
    generate_particles()

    
    # Read generated data file
    img = pi.read(output_file('particles'))

    # Analyze particles
    analyzers = 'volume coordinates pca boundingsphere'
    result = pi.newimage(ImageDataType.FLOAT32)
    pi.analyzeparticles(img, result, analyzers)

    # Get result data from the pi2 system
    pa = result.get_data()

    # Volume is in the first column (see the output of the 'headers' command)
    volume = pa[:, 0]

    print(f"Before filtering the shape of the results matrix is {pa.shape}")

    # Filter the measurement results using suitable Python methods.
    # Here we take only the largest particles
    limit = 20000
    pa = pa[volume > limit, :]

    print(f"After filtering the shape of the results matrix is {pa.shape}")

    # We could use set_data(...) function to push the filtered results
    # back to pi2, but we may as well use the NumPy array pa directly
    # in pi2 commands.
    #result.set_data(pa)

    # Fill the particles that we left into the pa array (i.e. the big ones)
    # Note that (for simple distributed processing support) the fillparticles
    # sets all non-filled particles to 1.
    img = pi.read(output_file('particles'))
    pi.fillparticles(img, analyzers, pa, 2)
    # Set 1 -> 255 and 2 -> 128 so that the colors of the filled particles
    # correspond to the original colors.
    pi.replace(img, 1, 255)
    pi.replace(img, 2, 128)
    pi.writeraw(img, output_file('particles_large_colored'))



    # Make another visualization by drawing an ellipsoid approximation of each particle
    # on top of the particles
    ellipsoid_vis = pi.read(output_file('particles'))
    
    # This draws the particles as ellipsoids whose volume equals particle volume, with color 128
    pi.drawellipsoids(ellipsoid_vis, analyzers, result, 128, EllipsoidType.VOLUME)
    
    # Save the result
    pi.writeraw(ellipsoid_vis, output_file('particles_ellipsoids'))




def histogram():
    """
    Demonstrates histogram calculation.
    """

    # Read image
    img = pi.read(input_file())

    # Calculate histogram
    hist_min = 0        # Start of the first bin
    hist_max = 1000     # End of the last bin
    bin_count = 100     # Count of bins
    hist = pi.newimage(ImageDataType.FLOAT32)
    bins = pi.newimage(ImageDataType.FLOAT32)
    pi.hist(img, hist, bins, hist_min, hist_max, bin_count)
    

    # Get histogram data as NumPy array
    data = hist.get_data()
    bin_starts = bins.get_data()

    # Most of the background pixels are in the first bin.
    # We set it to zero so that the output image becomes nicely scaled
    # automatically.
    data[0] = 0

    # Bin edges can be generated like this if they have not been grabbed
    # from hist(...) command
    #bin_edges = np.linspace(hist_min, hist_max, bin_count + 1)

    # Convert bin starts to bin centers
    bin_size = bin_starts[1] - bin_starts[0]
    bin_centers = bin_starts + bin_size / 2

    
    # Plot the histogram
    import matplotlib.pyplot as plt
    plt.bar(bin_centers, data, 0.8*bin_size)
    plt.xlabel('Gray value')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig(output_file('histogram.png'))





def bivariate_histogram():
    """
    Demonstrates determination of bivariate histogram.
    """

    # Read original
    img = pi.read(input_file())
    
    # Create a copy of the image and binarize it
    bin = pi.newlike(img)
    pi.copy(img, bin)
    pi.threshold(bin, 185)

    # Calculate thickness map.
    # Convert input image first to uint32 so that it can hold even large values.
    pi.convert(bin, ImageDataType.UINT32);
    tmap = pi.newlike(bin, ImageDataType.FLOAT32)
    pi.tmap(bin, tmap)

    # Calculate bivariate histogram of original and thickness map
    hist = pi.newimage(ImageDataType.FLOAT32)
    bins1 = pi.newimage(ImageDataType.FLOAT32)
    bins2 = pi.newimage(ImageDataType.FLOAT32)
    gray_max = 500   # Maximum gray value in the histogram bins
    thick_max = 30   # Maximum thickness in the histogram bins
    pi.hist2(img, 0, gray_max, 60, tmap, 0, thick_max, 15, hist, bins1, bins2) # 60 gray value bins, 15 thickness bins

    # Write output files to disk
    pi.writetif(img, output_file('original'))
    pi.writetif(tmap, output_file('thickness_map'))
    pi.writetif(hist, output_file('histogram'))
    pi.writetif(bins1, output_file('bins1'))
    pi.writetif(bins2, output_file('bins2'))

    # Get histogram data
    h = hist.get_data()

    # Make plot
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(4.5, 6))
    
    # Note that here we plot the transpose of the histogram as that fits better into an image.
    # We also flip it upside down so that both x- and y-values increase in the usual directions.
    pltimg = plt.imshow(np.flipud(h.transpose()), vmin=0, vmax=1.5e4, extent=(0, thick_max, 0, gray_max), aspect=1/8)
    cbar = fig.colorbar(pltimg)
    plt.xlabel('Thickness [pix]')
    plt.ylabel('Gray value [pix]')
    cbar.set_label('Count', rotation=90)
    plt.tight_layout()
    plt.show(block=False)
    plt.savefig(output_file('bivariate_histogram.png'))

    

def particle_segmentation():
    """
    Demonstrates segmentation of particles using watershed seeded
    by (filtered) local maxima of distance map.
    """

    # Generate particle image
    generate_particles(1000, 0)

    
    # Read generated data file
    img = pi.read(output_file('particles'))

    # Convert to wider data type so that we can label each particle
    pi.convert(img, ImageDataType.UINT16)

    # Calculate distance map
    dmap = pi.newimage(ImageDataType.FLOAT32)
    pi.dmap(img, dmap)

    # Find maxima
    # Note that in some cases the system is able to automatically
    # change the data type of output images, so we don't have to
    # specify any data type in the pi.newimage() command.
    # This does not work always, though, as there might be many
    # possible output data types.
    maxima = pi.newimage()
    pi.localmaxima(dmap, maxima)

    # Remove unnecessary maxima to avoid over-segmentation
    pi.cleanmaxima(dmap, maxima)

    # Create image with labeled maxima only.
	# First set all pixels to zero, then label maxima.
    pi.set(img, 0)
    pi.labelmaxima(img, maxima)

    # Grow labels back to original geometry, using distance map value
    # as filling priority.
    pi.grow(img, dmap)

    # Save result
    pi.writeraw(img, output_file('particles_watershed'))




def levelset_fill_cavity():
    """
    Demonstrates how to fill a non-convex cavity using a level-set method.
    """

    # Create an image that contains a rectangle with a cavity in its edge.
    # For this example we make a 2D image for easier visualization,
    # but everything should work the same for a 3D volume image, too.
    c = 30 # Radius of the image
    geom = pi.newimage(ImageDataType.UINT8, 2 * c, 2 * c)
    pi.box(geom, [c - 15, c - 15, 0], 30, 128)      # The box
    pi.sphere(geom, [c, c, 0], 8, 0)                # Part of cavity
    pi.sphere(geom, [c, c - 6, 0], 6, 0)            # Part of cavity
    pi.sphere(geom, [c, c - 12, 0], 6, 0)           # Part of cavity

    # Noise to make the image more realistic
    pi.noise(geom, 0, 20)
    
    # Save the geometry
    pi.writetif(geom, output_file('levelset_geom'));


    # Initialize an image that holds the level set function.
    # After iteration below, the segmentation is the part of phi
    # that has positive value.
    phi = pi.newlike(geom, ImageDataType.FLOAT32)

    # Update phi iteratively
    for n in np.arange(0, 200):
        print(f"Iteration {n}" )

        # Construct force field that acts on the surface defined by phi
        # The force will be the sum of three terms.
        F = pi.newlike(phi)



        # First term: the force is positive inside the object and negative everywhere else.
		# This makes the surface take the shape of the object.
        pi.copy(geom, F)
        pi.divide(F, 128)
        pi.subtract(F, 0.5)



        # Second term: Penalize high curvature by making curvy
        # regions less curved.
        kappa = pi.newlike(phi)
        pi.meancurvature(phi, kappa, 0.5)

        # Multiply the curvature values to scale them correctly.
        pi.multiply(kappa, -5)

        # Remove negative curvature values. They correspond to
        # convex shapes, and we want to zero those so that
        # they don't have any effect on the surface.
        pi.max(kappa, 0)

        # Add kappa to the force term
        pi.add(F, kappa)



        # Third term: Normal force
        # This term makes the surface move towards its normal.
        # This term is not required in this example.
        #L = pi.newlike(phi)
        #pi.gradientmagnitude(phi, L, 0.75, 0)
        #pi.multiply(L, 0.01)
        #pi.add(F, L)


        # Smooth the total force a little bit
        # This is not strictly by the book, but smoothing seems to
        # make phi converge faster to a smoother result.
        tmp = pi.newlike(F)
        pi.gaussfilter(F, tmp, 0.5, BoundaryCondition.NEAREST)
        pi.set(F, tmp)


        # Multiply by time step
        dt = 1
        pi.multiply(F, dt)

        # Add force*dt to the surface
        pi.add(phi, F)


        # Convert phi to segmentation and save it for visualization purposes
        vis = pi.newlike(phi)
        pi.copy(phi, vis)
        pi.threshold(vis, 0)
        pi.convert(vis, ImageDataType.UINT8)
        pi.writetif(vis, output_file(f"./levelset_result/{n}"))




def find_surface():
    """
    Demonstrates surface recognition using Edwards-Wilkinson model.
    """

    # Read image
    orig = pi.read(input_file())

    # Crop smaller piece so that the surface does not 'fall through'
    # as there's nothing in the edges of the image.
    # Usually surfaces are recognized from e.g. images of paper sheet
    # where the sheet spans the whole image.
    img = pi.newimage(orig.get_data_type())
    pi.crop(orig, img, [80, 90, 0], [85, 120, orig.get_depth()])

    pi.writeraw(img, output_file('findsurface_geometry'))

    # This image will hold the height map of the surface
    hmap = pi.newimage(ImageDataType.FLOAT32)

    # ... and this one will be a visualization of surface evolution
    vis = pi.newimage(orig.get_data_type())

    # Now find the surface.
    # We will stop at gray level 100 and perform 60 iterations
    # with surface tension 1.0.
    pi.findsurface(img, hmap, 100, Direction.DOWN, 1.0, 60, vis, img.get_height() / 2, 900)

    # Save the result surface
    pi.writeraw(hmap, output_file('findsurface_height_map'))

    # Save the visualization
    pi.writeraw(vis, output_file('findsurface_evolution'))



    # We can also draw the surface to the image
    surf_vis = pi.newimage()
    pi.set(surf_vis, img)
    pi.drawheightmap(surf_vis, hmap, 900)
    pi.writeraw(surf_vis, output_file('findsurface_full_vis'))

    # Or set pixels below or above the surface.
    # This is useful to, e.g., get rid of background noise above
    # or below the surface of the sample.
    before_vis = pi.newimage()
    pi.set(before_vis, img)
    pi.setbeforeheightmap(before_vis, hmap, 900)
    pi.writeraw(before_vis, output_file('findsurface_before_vis'))

    after_vis = pi.newimage()
    pi.set(after_vis, img)
    pi.setafterheightmap(after_vis, hmap, 900)
    pi.writeraw(after_vis, output_file('findsurface_after_vis'))


    # Or the image can be shifted according to the surface
    pi.negate(hmap)
    temp = pi.newlike(hmap)
    pi.gaussfilter(hmap, temp, 1) # This smooths the surface map a little bit
    pi.set(hmap, temp)
    pi.shiftz(img, hmap, True)
    pi.writeraw(img, output_file('findsurface_shiftz'))





def sample_orientations_from_vonmises_fisher(main_az, main_pol, kappa, n):
    """
    Sample directions from von Mises-Fisher distribution.
    main_az and main_pol give the azimuthal and polar angles of the main direction.
    kappa indicates the spread of the directions around the main direction.
    Large kappa means small spread.
    n gives the number of directions to generate.
    Returns n 3-component unit vectors.

    This code mostly from https://stats.stackexchange.com/questions/156729/sampling-from-von-mises-fisher-distribution-in-python
    but its correctness has not been checked. It seems to
    create plausible results, though.
    """

    import scipy as sc
    import scipy.stats
    import scipy.linalg as la

    def sample_tangent_unit(mu):
        mat = np.matrix(mu)

        if mat.shape[1]>mat.shape[0]:
            mat = mat.T

        U,_,_ = la.svd(mat)
        nu = np.matrix(np.random.randn(mat.shape[0])).T
        x = np.dot(U[:,1:],nu[1:,:])
        return x/la.norm(x)


    def rW(n, kappa, m):
        dim = m-1
        b = dim / (np.sqrt(4*kappa*kappa + dim*dim) + 2*kappa)
        x = (1-b) / (1+b)
        c = kappa*x + dim*np.log(1-x*x)

        y = []
        for i in range(0,n):
            done = False
            while not done:
                z = sc.stats.beta.rvs(dim/2,dim/2)
                w = (1 - (1+b)*z) / (1 - (1-b)*z)
                u = sc.stats.uniform.rvs()
                if kappa*w + dim*np.log(1-x*w) - c >= np.log(u):
                    done = True
            y.append(w)

        return np.array(y)


    def rvMF(n,theta):
        dim = len(theta)
        kappa = np.linalg.norm(theta)
        mu = theta / kappa

        w = rW(n, kappa, dim)

        result = []
        for sample in range(0,n):
        
            v = sample_tangent_unit(mu).transpose()
            v = np.asarray(v.transpose()).squeeze()

            result.append(np.sqrt(1-w[sample]**2)*v + w[sample]*mu)

        return result




    # First sample directions from the von Mises-Fisher distribution
    main_dir = np.array([np.cos(main_az) * np.sin(main_pol), 1 * np.sin(main_az) * np.sin(main_pol), 1 * np.cos(main_pol)])
    directions = rvMF(n, kappa * main_dir)

    # Convert to orientations (v and -v are the same)
    # by ensuring that all directions have positive x coordinate.
    # This is the same convention used in pi2.
    for dir in directions:
        if dir[0] < 0:
            dir *= -1

    return directions



def generate_orientation_image(main_az=np.pi/4, main_pol=np.pi/2, kappa=3, n=150):
    """
    Generates test image for orientation analysis demonstrations.
    By default, plots 100 capsules whose orientations are taken from from von Mises-
    Fisher distribution (main direction = x + 45 deg towards y, kappa = 3)

    Returns the resulting pi2 image and direction vectors.
    """

    directions = sample_orientations_from_vonmises_fisher(main_az, main_pol, kappa, n)

    size = 300
    img = pi.newimage(ImageDataType.UINT8, [size, size, size])
    for dir in directions:

        L = 100
        r = 5
        pos = np.random.rand(1, 3) * size;

        pi.capsule(img, pos - L * dir, pos + L * dir, r, 255)

    return img, directions



def orientation_analysis():
    """
    Demonstrates how to determine and visualize orientation of structures.
    """

    # Create test image with given main fibre direction
    main_azimuthal = np.pi/4
    main_polar = np.pi/2
    img, true_orientations = generate_orientation_image(main_azimuthal, main_polar)
    
    # Save it for later visualization
    pi.writetif(img, output_file('cylinders'))

    # Calculate orientation of cylinders
    # Note that we need to convert the input image to float32
    # format as it is used as an output image, too.
    # Additionally we need images to store the azimuthal and polar
    # orientation angles.
    pi.convert(img, ImageDataType.FLOAT32)
    azimuthal = pi.newimage(ImageDataType.FLOAT32)
    polar = pi.newimage(ImageDataType.FLOAT32)
    pi.cylinderorientation(img, azimuthal, polar, 1, 1)

    # Now img has been replaced with 'orientation energy'
    pi.writetif(img, output_file('cylinders_energy'))

    # Make a color-coded visualization of the orientations
    r = pi.newimage()
    g = pi.newimage()
    b = pi.newimage()
    pi.mainorientationcolor(img, azimuthal, polar, main_azimuthal, main_polar, r, g, b)
    #pi.axelssoncolor(img, azimuthal, polar, r, g, b) # This is another possibility if main orientation is not available.
    pi.writeraw(r, g, b, output_file('cylinders_main_orientation_color'))


    # Make orientation histogram.
    # Energy is used as weight
    hist = pi.newimage(ImageDataType.FLOAT32)
    bins1 = pi.newimage(ImageDataType.FLOAT32)
    bins2 = pi.newimage(ImageDataType.FLOAT32)
    pi.whist2(azimuthal, -np.pi, np.pi, 20, polar, 0, np.pi, 10, img, hist, bins1, bins2) # 20 azimuthal angle bins, 10 polar angle bins

    


    # Make a plot that compares the true orientation distribution to the estimated one
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(4.5, 6))

    # First plot the true orientations
    plt.subplot(2, 1, 1)

    # Convert directions to polar coordinates using the same convention that pi2 uses
    azs = []
    pols = []
    for dir in true_orientations:
        x = dir[0]
        y = dir[1]
        z = dir[2]
        r = np.sqrt(x * x + y * y + z * z)
        azimuthal = np.arctan2(y, x)
        polar = np.arccos(z / r)
        azs.append(azimuthal)
        pols.append(polar)

    # Calculate orientation histogram using the NumPy method
    hst, xedges, yedges = np.histogram2d(azs, pols, range=[[-np.pi, np.pi], [0, np.pi]], bins=[20, 10])

    # Plot the histogram
    pltimg = plt.imshow(hst.transpose(), extent=(xedges[0], xedges[-1], yedges[0], yedges[-1]))
    cbar = fig.colorbar(pltimg)
    plt.xlabel('Azimuthal angle [rad]')
    plt.ylabel('Polar angle [rad]')
    plt.title('True distribution of cylinder orientations')


    # Now plot the histogram estimated from the image
    plt.subplot(2, 1, 2)
    pltimg = plt.imshow(hist.get_data(), extent=(-np.pi, np.pi, 0, np.pi))
    cbar = fig.colorbar(pltimg)
    plt.xlabel('Azimuthal angle [rad]')
    plt.ylabel('Polar angle [rad]')
    plt.title('Distribution estimated from the image')

    # Show and print the figure
    plt.tight_layout()
    plt.show(block=False)

    plt.savefig(output_file('bivariate_histogram_comparison.png'))




def montage():
    """
    Demonstrates creation of a montage from a 3D image.
    """

    # Read the 3D image
    img = pi.read(input_file())

    # Create empty image for the montage
    montage = pi.newimage(img.get_data_type())

    # Make the montage, 4x3 slices, scaled to 0.5 of
    # original size
    pi.montage(img, montage, 4, 3, 0.5)

    # Save the output
    pi.writetif(montage, output_file('montage'))



def skeleton_types():
    """
    Demonstrates differences between skeletonization algorithms available in pi2.
    """

    # Create a 3D image that contains a thick T letter
    img = pi.newimage(ImageDataType.UINT8, 100, 100, 100)
    pi.box(img, [20, 20, 40], [61, 21, 26], 255)
    pi.box(img, [40, 20, 40], [21, 71, 26], 255)

    # Add a cavity (hole) to the T
    pi.sphere(img, [50, 50, 50], 5, 0)

    # Save the input geometry
    pi.writeraw(img, output_file('T'))

    # Create a skeleton that contains 1-pixel thick planes and 1-pixel thick lines.
    surface_skele_true = pi.newlike(img)
    pi.set(surface_skele_true, img)
    pi.surfaceskeleton(surface_skele_true, True)
    pi.writeraw(surface_skele_true, output_file('T_surface_skele_true'))

    # Create a skeleton that contains 1-pixel thick planes only around cavities
    # and 1-pixel thick lines elsewhere.
    surface_skele_false = pi.newlike(img)
    pi.set(surface_skele_false, img)
    pi.surfaceskeleton(surface_skele_false, False)
    pi.writeraw(surface_skele_false, output_file('T_surface_skele_false'))

    # Create a skeleton that contains only 1-pixel thick lines.
    line_skele = pi.newlike(img)
    pi.set(line_skele, img)
    pi.lineskeleton(line_skele)
    pi.writeraw(line_skele, output_file('T_line_skele'))




#create_and_access_images()
#read_and_write_image()
#help()
#math()
#convert_format()
#filtering()
#thickmap()
#watershed()
#skeleton_vtk()
#big_endian_and_little_endian()
#greedy_coloring()
#linefilters()
#binning_scaling()
#rotations()
#analyze_particles()
#fill_particles()
#histogram()
#bivariate_histogram()
#particle_segmentation()
#levelset_fill_cavity()
#find_surface()
#orientation_analysis()
#montage()
#skeleton_types()

