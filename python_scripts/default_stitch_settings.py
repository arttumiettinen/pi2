

def write_stitch_settings(sample_name, binning, positions, point_spacing=60, coarse_block_radius=120, coarse_binning=4, cluster_name='None', normalize_while_stitching=False, max_circle_diameter=-1):
    """
    Writes new stitch settings file to stitch_settings.txt.
    sample_name - prefix of output files.
    binning - desired value for binning.
    positions - Lines that will be written under heading [positions]
    """

    settings_text = f"""[stitch]

# Output name
sample_name = {sample_name}

# Output binning
binning = {binning}

# Space between points in the reference grid in the pairwise matching step.
# Specify such a value that the deformation is approximately linear between the grid points.
# The value should scale with image size.
point_spacing = {point_spacing}


# Size of image block in pairwise matching (coarse step).
# The block does not have to be cube, you can use syntax [a, b, c] to specify different dimension
# for each coordinate direction.
# The block radius is the maximal local shift that can be found between any two sub-images.
# The value should scale with image size.
coarse_block_radius = {coarse_block_radius}


# Amount of binning for coarse matching
# If set to 1, no fine matching is performed.
# Set to value greater than one to make a coarse match first and then a fine match using settings below.
# The value should scale with image size.
coarse_binning = {coarse_binning}


# Size of image block in pairwise matching (fine step).
# Indicates maximum shift on top of shift determined in coarse matching.
# The value should scale with image size.
fine_block_radius = 1


# Amount of binning in fine matching.
# If greater than one, the full resolution of the image is not used to find the shifts.
# Sometimes the results seem to be very good although this value is greater than one.
# If set to the same value than coarse_binning, the fine matching step is skipped.
# The value should scale with image size.
fine_binning = {coarse_binning}


# Set to true to try to optimize locations and orientations of the sub-images so that
# any discrepancies in overlapping regions are as small as possible before correcting
# them with non-rigid displacements.
global_optimization = True


# Allow rigid body rotation?
# Disallow if the rotations have not been determined accurately enough, or if
# errors in rotations propagate through the sample in unwanted way (e.g. in the case
# of long Nx1x1 mosaic)
allow_rotation = False


# Allow local deformations?
# Disallow to make a rigid stitch instead of non-rigid one. Use e.g. for testing
# the effect of the non-rigid transformation on the output.
# Note that the local displacement fields are calculated anyway for determination
# of optimal positions for the sub-images, even if local deformations are not
# allowed in the final mosaic.
allow_local_deformations = True


# Normalization of gray values in pairwise matching phase.
# Set to true to make means of the images in the overlapping regions the same before pairwise
# matching. The value does not have any effect on the gray levels of the final mosaic.
normalize_in_blockmatch = True


# Normalization of gray values in the final mosaic.
# Set to true to make mean of the images in the overlapping regions the same before
# assembly of the final mosaic. Enabling or disabling this option may cause global
# or local gray-value gradients, respectively.
normalize_while_stitching = {normalize_while_stitching}


# Set to True to treat pixels that have value 0 in the tiles as missing values.
# Affects the final mosaic building phase only, zeroes are still treated as zeroes in the matching phase.
zeroes_are_missing_values = True


# Selects tile pixel weighing method. Set to a negative value to use rect weight, where weight
# of a pixel is proportional to the distance from the edge of the tile. Set to 0 to use
# maximum circle weight, where the weight of a pixel is proportional to the distance from
# the edge of a maximal inscribed circle that fits inside the image in the xy-plane. Set
# to a positive value to use circle weight with user-specified diameter.
# This is useful if, e.g. tomographic slices contain bad values outside the
# well-reconstructed region.
# Note that currently masking is made after calculation of gray level normalization, see
# also normalize_while_stitching and normalize_in_blockmatch settings.
max_circle_diameter = {max_circle_diameter}


# Threshold value for displacement filtering.
# This is the T value in the filtering process that is done before determination of
# similarity transformations and before mosaic assembly.
# The value should scale with image size.
displacement_filter_threshold = 12


# Indicates which cluster manager should be used to submit jobs
cluster = {cluster_name}


# Maximum size of image block that is processed in one process is max_block_size^3.
# Set to zero to determine the value automatically based on available RAM.
# If create_goodness is false, set to such a value that
# (tile_pixel_size + 8) * max_block_size^3 < (available_memory_in_bytes).
# If create_goodness is true, set to such a value that 
# (tile_pixel_size + 12) * max_block_size^3 < (available_memory_in_bytes).
max_block_size = 0

# Defines the output format 
# candidates are: zarr, raw
output_format = raw


[positions]
# Enter image file name and approximate position of that image here.
# The positions should be given in the coordinate system of the first image, but the origin
# can be anything.
# You can skip the dimensions and extension as long as each file is uniquely identified.
{positions}
"""

    #settings_text = settings_text.format(binning=binning, sample_name=sample_name, positions=positions, coarse_binning=coarse_binning, coarse_block_radius=coarse_block_radius, point_spacing=point_spacing)

    with open("stitch_settings.txt", "w") as file:
        file.write(settings_text)
    


