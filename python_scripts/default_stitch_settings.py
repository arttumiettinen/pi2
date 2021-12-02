

def write_stitch_settings(sample_name, binning, positions, point_spacing=60, coarse_block_radius=120, coarse_binning=4, cluster_name='None', normalize_while_stitching=False):
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


# Threshold value for displacement filtering.
# This is the T value in the filtering process that is done before determination of
# similarity transformations and before mosaic assembly.
# The value should scale with image size.
displacement_filter_threshold = 12


# Set to true to re-do pairwise matching.
# If output files exist, the pairwise matching process is skipped for those images.
# Use this if parameters of pairwise matching process are changed and you don't want to
# manually remove the relevant output files.
redo = False


# Indicates which cluster manager should be used to submit jobs
cluster = {cluster_name}


# Maximum size of image block that is processed in one process is max_block_size^3.
# Set to such a value that (2 * pixel_size_in_bytes + 4) * max_block_size^3 < available memory in bytes.
max_block_size = 2500




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
    


