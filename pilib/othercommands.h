#pragma once

#include "command.h"
#include "commandsbase.h"
#include "overlapdistributable.h"
#include "parseexception.h"
#include "filtercommands.h"
#include "pointprocesscommands.h"
#include "structure.h"
#include "pilibutilities.h"
#include "io/io.h"
#include "registration.h"
#include "stitching.h"
#include "noise.h"
#include "misc.h"
#include "particleanalysis.h"
#include "histogram.h"

#include <vector>
#include <string>

using namespace std;
using namespace itl2;

namespace pilib
{


	template<typename pixel_t> class BlockMatchCommand : public Command
	{
	public:
		BlockMatchCommand() : Command("blockmatch", "Calculates displacement field between two images.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "reference image", "Reference image (non-moving image)."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "deformed image", "Deformed image (image to register to non-moving image)."),
				CommandArgument<coord_t>(ParameterDirection::In, "xmin", "X-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "xmax", "X-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "xstep", "Step between calculation points in x-direction."),
				CommandArgument<coord_t>(ParameterDirection::In, "ymin", "Y-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "ymax", "Y-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "ystep", "Step between calculation points in y-direction."),
				CommandArgument<coord_t>(ParameterDirection::In, "zmin", "Z-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "zmax", "Z-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "zstep", "Step between calculation points in z-direction."),
				CommandArgument<Vec3d>(ParameterDirection::In, "initial shift", "Initial shift between the images."),
				CommandArgument<string>(ParameterDirection::In, "file name prefix", "Prefix (and path) of files to write. The command will save point grid in the reference image, corresponding points in the deformed image, and goodness-of-fit. If the files exists, the current contents are erased."),
				CommandArgument<Vec3c>(ParameterDirection::In, "comparison radius", "Radius of comparison region.", Vec3c(25, 25, 25))
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& ref = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& def = *pop<Image<pixel_t>* >(args);
			coord_t xmin = pop<coord_t>(args);
			coord_t xmax = pop<coord_t>(args);
			coord_t xstep = pop<coord_t>(args);
			coord_t ymin = pop<coord_t>(args);
			coord_t ymax = pop<coord_t>(args);
			coord_t ystep = pop<coord_t>(args);
			coord_t zmin = pop<coord_t>(args);
			coord_t zmax = pop<coord_t>(args);
			coord_t zstep = pop<coord_t>(args);
			Vec3d initialShift = pop<Vec3d>(args);
			string fname = pop<string>(args);
			Vec3c compRadius = pop<Vec3c>(args);

			PointGrid3D<coord_t> refPoints(PointGrid1D<coord_t>(xmin, xmax, xstep), PointGrid1D<coord_t>(ymin, ymax, ystep), PointGrid1D<coord_t>(zmin, zmax, zstep));
			Image<Vec3d> defPoints(refPoints.pointCounts());
			Image<float32_t> fitGoodness(defPoints.dimensions());

			// Construct initial guess of the deformed points
			for (coord_t zi = 0; zi < defPoints.depth(); zi++)
			{
				for (coord_t yi = 0; yi < defPoints.height(); yi++)
				{
					for (coord_t xi = 0; xi < defPoints.width(); xi++)
					{
						defPoints(xi, yi, zi) = Vec3d(refPoints(xi, yi, zi)) + initialShift;
					}
				}
			}

			blockMatch(ref, def, refPoints, defPoints, fitGoodness, compRadius);

			//filterDisplacements(refPoints, defPoints, fitGoodness);
			writeBlockMatchResult(fname, refPoints, defPoints, fitGoodness);
		}
	};



	class BlockMatchPartialLoadCommand : public Command
	{
	public:
		BlockMatchPartialLoadCommand() : Command("blockmatchmemsave", "Calculates displacement between two images, loads only overlapping region from disk.",
			{
				CommandArgument<string>(ParameterDirection::In, "reference image file", "Name of reference image file (non-moving image)."),
				CommandArgument<string>(ParameterDirection::In, "deformed image file", "Name of deformed image file (image to register to non-moving image)."),
				CommandArgument<coord_t>(ParameterDirection::In, "xmin", "X-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "xmax", "X-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "xstep", "Step between calculation points in x-direction."),
				CommandArgument<coord_t>(ParameterDirection::In, "ymin", "Y-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "ymax", "Y-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "ystep", "Step between calculation points in y-direction."),
				CommandArgument<coord_t>(ParameterDirection::In, "zmin", "Z-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "zmax", "Z-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(ParameterDirection::In, "zstep", "Step between calculation points in z-direction."),
				CommandArgument<Vec3d>(ParameterDirection::In, "initial shift", "Initial shift between the images."),
				CommandArgument<string>(ParameterDirection::In, "file name prefix", "Prefix (and path) of files to write. The command will save point grid in the reference image, corresponding points in the deformed image, and goodness-of-fit. If the files exists, the current contents are erased."),
				CommandArgument<Vec3c>(ParameterDirection::In, "coarse comparison radius", "Radius of comparison region for coarse matching.", Vec3c(25, 25, 25)),
				CommandArgument<coord_t>(ParameterDirection::In, "coarse binning", "Amount of resolution reduction in coarse matching phase.", 2),
				CommandArgument<Vec3c>(ParameterDirection::In, "fine comparison radius", "Radius of comparison region for fine (full-resolution) matching.", Vec3c(10, 10, 10)),
				CommandArgument<coord_t>(ParameterDirection::In, "fine binning", "Amount of resolution reduction in fine matching phase. Set to same value than coarse binning to skip fine matching phase.", 2),
				CommandArgument<bool>(ParameterDirection::In, "normalize", "Indicates if the mean gray values of the two images should be made same in the overlapping region before matching.", true)
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			string refFile = pop<string>(args);
			string defFile = pop<string>(args);
			coord_t xmin = pop<coord_t>(args);
			coord_t xmax = pop<coord_t>(args);
			coord_t xstep = pop<coord_t>(args);
			coord_t ymin = pop<coord_t>(args);
			coord_t ymax = pop<coord_t>(args);
			coord_t ystep = pop<coord_t>(args);
			coord_t zmin = pop<coord_t>(args);
			coord_t zmax = pop<coord_t>(args);
			coord_t zstep = pop<coord_t>(args);
			Vec3d initialShift = pop<Vec3d>(args);
			string fname = pop<string>(args);
			Vec3c coarseCompRadius = pop<Vec3c>(args);
			coord_t coarseBinning = pop<coord_t>(args);
			Vec3c fineCompRadius = pop<Vec3c>(args);
			coord_t fineBinning = pop<coord_t>(args);
			bool normalize = pop<bool>(args);

			Vec3c refDimensions;
			ImageDataType refDT;
			if (!io::getInfo(refFile, refDimensions, refDT))
				throw ITLException("Unable to find dimensions and data type of reference image file.");

			//raw::expandRawFilename(refFile);
			//if (!raw::internals::parseDimensions(refFile, refDimensions, refDT))
			//	throw ITLException("Unable to parse dimensions of reference image from file name.");

			Vec3c defDimensions;
			ImageDataType defDT;
			if (!io::getInfo(defFile, defDimensions, defDT))
				throw ITLException("Unable to find dimensions and data type of deformed image file.");

			//raw::expandRawFilename(defFile);
			//if (!raw::internals::parseDimensions(defFile, defDimensions, defDT))
			//	throw ITLException("Unable to parse dimensions of deformed image from file name.");

			if (refDT != defDT)
				throw ITLException("Data types of reference and deformed images must be the same.");


			if (coarseBinning < 1)
				throw ITLException("Coarse binning must be greater than or equal to 1.");

			if (fineBinning < 1)
				throw ITLException("Fine binning must be greater than or equal to 1.");

			if(xmin > xmax || ymin > ymax || zmin > zmax)
				throw ITLException("Invalid reference grid definition.");


			PointGrid3D<coord_t> refPoints(PointGrid1D<coord_t>(xmin, xmax, xstep), PointGrid1D<coord_t>(ymin, ymax, ystep), PointGrid1D<coord_t>(zmin, zmax, zstep));
			Image<Vec3d> defPoints(refPoints.pointCounts());
			Image<float32_t> fitGoodness(defPoints.dimensions());
			double normFact;

			// Construct initial guess of the deformed points
			for (coord_t zi = 0; zi < defPoints.depth(); zi++)
			{
				for (coord_t yi = 0; yi < defPoints.height(); yi++)
				{
					for (coord_t xi = 0; xi < defPoints.width(); xi++)
					{
						defPoints(xi, yi, zi) = Vec3d(refPoints(xi, yi, zi)) + initialShift;
					}
				}
			}

			if (refDT == ImageDataType::UInt8)
			{
				blockMatchPartialLoad<uint8_t, uint8_t>(refFile, defFile, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, normalize, normFact);
			}
			else if (refDT == ImageDataType::UInt16)
			{
				blockMatchPartialLoad<uint16_t, uint16_t>(refFile, defFile, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, normalize, normFact);
			}
			else if (refDT == ImageDataType::Float32)
			{
				blockMatchPartialLoad<float32_t, float32_t>(refFile, defFile, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, normalize, normFact);
			}
			else
				throw ParseException("Unsupported image data type (Please add the data type to BlockMatchPartialLoadCommand in commands.h file).");

			writeBlockMatchResult(fname, refPoints, defPoints, fitGoodness, normFact);
		}
	};



	/**
	This version of pullback command reads arguments from disk.
	*/
	template<typename pixel_t> class PullbackCommand : public Command
	{
	public:
		PullbackCommand() : Command("pullback", "Applies reverse of a deformation (calculated using blockmatch command) to image. In other words, performs pull-back operation. Makes output image the same size than input image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "image", "Image that will be pulled back."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "pullback image", "Will store the result of the pullback operation."),
				CommandArgument<string>(ParameterDirection::In, "file name prefix", "File name prefix (and path) passed to blockmatch command."),
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& deformed = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& pullback = *pop<Image<pixel_t>* >(args);
			string fname = pop<string>(args);
			
			PointGrid3D<coord_t> refPoints;
			Image<Vec3d> defPoints;
			Image<float32_t> fitGoodness;
			double normFact;
			readBlockMatchResult(fname, refPoints, defPoints, fitGoodness, normFact);

			pullback.ensureSize(deformed);

			reverseDeformation(deformed, pullback, refPoints, defPoints, CubicInterpolator<pixel_t, pixel_t, double, double>(BoundaryCondition::Zero));
		}
	};

	/**
	This version of pullback command takes arguments as images.
	*/
	template<typename pixel_t> class PullbackNoDiskCommand : public Command
	{
	public:
		PullbackNoDiskCommand() : Command("pullback", "Applies reverse of a deformation (calculated using blockmatch command) to image. In other words, performs pull-back operation. Makes output image the same size than input image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "image", "Image that will be pulled back, i.e. the deformed image."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "pullback image", "Will store the result of the pullback operation, i.e. the deformed image transformed to coordinates of the reference image."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid start", "Start of reference point grid in the coordinates of the reference image."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid step", "Grid step in each coordinate direction."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid max", "End of reference point grid in the coordinates of the reference image. The grid will contain floor((max - start) / step) + 1 points in each coordinate direction. Difference between maximum and minimum does not need to be divisible by step."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "x", "X-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image must equal point counts in the reference grid."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "y", "Y-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image must equal point counts in the reference grid."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "z", "Z-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image must equal point counts in the reference grid.")
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& deformed = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& pullback = *pop<Image<pixel_t>* >(args);
			Vec3c gridStart = pop<Vec3c>(args);
			Vec3c gridStep = pop<Vec3c>(args);
			Vec3c gridEnd = pop<Vec3c>(args);
			Image<float32_t>& x = *pop<Image<float32_t>*>(args);
			Image<float32_t>& y = *pop<Image<float32_t>*>(args);
			Image<float32_t>& z = *pop<Image<float32_t>*>(args);


			PointGrid3D<coord_t> refPoints(
				PointGrid1D<coord_t>(gridStart.x, gridEnd.x, gridStep.x),
				PointGrid1D<coord_t>(gridStart.y, gridEnd.y, gridStep.y),
				PointGrid1D<coord_t>(gridStart.z, gridEnd.z, gridStep.z));

			if (x.dimensions() != refPoints.pointCounts())
				throw ITLException("Point counts in the reference grid must match sizes of x, y, and z images that contain reference grid points in deformed coordinates.");
			x.checkSize(y);
			x.checkSize(z);

			Image<Vec3d> defPoints;
			defPoints.ensureSize(x);

			for (coord_t n = 0; n < defPoints.pixelCount(); n++)
				defPoints(n) = Vec3d(x(n), y(n), z(n));

			pullback.ensureSize(deformed);

			reverseDeformation(deformed, pullback, refPoints, defPoints, CubicInterpolator<pixel_t, pixel_t, double, double>(BoundaryCondition::Zero));
		}
	};


	class FilterDisplacementsCommand : public Command
	{
	public:
		FilterDisplacementsCommand() : Command("filterdisplacements", "Helper for non-rigid stitching script. Filters displacements calculated by blockmatch commands.",
			{
				CommandArgument<string>(ParameterDirection::In, "file name prefix", "Value passed as file name prefix argument to blockmatch."),
				CommandArgument<double>(ParameterDirection::In, "threshold", "Threshold value for filtering. Displacements whose some component differs more than this value from median filtered displacements are considered to be bad.", 3.0)
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			string fname = pop<string>(args);
			double threshold = pop<double>(args);

			PointGrid3D<coord_t> refPoints;
			Image<Vec3d> defPoints;
			Image<float32_t> gof;
			double normFact;

			readBlockMatchResult(fname, refPoints, defPoints, gof, normFact);

			filterDisplacements(refPoints, defPoints, gof, 5, (float32_t)threshold);

			writeBlockMatchResult(fname + string("_filtered"), refPoints, defPoints, gof, normFact);
		}
	};


	//template<typename pixel_t> class StitchCommand : public Command
	//{
	//public:
	//	StitchCommand() : Command("stitch", "Stitches subimages to one big image, given geometric transformation for each subimage. NOTE: The size of the output image does not need to be correct. Pass in image of size (1, 1, 1) to save memory during the process and to let the command allocate the image after processing.",
	//		{
	//			CommandArgument<Image<pixel_t> >(ParameterDirection::In, "output image", "Output image."),
	//			//CommandArgument(ParameterDirection::In, ipt<uint8_t>(), "output weight image", "Output weight image."),
	//			CommandArgument<string>(ParameterDirection::In, "file list", "File name of index file that lists the files to be stitched."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "x", "X-coordinate of the region of the stitched image that will be output."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-coordinate of the region of the stitched image that will be output."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-coordinate of the region of the stitched image that will be output."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the output region."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the output region."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the output region."),
	//			CommandArgument(ParameterDirection::In, Bool, "normalize", "Set to true to make mean gray value of images the same in the overlapping region.", true)
	//		})
	//	{
	//	}

	//	virtual void run(vector<ParamVariant>& args) const
	//	{
	//		Image<pixel_t>* output = pop<Image<pixel_t>* >(args);
	//		//Image<uint8_t>* outputWeight = pop<Image<uint8_t>* >(args);
	//		string indexFile = pop<string>(args);
	//		coord_t x = pop<coord_t>(args);
	//		coord_t y = pop<coord_t>(args);
	//		coord_t z = pop<coord_t>(args);
	//		coord_t w = pop<coord_t>(args);
	//		coord_t h = pop<coord_t>(args);
	//		coord_t d = pop<coord_t>(args);
	//		bool normalize = pop<bool>(args);

	//		Vec3c pos(x, y, z);
	//		Vec3c size(w, h, d);

	//		// Read index file
	//		ifstream in(indexFile);
	//		vector<string> images;
	//		vector<string> transformations;
	//		//vector<string> srcWeights;
	//		while (in.good())
	//		{
	//			string imgFile;
	//			string transfFile;
	//			string weightFile;
	//			getline(in, imgFile);
	//			getline(in, transfFile);
	//			//getline(in, weightFile);
	//			if (imgFile.length() > 0 && transfFile.length() > 0)
	//			{
	//				images.push_back(imgFile);
	//				transformations.push_back(transfFile);
	//				//srcWeights.push_back(weightFile);
	//			}
	//		}
	//		
	//		
	//		//stitch<pixel_t>(images, transformations, srcWeights, pos, size, output, outputWeight, normalize);
	//		stitch<pixel_t>(images, transformations, pos, size, output, normalize);
	//	}
	//};


	class DetermineWorldToLocalCommand : public Command
	{
	public:
		DetermineWorldToLocalCommand() : Command("determine_world_to_local", "Helper for non-rigid stitching script. Used to run preprocessing for stitch_ver2 command by non-rigid stitcher script.",
			{
				CommandArgument<string>(ParameterDirection::In, "transformation file", "Name of transformation file to process."),
				CommandArgument<Vec3c>(ParameterDirection::In, "image size", "Size of the source image."),
				CommandArgument<string>(ParameterDirection::In, "world to local prefix", "Prefix for output files."),
				CommandArgument<bool>(ParameterDirection::In, "allow local shifts", "Set to true to allow non-rigid local deformations. Set to false to see the result without local deformations.")
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			string transfFile = pop<string>(args);
			Vec3c imageSize = pop<Vec3c>(args);
			string prefix = pop<string>(args);
			bool allowLocalShifts = pop<bool>(args);

			determineWorldToLocal(transfFile, imageSize, prefix, allowLocalShifts);
		}
	};

	template<typename pixel_t> class StitchVer2Command : public Command
	{
	public:
		StitchVer2Command() : Command("stitch_ver2", "Helper for non-rigid stitching script. Stitches subimages to one big image, given geometric transformation for each subimage. NOTE: The size of the output image does not need to be correct. Pass in image of size (1, 1, 1) to save memory during the process and to let the command allocate the image after processing.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "output image", "Output image."),
				//CommandArgument(ParameterDirection::In, ipt<uint8_t>(), "output weight image", "Output weight image."),
				CommandArgument<string>(ParameterDirection::In, "file list", "File name of index file that lists the files to be stitched."),
				CommandArgument<coord_t>(ParameterDirection::In, "x", "X-coordinate of the region of the stitched image that will be output."),
				CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-coordinate of the region of the stitched image that will be output."),
				CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-coordinate of the region of the stitched image that will be output."),
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the output region."),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the output region."),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the output region."),
				CommandArgument<bool>(ParameterDirection::In, "normalize", "Set to true to make mean gray value of images the same in the overlapping region.", true)
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& output = *pop<Image<pixel_t>* >(args);
			//Image<uint8_t>* outputWeight = pop<Image<uint8_t>* >(args);
			string indexFile = pop<string>(args);
			coord_t x = pop<coord_t>(args);
			coord_t y = pop<coord_t>(args);
			coord_t z = pop<coord_t>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);
			bool normalize = pop<bool>(args);

			Vec3c pos(x, y, z);
			Vec3c size(w, h, d);

			stitchVer2<pixel_t>(indexFile, pos, size, output, normalize);
		}
	};

	template<typename pixel_t> class DistanceMapCommand : public Command
	{
	public:
		DistanceMapCommand() : Command("dmap", "Calculates distance transform.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Input image where background is marked with background value given by the third argument."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "output image", "Output image (distance map) where pixel value equals Euclidean distance to the nearest background pixel."),
				CommandArgument<double>(ParameterDirection::In, "background value", "Pixels belonging to the background are marked with this value in the input image.", 0)
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& input = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& output = *pop<Image<float32_t>* >(args);
			double bgval = pop<double>(args);
			
			distanceTransform(input, output, NULL, pixelRound<pixel_t>(bgval));
		}
	};


	template<typename pixel_t> class HistogramCommand : public Command
	{
	public:
		HistogramCommand() : Command("hist", "Calculate histogram of an image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image whose histogram will be calculated."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "histogram", "Image that will contain the histogram on output."),
				CommandArgument<double>(ParameterDirection::In, "histogram minimum", "Minimum value to be included in the histogram.", 0),
				CommandArgument<double>(ParameterDirection::In, "histogram maximum", "The end of the last bin. This is the smallest value above minimum that is not included in the histogram.", 256),
				CommandArgument<coord_t>(ParameterDirection::In, "bin count", "Count of bins. The output image will be resized to contain this many pixels.", 256),
				CommandArgument<string>(ParameterDirection::In, "output file", "Name of file where the histogram data is to be saved in CSV format. Specify empty string to disable saving as text file.", "")
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& hist = *pop<Image<float32_t>* >(args);
			double hmin = pop<double>(args);
			double hmax = pop<double>(args);
			coord_t bincount = pop<coord_t>(args);
			string fname = pop<string>(args);

			hist.ensureSize(bincount);
			histogram(in, hist, Vec2d(hmin, hmax));

			if (fname.length() > 0)
			{
				ofstream out(fname);
				out << "Bin start [pixel]\tCount" << endl;
				for (coord_t n = 0; n < hist.pixelCount(); n++)
				{
					pixel_t binStart = pixelRound<pixel_t>(hmin + (double)n / (double)hist.pixelCount() * (hmax - hmin));
					out << binStart << "\t" << hist(n) << endl;
				}
			}
		}
	};

	template<typename pixel_t> class FloodFillCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		FloodFillCommand() : OneImageInPlaceCommand<pixel_t>("floodfill", "Performs flood fill. Fills start point and all its neighbours and their neighbours etc. recursively as long as the color of the pixel to be filled equals color of the start point.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "start point", "Starting point for the fill."),
				CommandArgument<double>(ParameterDirection::In, "fill value", "Fill color.")
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			Vec3c startPoint = pop<Vec3c>(args);
			pixel_t color = pixelRound<pixel_t>(pop<double>(args));


			floodfill(in, startPoint, color, color, Connectivity::AllNeighbours);
		}
	};

	template<typename pixel_t> class LabelCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		LabelCommand() : OneImageInPlaceCommand<pixel_t>("label", "Labels distinct regions in the image with individual colors.",
			{
				CommandArgument<double>(ParameterDirection::In, "region color", "Current color of regions that should be labeled.", 0)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			pixel_t color = pixelRound<pixel_t>(pop<double>(args));

			labelParticles(in, color);
		}
	};









	






	template<typename pixel_t> class NormalizeZCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		NormalizeZCommand() : OneImageInPlaceCommand<pixel_t>("normalizez", "Makes sure that all z-slices of the image have the same mean value.",
			{
				CommandArgument<double>(ParameterDirection::In, "target mean", "Global mean that the image should have after normalization. Specify nothing or nan to retain global mean of the image.", numeric_limits<double>::signaling_NaN())
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			double globalMean = pop<double>(args);
			normalizeZ(in, (float32_t)globalMean);
		}
	};

	




	template<typename pixel_t> class CannyPart1Command : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		CannyPart1Command() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("cannyPart1", "Performs first part of Canny edge detection; to be more precise, everything except edge tracking and final thresholding. This command is used in the distributed implementation of Canny edge detection. You probably should use 'canny' command instead of this one. Skips the initial Gaussian blurring step, please perform it separately if you want to do it. Calculates image derivatives using convolution with derivative of Gaussian.",
			{
				CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter for derivative calculation. Set to the preferred scale of edges that should be detected. Derivatives are calculated using convolutions with derivative of Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "lower threshold", "Edges that have gradient magnitude below lower threshold value are discarded. Edges that have gradient magnitude between lower and upper thresholds are included in the result only if they touch some edge that has gradient magnitude above upper threshold."),
				CommandArgument<double>(ParameterDirection::In, "upper threshold", "Edges that have gradient magnitude above upper threshold value are always included in the result. Edges that have gradient magnitude between lower and upper thresholds are included in the result only if they touch some edge that has gradient magnitude above upper threshold.")
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			double derSigma = pop<double>(args);
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			itl2::internals::cannyPart1(in, derSigma, loThreshold, hiThreshold);
		}

		virtual double calculateExtraMemory(vector<ParamVariant>& args) const
		{
			return 3.0 * sizeof(float32_t) / sizeof(pixel_t);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			double derSigma = get<double>(args[1]);

			coord_t margin = math::round(3 * derSigma) + 4;

			return Vec3c(margin, margin, margin);
		}
	};

	template<typename pixel_t> class CannyPart2Command : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		CannyPart2Command() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("cannyPart2", "Performs edge tracking part of Canny edge detection. This command is used in the distributed implementation of Canny edge detection. You probably should use 'canny' command instead of this one.",
			{})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			itl2::internals::cannyPart2(in);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			coord_t margin = 3;

			return Vec3c(margin, margin, margin);
		}
	};

	template<typename pixel_t> class CannyCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	public:
		CannyCommand() : OneImageInPlaceCommand<pixel_t>("canny", "Performs Canny edge detection. Skips the initial Gaussian blurring step, please perform it separately if you want to do it. Calculates image derivatives using convolution with derivative of Gaussian.",
			{
				CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter for derivative calculation. Set to the preferred scale of edges that should be detected. Derivatives are calculated using convolutions with derivative of Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "lower threshold", "Edges that have gradient magnitude below lower threshold value are discarded. Edges that have gradient magnitude between lower and upper thresholds are included in the result only if they touch some edge that has gradient magnitude above upper threshold."),
				CommandArgument<double>(ParameterDirection::In, "upper threshold", "Edges that have gradient magnitude above upper threshold value are always included in the result. Edges that have gradient magnitude between lower and upper thresholds are included in the result only if they touch some edge that has gradient magnitude above upper threshold.")
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			double derSigma = pop<double>(args);
			double loThreshold= pop<double>(args);
			double hiThreshold = pop<double>(args);

			canny(in, derSigma, loThreshold, hiThreshold);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			double derSigma = pop<double>(args);
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			CannyPart1Command<pixel_t> part1;
			part1.runDistributed(distributor, { &in, derSigma, loThreshold, hiThreshold });
			
			// Iterate edge tracking until there are no changes
			CannyPart2Command<pixel_t> part2;
			coord_t changed;
			do
			{
				vector<string> output = part2.runDistributed(distributor, { &in });

				changed = parseTotalCount(output, "pixels changed");
			} while (changed > 0);

			// Final thresholding to get rid of weak edges
			ThresholdConstantCommand<pixel_t> th;
			th.runDistributed(distributor, { &in, 1.0 });

			return vector<string>();
		}
	};


	template<typename pixel_t> class GrowCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	public:
		using Distributable::runDistributed;

		GrowCommand() : OneImageInPlaceCommand<pixel_t>("grow", "Grows regions with source color to regions with target color as much as possible.",
			{
				CommandArgument<double>(ParameterDirection::In, "source color", "Color that defines regions that are going to be grown."),
				CommandArgument<double>(ParameterDirection::In, "target color", "Color where the regions will grow.")
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			double src = pop<double>(args);
			double target = pop<double>(args);
			size_t changed = grow(in, pixelRound<pixel_t>(src), pixelRound<pixel_t>(target));
			cout << changed << " pixels changed." << endl;
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& img = *get<DistributedImage<pixel_t>* >(args[0]);
			double src = get<double>(args[1]);
			double target = get<double>(args[2]);

			coord_t changed;
			do
			{
				vector<string> output = distributor.distribute(this, args, 2, Vec3c(3, 3, 3));

				changed = parseTotalCount(output, "pixels changed.");
				cout << changed << " pixels changed." << endl;
			} while (changed > 0);

			return vector<string>();
		}
	};


	template<typename pixel_t> class DualThresholdCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	public:
		DualThresholdCommand() : OneImageInPlaceCommand<pixel_t>("dualthreshold", "First sets all pixels with value over upper threshold to 1. Then sets all regions to 1 that have value over lower threshold and that are connected to region that has value over upperThreshold.",
			{
				CommandArgument<double>(ParameterDirection::In, "lower threshold", "Regions that have value below lower threshold value are discarded. Regions that have value between lower and upper thresholds are included in the result only if they touch some region that has value above upper threshold."),
				CommandArgument<double>(ParameterDirection::In, "upper threshold", "Regions that have value above upper threshold value are always included in the result. Regions that have value between lower and upper thresholds are included in the result only if they touch some regoin that has value above upper threshold.")
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			dualThreshold(in, pixelRound<pixel_t>(loThreshold), pixelRound<pixel_t>(hiThreshold));
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>* >(args);
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			// Multi-threshold to two classes.
			DoubleThresholdCommand<pixel_t> mt;
			mt.runDistributed(distributor, { &img, loThreshold, hiThreshold });

			// Convert all those structures to "sure" that touch a "sure" structure.
			GrowCommand<pixel_t> grow;
			grow.runDistributed(distributor, { &img, 2.0, 1.0 });

			// Threshold so that only "sure" structures are left.
			ThresholdConstantCommand<pixel_t> th;
			th.runDistributed(distributor, { &img, 1.0 });

			return vector<string>();
		}
	};


	class CylindricalityCommand : public OverlapDistributable<OneImageInPlaceCommand<float32_t> >
	{
	public:
		CylindricalityCommand() : OverlapDistributable<OneImageInPlaceCommand<float32_t> >("cylindricality", "Estimates likelihood of structures being cylinders, based on eigenvalues of the local structure tensor.",
			{
				CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter for derivative calculation. Set to the preferred scale of edges that should be detected. Derivatives are calculated using convolutions with derivative of Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "The structure tensor is smoothed by convolution with Gaussian. This parameter defines the standard deviation of the smoothing Gaussian.", 1.0)
			})
		{
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			double derSigma = get<double>(args[1]);
			double smoothSigma = get<double>(args[2]);

			coord_t margin = math::round(3 * (derSigma + smoothSigma)) + 4;

			return Vec3c(margin, margin, margin);
		}

		virtual double calculateExtraMemory(vector<ParamVariant>& args) const
		{
			return 6.0;
		}

		virtual void run(Image<float32_t>& in, vector<ParamVariant>& args) const
		{
			double derSigma = pop<double>(args);
			double smoothSigma = pop<double>(args);

			structureTensor<float32_t>(in, derSigma, smoothSigma, 0, 0, 0, 0, 0, 0, 0, 0, 0, &in, 0, 0, 0.0);
		}

		virtual JobType getJobType() const
		{
			return JobType::Slow;
		}
	};

	// Structure tensor planarity does not work very well so the command is disabled for now.
	//class PlanarityCommand : public OneImageInPlaceCommand<float32_t>
	//{
	//public:
	//	PlanarityCommand() : OneImageInPlaceCommand<float32_t>("planarity", "Estimates likelihood of structures being planar, based on eigenvalues of the local structure tensor.",
	//		{
	//			CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter for derivative calculation. Set to the preferred scale of edges that should be detected. Derivatives are calculated using convolutions with derivative of Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
	//			CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "The structure tensor is smoothed by convolution with Gaussian. This parameter defines the standard deviation of the smoothing Gaussian.", 1.0)
	//		})
	//	{
	//	}

	//	virtual void run(Image<float32_t>& in, vector<ParamVariant>& args) const
	//	{
	//		double derSigma = pop<double>(args);
	//		double smoothSigma = pop<double>(args);

	//		structureTensor<float32_t>(in, derSigma, smoothSigma, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &in, 0, 0.0);
	//	}
	//};


	template<typename input_t> class NoiseCommand : public OneImageInPlaceCommand<input_t>
	{
	public:
		NoiseCommand() : OneImageInPlaceCommand<input_t>("noise", "Adds Gaussian noise to the image.",
			{
				CommandArgument<double>(ParameterDirection::In, "mean", "Mean value of the noise to add.", 0),
				CommandArgument<double>(ParameterDirection::In, "standard deviation", "Standard deviation of the noise to add.", 0)
			})
		{
		}

		virtual void run(Image<input_t>& in, vector<ParamVariant>& args) const
		{
			double mean = pop<double>(args);
			double std = pop<double>(args);
			noise(in, mean, std);
		}
	};
}

