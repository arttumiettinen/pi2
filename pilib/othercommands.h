#pragma once

#include "command.h"
#include "commandsbase.h"
#include "overlapdistributable.h"
#include "parseexception.h"
#include "filtercommands.h"
#include "pointprocesscommands.h"
#include "pilibutilities.h"
#include "registration.h"
#include "stitching.h"
#include "noise.h"
#include "misc.h"
#include "commandlist.h"
#include "standardhelp.h"
#include "montage.h"

#include <vector>
#include <string>

using namespace itl2;

namespace pilib
{


	template<typename pixel_t> class BlockMatchCommand : public Command
	{
	protected:
		friend class CommandList;
	
		BlockMatchCommand() : Command("blockmatch", "Calculates displacement field between two images. NOTE: This command is currently implemented in very old format, and thus it forcibly saves the results to a file.",
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
				CommandArgument<std::string>(ParameterDirection::In, "file name prefix", "Prefix (and path) of files to write. The command will save point grid in the reference image, corresponding points in the deformed image, and goodness-of-fit. If the files exists, the current contents are erased."),
				CommandArgument<Vec3c>(ParameterDirection::In, "comparison radius", "Radius of comparison region.", Vec3c(25, 25, 25))
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
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
			std::string fname = pop<std::string>(args);
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
	protected:
		friend class CommandList;
	
		BlockMatchPartialLoadCommand() : Command("blockmatchmemsave", "Calculates displacement between two images, loads only overlapping region from disk. NOTE: This command is currently implemented in very old format, and thus it forcibly saves the results to a file.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "reference image file", "Name of reference image file (non-moving image)."),
				CommandArgument<std::string>(ParameterDirection::In, "deformed image file", "Name of deformed image file (image to register to non-moving image)."),
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
				CommandArgument<std::string>(ParameterDirection::In, "file name prefix", "Prefix (and path) of files to write. The command will save point grid in the reference image, corresponding points in the deformed image, and goodness-of-fit. If the files exists, the current contents are erased."),
				CommandArgument<Vec3c>(ParameterDirection::In, "coarse comparison radius", "Radius of comparison region for coarse matching.", Vec3c(25, 25, 25)),
				CommandArgument<coord_t>(ParameterDirection::In, "coarse binning", "Amount of resolution reduction in coarse matching phase.", 2),
				CommandArgument<Vec3c>(ParameterDirection::In, "fine comparison radius", "Radius of comparison region for fine (full-resolution) matching.", Vec3c(10, 10, 10)),
				CommandArgument<coord_t>(ParameterDirection::In, "fine binning", "Amount of resolution reduction in fine matching phase. Set to same value than coarse binning to skip fine matching phase.", 2),
				CommandArgument<bool>(ParameterDirection::In, "normalize", "Indicates if the mean gray values of the two images should be made same in the overlapping region before matching.", true)
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			std::string refFile = pop<std::string>(args);
			std::string defFile = pop<std::string>(args);
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
			std::string fname = pop<std::string>(args);
			Vec3c coarseCompRadius = pop<Vec3c>(args);
			coord_t coarseBinning = pop<coord_t>(args);
			Vec3c fineCompRadius = pop<Vec3c>(args);
			coord_t fineBinning = pop<coord_t>(args);
			bool normalize = pop<bool>(args);

			Vec3c refDimensions;
			ImageDataType refDT;
			std::string reason;
			if (!io::getInfo(refFile, refDimensions, refDT, reason))
				throw ITLException(std::string("Unable to find dimensions and data type of reference image file. ") + reason);

			//raw::expandRawFilename(refFile);
			//if (!raw::internals::parseDimensions(refFile, refDimensions, refDT))
			//	throw ITLException("Unable to parse dimensions of reference image from file name.");

			Vec3c defDimensions;
			ImageDataType defDT;
			reason = "";
			if (!io::getInfo(defFile, defDimensions, defDT, reason))
				throw ITLException(std::string("Unable to find dimensions and data type of deformed image file. ") + reason);

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
	protected:
		friend class CommandList;
	
		PullbackCommand() : Command("pullback", "Applies reverse of a deformation (calculated using blockmatch command) to image. In other words, performs pull-back operation. Makes output image the same size than the input image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "image", "Image that will be pulled back."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "pullback image", "Will store the result of the pullback operation."),
				CommandArgument<std::string>(ParameterDirection::In, "file name prefix", "File name prefix (and path) passed to blockmatch command."),
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& deformed = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& pullback = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);
			
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
	protected:
		friend class CommandList;
	
		PullbackNoDiskCommand() : Command("pullback", "Applies reverse of a deformation (calculated using blockmatch command) to image. In other words, performs pull-back operation. Makes output image the same size than the input image.",
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

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
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
	protected:
		friend class CommandList;

		FilterDisplacementsCommand() : Command("filterdisplacements", "Helper for non-rigid stitching script. Filters displacements calculated by blockmatch commands.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "file name prefix", "Value passed as file name prefix argument to blockmatch."),
				CommandArgument<double>(ParameterDirection::In, "threshold", "Threshold value for filtering. Displacements whose some component differs more than this value from median filtered displacements are considered to be bad.", 3.0)
			})
		{
		}

	public:

		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			std::string fname = pop<std::string>(args);
			double threshold = pop<double>(args);

			PointGrid3D<coord_t> refPoints;
			Image<Vec3d> defPoints;
			Image<float32_t> gof;
			double normFact;

			readBlockMatchResult(fname, refPoints, defPoints, gof, normFact);

			filterDisplacements(refPoints, defPoints, gof, 5, (float32_t)threshold);

			writeBlockMatchResult(fname + std::string("_filtered"), refPoints, defPoints, gof, normFact);
		}
	};


	//template<typename pixel_t> class StitchCommand : public Command
	//{
	//public:
	//	StitchCommand() : Command("stitch", "Stitches subimages to one big image, given geometric transformation for each subimage. NOTE: The size of the output image does not need to be correct. Pass in image of size (1, 1, 1) to save memory during the process and to let the command allocate the image after processing.",
	//		{
	//			CommandArgument<Image<pixel_t> >(ParameterDirection::In, "output image", "Output image."),
	//			//CommandArgument(ParameterDirection::In, ipt<uint8_t>(), "output weight image", "Output weight image."),
	//			CommandArgument<std::string>(ParameterDirection::In, "file list", "File name of index file that lists the files to be stitched."),
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

	//	virtual void run(std::vector<ParamVariant>& args) const override
	//	{
	//		Image<pixel_t>* output = pop<Image<pixel_t>* >(args);
	//		//Image<uint8_t>* outputWeight = pop<Image<uint8_t>* >(args);
	//		std::string indexFile = pop<std::string>(args);
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
	//		std::vector<std::string> images;
	//		std::vector<std::string> transformations;
	//		//std::vector<std::string> srcWeights;
	//		while (in.good())
	//		{
	//			std::string imgFile;
	//			std::string transfFile;
	//			std::string weightFile;
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
	protected:
		friend class CommandList;

		DetermineWorldToLocalCommand() : Command("determine_world_to_local", "Helper for non-rigid stitching script. Used to run preprocessing for stitch_ver2 command by non-rigid stitcher script.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "transformation file", "Name of transformation file to process."),
				CommandArgument<Vec3c>(ParameterDirection::In, "image size", "Size of the source image."),
				CommandArgument<std::string>(ParameterDirection::In, "world to local prefix", "Prefix for output files."),
				CommandArgument<bool>(ParameterDirection::In, "allow local shifts", "Set to true to allow non-rigid local deformations. Set to false to see the result without local deformations.")
			})
		{
		}

	public:

		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			std::string transfFile = pop<std::string>(args);
			Vec3c imageSize = pop<Vec3c>(args);
			std::string prefix = pop<std::string>(args);
			bool allowLocalShifts = pop<bool>(args);

			determineWorldToLocal(transfFile, imageSize, prefix, allowLocalShifts);
		}
	};

	template<typename pixel_t> class StitchVer2Command : public Command
	{
	protected:
		friend class CommandList;

		StitchVer2Command() : Command("stitch_ver2", "Helper for non-rigid stitching script. Stitches subimages to one big image, given geometric transformation for each subimage. NOTE: The size of the output image does not need to be correct. Pass in image of size (1, 1, 1) to save memory during the process and to let the command allocate the image after processing.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::InOut, "output image", "Output image."),
				//CommandArgument(ParameterDirection::In, ipt<uint8_t>(), "output weight image", "Output weight image."),
				CommandArgument<std::string>(ParameterDirection::In, "file list", "File name of index file that lists the files to be stitched."),
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

	public:

		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& output = *pop<Image<pixel_t>* >(args);
			//Image<uint8_t>* outputWeight = pop<Image<uint8_t>* >(args);
			std::string indexFile = pop<std::string>(args);
			coord_t x = pop<coord_t>(args);
			coord_t y = pop<coord_t>(args);
			coord_t z = pop<coord_t>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);
			bool normalize = pop<bool>(args);

			Vec3c pos(x, y, z);
			Vec3c size(w, h, d);

			//stitchVer2<pixel_t>(indexFile, pos, size, output, normalize);
			stitchVer3<pixel_t>(indexFile, pos, size, output, nullptr, normalize);
		}
	};

	template<typename pixel_t> class StitchVer3Command : public Command
	{
	protected:
		friend class CommandList;

		StitchVer3Command() : Command("stitch_ver3", "Helper for non-rigid stitching script. Stitches subimages to one big image, given geometric transformation for each subimage. NOTE: The size of the output image does not need to be correct. Pass in image of size (1, 1, 1) to save memory during the process and to let the command allocate the image after processing. This is the same than stitch_ver2 command but creates additional goodness of stitching-output image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::InOut, "output image", "Output image."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::InOut, "goodness image", "This image will store indicator local goodness of match between the images. The indicator is standard deviation of all the overlapping images at each pixel."),
				CommandArgument<std::string>(ParameterDirection::In, "file list", "File name of index file that lists the files to be stitched."),
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

	public:

		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& output = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& goodness = *pop<Image<pixel_t>* >(args);
			std::string indexFile = pop<std::string>(args);
			coord_t x = pop<coord_t>(args);
			coord_t y = pop<coord_t>(args);
			coord_t z = pop<coord_t>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);
			bool normalize = pop<bool>(args);

			Vec3c pos(x, y, z);
			Vec3c size(w, h, d);

			stitchVer3<pixel_t>(indexFile, pos, size, output, &goodness, normalize);
		}
	};





	template<typename pixel_t> class FloodFillCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		FloodFillCommand() : OneImageInPlaceCommand<pixel_t>("floodfill", "Performs flood fill. Fills start point and all its neighbours and their neighbours etc. recursively as long as the color of the pixel to be filled equals color of the start point.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "start point", "Starting point for the fill."),
				CommandArgument<double>(ParameterDirection::In, "fill value", "Fill color."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the region to fill. ") + connectivityHelp(), Connectivity::AllNeighbours),
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			Vec3c startPoint = pop<Vec3c>(args);
			pixel_t color = pixelRound<pixel_t>(pop<double>(args));
			Connectivity connectivity = pop<Connectivity>(args);

			floodfill(in, startPoint, color, color, connectivity);
		}
	};









	






	template<typename pixel_t> class NormalizeZCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		NormalizeZCommand() : OneImageInPlaceCommand<pixel_t>("normalizez", "Makes sure that all z-slices of the image have the same mean value.",
			{
				CommandArgument<double>(ParameterDirection::In, "target mean", "Global mean that the image should have after normalization. Specify nothing or nan to retain global mean of the image.", std::numeric_limits<double>::signaling_NaN())
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			double globalMean = pop<double>(args);
			normalizeZ(in, (float32_t)globalMean);
		}
	};

	




	template<typename pixel_t> class CannyPart1Command : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		CannyPart1Command() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("cannyPart1", "Performs first part of Canny edge detection; to be more precise, everything except edge tracking and final thresholding. This command is used in the distributed implementation of Canny edge detection. You probably should use `canny` command instead of this one. Skips the initial Gaussian blurring step, please perform it separately if you want to do it. Calculates image derivatives using convolution with derivative of Gaussian.",
			{
				CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter for derivative calculation. Set to the preferred scale of edges that should be detected. Derivatives are calculated using convolutions with derivative of Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "lower threshold", "Edges that have gradient magnitude below lower threshold value are discarded. Edges that have gradient magnitude between lower and upper thresholds are included in the result only if they touch some edge that has gradient magnitude above upper threshold."),
				CommandArgument<double>(ParameterDirection::In, "upper threshold", "Edges that have gradient magnitude above upper threshold value are always included in the result. Edges that have gradient magnitude between lower and upper thresholds are included in the result only if they touch some edge that has gradient magnitude above upper threshold.")
			})
		{
		}

	public:
		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			double derSigma = pop<double>(args);
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			itl2::internals::cannyPart1(in, derSigma, loThreshold, hiThreshold);
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
		{
			return 3.0 * sizeof(float32_t) / sizeof(pixel_t);
		}

		virtual Vec3c calculateOverlap(const std::vector<ParamVariant>& args) const override
		{
			double derSigma = std::get<double>(args[1]);

			coord_t margin = itl2::round(3 * derSigma) + 4;

			return Vec3c(margin, margin, margin);
		}
	};

	template<typename pixel_t> class CannyPart2Command : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		CannyPart2Command() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("cannyPart2", "Performs edge tracking part of Canny edge detection. This command is used in the distributed implementation of Canny edge detection. You probably should use `canny` command instead of this one.",
			{})
		{
		}

	public:
		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			itl2::internals::cannyPart2(in);
		}

		virtual Vec3c calculateOverlap(const std::vector<ParamVariant>& args) const override
		{
			coord_t margin = 3;

			return Vec3c(margin, margin, margin);
		}

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			// This command relies on output, so we cannot delay it.
			return false;
		}
	};

	template<typename pixel_t> class CannyCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		CannyCommand() : OneImageInPlaceCommand<pixel_t>("canny", "Performs Canny edge detection. Skips the initial Gaussian blurring step, please perform it separately if you want to do it. Calculates image derivatives using convolution with derivative of Gaussian.",
			{
				CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter for derivative calculation. Set to the preferred scale of edges that should be detected. Derivatives are calculated using convolutions with derivative of Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "lower threshold", "Edges that have gradient magnitude below lower threshold value are discarded. Edges that have gradient magnitude between lower and upper thresholds are included in the result only if they touch some edge that has gradient magnitude above upper threshold."),
				CommandArgument<double>(ParameterDirection::In, "upper threshold", "Edges that have gradient magnitude above upper threshold value are always included in the result. Edges that have gradient magnitude between lower and upper thresholds are included in the result only if they touch some edge that has gradient magnitude above upper threshold.")
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			double derSigma = pop<double>(args);
			double loThreshold= pop<double>(args);
			double hiThreshold = pop<double>(args);

			canny(in, derSigma, loThreshold, hiThreshold);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			double derSigma = pop<double>(args);
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			auto& part1 = CommandList::get<CannyPart1Command<pixel_t> >();
			part1.runDistributed(distributor, { &in, derSigma, loThreshold, hiThreshold });
			
			// Iterate edge tracking until there are no changes
			auto& part2 = CommandList::get<CannyPart2Command<pixel_t> >();
			coord_t changed;
			do
			{
				std::vector<std::string> output = part2.runDistributed(distributor, { &in });

				changed = parseTotalCount(output, "pixels changed");
			} while (changed > 0);

			// Final thresholding to get rid of weak edges
			auto& th = CommandList::get<ThresholdConstantCommand<pixel_t> >();
			th.runDistributed(distributor, { &in, 1.0 });

			return std::vector<std::string>();
		}
	};


	template<typename pixel_t> class GrowCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		GrowCommand() : OneImageInPlaceCommand<pixel_t>("grow", "Grows regions with source color to regions with target color as much as possible.",
			{
				CommandArgument<double>(ParameterDirection::In, "source color", "Color that defines regions that are going to be grown."),
				CommandArgument<double>(ParameterDirection::In, "target color", "Color where the regions will grow."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the regions to grow. ") + connectivityHelp(), Connectivity::NearestNeighbours),
			},
			"growlabels, floodfill, regionremoval")
		{
		}

	public:
		using Distributable::runDistributed;

		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			double src = pop<double>(args);
			double target = pop<double>(args);
			Connectivity connectivity = pop<Connectivity>(args);

			size_t changed = grow(in, pixelRound<pixel_t>(src), pixelRound<pixel_t>(target), connectivity);
			std::cout << std::endl << changed << " pixels changed." << std::endl;
		}

		virtual Vec3c getMargin(const std::vector<ParamVariant>& args) const override
		{
			return Vec3c(3, 3, 3);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Normal;
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const
		{
			// Allocate some extra memory for priority queue
			return 1.0;
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			// TODO: This could be made with IterableDistributable
			coord_t changed;
			do
			{
				std::vector<std::string> output = distributor.distribute(this, args);

				changed = parseTotalCount(output, "pixels changed.");
				std::cout << std::endl << changed << " pixels changed." << std::endl;
			} while (changed > 0);

			return std::vector<std::string>();
		}
	};


	template<typename label_t, typename weight_t> class GrowPriorityCommand : public TwoImageInputParamCommand<label_t, weight_t>
	{
	protected:
		friend class CommandList;

		GrowPriorityCommand() : TwoImageInputParamCommand<label_t, weight_t>("grow",
			"Grows regions from seed points outwards. Seeds points are all nonzero pixels in the input image, pixel value defining region label. Each seed is grown towards surrounding zero pixels. Fill priority for each pixel is read from the corresponding pixel in the parameter image. Pixels for which priority is zero or negative are never filled. This process is equal to Meyer's watershed algorithm for given set of seeds, and watershed cuts are borders between filled regions in the output image.",
			{},
			"grow, growlabels, floodfill, regionremoval")
		{
		}

	public:
		virtual void run(Image<label_t>& labels, Image<weight_t>& weights, std::vector<ParamVariant>& args) const override
		{
			grow(labels, weights);
		}
	};



	template<typename pixel_t> class GrowLabelsCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		GrowLabelsCommand() : OneImageInPlaceCommand<pixel_t>("growlabels",
			"Grows all colored regions as much as possible into pixels that have a specific color. "
			"In practice, this command first finds all unique colors in the image, and uses each set of "
			"pixels having the same color as seed points for a flood fill that proceeds to pixels whose value is given in the 'allowed color' argument. "
			"\n\n"
			"This growing method is suited only for situations where separate parts of the original structure are labelled and "
			"the labels must be grown back to the original structure. **If there are multiple labels in "
			"a connected component, non-labeled pixels are assigned the smallest label in the non-distributed version "
			"and (mostly) random label among all the possibilities in the distributed version.** "
			"Therefore, **this function is suited only for images containing separate blobs or particles**, where each "
			"particle contains seed point(s) of only single value. "
			"\n\n"
			"An alternative to this command is `morphorec`. "
			"It works such that each pixel will get the label of the nearest labeled pixel.",
			{
				CommandArgument<double>(ParameterDirection::In, "allowed color", "Color where other colors will be grown into."),
				CommandArgument<double>(ParameterDirection::In, "background color", "Background color. Values of pixels having this color are not changed. Set to the same value than allowed color to fill to all pixels."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the regions to grow. ") + connectivityHelp(), Connectivity::NearestNeighbours),
			},
			"grow, growlabels, floodfill, regionremoval, morphorec")
		{
		}

	public:
		using Distributable::runDistributed;

		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			double allowed = pop<double>(args);
			double bg = pop<double>(args);
			Connectivity connectivity = pop<Connectivity>(args);

			size_t changed = growAll(in, pixelRound<pixel_t>(allowed), pixelRound<pixel_t>(bg), connectivity);
			std::cout << std::endl << changed << " pixels changed." << std::endl;
		}

		virtual Vec3c getMargin(const std::vector<ParamVariant>& args) const override
		{
			return Vec3c(3, 3, 3);
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const
		{
			// Allocate some extra memory for priority queue in filling
			return 2.0;
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			// TODO: This could be made with IterableDistributable
			coord_t changed;
			do
			{
				std::vector<std::string> output = distributor.distribute(this, args);

				changed = parseTotalCount(output, "pixels changed.");
				std::cout << std::endl << changed << " pixels changed." << std::endl;
			} while (changed > 0);

			return std::vector<std::string>();
		}
	};








	template<typename pixel_t> class DualThresholdCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		DualThresholdCommand() : OneImageInPlaceCommand<pixel_t>("dualthreshold", "First sets all pixels with value over upper threshold to 1. Then sets all regions to 1 that have value over lower threshold and that are connected to region that has value over upper threshold.",
			{
				CommandArgument<double>(ParameterDirection::In, "lower threshold", "Regions that have value below lower threshold value are discarded. Regions that have value between lower and upper thresholds are included in the result only if they touch some region that has value above upper threshold."),
				CommandArgument<double>(ParameterDirection::In, "upper threshold", "Regions that have value above upper threshold value are always included in the result. Regions that have value between lower and upper thresholds are included in the result only if they touch some regoin that has value above upper threshold.")
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			dualThreshold(in, pixelRound<pixel_t>(loThreshold), pixelRound<pixel_t>(hiThreshold));
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>* >(args);
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			// Multi-threshold to two classes.
			auto& mt = CommandList::get<DoubleThresholdCommand<pixel_t> >();
			mt.runDistributed(distributor, { &img, loThreshold, hiThreshold });

			// Convert all those structures to "sure" that touch a "sure" structure.
			auto& grow = CommandList::get<GrowCommand<pixel_t> >();
			grow.runDistributed(distributor, { &img, 2.0, 1.0, Connectivity::NearestNeighbours });

			// Threshold so that only "sure" structures are left.
			auto& th = CommandList::get<ThresholdConstantCommand<pixel_t> >();
			th.runDistributed(distributor, { &img, 1.0 });

			return std::vector<std::string>();
		}
	};





	template<typename input_t> class NoiseCommand : public OneImageInPlaceCommand<input_t>
	{
	protected:
		friend class CommandList;

		NoiseCommand() : OneImageInPlaceCommand<input_t>("noise", "Adds Gaussian noise to the image.",
			{
				CommandArgument<double>(ParameterDirection::In, "mean", "Mean value of the noise to add.", 0),
				CommandArgument<double>(ParameterDirection::In, "standard deviation", "Standard deviation of the noise to add.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "seed", "Seed value. Set to zero to use time-based seed.", 0)
			})
		{
		}

	public:
		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			double mean = pop<double>(args);
			double std = pop<double>(args);
			coord_t seed = pop<coord_t>(args);

			noise(in, mean, std, (unsigned int)seed);
		}
	};




	template<typename pixel_t> class MontageCommand : public TwoImageInputOutputCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		MontageCommand() : TwoImageInputOutputCommand<pixel_t>("montage", "Makes a 2D montage of a 3D image.",
			{
				CommandArgument<size_t>(ParameterDirection::In, "columns", "Number of 2D slices in the montage in the horizontal direction."),
				CommandArgument<size_t>(ParameterDirection::In, "rows", "Number of 2D slices in the montage in the vertical direction."),
				CommandArgument<double>(ParameterDirection::In, "scale", "Scaling factor between slices in the original image and the slices in the montage.", 1.0),
				CommandArgument<size_t>(ParameterDirection::In, "first slice", "First slice to include in the montage.", 0),
				CommandArgument<size_t>(ParameterDirection::In, "last slice", "Last slice to include in the montage. Note that the columns and rows parameters define the maximum number of slices that will fit into the montage.", std::numeric_limits<size_t>::max()),
				CommandArgument<size_t>(ParameterDirection::In, "step", "Step between slices to include in the montage. Specify zero to set the step to a value that accommodates approximately all the stack slices in the montage.", 0),
				CommandArgument<size_t>(ParameterDirection::In, "border width", "Width of borders between the slices in the montage.", 0),
				CommandArgument<double>(ParameterDirection::In, "border color", "Color of borders between slices in the montage.", 0)
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			size_t columns = pop<size_t>(args);
			size_t rows = pop<size_t>(args);
			double scale = pop<double>(args);
			size_t firstSlice = pop<size_t>(args);
			size_t lastSlice = pop<size_t>(args);
			size_t step = pop<size_t>(args);
			size_t borderWidth = pop<size_t>(args);
			double borderColor = pop<double>(args);

			montage(in, out, columns, rows, scale, firstSlice, lastSlice, step, borderWidth, pixelRound<pixel_t>(borderColor));
		}
	};
}

