#pragma once

#include "command.h"
#include "commandsbase.h"
#include "overlapdistributable.h"
#include "parseexception.h"
#include "filtercommands.h"
#include "pointprocesscommands.h"
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

	inline std::string blockMatchSeeAlso()
	{
		return "blockmatch, blockmatchmemsave, pullback, pointstodeformed";
	}



	template<typename pixel_t> class BlockMatchNoDiskCommand : public Command
	{
	protected:
		friend class CommandList;

		BlockMatchNoDiskCommand() : Command("blockmatch", "Calculates displacement field between two images.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "reference image", "Reference image (non-moving image)."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "deformed image", "Deformed image (image to register to non-moving image)."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid start", "Start of reference point grid in the coordinates of the reference image."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid step", "Grid step in each coordinate direction."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid max", "End of reference point grid in the coordinates of the reference image. The grid will contain floor((max - start) / step) + 1 points in each coordinate direction. Difference between maximum and minimum does not need to be divisible by step."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "x", "At output, contains the estimated X-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image are set to point counts in the reference grid."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "y", "At output, contains the estimated Y-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image are set to point counts in the reference grid."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "z", "At output, contains the estimated Z-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image are set to point counts in the reference grid."),
				CommandArgument<Vec3d>(ParameterDirection::In, "initial shift", "Initial shift between the images."),
				CommandArgument<Vec3c>(ParameterDirection::In, "comparison radius", "Radius of comparison region.", Vec3c(25, 25, 25)),
				CommandArgument<std::string>(ParameterDirection::In, "subpixel accuracy", "Subpixel accuracy mode. Can be 'none', 'quadratic', or 'centroid'.", "centroid")
			},
			blockMatchSeeAlso())
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& ref = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& def = *pop<Image<pixel_t>* >(args);

			Vec3c gridStart = pop<Vec3c>(args);
			Vec3c gridStep = pop<Vec3c>(args);
			Vec3c gridEnd = pop<Vec3c>(args);
			Image<float32_t>& x = *pop<Image<float32_t>*>(args);
			Image<float32_t>& y = *pop<Image<float32_t>*>(args);
			Image<float32_t>& z = *pop<Image<float32_t>*>(args);

			Vec3d initialShift = pop<Vec3d>(args);
			Vec3c compRadius = pop<Vec3c>(args);
			SubpixelAccuracy mode = fromString<SubpixelAccuracy>(pop<string>(args));

			PointGrid3D<coord_t> refPoints(
				PointGrid1D<coord_t>(gridStart.x, gridEnd.x, gridStep.x),
				PointGrid1D<coord_t>(gridStart.y, gridEnd.y, gridStep.y),
				PointGrid1D<coord_t>(gridStart.z, gridEnd.z, gridStep.z));

			x.ensureSize(refPoints.pointCounts());
			y.ensureSize(refPoints.pointCounts());
			z.ensureSize(refPoints.pointCounts());

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

			blockMatch(ref, def, refPoints, defPoints, fitGoodness, compRadius, mode);

			// Copy blockmatch output to output images.
			for (coord_t zi = 0; zi < defPoints.depth(); zi++)
			{
				for (coord_t yi = 0; yi < defPoints.height(); yi++)
				{
					for (coord_t xi = 0; xi < defPoints.width(); xi++)
					{
						Vec3d dp = defPoints(xi, yi, zi);
						x(xi, yi, zi) = (float32_t)dp.x;
						y(xi, yi, zi) = (float32_t)dp.y;
						z(xi, yi, zi) = (float32_t)dp.z;
					}
				}
			}
		}
	};


	template<typename pixel_t> class BlockMatchNoDiskMultiCommand : public Command
	{
	protected:
		friend class CommandList;

		BlockMatchNoDiskMultiCommand() : Command("blockmatch", "Calculates displacement field between two images with two-step multi-resolution approach, where coarse displacement is first calculated with larger block size (and binning) and the result is refined in second phase with smaller block size (and binning).",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "reference image", "Reference image (non-moving image)."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "deformed image", "Deformed image (image to register to non-moving image)."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid start", "Start of reference point grid in the coordinates of the reference image."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid step", "Grid step in each coordinate direction."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid max", "End of reference point grid in the coordinates of the reference image. The grid will contain floor((max - start) / step) + 1 points in each coordinate direction. Difference between maximum and minimum does not need to be divisible by step."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "x", "At output, contains the estimated X-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image are set to point counts in the reference grid."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "y", "At output, contains the estimated Y-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image are set to point counts in the reference grid."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "z", "At output, contains the estimated Z-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image are set to point counts in the reference grid."),
				CommandArgument<Vec3d>(ParameterDirection::In, "initial shift", "Initial shift between the images."),
				CommandArgument<Vec3c>(ParameterDirection::In, "coarse comparison radius", "Radius of comparison region for coarse matching.", Vec3c(25, 25, 25)),
				CommandArgument<size_t>(ParameterDirection::In, "coarse binning", "Amount of resolution reduction in coarse matching phase.", 2),
				CommandArgument<Vec3c>(ParameterDirection::In, "fine comparison radius", "Radius of comparison region for fine (full-resolution) matching.", Vec3c(10, 10, 10)),
				CommandArgument<size_t>(ParameterDirection::In, "fine binning", "Amount of resolution reduction in fine matching phase. Set to same value than coarse binning to skip fine matching phase.", 1),
				CommandArgument<std::string>(ParameterDirection::In, "subpixel accuracy", "Subpixel accuracy mode. Can be 'none', 'quadratic', or 'centroid'.", "centroid")
			},
			blockMatchSeeAlso())
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& ref = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& def = *pop<Image<pixel_t>* >(args);

			Vec3c gridStart = pop<Vec3c>(args);
			Vec3c gridStep = pop<Vec3c>(args);
			Vec3c gridEnd = pop<Vec3c>(args);
			Image<float32_t>& x = *pop<Image<float32_t>*>(args);
			Image<float32_t>& y = *pop<Image<float32_t>*>(args);
			Image<float32_t>& z = *pop<Image<float32_t>*>(args);

			Vec3d initialShift = pop<Vec3d>(args);

			Vec3c coarseCompRadius = pop<Vec3c>(args);
			size_t coarseBinning = pop<size_t>(args);
			Vec3c fineCompRadius = pop<Vec3c>(args);
			size_t fineBinning = pop<size_t>(args);

			SubpixelAccuracy mode = fromString<SubpixelAccuracy>(pop<string>(args));

			PointGrid3D<coord_t> refPoints(
				PointGrid1D<coord_t>(gridStart.x, gridEnd.x, gridStep.x),
				PointGrid1D<coord_t>(gridStart.y, gridEnd.y, gridStep.y),
				PointGrid1D<coord_t>(gridStart.z, gridEnd.z, gridStep.z));

			x.ensureSize(refPoints.pointCounts());
			y.ensureSize(refPoints.pointCounts());
			z.ensureSize(refPoints.pointCounts());

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

			blockMatchMulti(ref, def, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, mode);

			// Copy blockmatch output to output images.
			for (coord_t zi = 0; zi < defPoints.depth(); zi++)
			{
				for (coord_t yi = 0; yi < defPoints.height(); yi++)
				{
					for (coord_t xi = 0; xi < defPoints.width(); xi++)
					{
						Vec3d dp = defPoints(xi, yi, zi);
						x(xi, yi, zi) = (float32_t)dp.x;
						y(xi, yi, zi) = (float32_t)dp.y;
						z(xi, yi, zi) = (float32_t)dp.z;
					}
				}
			}
		}
	};




	template<typename pixel_t> class BlockMatchCommand : public Command
	{
	protected:
		friend class CommandList;

		BlockMatchCommand() : Command("blockmatch", "Calculates displacement field between two images. NOTE: This command is deprecated as it forcibly saves the results to disk. Consider using the version with output to variables.",
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
				CommandArgument<Vec3c>(ParameterDirection::In, "comparison radius", "Radius of comparison region.", Vec3c(25, 25, 25)),
				CommandArgument<std::string>(ParameterDirection::In, "subpixel accuracy", "Subpixel accuracy mode. Can be 'none', 'quadratic', or 'centroid'.", "centroid")
			},
			blockMatchSeeAlso())
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
			SubpixelAccuracy mode = fromString<SubpixelAccuracy>(pop<string>(args));

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

			blockMatch(ref, def, refPoints, defPoints, fitGoodness, compRadius, mode);

			writeBlockMatchResult(fname, refPoints, defPoints, fitGoodness, 0, 1, 0);
		}
	};

	template<typename pixel_t> class BlockMatchMultiCommand : public Command
	{
	protected:
		friend class CommandList;

		BlockMatchMultiCommand() : Command("blockmatch", "Calculates displacement field between two images with two-step multi-resolution approach, where coarse displacement is first calculated with larger block size (and binning) and the result is refined in second phase with smaller block size (and binning). NOTE: This command is deprecated as it forcibly saves the results to disk. Consider using the version with output to variables.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "reference image", "Reference image (non-moving image)."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "deformed image", "Deformed image (image to register to non-moving image)."),
				CommandArgument<Vec3c>(ParameterDirection::In, "x grid", "Calculation point grid definition in X-direction. The format is [start coordinate, end coordinate, step]."),
				CommandArgument<Vec3c>(ParameterDirection::In, "y grid", "Calculation point grid definition in Y-direction. The format is [start coordinate, end coordinate, step]."),
				CommandArgument<Vec3c>(ParameterDirection::In, "z grid", "Calculation point grid definition in Z-direction. The format is [start coordinate, end coordinate, step]."),
				CommandArgument<Vec3d>(ParameterDirection::In, "initial shift", "Initial shift between the images."),
				CommandArgument<std::string>(ParameterDirection::In, "file name prefix", "Prefix (and path) of files to write. The command will save point grid in the reference image, corresponding points in the deformed image, and goodness-of-fit. If the files exists, the current contents are erased."),
				CommandArgument<Vec3c>(ParameterDirection::In, "coarse comparison radius", "Radius of comparison region for coarse matching.", Vec3c(25, 25, 25)),
				CommandArgument<size_t>(ParameterDirection::In, "coarse binning", "Amount of resolution reduction in coarse matching phase.", 2),
				CommandArgument<Vec3c>(ParameterDirection::In, "fine comparison radius", "Radius of comparison region for fine (full-resolution) matching.", Vec3c(10, 10, 10)),
				CommandArgument<size_t>(ParameterDirection::In, "fine binning", "Amount of resolution reduction in fine matching phase. Set to same value than coarse binning to skip fine matching phase.", 1),
				CommandArgument<std::string>(ParameterDirection::In, "subpixel accuracy", "Subpixel accuracy mode. Can be 'none', 'quadratic', or 'centroid'.", "centroid")
			},
			blockMatchSeeAlso())
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& ref = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& def = *pop<Image<pixel_t>* >(args);
			Vec3c xGrid = pop<Vec3c>(args);
			Vec3c yGrid = pop<Vec3c>(args);
			Vec3c zGrid = pop<Vec3c>(args);
			Vec3d initialShift = pop<Vec3d>(args);
			std::string fname = pop<std::string>(args);
			Vec3c coarseCompRadius = pop<Vec3c>(args);
			size_t coarseBinning = pop<size_t>(args);
			Vec3c fineCompRadius = pop<Vec3c>(args);
			size_t fineBinning = pop<size_t>(args);
			SubpixelAccuracy mode = fromString<SubpixelAccuracy>(pop<string>(args));

			coord_t xmin = xGrid.x;
			coord_t xmax = xGrid.y;
			coord_t xstep = xGrid.z;
			coord_t ymin = yGrid.x;
			coord_t ymax = yGrid.y;
			coord_t ystep = yGrid.z;
			coord_t zmin = zGrid.x;
			coord_t zmax = zGrid.y;
			coord_t zstep = zGrid.z;

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

			blockMatchMulti(ref, def, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, mode);

			writeBlockMatchResult(fname, refPoints, defPoints, fitGoodness, 0, 1, 0);
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
				CommandArgument<bool>(ParameterDirection::In, "normalize", "Indicates if the mean gray values of the two images should be made same in the overlapping region before matching.", true),
				CommandArgument<string>(ParameterDirection::In, "subpixel accuracy", "Subpixel accuracy mode. Can be 'none', 'quadratic', or 'centroid'.", "centroid")
			},
			blockMatchSeeAlso())
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
			SubpixelAccuracy mode = fromString<SubpixelAccuracy>(pop<string>(args));

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
			double normFact, normFactStd, meanDef;

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
				blockMatchPartialLoad<uint8_t, uint8_t>(refFile, defFile, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, normalize, normFact, normFactStd, meanDef, mode);
			}
			else if (refDT == ImageDataType::UInt16)
			{
				blockMatchPartialLoad<uint16_t, uint16_t>(refFile, defFile, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, normalize, normFact, normFactStd, meanDef, mode);
			}
			else if (refDT == ImageDataType::Float32)
			{
				blockMatchPartialLoad<float32_t, float32_t>(refFile, defFile, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, normalize, normFact, normFactStd, meanDef, mode);
			}
			else
				throw ParseException("Unsupported image data type (Please add the data type to BlockMatchPartialLoadCommand in commands.h file).");

			writeBlockMatchResult(fname, refPoints, defPoints, fitGoodness, normFact, normFactStd, meanDef);
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
				CommandArgument<InterpolationMode>(ParameterDirection::In, "interpolation mode", string("Interpolation mode. ") + interpolationHelp(), InterpolationMode::Cubic),
			},
			blockMatchSeeAlso())
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& deformed = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& pullback = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);
			InterpolationMode imode = pop<InterpolationMode>(args);
			auto interp = createInterpolator<pixel_t, pixel_t, double, double>(imode, BoundaryCondition::Zero);
			
			PointGrid3D<coord_t> refPoints;
			Image<Vec3d> defPoints;
			Image<float32_t> fitGoodness;
			double normFact, normFactStd, meanDef;
			readBlockMatchResult(fname, refPoints, defPoints, fitGoodness, normFact, normFactStd, meanDef);

			pullback.ensureSize(deformed);

			reverseDeformation(deformed, pullback, refPoints, defPoints, *interp, false);
		}
	};




	class PointsToDeformedCommand : public Command
	{
	protected:
		friend class CommandList;

		PointsToDeformedCommand() : Command("pointstodeformed", "Projects points from reference configuration to deformed configuration, using a transformation determined with the `blockmatch` command.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "points", "Image that contains the points that will be transformed. The size of the image must be 3xN where N is the number of points to transform."),
				CommandArgument<std::string>(ParameterDirection::In, "file name prefix", "File name prefix (and path) passed to blockmatch command."),
			},
			blockMatchSeeAlso())
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<float32_t>& points = *pop<Image<float32_t>* >(args);
			std::string fname = pop<std::string>(args);

			PointGrid3D<coord_t> refPoints;
			Image<Vec3d> defPoints;
			Image<float32_t> fitGoodness;
			double normFact, normFactStd, meanDef;
			readBlockMatchResult(fname, refPoints, defPoints, fitGoodness, normFact, normFactStd, meanDef);

			vector<Vec3d> pointsv;
			pointsv.reserve(points.height());
			for (coord_t n = 0; n < points.height(); n++)
			{
				Vec3d v(points(0, n), points(1, n), points(2, n));
				pointsv.push_back(v);
			}

			pointsToDeformed(pointsv, refPoints, defPoints);

			for (size_t n = 0; n < pointsv.size(); n++)
			{
				const auto& v = pointsv[n];
				points(0, n) = (float32_t)v.x;
				points(1, n) = (float32_t)v.y;
				points(2, n) = (float32_t)v.z;
			}
		}
	};




	/**
	This version of pullback command takes arguments as images.
	*/
	template<typename pixel_t> class PullbackNoDiskCommand : public Command
	{
	protected:
		friend class CommandList;
	
		PullbackNoDiskCommand() : Command("pullback", "Applies reverse of a deformation (calculated, e.g., using the blockmatch command) to image. In other words, performs pull-back operation.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "image", "Image that will be pulled back, i.e. the deformed image."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "pullback image", "Will store the result of the pullback operation, i.e. the deformed image transformed to coordinates of the reference image."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid start", "Start of reference point grid in the coordinates of the reference image."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid step", "Grid step in each coordinate direction."),
				CommandArgument<Vec3c>(ParameterDirection::In, "grid max", "End of reference point grid in the coordinates of the reference image. The grid will contain floor((max - start) / step) + 1 points in each coordinate direction. Difference between maximum and minimum does not need to be divisible by step."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "x", "X-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image must equal point counts in the reference grid."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "y", "Y-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image must equal point counts in the reference grid."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "z", "Z-coordinate of each reference grid point in the coordinates of the deformed image. Dimensions of this image must equal point counts in the reference grid."),
				CommandArgument<InterpolationMode>(ParameterDirection::In, "interpolation mode", string("Interpolation mode. ") + interpolationHelp(), InterpolationMode::Cubic),
				CommandArgument<Vec3d>(ParameterDirection::In, "pullback position", "Position of region to be pulled back in reference image coordinates.", Vec3d(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "pullback size", "Size of the region to be pulled back. Specify zeroes to default to the size of the deformed image.", Vec3c(0, 0, 0))
			},
			blockMatchSeeAlso())
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
			InterpolationMode imode = pop<InterpolationMode>(args);
			auto interp = createInterpolator<pixel_t, pixel_t, double, double>(imode, BoundaryCondition::Zero);
			Vec3d pullbackPos = pop<Vec3d>(args);
			Vec3c pullbackSize = pop<Vec3c>(args);

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

			if (pullbackSize == Vec3c(0, 0, 0))
				pullbackSize = deformed.dimensions();

			pullback.ensureSize(pullbackSize);

			reverseDeformation(deformed, pullback, refPoints, defPoints, pullbackPos, *interp, false);
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
			},
			blockMatchSeeAlso())
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
			double normFact, normFactStd, meanDef;

			readBlockMatchResult(fname, refPoints, defPoints, gof, normFact, normFactStd, meanDef);

			filterDisplacements(refPoints, defPoints, gof, 5, (float32_t)threshold);

			writeBlockMatchResult(fname + std::string("_filtered"), refPoints, defPoints, gof, normFact, normFactStd, meanDef);
		}
	};



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
			},
			blockMatchSeeAlso())
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
				CommandArgument<std::string>(ParameterDirection::In, "file list", "File name of index file that lists the files to be stitched."),
				CommandArgument<coord_t>(ParameterDirection::In, "x", "X-coordinate of the region of the stitched image that will be output."),
				CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-coordinate of the region of the stitched image that will be output."),
				CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-coordinate of the region of the stitched image that will be output."),
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the output region."),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the output region."),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the output region."),
				CommandArgument<bool>(ParameterDirection::In, "normalize", "Set to true to make mean gray value of images the same in the overlapping region.", true),
				CommandArgument<bool>(ParameterDirection::In, "mask max circle", "Set to true to use only data in the maximum inscribed circle in each xy-slice of the input image.", false)
			},
			blockMatchSeeAlso())
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
			std::string indexFile = pop<std::string>(args);
			coord_t x = pop<coord_t>(args);
			coord_t y = pop<coord_t>(args);
			coord_t z = pop<coord_t>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);
			bool normalize = pop<bool>(args);
			bool maxCircle = pop<bool>(args);

			Vec3c pos(x, y, z);
			Vec3c size(w, h, d);

			stitchVer3<pixel_t>(indexFile, pos, size, output, nullptr, normalize, maxCircle);
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
				CommandArgument<bool>(ParameterDirection::In, "normalize", "Set to true to make mean gray value of images the same in the overlapping region.", true),
				CommandArgument<bool>(ParameterDirection::In, "mask max circle", "Set to true to use only data in the maximum inscribed circle in each xy-slice of the input image.", false)
			},
			blockMatchSeeAlso())
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
			bool maxCircle = pop<bool>(args);

			Vec3c pos(x, y, z);
			Vec3c size(w, h, d);

			stitchVer3<pixel_t>(indexFile, pos, size, output, &goodness, normalize, maxCircle);
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







	template<typename input_t> class NoiseCommand : public OneImageInPlaceCommand<input_t>
	{
	protected:
		friend class CommandList;

		NoiseCommand() : OneImageInPlaceCommand<input_t>("noise", "Adds additive Gaussian noise to the image.",
			{
				CommandArgument<double>(ParameterDirection::In, "mean", "Mean value of the noise to add.", 0),
				CommandArgument<double>(ParameterDirection::In, "standard deviation", "Standard deviation of the noise to add. Specify zero to select standard deviation based on typical maximum value range of the pixel data type.", 0),
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

