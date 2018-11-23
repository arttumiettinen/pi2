#pragma once

#include "command.h"
#include "commandsbase.h"
#include "overlapdistributable.h"
#include "parseexception.h"

#include <vector>
#include <string>

#include "itl2.h"

using namespace std;
using namespace itl2;

namespace pilib
{

	namespace internals
	{
		/*
		Get value at the first line that reads "key = value".
		*/
		double getValue(const string& s, const string& key)
		{
			size_t startPos = s.find(key);
			if (startPos != string::npos)
			{
				startPos += key.length();
				size_t lineEnd = s.find('\n', startPos);
				string number = s.substr(startPos, lineEnd - startPos - 1);
				return fromString<double>(number);
			}
			else
			{
				throw ITLException("Key not found.");
			}
		}

		double sumReducer(const vector<double> vals, const vector<double> counts)
		{
			return sum(vals);
		}

		double minReducer(const vector<double> vals, const vector<double> counts)
		{
			return min(vals);
		}

		double maxReducer(const vector<double> vals, const vector<double> counts)
		{
			return max(vals);
		}

		double meanReducer(const vector<double> vals, const vector<double> counts)
		{
			return mean(vals, counts);
		}
	}

#define CONCAT(str1, str2) #str1 #str2
#define DEF_PROJECT(classname, funcname, cmdname) \
	template<typename in_t> class classname##AllPixelsCommand : public TwoImageInputOutputCommand<in_t, float32_t>, public Distributable	\
	{																										\
	public:																									\
		classname##AllPixelsCommand() : TwoImageInputOutputCommand<in_t, float32_t>(#cmdname, "Calculates " #funcname " of all pixels in the input image. The output is a 1x1x1 image.",	\
			{																								\
				CommandArgument<bool>(In, "print to log", "Set to true to print the results to the log.", false)	\
			}) {}																							\
																											\
		virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const			\
		{																									\
			bool print = pop<bool>(args);																	\
			double res = (double)funcname (in);																\
			out.ensureSize(1, 1, 1);																		\
			out(0) = (float32_t)res;																		\
			if(print)																						\
			{																								\
				cout << #funcname << " = " << out(0) << endl;												\
				cout << "count = " << in.pixelCount() << endl;												\
			}																								\
		}																									\
																											\
		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const				\
		{																									\
			DistributedImage<float32_t>& out = *get<DistributedImage<float32_t>* >(args[1]);				\
			bool print = get<bool>(args[2]);																\
			out.ensureSize(1, 1, 1);																		\
																											\
			vector<ParamVariant> args2;																		\
			args2.push_back(args[0]);																		\
			args2.push_back(args[1]);																		\
			ParamVariant p;																					\
			p.bval = true;																					\
			args2.push_back(p);																				\
																											\
			vector<string> results = distributor.distribute(this, args2, 2, Vec3c(0, 0, 0), 0, 0);			\
			vector<double> vals, counts;																	\
			for (const string& s : results)																	\
			{																								\
				vals.push_back(internals::getValue(s, #funcname " = "));									\
				counts.push_back(internals::getValue(s, "count = "));										\
			}																								\
																											\
			Image<float32_t> tmp;																			\
			out.readTo(tmp);																				\
			tmp(0) = (float32_t)internals::funcname##Reducer(vals, counts);									\
			out.setData(tmp);																				\
																											\
			if (print)																						\
			{																								\
				cout << #funcname << " = " << tmp(0) << endl;												\
				cout << "count = " << sum(counts) << endl;													\
			}																								\
		}																									\
																											\
		virtual void getCorrespondingBlock(vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const	\
		{																									\
			if (argIndex == 1)																				\
			{																								\
				readStart = Vec3c(0, 0, 0);																	\
				readSize = Vec3c(1, 1, 1);																	\
				writeFilePos = Vec3c(0, 0, 0);																\
				writeImPos = Vec3c(0, 0, 0);																\
				writeSize = Vec3c(1, 1, 1);																	\
			}																								\
		}																									\
	};																										\
																											\
	template<typename in_t> class classname##ProjectCommand : public TwoImageInputOutputCommand<in_t, float32_t>, public Distributable	\
	{																										\
	public:																									\
		classname##ProjectCommand() : TwoImageInputOutputCommand<in_t, float32_t>(CONCAT(funcname, project), "Calculates projection of the input image. The dimensionality of the output image is the dimensionality of the input image subtracted by one.",	\
		{ CommandArgument<size_t>(In, "dimension", "Dimension to project over, zero corresponding to x, one corresponding to y, and 2 corresponding to z.", 2) }) {}	\
																											\
		virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const			\
		{																									\
			size_t dim = pop<size_t>(args);																	\
			if(dim < 0 || dim > 2)																			\
				throw ITLException("Invalid dimension specification.");										\
																											\
			funcname (in, dim, out);																		\
		}																									\
																											\
		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const				\
		{																									\
			DistributedImage<in_t>& in = *get<DistributedImage<in_t>* >(args[0]);							\
			DistributedImage<float32_t>& out = *get<DistributedImage<float32_t>* >(args[1]);				\
			size_t dim = get<size_t>(args[2]);																\
																											\
			size_t distrDim;																				\
			if (dim == 2)																					\
			{																								\
				/* Z project */																				\
				out.ensureSize(in.width(), in.height());													\
				distrDim = 1;																				\
			}																								\
			else if (dim == 1)																				\
			{																								\
				/* Y project */																				\
				out.ensureSize(in.width(), in.depth());														\
				distrDim = 1;																				\
			}																								\
			else if (dim == 0)																				\
			{																								\
				/* X project. Swap z and y in output to make the image more logical (in some sense...) */	\
				out.ensureSize(in.depth(), in.height());													\
				distrDim = 1;																				\
			}																								\
			else																							\
			{																								\
				throw ITLException("Invalid dimension specification.");										\
			}																								\
																											\
			distributor.distribute(this, args, distrDim, Vec3c(0, 0, 0));									\
		}																									\
																											\
		virtual void getCorrespondingBlock(vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const	\
		{																									\
			if (argIndex == 0)																				\
			{																								\
				DistributedImage<in_t>& in = *get<DistributedImage<in_t>* >(args[0]);						\
				size_t dim = get<size_t>(args[2]);															\
				if(dim == 2)																				\
				{																							\
					readSize = Vec3c(readSize.x, readSize.y, in.depth());									\
					readStart = Vec3c(readStart.x, readStart.y, 0);											\
				}																							\
				else if (dim == 1)																			\
				{																							\
					readSize = Vec3c(readSize.x, in.height(), readSize.y);									\
					readStart = Vec3c(readStart.x, 0, readStart.y);											\
				}																							\
				else /* dim == 0 */																			\
				{																							\
					readSize = Vec3c(in.width(), readSize.y, readSize.x);									\
					readStart = Vec3c(0, readStart.y, readStart.x);											\
				}																							\
			}																								\
		}																									\
	};

#undef min
#undef max
	DEF_PROJECT(Sum, sum, sum)
	DEF_PROJECT(Min, min, minval)
	DEF_PROJECT(Max, max, maxval)
	DEF_PROJECT(Mean, mean, mean)



	template<typename pixel_t> class BlockMatchCommand : public Command
	{
	public:
		BlockMatchCommand() : Command("blockmatch", "Calculates displacement field between two images.",
			{
				CommandArgument<Image<pixel_t> >(In, "reference image", "Reference image (non-moving image)."),
				CommandArgument<Image<pixel_t> >(In, "deformed image", "Deformed image (image to register to non-moving image)."),
				CommandArgument<coord_t>(In, "xmin", "X-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(In, "xmax", "X-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(In, "xstep", "Step between calculation points in x-direction."),
				CommandArgument<coord_t>(In, "ymin", "Y-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(In, "ymax", "Y-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(In, "ystep", "Step between calculation points in y-direction."),
				CommandArgument<coord_t>(In, "zmin", "Z-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(In, "zmax", "Z-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(In, "zstep", "Step between calculation points in z-direction."),
				CommandArgument<Vec3d>(In, "initial shift", "Initial shift between the images."),
				CommandArgument<string>(In, "file name prefix", "Prefix (and path) of files to write. The command will save point grid in the reference image, corresponding points in the deformed image, and goodness-of-fit. If the files exists, the current contents are erased."),
				CommandArgument<Vec3c>(In, "comparison radius", "Radius of comparison region.", Vec3c(25, 25, 25))
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
				CommandArgument<string>(In, "reference image file", "Name of reference image file (non-moving image)."),
				CommandArgument<string>(In, "deformed image file", "Name of deformed image file (image to register to non-moving image)."),
				CommandArgument<coord_t>(In, "xmin", "X-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(In, "xmax", "X-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(In, "xstep", "Step between calculation points in x-direction."),
				CommandArgument<coord_t>(In, "ymin", "Y-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(In, "ymax", "Y-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(In, "ystep", "Step between calculation points in y-direction."),
				CommandArgument<coord_t>(In, "zmin", "Z-coordinate of the first calculation point in the reference image."),
				CommandArgument<coord_t>(In, "zmax", "Z-coordinate of the last calculation point in the reference image."),
				CommandArgument<coord_t>(In, "zstep", "Step between calculation points in z-direction."),
				CommandArgument<Vec3d>(In, "initial shift", "Initial shift between the images."),
				CommandArgument<string>(In, "file name prefix", "Prefix (and path) of files to write. The command will save point grid in the reference image, corresponding points in the deformed image, and goodness-of-fit. If the files exists, the current contents are erased."),
				CommandArgument<Vec3c>(In, "coarse comparison radius", "Radius of comparison region for coarse matching.", Vec3c(25, 25, 25)),
				CommandArgument<coord_t>(In, "coarse binning", "Amount of resolution reduction in coarse matching phase.", 2),
				CommandArgument<Vec3c>(In, "fine comparison radius", "Radius of comparison region for fine (full-resolution) matching.", Vec3c(10, 10, 10)),
				CommandArgument<coord_t>(In, "fine binning", "Amount of resolution reduction in fine matching phase. Set to same value than coarse binning to skip fine matching phase.", 2),
				CommandArgument<bool>(In, "normalize", "Indicates if the mean gray values of the two images should be made same in the overlapping region before matching.", true)
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
			if (!raw::internals::parseDimensions(refFile, refDimensions, refDT))
				throw ITLException("Unable to parse dimensions of reference image from file name.");

			Vec3c defDimensions;
			ImageDataType defDT;
			if (!raw::internals::parseDimensions(defFile, defDimensions, defDT))
				throw ITLException("Unable to parse dimensions of deformed image from file name.");

			if (refDT != defDT)
				throw ITLException("Data types of reference and deformed images must be the same.");

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

			if (refDT == UInt8)
			{
				blockMatchPartialLoad<uint8_t, uint8_t>(refFile, refDimensions, defFile, defDimensions, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, normalize, normFact);
			}
			else if (refDT == UInt16)
			{
				blockMatchPartialLoad<uint16_t, uint16_t>(refFile, refDimensions, defFile, defDimensions, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, normalize, normFact);
			}
			else if (refDT == Float32)
			{
				blockMatchPartialLoad<float32_t, float32_t>(refFile, refDimensions, defFile, defDimensions, refPoints, defPoints, fitGoodness, coarseCompRadius, coarseBinning, fineCompRadius, fineBinning, normalize, normFact);
			}
			else
				throw ParseException("Unsupported image data type (Please add the data type to BlockMatchPartialLoadCommand in commands.h file).");

			writeBlockMatchResult(fname, refPoints, defPoints, fitGoodness, normFact);
		}
	};


	template<typename pixel_t> class PullbackCommand : public Command
	{
	public:
		PullbackCommand() : Command("pullback", "Applies reverse of a deformation (calculated using blockmatch command) to image. In other words, performs pull-back operation. Makes output image the same size than input image.",
			{
				CommandArgument<Image<pixel_t> >(In, "image", "Image that will be pulled back."),
				CommandArgument<Image<pixel_t> >(Out, "pullback image", "Will store the result of the pullback operation."),
				CommandArgument<string>(In, "file name prefix", "File name prefix (and path) passed to blockmatch command."),
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

			reverseDeformation(deformed, pullback, refPoints, defPoints, CubicInterpolator<pixel_t, pixel_t, double, double>(Zero));
		}
	};


	class FilterDisplacementsCommand : public Command
	{
	public:
		FilterDisplacementsCommand() : Command("filterdisplacements", "Helper for non-rigid stitching script. Filters displacements calculated by blockmatch commands.",
			{
				CommandArgument<string>(In, "file name prefix", "Value passed as file name prefix argument to blockmatch."),
				CommandArgument<double>(In, "threshold", "Threshold value for filtering. Displacements whose some component differs more than this value from median filtered displacements are considered to be bad.", 3.0)
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
	//			CommandArgument<Image<pixel_t> >(In, "output image", "Output image."),
	//			//CommandArgument(In, ipt<uint8_t>(), "output weight image", "Output weight image."),
	//			CommandArgument<string>(In, "file list", "File name of index file that lists the files to be stitched."),
	//			CommandArgument<coord_t>(In, "x", "X-coordinate of the region of the stitched image that will be output."),
	//			CommandArgument<coord_t>(In, "y", "Y-coordinate of the region of the stitched image that will be output."),
	//			CommandArgument<coord_t>(In, "z", "Z-coordinate of the region of the stitched image that will be output."),
	//			CommandArgument<coord_t>(In, "width", "Width of the output region."),
	//			CommandArgument<coord_t>(In, "height", "Height of the output region."),
	//			CommandArgument<coord_t>(In, "depth", "Depth of the output region."),
	//			CommandArgument(In, Bool, "normalize", "Set to true to make mean gray value of images the same in the overlapping region.", "true")
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
				CommandArgument<string>(In, "transformation file", "Name of transformation file to process."),
				CommandArgument<Vec3c>(In, "image size", "Size of the source image."),
				CommandArgument<string>(In, "world to local prefix", "Prefix for output files."),
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			string transfFile = pop<string>(args);
			Vec3c imageSize = pop<Vec3c>(args);
			string prefix = pop<string>(args);

			determineWorldToLocal(transfFile, imageSize, prefix);
		}
	};

	template<typename pixel_t> class StitchVer2Command : public Command
	{
	public:
		StitchVer2Command() : Command("stitch_ver2", "Helper for non-rigid stitching script. Stitches subimages to one big image, given geometric transformation for each subimage. NOTE: The size of the output image does not need to be correct. Pass in image of size (1, 1, 1) to save memory during the process and to let the command allocate the image after processing.",
			{
				CommandArgument<Image<pixel_t> >(In, "output image", "Output image."),
				//CommandArgument(In, ipt<uint8_t>(), "output weight image", "Output weight image."),
				CommandArgument<string>(In, "file list", "File name of index file that lists the files to be stitched."),
				CommandArgument<coord_t>(In, "x", "X-coordinate of the region of the stitched image that will be output."),
				CommandArgument<coord_t>(In, "y", "Y-coordinate of the region of the stitched image that will be output."),
				CommandArgument<coord_t>(In, "z", "Z-coordinate of the region of the stitched image that will be output."),
				CommandArgument<coord_t>(In, "width", "Width of the output region."),
				CommandArgument<coord_t>(In, "height", "Height of the output region."),
				CommandArgument<coord_t>(In, "depth", "Depth of the output region."),
				CommandArgument<bool>(In, "normalize", "Set to true to make mean gray value of images the same in the overlapping region.", true)
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
				CommandArgument<Image<pixel_t> >(In, "input image", "Input image where background is marked with background value given by the third argument."),
				CommandArgument<Image<float32_t> >(Out, "output image", "Output image (distance map) where pixel value equals Euclidean distance to the nearest background pixel."),
				CommandArgument<double>(In, "background value", "Pixels belonging to the background are marked with this value in the input image.", 0)
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
				CommandArgument<Image<pixel_t> >(In, "input image", "Image whose histogram will be calculated."),
				CommandArgument<Image<float32_t> >(Out, "histogram", "Image that will contain the histogram on output."),
				CommandArgument<double>(In, "histogram minimum", "Minimum value to be included in the histogram.", 0),
				CommandArgument<double>(In, "histogram maximum", "The end of the last bin. This is the smallest value above minimum that is not included in the histogram.", 256),
				CommandArgument<coord_t>(In, "bin count", "Count of bins. The output image will be resized to contain this many pixels.", 256),
				CommandArgument<string>(In, "output file", "Name of file where the histogram data is to be saved in CSV format. Specify empty string to disable saving as text file.", "")
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
				CommandArgument<Vec3c>(In, "start point", "Starting point for the fill."),
				CommandArgument<double>(In, "fill value", "Fill color.")
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			Vec3c startPoint = pop<Vec3c>(args);
			pixel_t color = pixelRound<pixel_t>(pop<double>(args));


			floodfill(in, startPoint, color, color, AllNeighbours);
		}
	};









	






	template<typename pixel_t> class NormalizeZCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		NormalizeZCommand() : OneImageInPlaceCommand<pixel_t>("normalizez", "Makes sure that all z-slices of the image have the same mean value.",
			{
				CommandArgument<double>(In, "target mean", "Global mean that the image should have after normalization. Specify nothing or nan to retain global mean of the image.", numeric_limits<double>::signaling_NaN())
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			double globalMean = pop<double>(args);
			normalizeZ(in, (float32_t)globalMean);
		}
	};

	template<typename pixel_t> class RegionRemovalCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		RegionRemovalCommand() : OneImageInPlaceCommand<pixel_t>("regionremoval", "Removes nonzero regions smaller than given threshold.",
			{
				CommandArgument<coord_t>(In, "volume threshold", "All nonzero regions smaller than this threshold are removed.", 600)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			size_t volumeLimit = pop<size_t>(args);
			regionRemoval(in, volumeLimit);
		}
	};


	template<typename pixel_t> class RampCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		RampCommand() : OneImageInPlaceCommand<pixel_t>("ramp", "Fills image with ramp in given dimension, i.e. performs img[r] = r[dimension].",
			{
				CommandArgument<coord_t>(In, "dimension", "Dimension of the ramp.", 0)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			size_t dim = pop<size_t>(args);

			ramp(in, dim);
		}
	};
}

