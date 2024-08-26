#pragma once

#include <fstream>

#include "commandsbase.h"
#include "distributable.h"
#include "trivialdistributable.h"
#include "distributor.h"
#include "particleanalysis.h"
#include "regionremoval.h"
#include "csa.h"

#include "othercommands.h"
#include "pointprocesscommands.h"
#include "specialcommands.h"
#include "commandlist.h"

#include "pilibutilities.h"
#include "io/vectorio.h"

#include "standardhelp.h"

#include "pisystem.h"

namespace pilib
{

	inline std::string particleSeeAlso()
	{
		return "analyzeparticles, listanalyzers, headers, fillparticles, drawellipsoids, label, analyzelabels, regionremoval, greedycoloring, csa";
	}


	template<typename pixel_t> class LabelCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		LabelCommand() : OneImageInPlaceCommand<pixel_t>("label", "Labels distinct regions in the image with individual colors. The regions do not need to be separated by background. Execution is terminated in an error if the pixel data type does not support large enough values to label all the particles.",
			{
				CommandArgument<double>(ParameterDirection::In, "region color", "Current color of regions that should be labeled. Set to zero to label all non-zero regions.", 0),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the particles. ") + connectivityHelp(), Connectivity::NearestNeighbours),
			},
			particleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			pixel_t color = pixelRound<pixel_t>(pop<double>(args));
			Connectivity connectivity = pop<Connectivity>(args);

			labelParticles(in, color, (pixel_t)1, connectivity);
		}
	};



	class ListAnalyzersCommand : public TrivialDistributable
	{
	protected:
		friend class CommandList;

		ListAnalyzersCommand() : Command("listanalyzers", "Shows names of all available analysis methods that can be used in conjunction with `analyzeparticles` or `csa` commands. See also `headers` command.", {}, particleSeeAlso())
		{

		}

	private:

		static void print(const AnalyzerSet<Vec3sc, uint8_t>& analyzers)
		{
			for (size_t n = 0; n < analyzers.size(); n++)
			{
				std::cout << analyzers[n]->name() << std::endl;
				std::cout << analyzers[n]->description() << std::endl;
				if (n < analyzers.size() - 1)
					std::cout << std::endl;
			}
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			std::cout << "For generic particle analysis:" << std::endl;
			std::cout << "------------------------------" << std::endl;
			std::cout << std::endl;
			print(allAnalyzers<uint8_t>(Vec3c(1, 1, 1)));
			std::cout << std::endl;
			std::cout << "For fibre cross-section analysis:" << std::endl;
			std::cout << "--------------------------------" << std::endl;
			std::cout << std::endl;
			print(allCrossSectionAnalyzers<uint8_t>());
		}
	};

	class HeadersCommand : public TrivialDistributable
	{
	protected:
		friend class CommandList;

		HeadersCommand() : Command("headers", "Shows the column headers of a particle analysis result table.",
			{
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers that were used. Use the same value that was passed to `analyzeparticles` command. Separate the analyzer names with any non-alphanumeric character sequence."),
			},
			particleSeeAlso())
		{

		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			string titles = pop<string>(args);
			AnalyzerSet<Vec3sc, uint8_t> analyzers = createAnalyzers<uint8_t>(titles, Vec3c(1, 1, 1));
			//AnalyzerSet<Vec3sc, uint8_t> analyzers;
			//if (contains(titles, "2d"))
			//{
			//	// Probably cross-section analyzers
			//	analyzers = createCrossSectionAnalyzers<uint8_t>(titles);
			//}
			//else
			//{
			//	// Probably 3D analyzers
			//	analyzers = createAnalyzers<uint8_t>(titles, Vec3c(1, 1, 1));
			//}

			std::cout << analyzers.headers() << std::endl;
		}
	};

	class Headers2Command : public TrivialDistributable
	{
	protected:
		friend class CommandList;

		Headers2Command() : Command("headers", "Gets a comma-separated list of the column headers of particle analysis result table.",
			{
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers that were used. Use the same value that was passed to `analyzeparticles` command. Separate the analyzer names with any non-alphanumeric character sequence."),
				CommandArgument<string>(ParameterDirection::Out, "value", "The headers are placed into this string as a comma-separated list."),
			},
			particleSeeAlso())
		{

		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			string titles = pop<string>(args);
			string* value = std::get<string*>(args[1]);
			AnalyzerSet<Vec3sc, uint8_t> analyzers = createAnalyzers<uint8_t>(titles, Vec3c(1, 1, 1));
			std::stringstream s;
			s << analyzers.headers();
			*value = s.str();
		}
	};


	template<typename pixel_t> class AnalyzeParticlesBlockCommand : public Command, public Distributable
	{
	private:

		//static void writeList(const string& filename, const vector<coord_t>& v)
		//{
		//	std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
		//	if (!out)
		//		throw ITLException(string("Unable to write to: ") + filename);

		//	size_t s = v.size();
		//	out.write((const char*)&s, sizeof(size_t));

		//	for (size_t n = 0; n < v.size(); n++)
		//	{
		//		out.write((char*)&v[n], sizeof(coord_t));
		//	}
		//}

		//static void writeList(std::ofstream& out, const vector<Vec3sc>& v)
		//{
		//	size_t s = v.size();
		//	out.write((const char*)&s, sizeof(size_t));

		//	for (size_t n = 0; n < v.size(); n++)
		//	{
		//		out.write((char*)&v[n].x, sizeof(int32_t));
		//		out.write((char*)&v[n].y, sizeof(int32_t));
		//		out.write((char*)&v[n].z, sizeof(int32_t));
		//	}
		//}

		//static void writeList(const string& filename, const vector<Vec3sc>& v)
		//{
		//	std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
		//	if (!out)
		//		throw ITLException(string("Unable to write to: ") + filename);

		//	writeList(out, v);
		//}

		//static void writeList(const string& filename, const vector<vector<Vec3sc> >& v)
		//{
		//	std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
		//	if (!out)
		//		throw ITLException(string("Unable to write to: ") + filename);

		//	size_t s = v.size();
		//	out.write((const char*)&s, sizeof(size_t));

		//	for (size_t n = 0; n < v.size(); n++)
		//	{
		//		writeList(out, v[n]);
		//	}
		//}

		static void writeList(const string& filename, const vector<vector<Vec3sc> >& v)
		{
			itl2::writeListFile(filename, v, [=](std::ofstream& out, const std::vector<Vec3sc>& v) { itl2::writeList<Vec3sc>(out, v); });
		}
		
	protected:
		friend class CommandList;

		AnalyzeParticlesBlockCommand() : Command("analyzeparticlesblock", "This is an internal command used by the `analyzeparticles` command to analyze a block of the source image when distributed processing is enabled.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "img", "Image containing the particles. The pixels of the image will be set to various marker values; background pixels remain zero."),
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of analyzers to use. Separate the analyzer names with any non-alphanumeric character sequence."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the particles. ") + connectivityHelp(), Connectivity::NearestNeighbours),
				CommandArgument<size_t>(ParameterDirection::In, "volume limit", "Maximum size of particles to consider. Specify zero to consider all particles.", 0),
				CommandArgument<Vec3c>(ParameterDirection::In, "original dimensions", "Dimensions of the full image."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name template for files where the resulting data will be saved."),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing."),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current block in coordinates of the full image.")
			})
		{
		}

	public:
		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& img = *pop<Image<pixel_t>* >(args);
			string analyzerNames = pop<string>(args);
			Connectivity connectivity = pop<Connectivity>(args);
			size_t volumeLimit = pop<size_t>(args);
			Vec3c originalDimensions = pop<Vec3c>(args);
			string filename = pop<string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE index = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

			// Analyze particles but do not combine incomplete ones on the borders of the calculation blocks.
			// Combination is done afterwards after all nodes have finished processing.
			pixel_t fillColor = itl2::internals::SpecialColors<pixel_t>::fillColor();
			pixel_t largeColor = itl2::internals::SpecialColors<pixel_t>::largeColor();

			auto analyzers = createAnalyzers<pixel_t>(analyzerNames, originalDimensions);

			Results results;
			vector<vector<Vec3sc> > incompleteParticles;
			vector<vector<Vec3sc> > largeEdgePoints;
			vector<coord_t> edgeZ;

			results.headers() = analyzers.headers();

			itl2::internals::analyzeParticlesBlocks(img, analyzers, results, largeEdgePoints, incompleteParticles, connectivity, volumeLimit, fillColor, largeColor, Vec3sc(origin), edgeZ);

			// Write outputs to the output file
			filename += "_" + itl2::toString(index);

			writeText(filename + "_results.txt", results.str());
			writeList(filename + "_incomplete_particles.dat", incompleteParticles);
			writeList(filename + "_large_edge_points.dat", largeEdgePoints);
			itl2::writeListFile(filename + "_edge_z.dat", edgeZ);
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}

		virtual double calculateExtraMemory(const vector<ParamVariant>& args) const override
		{
			// Reserve some extra memory for the processing queues.
			return 0.75;
		}
	};

	template<typename pixel_t> class PrepareAnalyzeParticlesCommand : public InPlacePointProcess<pixel_t>
	{
	protected:
		friend class CommandList;

		PrepareAnalyzeParticlesCommand() : InPlacePointProcess<pixel_t>("prepareanalyzeparticles", "This is an internal command used to prepare image for particle analysis.",
			{
				CommandArgument<double>(ParameterDirection::In, "fill color", "Fill color."),
				CommandArgument<double>(ParameterDirection::In, "large color", "Large color."),
			})
		{
		}

	public:

		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			double fillColor = pop<double>(args);
			double largeColor = pop<double>(args);

			prepareParticleAnalysis(in, pixelRound<pixel_t>(fillColor), pixelRound<pixel_t>(largeColor));
		}

	};

	template<typename pixel_t> class AnalyzeParticlesCommand : public Command, public Distributable
	{
	private:

		//static void readList(const string& filename, vector<coord_t>& v)
		//{
		//	std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);
		//	if (!in)
		//		throw ITLException(string("Unable to open file: ") + filename);

		//	size_t s = 0;
		//	in.read((char*)&s, sizeof(size_t));

		//	v.reserve(v.size() + s);

		//	for (size_t n = 0; n < s; n++)
		//	{
		//		coord_t val;
		//		in.read((char*)&val, sizeof(coord_t));
		//		v.push_back(val);
		//	}
		//}

		//static void readList(std::ifstream& in, vector<Vec3sc>& v)
		//{
		//	size_t s = 0;
		//	in.read((char*)&s, sizeof(size_t));

		//	v.reserve(v.size() + s);

		//	for (size_t n = 0; n < s; n++)
		//	{
		//		Vec3sc val;
		//		in.read((char*)&val.x, sizeof(uint32_t));
		//		in.read((char*)&val.y, sizeof(uint32_t));
		//		in.read((char*)&val.z, sizeof(uint32_t));
		//		v.push_back(val);
		//	}
		//}

		//static void readList(const string& filename, vector<Vec3sc>& v)
		//{
		//	std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);
		//	if (!in)
		//		throw ITLException(string("Unable to open file: ") + filename);

		//	readList(in, v);
		//}

		//static void readList(const string& filename, vector<vector<Vec3sc> >& v)
		//{
		//	std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);
		//	if (!in)
		//		throw ITLException(string("Unable to open file: ") + filename);

		//	size_t s = 0;
		//	in.read((char*)&s, sizeof(size_t));

		//	v.reserve(v.size() + s);

		//	vector<Vec3sc> val;
		//	for (size_t n = 0; n < s; n++)
		//	{
		//		val.clear();
		//		readList(in, val);
		//		v.push_back(val);
		//	}
		//}

		static void readList(const string& filename, vector<vector<Vec3sc> >& v)
		{
			itl2::readListFile(filename, v, [=](std::ifstream& in, std::vector<Vec3sc>& v) { itl2::readList<Vec3sc>(in, v); });
		}

	protected:
		friend class CommandList;

		AnalyzeParticlesCommand() : Command("analyzeparticles", "Analyzes shape of blobs or other particles (separate nonzero regions) in the input image. Assumes all the particles have the same color. All the nonzero pixels in the input image will be set to same value. Output image will contain results of the measurements. There will be one row for each particle found in the input image. Use command `headers` to get interpretation of the columns. The order of the particles in the results may be different in normal and distributed processing modes. If you wish to analyze labeled particles, see `analyzelabels`.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::InOut, "input image", "Input image. The particles in this image will be filled with temporary color."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "results", "Image where analysis results are placed. This image will contain one row for each particle found in the input image. Use command `headers` to retrieve meanings of columns."),
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers to use. Use command `listanalyzers` to see all the names that can be specified. Separate the analyzer names with any non-alphanumeric character sequence.", "coordinates, volume"),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the particles. ") + connectivityHelp(), Connectivity::NearestNeighbours),
				CommandArgument<size_t>(ParameterDirection::In, "volume limit", "Maximum size of particles to consider, in pixels. Specify zero to consider all particles.", 0),
				CommandArgument<bool>(ParameterDirection::In, "single-threaded", "Set to true to use single-threaded processing. That might be faster for some cases. Does not have any effect if distributed processing is enabled.", false)
			},
			particleSeeAlso())
		{
		}

	public:

		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& out = *pop<Image<float32_t>* >(args);
			string analyzerNames = pop<string>(args);
			Connectivity connectivity = pop<Connectivity>(args);
			size_t volumeLimit = pop<size_t>(args);
			bool singleThreaded = pop<bool>(args);


			auto analyzers = createAnalyzers<pixel_t>(analyzerNames, in.dimensions());
			Results results;
			if(!singleThreaded)
				analyzeParticles(in, analyzers, results, connectivity, volumeLimit);
			else
				analyzeParticlesSingleThreaded(in, analyzers, results, connectivity, volumeLimit);

			// Convert results to output image.
			results.toImage(out);
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<float32_t>& out = *pop<DistributedImage<float32_t>* >(args);
			string analyzerNames = pop<string>(args);
			Connectivity connectivity = pop<Connectivity>(args);
			size_t volumeLimit = pop<size_t>(args);
			bool singleThreaded = pop<bool>(args);

			// This is not needed in this phase but it checks analyzer names, so run it before doing anything else.
			auto analyzers = createAnalyzers<pixel_t>(analyzerNames, in.dimensions());
			Results results;
			results.headers() = analyzers.headers();

			pixel_t fillColor = itl2::internals::SpecialColors<pixel_t>::fillColor();
			pixel_t largeColor = itl2::internals::SpecialColors<pixel_t>::largeColor();
			CommandList::get<PrepareAnalyzeParticlesCommand<pixel_t> >().runDistributed(distributor, { &in, (double)fillColor, (double)largeColor });


			// Create temp file path
			//unsigned int seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
			//std::mt19937 gen(seed);
			//string tempFilename = string("./tmp_images/particle_analysis_data_") + itl2::toString(gen());
			//fs::remove(tempFilename);
			string tempFilename = createTempFilename("particle_analysis_data");

			vector<ParamVariant> distributedArgs = { &in, analyzerNames, connectivity, volumeLimit, in.dimensions(), tempFilename, Distributor::BLOCK_INDEX_ARG_TYPE(), Distributor::BLOCK_ORIGIN_ARG_TYPE() };
			vector<string> output = CommandList::get<AnalyzeParticlesBlockCommand<pixel_t> >().runDistributed(distributor, distributedArgs);

			std::cout << "Loading data..." << std::endl;

			// Load data files and combine results (local processing, loads files one at a time)
			vector<vector<Vec3sc> > incompleteParticles;
			vector<vector<Vec3sc> > largeEdgePoints;
			vector<coord_t> edgeZ;
			for (size_t n = 0; n < output.size(); n++)
			{
				string resultsName = tempFilename + "_" + itl2::toString(n) + "_results.txt";
				string incompleteName = tempFilename + "_" + itl2::toString(n) + "_incomplete_particles.dat";
				string largeName = tempFilename + "_" + itl2::toString(n) + "_large_edge_points.dat";
				string edgezName = tempFilename + "_" + itl2::toString(n) + "_edge_z.dat";

				std::cout << "Reading results of job " << n << std::endl;

				results.readText(resultsName);
				readList(incompleteName, incompleteParticles);
				readList(largeName, largeEdgePoints);
				itl2::readListFile(edgezName, edgeZ);

				itl2::internals::combineParticleAnalysisResults(analyzers, results, largeEdgePoints, incompleteParticles, volumeLimit, connectivity, edgeZ, n < output.size() - 1);
			}

			// Load data files and combine results (local processing, loads all data files at once -> memory problem for large images)
			//vector<vector<Vec3sc> > incompleteParticles;
			//vector<vector<Vec3sc> > largeEdgePoints;
			//vector<coord_t> edgeZ;
			//for (size_t n = 0; n < output.size(); n++)
			//{
			//	string resultsName = tempFilename + "_" + itl2::toString(n) + "_results.txt";
			//	string incompleteName = tempFilename + "_" + itl2::toString(n) + "_incomplete_particles.dat";
			//	string largeName = tempFilename + "_" + itl2::toString(n) + "_large_edge_points.dat";
			//	string edgezName = tempFilename + "_" + itl2::toString(n) + "_edge_z.dat";

			//	std::cout << "Reading results of job " << n << std::endl;

			//	results.readText(resultsName);
			//	readList(incompleteName, incompleteParticles);
			//	readList(largeName, largeEdgePoints);
			//	readList(edgezName, edgeZ);
			//}
			//
			//itl2::internals::combineParticleAnalysisResults(analyzers, results, largeEdgePoints, incompleteParticles, volumeLimit, connectivity, edgeZ);
			
			// Do not remove the files until here so that if something goes wrong we can use them for debugging
			for (size_t n = 0; n < output.size(); n++)
			{
				string resultsName = tempFilename + "_" + itl2::toString(n) + "_results.txt";
				string incompleteName = tempFilename + "_" + itl2::toString(n) + "_incomplete_particles.dat";
				string largeName = tempFilename + "_" + itl2::toString(n) + "_large_edge_points.dat";
				string edgezName = tempFilename + "_" + itl2::toString(n) + "_edge_z.dat";
				
				fs::remove(resultsName);
				fs::remove(incompleteName);
				fs::remove(largeName);
				fs::remove(edgezName);
			}

			// Convert results to output image.
			Image<float32_t> outLocal;
			results.toImage(outLocal);
			out.setData(outLocal);

			return vector<string>();
		}
	};

	template<typename pixel_t> class FillParticlesCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		FillParticlesCommand() : OneImageInPlaceCommand<pixel_t>("fillparticles", "Fills particles that correspond to an entry in a list of particles with specified value. All other particles will be set to value 1. This command does not support cases where the particles have different colors.",
			{
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers that have been used to analyze the particles in the `analyzeparticles` command. The analyzers must contain 'coordinates' analyzer; otherwise this command does not know where the particles are."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "results", "Analysis results image."),
				CommandArgument<double>(ParameterDirection::In, "fill color", "Fill color."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the particles. ") + connectivityHelp(), Connectivity::NearestNeighbours),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Shift that is to be applied to the image before filling the particles. This argument is used internally in distributed processing.", Distributor::BLOCK_ORIGIN_ARG_TYPE())
			},
			particleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const override
		{
			string analyzerNames = pop<string>(args);
			Image<float32_t>& resultsImg = *pop<Image<float32_t>*>(args);
			double color = pop<double>(args);
			Connectivity connectivity = pop<Connectivity>(args);
			Vec3c blockPos = pop<Vec3c>(args);

			// This is used to make the image binary so that distributed processing is guaranteed to succeed.
			// There we assume that particles have color 1 and background is zero.
			// NOTE: To remove that assumption (and make the system work with multicolored particles) we need to do the whole
			// region growing process in some other way. The problem is that there might be two particles with different color
			// that touch each other and the region where they touch is located at the overlapping portion of the calculation blocks,
			// such that only one of the particles is visible in each block.
			// If we fill one of the particles beginning from the first block, we don't know in the growing phase whether we should
			// grow into the second particle or not!
			// By making sure that we support only binary particles removes this ambiguity, but it could be supported
			// in a more complex implementation.
			threshold(img, 0.0);

			auto analyzers = createAnalyzers<pixel_t>(analyzerNames, img.dimensions());
			Results results;
			results.fromImage(analyzers.headers(), resultsImg);
			fillParticles(img, results, pixelRound<pixel_t>(color), connectivity, -blockPos);
		}

		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 2)
			{
				DistributedImage<float32_t>& resultsImg = *std::get<DistributedImage<float32_t>*>(args[2]);

				// Always read the whole particle analysis results image.
				readStart = Vec3c(0, 0, 0);
				readSize = resultsImg.dimensions();
			}
		}

		using Distributable::runDistributed;

		virtual Vec3c getMargin(const vector<ParamVariant>& args) const override
		{
			return Vec3c(0, 0, 1);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& img = *std::get<DistributedImage<pixel_t>*>(args[0]);
			double fillColor = std::get<double>(args[3]);
			Connectivity connectivity = std::get<Connectivity>(args[4]);

			// NOTE: Here we cannot fill with any existing color, as the grow phase would erroneously grow
			// the existing (non-filled) regions, too.
			// The non-distributed version makes the image binary, so we first fill with color 2
			// and then replace 2 with fillColor.

			// Run this command
			//distributor.distribute(this, args, 2, Vec3c(0, 0, 1));
			vector<ParamVariant> tmpArgs = args;
			tmpArgs[3] = 2.0;
			distributor.distribute(this, tmpArgs);
			
			// Grow all regions whose color is fillColor towards regions whose color is 1
			auto& grow = CommandList::get<GrowCommand<pixel_t> >();
			grow.runDistributed(distributor, { &img, 2.0, 1.0, connectivity });

			// Now replace color 2 by fillColor
			auto& rc = CommandList::get<ReplaceCommand<pixel_t> >();
			rc.runDistributed(distributor, { &img, 2.0, fillColor});

			return vector<string>();
		}
	};



	template<typename pixel_t> class DrawEllipsoidsCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		DrawEllipsoidsCommand() : OneImageInPlaceCommand<pixel_t>("drawellipsoids",
			"Visualizes particles that have been analyzed with the `analyzeparticles` command by drawing "
			"a (scaled) principal component ellipsoid of each particle. The semi-axis directions of the "
			"ellipsoid are the principal directions of the particle, and semi-axis lengths are derived from the "
			"size of the particle in the principal directions as detailed in the 'ellipsoid type' parameter. "
			"If ellipsoid type parameter is set to 'BoundingSphere', the bounding sphere of the particle is drawn instead of ellipsoid.",
			{
				CommandArgument<string>(ParameterDirection::In, "analyzers",
				"List of names of analyzers that have been used to analyze the particles in the `analyzeparticles` command. "
				"The analyzers argument must contain the 'pca' analyzer if the ellipsoid type argument is set to any type of ellipsoid. "
				"Additionally, 'volume' analyzer is required if the 'Volume' ellipsoid type is selected, and 'boundingsphere' analyzer "
				"is required if the 'BoundingSphere' ellipsoid type is selected."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "results", "Analysis results image."),
				CommandArgument<double>(ParameterDirection::In, "fill color", "Fill color."),
				CommandArgument<string>(ParameterDirection::In, "ellipsoid type", "Type of ellipsoid to draw. 'Principal' denotes the principal axis ellipsoid without scaling, 'Bounding' results in an ellipsoid that covers all the points of the particle, and 'Volume' results in an ellipsoid whose volume equals the volume of the particle. 'BoundingSphere' denotes bounding sphere of the particle.", "principal"),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Shift that is to be applied to the image before filling the particles. This argument is used internally in distributed processing.", Distributor::BLOCK_ORIGIN_ARG_TYPE())
			},
			particleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const override
		{
			string analyzerNames = pop<string>(args);
			Image<float32_t>& resultsImg = *pop<Image<float32_t>*>(args);
			double color = pop<double>(args);
			string types = pop<string>(args);
			Vec3c blockPos = pop<Vec3c>(args);

			EllipsoidType type = fromString<EllipsoidType>(types);

			auto analyzers = createAnalyzers<pixel_t>(analyzerNames, img.dimensions());
			Results results;
			results.fromImage(analyzers.headers(), resultsImg);

			if (blockPos != Vec3c(0, 0, 0))
			{
				// Adjust ellipsoid locations to account for the new origin.
				size_t cxi = results.getColumnIndex("CX [pixel]");
				size_t cyi = results.getColumnIndex("CY [pixel]");
				size_t czi = results.getColumnIndex("CZ [pixel]");
				for (size_t n = 0; n < results.size(); n++)
				{
					results[n][cxi] -= blockPos.x;
					results[n][cyi] -= blockPos.y;
					results[n][czi] -= blockPos.z;
				}
			}

			drawEllipsoids(img, results, pixelRound<pixel_t>(color), type);
		}

		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 2)
			{
				DistributedImage<float32_t>& resultsImg = *std::get<DistributedImage<float32_t>*>(args[2]);

				// Always read the whole particle analysis results image.
				readStart = Vec3c(0, 0, 0);
				readSize = resultsImg.dimensions();
			}
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			return 1;
		}
	};




	template<typename pixel_t> class RegionRemovalCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		RegionRemovalCommand() : OneImageInPlaceCommand<pixel_t>("regionremoval", "Removes nonzero foreground regions smaller than given threshold. This command supports only binary images (where regions are separated by backround).",
			{
				CommandArgument<size_t>(ParameterDirection::In, "volume threshold", "All nonzero regions consisting of less than this many pixels are removed.", 600),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the particles. ") + connectivityHelp(), Connectivity::NearestNeighbours),
				CommandArgument<bool>(ParameterDirection::In, "allow multi-threading", "Set to true to allow multi-threaded processing. Set to false to use single-threaded processing. Single-threaded processing is often faster if it is known in advance that there are only a few particles or if the image is small. This argument has no effect in the distributed processing mode. There, the processing is always multi-threaded.", true),
			},
			"openingfilter, closingfilter, analyzeparticles")
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			size_t volumeLimit = pop<size_t>(args);
			Connectivity connectivity = pop<Connectivity>(args);
			bool multiThreaded = pop<bool>(args);

			// TODO: Add preserveEdges if that is required
			regionRemoval(in, volumeLimit, false, connectivity, multiThreaded);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>* >(args);
			size_t volumeLimit = pop<size_t>(args);
			Connectivity connectivity = pop<Connectivity>(args);
			bool multiThreaded = pop<bool>(args);
			//coord_t threadCount = (coord_t)pop<size_t>(args);

			// TODO: Add preserveEdges if that is required

			// Distributed analyze particles
			const string analyzers = "volume coordinates";


			// Create (temporary) image for analysis result
			// NOTE: This is image name, not file name.
			//string tempName = "regionremoval_particle_analysis_results_" + itl2::toString(randc(10000));
			//string datatype = "float32";
			//CommandList::get<NewImageCommand>().runDistributed(distributor, { tempName, datatype, (coord_t)1, (coord_t)1, (coord_t)1 });
			//DistributedImage<float32_t>& results = *(DistributedImage<float32_t>*)distributor.getSystem()->getDistributedImage(tempName);
			DistributedTempImage<float32_t> resultsImg(distributor, "regionremoval_particle_analysis_results", 1, DistributedImageStorageType::Raw);
			DistributedImage<float32_t>& results = resultsImg.get();

			// Find small particles using particle analysis
			//input image pixel_t,
			//results image float32_t
			//analyzers string
			//connectivity Connectivity
			//volume limit coord_t
			CommandList::get<AnalyzeParticlesCommand<pixel_t> >().runDistributed(distributor, { &img, &results, analyzers, connectivity, volumeLimit, false });

			// Fill small particles using distributed fillParticles

			// Parameters
			//img image pixel_t
			//analyzers string
			//results image float32_t
			//fill color double
			//connectivity Connectivity
			//placeholder Distributor::BLOCK_ORIGIN_ARG_TYPE
			CommandList::get<FillParticlesCommand<pixel_t> >().runDistributed(distributor, { &img, analyzers, &results, 0.0, connectivity, Distributor::BLOCK_ORIGIN_ARG_TYPE() });

			// Clear temp image
			//CommandList::get<ClearCommand>().runDistributed(distributor, { tempName });


			return vector<string>();
		}
	};

	template<typename pixel_t> class AnalyzeLabelsCommand : public Command
	{
	protected:
		friend class CommandList;

		AnalyzeLabelsCommand() : Command("analyzelabels", "Analyzes labeled regions of the input image. The regions do not need to be connected. Region having value zero is skipped.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "image", "The input image where each particle is labeled with different color."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "results", "Analysis results image."),
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers to use. Use command `listanalyzers` to see all the names that can be specified. Separate the analyzer names with any non-alphanumeric character sequence.", "coordinates, volume"),
			},
			particleSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& out = *pop<Image<float32_t>* >(args);
			string analyzerNames = pop<string>(args);
			
			Results results;
			analyzeLabels(in, analyzerNames, results);

			results.toImage(out);
		}
	};


	template<typename pixel_t> class GreedyColoringCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		GreedyColoringCommand() : OneImageInPlaceCommand<pixel_t>("greedycoloring", "Perform greedy coloring of regions. Colors each region in image such that its neighbours are colored with different colors, and uses as little colors as possible. Uses greedy algorithm so the count of colors used might not be minimal. Assumes background to have value 0 and colors all non-zero regions.",
			{
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the regions. ") + connectivityHelp(), Connectivity::AllNeighbours),
			},
			particleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const override
		{
			Connectivity connectivity = pop<Connectivity>(args);

			greedyColoring(img, connectivity);
		}

	};


	namespace internals
	{
		template<typename pixel_t> AnalyzerSet<Vec3sc, pixel_t> csaAnalyzers()
		{
			return createCrossSectionAnalyzers<pixel_t>("coordinates2d, volume, pca2d, convexhull2d, bounds2d");
		}
	}


	class CSAHeaders : public TrivialDistributable
	{
	protected:
		friend class CommandList;

		CSAHeaders() : Command("csaheaders", "Retrieves information on the meaning of the columns in the result image of `csa` command.", {}, "csa")
		{

		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			auto analyzers = internals::csaAnalyzers<uint8_t>();
			std::cout << analyzers.headers();

			// Additional columns that do not come from the analyzers.
			std::cout << ", X3 [pixel], Y3 [pixel], Z3 [pixel]";
			std::cout << ", CX3 [pixel], CY3 [pixel], CZ3 [pixel]";
			std::cout << ", length min, length max, length mean, length mode";

			std::cout << std::endl;

			std::cout << "The (X3, Y3, Z3) columns give the position of the center of the slice image in the coordinates of the original image." << std::endl;
			std::cout << "The (CX3, CY3, CZ3) columns give the position of the centroid of the cross-section in the coordinates of the original image. These columns are available only if CX and CY columns are available." << std::endl;
			std::cout << "The length columns are available only if length input image is given to the csa command." << std::endl;
		}
	};


	template<typename pixel_t> class CSACommand : public Command
	{
	protected:
		friend class CommandList;

		CSACommand() : Command("csa", "Analyzes cross-sections of cylindrical, tubular, or fibre-like objects. Use to e.g. measure cross-sectional area of fibres. Requires that the local orientation of the fibres has been determined using, e.g., `cylinderorientation` command.",
			{
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "original", "Binary image containing the fibres as foreground."),
				CommandArgument<Image<float32_t>>(ParameterDirection::In, "energy", "Image corresponding to the energy output image of cylinderorientation command. Non-zero energy signifies that local orientation is available at that location."),
				CommandArgument<Image<float32_t>>(ParameterDirection::In, "phi", R"(The azimuthal angle of the local fibre orientation direction. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.)"),
				CommandArgument<Image<float32_t>>(ParameterDirection::In, "theta", R"(The polar angle of the local fibre orientation direction. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.)"),
				CommandArgument<Image<float32_t>>(ParameterDirection::Out, "results", "Image where analysis results are placed. This image will contain one row for each fibre cross-section analyzed. Use command `csaheaders` to retrieve meanings of columns."),
				//CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers to use. Use command `listanalyzers` to see all the names that can be specified. Separate the analyzer names with any non-alphanumeric character sequence.", "coordinates2d, volume, pca2d, convexhull2d, bounds2d"),
				CommandArgument<size_t>(ParameterDirection::In, "slice radius", "Half width of one cross-sectional slice", 40),
				CommandArgument<size_t>(ParameterDirection::In, "slice count", "Count of cross-sectional slices to analyze", 1000),
				CommandArgument<size_t>(ParameterDirection::In, "random seed", "Seed for random number generator.", 123),
				CommandArgument<Image<pixel_t>>(ParameterDirection::Out, "slices", "The extracted slices are placed into this image."),
				CommandArgument<Image<pixel_t>>(ParameterDirection::Out, "visualization", "The locations where the slices are extracted from are drawn into this image"),
				CommandArgument<Image<float32_t>>(ParameterDirection::In, "length", "Image containing local fibre length."),
				CommandArgument<Image<float32_t>>(ParameterDirection::Out, "length slices", "The extracted slices from the length image are placed into this image."),
			},
			particleSeeAlso() + ", cylinderorientation")
		{

		}
	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>*>(args);
			Image<float32_t>& energy = *pop<Image<float32_t>*>(args);
			Image<float32_t>& phi = *pop<Image<float32_t>*>(args);
			Image<float32_t>& theta = *pop<Image<float32_t>*>(args);
			Image<float32_t>& out = *pop<Image<float32_t>*>(args);
			//string analyzerNames = pop<string>(args);
			size_t sliceRadius = pop<size_t>(args);
			size_t sliceCount = pop<size_t>(args);
			size_t randSeed = pop<size_t>(args);

			Image<pixel_t>& slices = *pop<Image<pixel_t>*>(args);
			Image<pixel_t>& vis = *pop<Image<pixel_t>*>(args);
			Image<float32_t>& length = *pop<Image<float32_t>*>(args);
			Image<float32_t>& lengthSlices = *pop<Image<float32_t>*>(args);

			auto analyzers = internals::csaAnalyzers<pixel_t>();

			Results results;
			
			csa<pixel_t>(in, energy, phi, theta, &length, analyzers, results, sliceRadius, sliceCount, randSeed, &slices, &lengthSlices, &vis);

			results.toImage(out);
		}
	};

	template<typename pixel_t> class CSA2Command : public Command
	{
	protected:
		friend class CommandList;

		CSA2Command() : Command("csa", "Analyzes cross-sections of cylindrical, tubular, or fibre-like objects. Use to e.g. measure cross-sectional area of fibres. Requires that the local orientation of the fibres has been determined using, e.g., `cylinderorientation` command.",
			{
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "original", "Binary image containing the fibres as foreground."),
				CommandArgument<Image<float32_t>>(ParameterDirection::In, "energy", "Image corresponding to the energy output image of cylinderorientation command. Non-zero energy signifies that local orientation is available at that location."),
				CommandArgument<Image<float32_t>>(ParameterDirection::In, "phi", R"(The azimuthal angle of the local fibre orientation direction. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.)"),
				CommandArgument<Image<float32_t>>(ParameterDirection::In, "theta", R"(The polar angle of the local fibre orientation direction. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.)"),
				CommandArgument<Image<float32_t>>(ParameterDirection::Out, "results", "Image where analysis results are placed. This image will contain one row for each fibre cross-section analyzed. Use command `csaheaders` to retrieve meanings of columns."),
				//CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers to use. Use command `listanalyzers` to see all the names that can be specified. Separate the analyzer names with any non-alphanumeric character sequence.", "coordinates2d, volume, pca2d, convexhull2d, bounds2d"),
				CommandArgument<size_t>(ParameterDirection::In, "slice radius", "Half width of one cross-sectional slice", 40),
				CommandArgument<size_t>(ParameterDirection::In, "slice count", "Count of cross-sectional slices to analyze", 1000),
				CommandArgument<size_t>(ParameterDirection::In, "random seed", "Seed for random number generator.", 123),
				CommandArgument<Image<pixel_t>>(ParameterDirection::Out, "slices", "The extracted slices are placed into this image."),
				CommandArgument<Image<pixel_t>>(ParameterDirection::Out, "visualization", "The locations where the slices are extracted from are drawn into this image"),
			},
			particleSeeAlso() + ", cylinderorientation")
		{

		}
	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>*>(args);
			Image<float32_t>& energy = *pop<Image<float32_t>*>(args);
			Image<float32_t>& phi = *pop<Image<float32_t>*>(args);
			Image<float32_t>& theta = *pop<Image<float32_t>*>(args);
			Image<float32_t>& out = *pop<Image<float32_t>*>(args);
			//string analyzerNames = pop<string>(args);
			size_t sliceRadius = pop<size_t>(args);
			size_t sliceCount = pop<size_t>(args);
			size_t randSeed = pop<size_t>(args);

			Image<pixel_t>& slices = *pop<Image<pixel_t>*>(args);
			Image<pixel_t>& vis = *pop<Image<pixel_t>*>(args);

			auto analyzers = internals::csaAnalyzers<pixel_t>();

			Results results;

			csa<pixel_t>(in, energy, phi, theta, nullptr, analyzers, results, sliceRadius, sliceCount, randSeed, &slices, nullptr, &vis);

			results.toImage(out);
		}
	};
}
