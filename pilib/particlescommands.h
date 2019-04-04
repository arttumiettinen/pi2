#pragma once

#include <fstream>

#include "commandsbase.h"
#include "distributable.h"
#include "trivialdistributable.h"
#include "distributor.h"
#include "particleanalysis.h"
#include "regionremoval.h"

#include "othercommands.h"
#include "pointprocesscommands.h"
#include "specialcommands.h"
//#include "iocommands.h"

#include "pisystem.h"

namespace pilib
{

	class ListAnalyzersCommand : public TrivialDistributable
	{
	public:
		ListAnalyzersCommand() : Command("listanalyzers", "Shows names of all available analysis methods that can be used in conjunction with particleanalysis command. See also headers command.")
		{

		}

		virtual void run(vector<ParamVariant>& args) const
		{
			auto analyzers = allAnalyzers<uint8_t>(Vec3c(1, 1, 1));
			for(size_t n = 0; n < analyzers.size(); n++)
			{
				cout << analyzers[n]->name() << endl;
				cout << analyzers[n]->description() << endl;
				if (n < analyzers.size() - 1)
					cout << endl;
			}
		}
	};

	class HeadersCommand : public TrivialDistributable
	{
	public:
		HeadersCommand() : Command("headers", "Shows headers of particle analysis result table.",
			{
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers that were used. Use the same value that was passed to analyzeparticles command. Separate the analyzer names with any non-alphanumeric character sequence."),
			})
		{

		}

		virtual void run(vector<ParamVariant>& args) const
		{
			string titles = pop<string>(args);

			auto analyzers = createAnalyzers<uint8_t>(titles, Vec3c(1, 1, 1));

			cout << analyzers.headers() << endl;
		}
	};


	template<typename pixel_t> class AnalyzeParticlesBlockCommand : public Command, public Distributable
	{
	private:

		static void writeList(const string& filename, const vector<coord_t>& v)
		{
			std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
			if (!out)
				throw ITLException(string("Unable to write to: ") + filename);

			size_t s = v.size();
			out.write((const char*)&s, sizeof(size_t));

			for (size_t n = 0; n < v.size(); n++)
			{
				out.write((char*)&v[n], sizeof(coord_t));
			}
		}

		static void writeList(std::ofstream& out, const vector<Vec3sc>& v)
		{
			size_t s = v.size();
			out.write((const char*)&s, sizeof(size_t));

			for (size_t n = 0; n < v.size(); n++)
			{
				out.write((char*)&v[n].x, sizeof(int32_t));
				out.write((char*)&v[n].y, sizeof(int32_t));
				out.write((char*)&v[n].z, sizeof(int32_t));
			}
		}

		static void writeList(const string& filename, const vector<Vec3sc>& v)
		{
			std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
			if (!out)
				throw ITLException(string("Unable to write to: ") + filename);

			writeList(out, v);
		}

		static void writeList(const string& filename, const vector<vector<Vec3sc> >& v)
		{
			std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
			if (!out)
				throw ITLException(string("Unable to write to: ") + filename);

			size_t s = v.size();
			out.write((const char*)&s, sizeof(size_t));

			for (size_t n = 0; n < v.size(); n++)
			{
				writeList(out, v[n]);
			}
		}
		
	public:
		AnalyzeParticlesBlockCommand() : Command("analyzeparticlesblock", "This is an internal command used by the analyzeparticles command to analyze a block of the source image when distributed processing is enabled.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "img", "Image containing the particles. The pixels of the image will be set to various marker values; background pixels remain zero."),
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of analyzers to use. Separate the analyzer names with any non-alphanumeric character sequence."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", "Connectivity of the particles.", Connectivity::NearestNeighbours),
				CommandArgument<size_t>(ParameterDirection::In, "volume limit", "Maximum size of particles to consider. Specify zero to consider all particles.", 0),
				CommandArgument<Vec3c>(ParameterDirection::In, "original dimensions", "Dimensions of the full image."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name template for files where the resulting data will be saved."),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing."),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current block in coordinates of the full image.")
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
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
			writeList(filename + "_edge_z.dat", edgeZ);
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			return distributor.distribute(this, args, 2, Vec3c(0, 0, 0));
		}
	};

	template<typename pixel_t> class AnalyzeParticlesCommand : public TwoImageInputOutputCommand<pixel_t, float32_t>, public Distributable
	{
	private:

		static void readList(const string& filename, vector<coord_t>& v)
		{
			std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);
			if (!in)
				throw ITLException(string("Unable to open file: ") + filename);

			size_t s = 0;
			in.read((char*)&s, sizeof(size_t));

			v.reserve(v.size() + s);

			for (size_t n = 0; n < s; n++)
			{
				coord_t val;
				in.read((char*)&val, sizeof(coord_t));
				v.push_back(val);
			}
		}

		static void readList(ifstream& in, vector<Vec3sc>& v)
		{
			size_t s = 0;
			in.read((char*)&s, sizeof(size_t));

			v.reserve(v.size() + s);

			for (size_t n = 0; n < s; n++)
			{
				Vec3sc val;
				in.read((char*)&val.x, sizeof(uint32_t));
				in.read((char*)&val.y, sizeof(uint32_t));
				in.read((char*)&val.z, sizeof(uint32_t));
				v.push_back(val);
			}
		}

		static void readList(const string& filename, vector<Vec3sc>& v)
		{
			std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);
			if (!in)
				throw ITLException(string("Unable to open file: ") + filename);

			readList(in, v);
		}

		static void readList(const string& filename, vector<vector<Vec3sc> >& v)
		{
			std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);
			if (!in)
				throw ITLException(string("Unable to open file: ") + filename);

			size_t s = 0;
			in.read((char*)&s, sizeof(size_t));

			v.reserve(v.size() + s);

			vector<Vec3sc> val;
			for (size_t n = 0; n < s; n++)
			{
				val.clear();
				readList(in, val);
				v.push_back(val);
			}
		}

	public:
		AnalyzeParticlesCommand() : TwoImageInputOutputCommand<pixel_t, float32_t>("analyzeparticles", "Analyzes shape of blobs or other particles (separate nonzero regions) in the input image. All the nonzero pixels in the input image will be set to same value. Output image will contain results of the measurements. There will be one row for each particle found in the input image. Use command headers to get interpretation of the columns.",
			{
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers to use. Use command listanalyzers to see all the names that can be specified. Separate the analyzer names with any non-alphanumeric character sequence.", "coordinates, volume"),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", "Connectivity of the particles.", Connectivity::NearestNeighbours),
				CommandArgument<size_t>(ParameterDirection::In, "volume limit", "Maximum size of particles to consider. Specify zero to consider all particles.", 0),
			})
		{
		}

		virtual void run(Image<pixel_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const
		{
			string analyzerNames = pop<string>(args);
			Connectivity connectivity = pop<Connectivity>(args);
			size_t volumeLimit = pop<size_t>(args);


			auto analyzers = createAnalyzers<pixel_t>(analyzerNames, in.dimensions());
			Results results;
			analyzeParticles(in, analyzers, results, connectivity, volumeLimit);

			// Convert results to output image.
			results.toImage(out);
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<float32_t>& out = *pop<DistributedImage<float32_t>* >(args);
			string analyzerNames = pop<string>(args);
			Connectivity connectivity = pop<Connectivity>(args);
			size_t volumeLimit = pop<size_t>(args);

			// This is not needed in this phase but it checks analyzer names, so run it before doing anything else.
			auto analyzers = createAnalyzers<pixel_t>(analyzerNames, in.dimensions());
			Results results;
			results.headers() = analyzers.headers();

			// Create temp file path
			unsigned int seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
			std::mt19937 gen(seed);
			string tempFilename = string("./tmp_images/particle_analysis_data_") + itl2::toString(gen());
			fs::remove(tempFilename);

			AnalyzeParticlesBlockCommand<pixel_t> cmd;
			vector<string> output = cmd.runDistributed(distributor, { &in, analyzerNames, connectivity, volumeLimit, in.dimensions(), tempFilename, Distributor::BLOCK_INDEX_ARG_TYPE(), Distributor::BLOCK_ORIGIN_ARG_TYPE() });

			// Load data files and combine results (local processing)
			cout << "Loading data..." << endl;

			vector<vector<Vec3sc> > incompleteParticles;
			vector<vector<Vec3sc> > largeEdgePoints;
			vector<coord_t> edgeZ;
			for (size_t n = 0; n < output.size(); n++)
			{
				string resultsName = tempFilename + "_" + itl2::toString(n) + "_results.txt";
				string incompleteName = tempFilename + "_" + itl2::toString(n) + "_incomplete_particles.dat";
				string largeName = tempFilename + "_" + itl2::toString(n) + "_large_edge_points.dat";
				string edgezName = tempFilename + "_" + itl2::toString(n) + "_edge_z.dat";

				cout << "Reading results of job " << n << endl;

				results.readText(resultsName);
				readList(incompleteName, incompleteParticles);
				readList(largeName, largeEdgePoints);
				readList(edgezName, edgeZ);
			}
			
			
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

			itl2::internals::combineParticleAnalysisResults(analyzers, results, largeEdgePoints, incompleteParticles, volumeLimit, connectivity, edgeZ);

			// Convert results to output image.
			Image<float32_t> outLocal;
			results.toImage(outLocal);
			out.setData(outLocal);

			return vector<string>();
		}
	};

	template<typename pixel_t> class FillParticlesCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	public:
		FillParticlesCommand() : OneImageInPlaceCommand<pixel_t>("fillparticles", "Fills particles that correspond to an entry in a list of particles. This command is not guaranteed to work properly if all the particles do not have the same color, i.e., if the input image is not a binary image.",
			{
				CommandArgument<string>(ParameterDirection::In, "analyzers", "List of names of analyzers that have been used to analyze the particles in analyzeparticles command. The analyzers must contain 'coordinates' analyzer; otherwise this command does not know where the particles are."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "results", "Analysis results image."),
				CommandArgument<double>(ParameterDirection::In, "fill color", "Fill color."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", "Connectivity of the particles.", Connectivity::NearestNeighbours),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Shift that is to be applied to the image before filling the particles. This argument is used internally in distributed processing.", Distributor::BLOCK_ORIGIN_ARG_TYPE())
			})
		{
		}

		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const
		{
			string analyzerNames = pop<string>(args);
			Image<float32_t>& resultsImg = *pop<Image<float32_t>*>(args);
			double color = pop<double>(args);
			Connectivity connectivity = pop<Connectivity>(args);
			Vec3c blockPos = pop<Vec3c>(args);

			// This is used to make the image binary so that distributed processing is guaranteed to succeed.
			// There we assume that particles have color 1 and background is zero.
			// NOTE: To remove that assumption (and make the system work with multicolored particles) we need to do the whole
			// region growing process in some other way. The problem is that there might be two particles that touch each other
			// and the region where they touch is located at the overlapping portion of the calculation blocks.
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

		virtual void getCorrespondingBlock(vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 2)
			{
				DistributedImage<float32_t>& resultsImg = *get<DistributedImage<float32_t>*>(args[2]);

				// Always read the whole particle analysis results image.
				readStart = Vec3c(0, 0, 0);
				readSize = resultsImg.dimensions();
			}
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& img = *get<DistributedImage<pixel_t>*>(args[0]);
			double fillColor = get<double>(args[3]);

			// NOTE: Here we cannot fill with any existing color, as the grow phase would erroneously grow
			// the existing (non-filled) regions, too.
			// The non-distributed version makes the image binary, so we first fill with color 2
			// and then replace 2 with fillColor.

			// Run this command
			//distributor.distribute(this, args, 2, Vec3c(0, 0, 1));
			vector<ParamVariant> tmpArgs = args;
			tmpArgs[3] = 2.0;
			distributor.distribute(this, tmpArgs, 2, Vec3c(0, 0, 1));
			
			// Grow all regions whose color is fillColor towards regions whose color is 1
			GrowCommand<pixel_t> grow;
			grow.runDistributed(distributor, { &img, 2.0, 1.0 });

			// Now replace color 2 by fillColor
			ReplaceCommand<pixel_t> rc;
			rc.runDistributed(distributor, { &img, 2.0, fillColor});

			return vector<string>();
		}
	};

	template<typename pixel_t> class RegionRemovalCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	public:
		RegionRemovalCommand() : OneImageInPlaceCommand<pixel_t>("regionremoval", "Removes nonzero foreground regions smaller than given threshold.",
			{
				CommandArgument<size_t>(ParameterDirection::In, "volume threshold", "All nonzero regions smaller than this threshold are removed.", 600),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", "Connectivity of the particles.", Connectivity::NearestNeighbours)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			size_t volumeLimit = pop<size_t>(args);
			Connectivity connectivity = pop<Connectivity>(args);

			// TODO: Add preserveEdges if that is required
			regionRemoval(in, volumeLimit, false, connectivity);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>* >(args);
			size_t volumeLimit = pop<size_t>(args);
			Connectivity connectivity = pop<Connectivity>(args);

			// TODO: Add preserveEdges if that is required

			// Distributed analyze particles
			const string analyzers = "volume coordinates";


			// Create (temporary) image for analysis results
			unsigned int seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
			std::mt19937 gen(seed);
			string tempName = string("regionremoval_particle_analysis") + itl2::toString(gen());
			fs::remove(tempName);

			NewImageCommand ni;
			string datatype = "float32";
			ni.runDistributed(distributor, { tempName, datatype, (coord_t)1, (coord_t)1, (coord_t)1 });

			DistributedImage<float32_t>& results = *(DistributedImage<float32_t>*)distributor.getSystem()->getDistributedImage(tempName);

			// Find small particles using particle analysis
			//input image pixel_t,
			//results image float32_t
			//analyzers string
			//connectivity Connectivity
			//volume limit coord_t
			AnalyzeParticlesCommand<pixel_t> ap;
			ap.runDistributed(distributor, { &img, &results, analyzers, connectivity, volumeLimit });


			//WriteRawCommand<float32_t> wr;
			//string tmp = "rr_pa";
			//wr.runDistributed(distributor, { &results, tmp});


			// Fill small particles using distributed fillParticles

			// Parameters
			//img image pixel_t
			//analyzers string
			//results image float32_t
			//fill color double
			//connectivity Connectivity
			//placeholder Distributor::BLOCK_ORIGIN_ARG_TYPE
			FillParticlesCommand<pixel_t> fp;
			fp.runDistributed(distributor, { &img, analyzers, &results, 0.0, connectivity, Distributor::BLOCK_ORIGIN_ARG_TYPE() });

			// Clear temp image
			ClearCommand cl;
			cl.runDistributed(distributor, { tempName });


			return vector<string>();
		}
	};
}
