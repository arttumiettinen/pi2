#pragma once

#include "commandsbase.h"
#include "overlapdistributable.h"

namespace pilib
{
	template<typename pixel_t> class HybridThinCommand : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		HybridThinCommand() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("hybridthin", "Thins one layer of pixels from the foreground of the image. Nonzero pixels are assumed to belong to the foreground. Run iteratively to calculate a hybrid skeleton.")
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			size_t changed = hybridThin(in);
			cout << changed << " pixels removed." << endl;
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			return Vec3c(10, 10, 10);
		}
	};

	template<typename pixel_t> class LineThinCommand : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		LineThinCommand() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("linethin", "Thins one layer of pixels from the foreground of the image. Nonzero pixels are assumed to belong to the foreground. Run iteratively to calculate a line skeleton.")
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			size_t changed = lineThin(in);
			cout << changed << " pixels removed." << endl;
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			return Vec3c(10, 10, 10);
		}
	};




	template<typename command_t, typename base_t> class IterableDistributable : public base_t, public Distributable
	{
	public:
		IterableDistributable(const string& name, const string& help, const vector<CommandArgumentBase>& extraArgs = {}) : base_t(name, help, extraArgs)
		{
		}

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			size_t lastTotalChanged = 0;
			size_t n = 0;
			while (true)
			{
				cout << "Iteration " << n << endl;

				// Run one iteration of thinning
				command_t cmd;
				vector<string> output = distributor.distribute(&cmd, args, 2, Vec3c(10, 10, 10));

				// Calculate total number of changed pixels
				size_t totalChanged = 0;
				for(size_t n = 0; n < output.size(); n++)
				{
    				string out = output[n];
    				
					size_t startPos = out.find("pixels removed.");
					if (startPos != string::npos)
					{
						size_t lineStart = out.rfind('\n', startPos);
						if (lineStart == string::npos)
							lineStart = 0;
						string number = out.substr(lineStart + 1, startPos - lineStart - 1);
						size_t val = fromString<size_t>(number);
						totalChanged += val;
					}
					else
					{
						throw ITLException(string("Invalid output from distributed thin command (job ") + itl2::toString(n) + string(", removed pixel count not found):\n") + out);
					}
				}

				cout << totalChanged << " pixels removed." << endl;

				if (totalChanged == lastTotalChanged)
				{
					// We have possibly reached the end of the iteration.
					// TODO: Subtract the current image from the previous one to make sure that everything's done;
					// or
					// adjust hybridThin command somehow so that it does not account for pixels in the overlapping regions when
					// calculating count of changed pixels;
					// or
					// make distributor.distribute(...) optionally insert a command that counts nonzero pixels in the non-overlapping region -
					// that number could be used in iteration end condition instead of number of changed pixels.
					break;
				}

				lastTotalChanged = totalChanged;
				n++;
			}

		}
	};

	template<typename pixel_t> class HybridSkeletonCommand : public IterableDistributable<HybridThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >
	{
	public:
		HybridSkeletonCommand() : IterableDistributable<HybridThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >("hybridskeleton", "Calculates skeleton of the foreground of the given image. Nonzero pixels are assumed to belong to the foreground. The skeleton contains both lines and plates.")
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			hybridSkeleton(in);
		}
	};

	template<typename pixel_t> class LineSkeletonCommand : public IterableDistributable<HybridThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >
	{
	public:
		LineSkeletonCommand() : IterableDistributable<HybridThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >("lineskeleton", "Calculates skeleton of the foreground of the given image. Nonzero pixels are assumed to belong to the foreground. The skeleton contains only lines (no plates). Note that if a line skeleton is required, it might be better idea to fill all holes in the structure and use hybridskeleton command as that seems to produce cleaner skeletons.")
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			lineSkeleton(in);
		}
	};



	template<typename pixel_t> class ClassifySkeletonCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		ClassifySkeletonCommand() : OneImageInPlaceCommand<pixel_t>("classifyskeleton", "Classifies line skeleton to end points, branch points, intersection points, and edge points.")
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			classifySkeleton(in, true, false, true);
		}
	};




	template<typename pixel_t> class TraceLineSkeletonBlockCommand : public Command, public Distributable
	{
	public:
		TraceLineSkeletonBlockCommand() : Command("tracelineskeletonblock", "This is an internal command used by the tracelineskeleton command to trace a block of a line skeleton when distributed processing is enabled.",
			{
				CommandArgument<Image<pixel_t> >(In, "skeleton", "Image containing the skeleton. The pixels of the image will be set to zero."),
				CommandArgument<Image<pixel_t> >(In, "original", "Original image from which the skeleton has been calculated. This image is used for branch shape measurements."),
				CommandArgument<string>(In, "filename", "Name template for file where the resulting network will be saved."),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing."),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current block in coordinates of the full image."),
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& orig = *pop<Image<pixel_t>* >(args);
			string filename = pop<string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE index = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

			filename += "_" + itl2::toString(index) + ".dat";

			// Trace (in multithreaded manner)
			vector<Network> nets;
			internals::traceLineSkeletonBlocks(in, orig, nets, origin);

			// Write all networks to the output file
			cout << "Writing " << nets.size() << " graphs to " << filename << endl;
			for (const Network& net : nets)
				net.write(filename, true);
		}

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			throw ITLException("Don't use this command directly - use tracelineskeleton instead.");
		}
	};


	template<typename pixel_t> class TraceLineSkeletonCommand : public Command, public Distributable
	{
	public:
		TraceLineSkeletonCommand() : Command("tracelineskeleton", "Traces a line skeleton into a graph structure. Each branch intersection point becomes a vertex in the graph and each branch becomes an edge.",
			{
				CommandArgument<Image<pixel_t> >(In, "skeleton", "Image containing the skeleton. The pixels of the image will be set to zero."),
				CommandArgument<Image<pixel_t> >(In, "original", "Original image from which the skeleton has been calculated. This image is used for branch shape measurements."),
				CommandArgument<Image<float32_t> >(Out, "vertices", "Image where vertex coordinates are stored. The size of the image is set to 3xN during processing, where N is the number of vertices in the graph."),
				CommandArgument<Image<uint64_t> >(Out, "edges", "Image where vertex indices corresponding to each edge will be set. The size of the image is set to 2xM where M is the number of edges. Each row of the image consists of a pair of indices to the vertex array."),
				CommandArgument<Image<float32_t> >(Out, "edge measurements", "Image that stores (pointCount, length, cross-sectional area, end distance, adjusted end distance) for each edge. The size of the image is set to 5xN during processing, where N is the number of edges in the graph. Each row contains properties of edge at corresponding row in the edges image."),
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& orig = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& vertices = *pop<Image<float32_t>* >(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>* >(args);
			Image<float32_t>& measurements = *pop<Image<float32_t>* >(args);
			
			Network net;
			traceLineSkeleton(in, orig, net);
			net.toImage(vertices, edges, &measurements);
		}

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<pixel_t>& orig = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<float32_t>& vertices = *pop<DistributedImage<float32_t>* >(args);
			DistributedImage<uint64_t>& edges = *pop<DistributedImage<uint64_t>* >(args);
			DistributedImage<float32_t>& measurements = *pop<DistributedImage<float32_t>* >(args);

			// Command that traces skeleton without combining incomplete vertices, and saves every subnetwork to given file
			TraceLineSkeletonBlockCommand<pixel_t> cmd;
			
			
			// Create temp file path
			unsigned int seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
			std::mt19937 gen(seed);
			string tempFilename = string("./tmp_images/skeleton_data_") + itl2::toString(gen());
			fs::remove(tempFilename);

			vector<ParamVariant> args2;
			ParamVariant p1;		// The source image
			p1.dimgval = &in;
			ParamVariant p11;		// The original image
			p11.dimgval = &orig;
			ParamVariant p2;		// Name of target file where all the graphs are saved
			p2.sval = &tempFilename;
			ParamVariant p3;		// Placeholder for block index (will be filled by the distributor)
			p3.ival = 0;
			ParamVariant p4;		// Placeholder for block origin (will be filled by the distributor)
			p4.vix = 0;
			p4.viy = 0;
			p4.viz = 0;
			args2.push_back(p1);
			args2.push_back(p11);
			args2.push_back(p2);
			args2.push_back(p3);
			args2.push_back(p4);

            cout << "Distributed tracing..." << endl;
			vector<string> output = distributor.distribute(&cmd, args2, 2, Vec3c(0, 0, 0));
            
            /*
            vector<string> output;
            output.push_back("abc1");
            output.push_back("abc2");
            output.push_back("abc3");
            output.push_back("abc4");
            output.push_back("abc5");
            output.push_back("abc6");
			string tempFilename = "./tmp_images/skeleton_data_2076255946";
			*/
            
            cout << "Loading data..." << endl;
			// Load the data files and combine all the graphs
			vector<Network> nets;
			for (size_t n = 0; n < output.size(); n++)
			{
				vector<Network> subnets;
				string fname = tempFilename + "_" + itl2::toString(n) + ".dat";
				
				cout << "Reading " << fname << endl;
				Network::read(fname, subnets);
				//cout << "Read " << subnets.size() << " subgraphs." << endl;
				nets.insert(nets.end(), subnets.begin(), subnets.end());
				fs::remove(fname);
			}
			
			//cout << "Reading done." << endl;

			Network fullnet;
			internals::combineTracedBlocks(nets, fullnet, true);

			// Now fullnet contains the whole network (and nets contains empty networks)
			// Convert it to images locally and set it to the outputs.
			Image<float32_t> verticesLocal;
			Image<uint64_t> edgesLocal;
			Image<float32_t> measurementsLocal;
			vertices.readTo(verticesLocal);
			edges.readTo(edgesLocal);
			measurements.readTo(measurementsLocal);

			fullnet.toImage(verticesLocal, edgesLocal, &measurementsLocal);

			vertices.setData(verticesLocal);
			edges.setData(edgesLocal);
			measurements.setData(measurementsLocal);
		}
	};

	class CleanSkeletonCommand : public Command
	{
	public:
		CleanSkeletonCommand() : Command("cleanskeleton", "Removes straight-through and isolated nodes from the network (i.e. all nodes that have either 0 or 2 neighbours, i.e. all nodes whose degree is 0 or 2).",
			{
				CommandArgument<Image<float32_t> >(InOut, "vertices", "Image where vertex coordinates are stored."),
				CommandArgument<Image<uint64_t> >(InOut, "edges", "Image where vertex indices corresponding to each edge are stored."),
				CommandArgument<Image<float32_t> >(InOut, "edge measurements", "Image where length and cross-sectional area of each edge is stored."),
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<float32_t>& vertices = *pop<Image<float32_t>* >(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>* >(args);
			Image<float32_t>& measurements = *pop<Image<float32_t>* >(args);

			Network net;
			net.fromImage(vertices, edges, &measurements);
			net.disconnectStraightThroughNodes(true);
			net.removeIsolatedNodes(true);
			net.toImage(vertices, edges, &measurements);
		}
	};

}
