#pragma once

#include "commandsbase.h"
#include "overlapdistributable.h"
#include "pilibutilities.h"
#include "network.h"
#include "traceskeleton.h"
#include "lineskeleton.h"
#include "surfaceskeleton.h"
#include "commandlist.h"
#include "fillskeleton.h"
#include "othercommands.h"
#include "traceskeletonpoints.h"

namespace pilib
{
	inline std::string skeleDistributionNote()
	{
		return "This command is not guaranteed to give the same result in both normal and distributed processing mode. Despite that, both modes should give a valid result.";
	}

	inline std::string skeleSeeAlso()
	{
		return "surfacethin, surfaceskeleton, linethin, lineskeleton, tracelineskeleton, classifyskeleton";
	}

	template<typename pixel_t> class SurfaceThinCommand : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		SurfaceThinCommand() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("surfacethin", "Thins one layer of pixels from the foreground of the image. Positive pixels are assumed to belong to the foreground. Run iteratively to calculate a surface skeleton. " + skeleDistributionNote(),
			{
				CommandArgument<bool>(ParameterDirection::In, "retain surfaces", "Set to false to allow thinning of surfaces to lines if the surface does not surround a cavity.", true)
			},
			skeleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			bool retainSurfaces = pop<bool>(args);
			size_t changed = thin(in, retainSurfaces);
			std::cout << std::endl << changed << " pixels removed." << std::endl;
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			return Vec3c(10, 10, 10);
		}

		virtual bool canDelay(const vector<ParamVariant>& args) const override
		{
			// We use output so no delaying is allowed.
			return false;
		}
	};

	template<typename pixel_t> class LineThinCommand : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		LineThinCommand() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("linethin", "Thins one layer of pixels from the foreground of the image. Positive pixels are assumed to belong to the foreground. Run iteratively to calculate a line skeleton. " + skeleDistributionNote(),
			{},
			skeleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			size_t changed = lineThin(in);
			std::cout << std::endl << changed << " pixels removed." << std::endl;
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			return Vec3c(10, 10, 10);
		}

		virtual bool canDelay(const vector<ParamVariant>& args) const override
		{
			// We use output so no delaying is allowed.
			return false;
		}
	};

	template<typename command_t, typename base_t> class IterableDistributable : public base_t, public Distributable
	{
	protected:
		friend class CommandList;

		IterableDistributable(const string& name, const string& help, const vector<CommandArgumentBase>& extraArgs = {}, const string& seeAlso = "") : base_t(name, help, extraArgs, seeAlso)
		{
		}

	public:
		virtual Vec3c getMargin(const vector<ParamVariant>& args) const override
		{
			return Vec3c(10, 10, 10);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			size_t lastTotalChanged = 0;
			size_t n = 0;
			while (true)
			{
				std::cout << "Iteration " << n << std::endl;

				// Run one iteration of thinning
				// TODO: It is not necessary to run all blocks if those and their neighbours did not change in the last iteration.
				command_t& cmd = CommandList::get<command_t>();
				vector<string> output = distributor.distribute(&cmd, args);

				// Calculate total number of changed pixels
				size_t totalChanged = parseTotalCount(output, "pixels removed");
				
				std::cout << std::endl << totalChanged << " pixels removed." << std::endl;

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

			return vector<string>();
		}
	};

	//template<typename pixel_t> class HybridSkeletonCommand : public IterableDistributable<HybridThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >
	//{
	//protected:
	//	friend class CommandList;

	//	HybridSkeletonCommand() : IterableDistributable<HybridThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >("hybridskeleton", "Calculates skeleton of the foreground of the given image. Positive pixels are assumed to belong to the foreground. The skeleton contains both lines and plates.")
	//	{
	//	}

	//public:
	//	virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
	//	{
	//		hybridSkeleton(in);
	//	}
	//};

	template<typename pixel_t> class SurfaceSkeletonCommand : public IterableDistributable<SurfaceThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		SurfaceSkeletonCommand() : IterableDistributable<SurfaceThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >("surfaceskeleton", "Calculates skeleton of the foreground of the given image. Positive pixels are assumed to belong to the foreground. The skeleton may contain both lines and plates. " + skeleDistributionNote(),
			{
				CommandArgument<bool>(ParameterDirection::In, "retain surfaces", "Set to false to allow thinning of surfaces to lines if the surface does not surround a cavity.", true)
			},
			skeleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			bool retainSurfaces = pop<bool>(args);
			surfaceSkeleton(in, retainSurfaces);
		}
	};

	template<typename pixel_t> class LineSkeletonCommand : public IterableDistributable<LineThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		LineSkeletonCommand() : IterableDistributable<LineThinCommand<pixel_t>, OneImageInPlaceCommand<pixel_t> >("lineskeleton", "Calculates skeleton of the foreground of the given image. Positive pixels are assumed to belong to the foreground. The skeleton contains only lines (no plates). " + skeleDistributionNote(),
			{},
			skeleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			lineSkeleton(in);
		}
	};



	template<typename pixel_t> class ClassifySkeletonCommand : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		ClassifySkeletonCommand() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("classifyskeleton",
"Classifies line skeleton to end points, curve points, branch points, junction points, internal points, and edge points according to Arcelli, From 3D Discrete Surface Skeletons to Curve Skeletons."
"End points are given value 2, curve points value 3, branch points value 4, junction points value 5, internal points value 6, and edge points value 7.",
			{},
			skeleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			classifySkeleton(in, true, false, true);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			return Vec3c(10, 10, 10);
		}
	};




	template<typename pixel_t> class ClassifyForTracingCommand : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		ClassifyForTracingCommand() : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >("classifyskeletonfortracing", "Classifies line skeleton to branch points and curve points. This command is used internally during processing of tracelineskeleton command in distributed mode.")
		{
		}

	public:
		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			itl2::internals::classifyForTracing(in);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			return Vec3c(10, 10, 10);
		}
	};


	inline std::string tracedSeeAlso()
	{
		return "cleanskeleton, pruneskeleton, removeedges, fillskeleton, getpointsandlines";
	}
	
	template<typename pixel_t> class TraceLineSkeletonBlockCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		TraceLineSkeletonBlockCommand() : Command("tracelineskeletonblock", "This is an internal command used by the tracelineskeleton command to trace a block of a line skeleton when distributed processing is enabled. The input image must be classified before a call to this command.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "skeleton", "Image containing the skeleton. The pixels of the image will be set to zero."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "original", "Original image from which the skeleton has been calculated. This image is used for branch shape measurements."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name template for file where the resulting network will be saved."),
				CommandArgument<bool>(ParameterDirection::In, "store all edge points", "Set to true to store all points of each edge to edge points image. If set to false, only single point on each edge is stored. This argument is required to be set to true if the graph will be converted to points and lines format later.", false),
				CommandArgument<size_t>(ParameterDirection::In, "thread count", "Count of threads to use in the tracing process. Set to zero to determine count of threads automatically. Set to one to use single-threaded processing.", 0),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "Standard deviation value to be used in smoothing of the edges (using anchored convolution) before length measurements.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "max displacement", "Maximum displacement of points in anchored convolution done before length measurements.", 0.5),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing."),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current block in coordinates of the full image."),
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
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>* pOrig = pop<Image<pixel_t>* >(args);
			string filename = pop<string>(args);
			bool storeAllEdgePoints = pop<bool>(args);
			coord_t threadCount = (coord_t)pop<size_t>(args);
			double smoothingSigma = pop<double>(args);
			double maxDisplacement = pop<double>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE index = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

			filename += "_" + itl2::toString(index) + ".dat";

			// Trace (in multithreaded manner)
			vector<Network> nets;
			itl2::internals::traceLineSkeletonBlocks(in, pOrig, storeAllEdgePoints, smoothingSigma, maxDisplacement, nets, Vec3sc(origin), threadCount);

			// Write all networks to the output file
			std::cout << "Writing " << nets.size() << " graphs to " << filename << std::endl;
			for (const Network& net : nets)
				net.write(filename, true);

			// Combine what we can combine now as here we have the original image still available for updated area measurements.
			// NOTE: This should not matter as traceLineSkeletonBlocks has the same pOrig data available during tracing.
			//Network finalNet;
			//internals::combineTracedBlocks(nets, finalNet, in.dimensions(), true, pOrig, false);
			//finalNet.write(filename, false);
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}
	};

	template<typename pixel_t> class TraceLineSkeletonBlock2Command : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		TraceLineSkeletonBlock2Command() : Command("tracelineskeletonblock", "This is an internal command used by the tracelineskeleton command to trace a block of a line skeleton when distributed processing is enabled. The input image must be classified before a call to this command.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "skeleton", "Image containing the skeleton. The pixels of the image will be set to zero."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name template for file where the resulting network will be saved."),
				CommandArgument<bool>(ParameterDirection::In, "store all edge points", "Set to true to store all points of each edge to edge points image. If set to false, only single point on each edge is stored. This argument is required to be set to true if the graph will be converted to points and lines format later.", false),
				CommandArgument<size_t>(ParameterDirection::In, "thread count", "Count of threads to use in the tracing process. Set to zero to determine count of threads automatically. Set to one to use single-threaded processing.", 0),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "Standard deviation value to be used in smoothing of the edges (using anchored convolution) before length measurements.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "max displacement", "Maximum displacement of points in anchored convolution done before length measurements.", 0.5),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing."),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current block in coordinates of the full image."),
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
			// Insert null pointer to arguments list instead of pointer to original image.
			args.insert(args.begin() + 1, (Image<pixel_t>*)0);
			CommandList::get<TraceLineSkeletonBlockCommand<pixel_t> >().run(args);
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}
	};



	class CombineTracedBlocksCommand : public Command//, public Distributable
	{
	protected:
		friend class CommandList;

		CombineTracedBlocksCommand() : Command("combinetracedblocks", "This is an internal command used by the tracelineskeleton command. It combines results of tracelineskeletonblock commands.",
			{
				CommandArgument<string>(ParameterDirection::In, "temp file prefix", ""),
				CommandArgument<size_t>(ParameterDirection::In, "temp file count", ""),
				CommandArgument<Vec3c>(ParameterDirection::In, "original dimensions", ""),
				CommandArgument<bool>(ParameterDirection::In, "store all edge points", ""),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", ""),
				CommandArgument<double>(ParameterDirection::In, "max displacement", ""),
				CommandArgument<string>(ParameterDirection::In, "vertices filename", ""),
				CommandArgument<string>(ParameterDirection::In, "edges filename", ""),
				CommandArgument<string>(ParameterDirection::In, "measurements filename", ""),
				CommandArgument<string>(ParameterDirection::In, "points filename", ""),
			},
			"tracelineskeleton")
		{
		}

	public:
		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(vector<ParamVariant>& args) const override
		{
			string tempFilename = pop<string>(args);	// tempFileName, prefix of {prefix}_{n}.dat files that contain the subnets
			size_t outputCount = pop<size_t>(args);			// output.size(), subnet count
			Vec3c inDimensions = pop<Vec3c>(args);
			bool storeAllEdgePoints = pop<bool>(args);
			double smoothingSigma = pop<double>(args);
			double maxDisplacement = pop<double>(args);
			string verticesFilename = pop<string>(args);
			string edgesFilename = pop<string>(args);
			string measurementsFilename = pop<string>(args);
			string pointsFilename = pop<string>(args);

			std::cout << "Loading data..." << std::endl;

			// Load the data files and combine all the graphs
			vector<Network> nets;
			for (size_t n = 0; n < outputCount; n++)
			{
				vector<Network> subnets;
				string fname = tempFilename + "_" + itl2::toString(n) + ".dat";

				std::cout << "Reading " << fname << std::endl;
				Network::read(fname, subnets);
				//std::cout << "Read " << subnets.size() << " subgraphs." << std::endl;
				nets.insert(nets.end(), subnets.begin(), subnets.end());
			}

			std::cout << "Reading done." << std::endl;

			Network fullnet;
			// NOTE: We use uint8_t specialization as we don't have any image.
			itl2::internals::combineTracedBlocks<uint8_t>(nets, fullnet, inDimensions, true, nullptr, true, storeAllEdgePoints, smoothingSigma, maxDisplacement);
			itl2::internals::finalize(fullnet);

			// Now fullnet contains the whole network (and nets contains empty networks)
			// Convert it to images locally and set it to the outputs.
			Image<float32_t> verticesLocal;
			Image<uint64_t> edgesLocal;
			Image<float32_t> measurementsLocal;
			Image<int32_t> pointsLocal;
			fullnet.toImage(verticesLocal, edgesLocal, &measurementsLocal, &pointsLocal);

			raw::writed(verticesLocal, verticesFilename);
			raw::writed(edgesLocal, edgesFilename);
			raw::writed(measurementsLocal, measurementsFilename);
			raw::writed(pointsLocal, pointsFilename);
		}

		//using Distributable::runDistributed;

		//virtual vector<string> runDistributed(Distributor & distributor, vector<ParamVariant> & args) const override
		//{
		//	

		//	return vector<string>();
		//}
	};



	template<typename pixel_t> class TraceLineSkeletonCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		TraceLineSkeletonCommand() : Command("tracelineskeleton", "Traces a line skeleton into a graph structure. Each branch intersection point becomes a vertex in the graph and each branch becomes an edge.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "skeleton", "Image containing the skeleton. The pixels of the image will be set to zero."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "original", "Original image from which the skeleton has been calculated. This image is used for branch shape measurements."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "vertices", "Image where vertex coordinates are stored. The size of the image is set to 3xN during processing, where N is the number of vertices in the graph."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::Out, "edges", "Image where vertex indices corresponding to each edge will be set. The size of the image is set to 2xM where M is the number of edges. Each row of the image consists of a pair of indices to the vertex array."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "edge measurements", "Image that stores (pointCount, length, cross-sectional area) for each edge. The size of the image is set to 3xN during processing, where N is the number of edges in the graph. Each row contains properties of edge at corresponding row in the edges image. Area measurements are approximations if distributed processing is allowed."),
				CommandArgument<Image<int32_t> >(ParameterDirection::Out, "edge points", "Image that stores some points on each edge. The points are required for skeleton filling commands to work correctly even if the skeleton is pruned. The points are stored such that the first edgeCount pixels of the image store count of points for each edge. The remaining pixels store (x, y, z) coordinates of each point and each edge sequentially. For example, the format for two edges that have 1 and 2 points is therefore '1 2 x11 y11 z11 x21 y21 z21 x22 y22 z22', where Aij is A-component of j:th point of edge i."),
				CommandArgument<bool>(ParameterDirection::In, "store all edge points", "Set to true to store all points of each edge to edge points image. If set to false, only single point on each edge is stored. This argument is required to be set to true if the graph will be converted to points and lines format later.", false),
				CommandArgument<size_t>(ParameterDirection::In, "thread count", "Count of threads to use in the tracing process. Set to zero to determine count of threads automatically. Set to one to use single-threaded processing.", 0),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "Standard deviation value to be used in smoothing of the edges (using anchored convolution) before length measurements.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "max displacement", "Maximum displacement of points in anchored convolution done before length measurements.", 0.5),
			},
			"surfaceskeleton, lineskeleton, " + tracedSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>* pOrig = pop<Image<pixel_t>* >(args);
			Image<float32_t>& vertices = *pop<Image<float32_t>* >(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>* >(args);
			Image<float32_t>& measurements = *pop<Image<float32_t>* >(args);
			Image<int32_t>& points = *pop<Image<int32_t>*>(args);
			bool storeAllEdgePoints = pop<bool>(args);
			coord_t threadCount = (coord_t)pop<size_t>(args);
			double smoothingSigma = pop<double>(args);
			double maxDisplacement = pop<double>(args);
			
			Network net;
			traceLineSkeleton(in, pOrig, storeAllEdgePoints, smoothingSigma, maxDisplacement, net, (int)threadCount);
			net.toImage(vertices, edges, &measurements, &points);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<pixel_t>* pOrig = pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<float32_t>& vertices = *pop<DistributedImage<float32_t>* >(args);
			DistributedImage<uint64_t>& edges = *pop<DistributedImage<uint64_t>* >(args);
			DistributedImage<float32_t>& measurements = *pop<DistributedImage<float32_t>* >(args);
			DistributedImage<int32_t>& points = *pop<DistributedImage<int32_t>*>(args);
			bool storeAllEdgePoints = pop<bool>(args);
			size_t threadCount = pop<size_t>(args);
			double smoothingSigma = pop<double>(args);
			double maxDisplacement = pop<double>(args);
			

			// Classify skeleton
			auto& cft = CommandList::get<ClassifyForTracingCommand<pixel_t> >();
			cft.runDistributed(distributor, { &in });

			
			// Create temp file path
			string tempFilename = createTempFilename("skeleton_data");

			// Command that traces skeleton without combining incomplete vertices, and saves every subnetwork to given file
			vector<string> output;
			if (pOrig)
			{
				in.checkSize(pOrig->dimensions());
				auto& cmd = CommandList::get<TraceLineSkeletonBlockCommand<pixel_t> >();
				output = cmd.runDistributed(distributor, { &in, pOrig, tempFilename, storeAllEdgePoints, threadCount, smoothingSigma, maxDisplacement, Distributor::BLOCK_INDEX_ARG_TYPE(), Distributor::BLOCK_ORIGIN_ARG_TYPE() });
			}
			else
			{
				auto& cmd = CommandList::get<TraceLineSkeletonBlock2Command<pixel_t>>();
				output = cmd.runDistributed(distributor, { &in, tempFilename, storeAllEdgePoints, threadCount, smoothingSigma, maxDisplacement, Distributor::BLOCK_INDEX_ARG_TYPE(), Distributor::BLOCK_ORIGIN_ARG_TYPE() });
			}


			combineBlocks(distributor, 
				tempFilename, output.size(), in.dimensions(),
				storeAllEdgePoints, smoothingSigma, maxDisplacement,
				vertices, edges, measurements, points);

			std::cout << "Deleting temporary files..." << std::endl;
			for (size_t n = 0; n < output.size(); n++)
			{
				string fname = tempFilename + "_" + itl2::toString(n) + ".dat";
				fs::remove(fname);
			}

			///////////////////
   //         
   //         std::cout << "Loading data..." << std::endl;
			//// Load the data files and combine all the graphs
			//vector<Network> nets;
			//for (size_t n = 0; n < output.size(); n++)
			//{
			//	vector<Network> subnets;
			//	string fname = tempFilename + "_" + itl2::toString(n) + ".dat";
			//	
			//	std::cout << "Reading " << fname << std::endl;
			//	Network::read(fname, subnets);
			//	//std::cout << "Read " << subnets.size() << " subgraphs." << std::endl;
			//	nets.insert(nets.end(), subnets.begin(), subnets.end());
			//	fs::remove(fname);
			//}
			//
			////std::cout << "Reading done." << std::endl;

			//Network fullnet;
			//itl2::internals::combineTracedBlocks<pixel_t>(nets, fullnet, in.dimensions(), true, nullptr, true, storeAllEdgePoints, smoothingSigma, maxDisplacement);
			//itl2::internals::finalize(fullnet);

			//// Now fullnet contains the whole network (and nets contains empty networks)
			//// Convert it to images locally and set it to the outputs.
			//Image<float32_t> verticesLocal;
			//Image<uint64_t> edgesLocal;
			//Image<float32_t> measurementsLocal;
			//Image<int32_t> pointsLocal;
			//vertices.readTo(verticesLocal);
			//edges.readTo(edgesLocal);
			//measurements.readTo(measurementsLocal);


			//fullnet.toImage(verticesLocal, edgesLocal, &measurementsLocal, &pointsLocal);

			//vertices.setData(verticesLocal);
			//edges.setData(edgesLocal);
			//measurements.setData(measurementsLocal);
			//points.setData(pointsLocal);

			///////////////////////////

			return vector<string>();
		}

	private:

		void combineBlocks(Distributor& distributor,
			const string& tempFilename, size_t outputCount, const Vec3c& inDimensions,
			bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement,
			DistributedImage<float32_t>& vertices,
			DistributedImage<uint64_t>& edges,
			DistributedImage<float32_t>& measurements,
			DistributedImage<int32_t>& points) const
		{
			// TODO: This could be automated in distributor such that all non-distributable commands can be run as jobs,
			// and reading and writing of the required images would be automatic.

			// Generate temp file prefixes
			string verticesFilename = createTempFilename("temp_vertices");
			string edgesFilename = createTempFilename("temp_edges");
			string measurementsFilename = createTempFilename("temp_measurements");
			string pointsFilename = createTempFilename("temp_points");

			// Build pi2 script. We make the script manually as we don't know the size of the output images yet.
			std::stringstream script;

			// Init so that we always print something (required at least in SLURM distributor)
			script << "echo(true, false);" << std::endl;

			// Save network to temp images
			script << "combinetracedblocks(" << tempFilename << ", "
											<< outputCount << ", "
											<< inDimensions << ", "
											<< storeAllEdgePoints << ", "
											<< smoothingSigma << ", "
											<< maxDisplacement << ", "
											<< verticesFilename << ", "
											<< edgesFilename << ", "
											<< measurementsFilename << ", "
											<< pointsFilename << ")" << std::endl;

			// Run the job
			distributor.flush();
			distributor.submitJob(script.str(), JobType::Normal);
			distributor.waitForJobs();

			// Get sizes of temp images
			Vec3c verticesDims, edgesDims, measurementsDims, pointsDims;
			ImageDataType dummy;
			string dummyString;
			raw::getInfo(verticesFilename, verticesDims, dummy, dummyString);
			raw::getInfo(edgesFilename, edgesDims, dummy, dummyString);
			raw::getInfo(measurementsFilename, measurementsDims, dummy, dummyString);
			raw::getInfo(pointsFilename, pointsDims, dummy, dummyString);

			// Set sizes of output DistributedImages
			vertices.ensureSize(verticesDims);
			edges.ensureSize(edgesDims);
			measurements.ensureSize(measurementsDims);
			points.ensureSize(pointsDims);

			// File system move temp images to CurrentReadSource of DistributedImages
			// => No temp images are in RAM + no temp images need to be read (file rename is enough)

			raw::internals::expandRawFilename(verticesFilename);
			raw::internals::expandRawFilename(edgesFilename);
			raw::internals::expandRawFilename(measurementsFilename);
			raw::internals::expandRawFilename(pointsFilename);

			std::cout << "Removing " << vertices.currentWriteTarget() << std::endl;
			
			fs::remove(vertices.currentWriteTarget());
			fs::remove(edges.currentWriteTarget());
			fs::remove(measurements.currentWriteTarget());
			fs::remove(points.currentWriteTarget());
			
			std::cout << "Renaming " << verticesFilename << " to " << vertices.currentWriteTarget() << std::endl;

			fs::rename(verticesFilename, vertices.currentWriteTarget());
			fs::rename(edgesFilename, edges.currentWriteTarget());
			fs::rename(measurementsFilename, measurements.currentWriteTarget());
			fs::rename(pointsFilename, points.currentWriteTarget());

			vertices.writeComplete();
			edges.writeComplete();
			measurements.writeComplete();
			points.writeComplete();
		}
	};

	template<typename pixel_t> class TraceLineSkeleton2Command : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		TraceLineSkeleton2Command() : Command("tracelineskeleton", "Traces a line skeleton into a graph structure. Each branch intersection point becomes a vertex in the graph and each branch becomes an edge.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "skeleton", "Image containing the skeleton. The pixels of the image will be set to zero."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "vertices", "Image where vertex coordinates are stored. The size of the image is set to 3xN during processing, where N is the number of vertices in the graph."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::Out, "edges", "Image where vertex indices corresponding to each edge will be set. The size of the image is set to 2xM where M is the number of edges. Each row of the image consists of a pair of indices to the vertex array."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "edge measurements", "Image that stores (pointCount, length, cross-sectional area) for each edge. The size of the image is set to 3xN during processing, where N is the number of edges in the graph. Each row contains properties of edge at corresponding row in the edges image."),
				CommandArgument<Image<int32_t> >(ParameterDirection::Out, "edge points", "Image that stores some points on each edge. The points are required for skeleton filling commands to work correctly even if the skeleton is pruned. The points are stored such that the first edgeCount pixels of the image store count of points for each edge. The remaining pixels store (x, y, z) coordinates of each point and each edge sequentially. For example, the format for two edges that have 1 and 2 points is therefore '1 2 x11 y11 z11 x21 y21 z21 x22 y22 z22', where Aij is A-component of j:th point of edge i."),
				CommandArgument<bool>(ParameterDirection::In, "store all edge points", "Set to true to store all points of each edge to edge points image. If set to false, only single point on each edge is stored. This argument is required to be set to true if the graph will be converted to points and lines format later.", false),
				CommandArgument<size_t>(ParameterDirection::In, "thread count", "Count of threads to use in the tracing process. Set to zero to determine count of threads automatically. Set to one to use single-threaded processing.", 0),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "Standard deviation value to be used in smoothing of the edges (using anchored convolution) before length measurements.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "max displacement", "Maximum displacement of points in anchored convolution done before length measurements.", 0.5),
			},
			"surfaceskeleton, lineskeleton, " + tracedSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			args.insert(args.begin() + 1, (Image<pixel_t>*)0);
			CommandList::get<TraceLineSkeletonCommand<pixel_t> >().run(args);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			args.insert(args.begin() + 1, (DistributedImage<pixel_t>*)0);
			return CommandList::get<TraceLineSkeletonCommand<pixel_t> >().runDistributed(distributor, args);
		}
	};

	class CleanSkeletonCommand : public Command
	{
	protected:
		friend class CommandList;

		CleanSkeletonCommand() : Command("cleanskeleton", "Removes straight-through and isolated nodes from a network traced from a skeleton (i.e. all nodes that have either 0 or 2 neighbours, i.e. all nodes whose degree is 0 or 2).",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "vertices", "Image where vertex coordinates are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::InOut, "edges", "Image where vertex indices corresponding to each edge are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "edge measurements", "Image where properties of each edge are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<int32_t> >(ParameterDirection::InOut, "edge points", "Image that stores some points on each edge. See `tracelineskeleton` command.")
			},
			tracedSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<float32_t>& vertices = *pop<Image<float32_t>* >(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>* >(args);
			Image<float32_t>& measurements = *pop<Image<float32_t>* >(args);
			Image<int32_t>& points = *pop<Image<int32_t>*>(args);

			Network net;
			net.fromImage(vertices, edges, &measurements, &points);
			net.removeStraightThroughNodes(true);
			net.toImage(vertices, edges, &measurements, &points);
		}
	};

	template<typename pixel_t> class RemoveEdgesCommand : public Command
	{
	protected:
		friend class CommandList;

		RemoveEdgesCommand() : Command("removeedges", "Prunes the network by removing user-selected edges from it. The edges to be removed are specified in an image, having one pixel for each edge, where non-zero value specifies that the corresponding edge is to be removed.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "vertices", "Image where vertex coordinates are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::InOut, "edges", "Image where vertex indices corresponding to each edge are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "edge measurements", "Image where properties of each edge are stored.  ee `tracelineskeleton` command."),
				CommandArgument<Image<int32_t> >(ParameterDirection::InOut, "edge points", "Image that stores some points on each edge. See `tracelineskeleton` command."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "flags", "Image that has as many pixels as there are edges. Edges corresponding to non-zero pixels are removed."),
				CommandArgument<bool>(ParameterDirection::In, "disconnect straight-through nodes", "If set to true, all straight-through nodes (nodes with degree = 2) are removed after pruning. This operation might change the network even if no edges are pruned."),
				CommandArgument<bool>(ParameterDirection::In, "remove isolated nodes", "If set to true, all isolated nodes (nodes with degree = 0) are removed after pruning. This operation might change the network even if no edges are pruned.")
			},
			tracedSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<float32_t>& vertices = *pop<Image<float32_t>* >(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>* >(args);
			Image<float32_t>& measurements = *pop<Image<float32_t>* >(args);
			Image<int32_t>& points = *pop<Image<int32_t>*>(args);
			Image<pixel_t>& flags = *pop<Image<pixel_t>*>(args);
			bool disconnectStraightThrough = pop<bool>(args);
			bool removeIsolated = pop<bool>(args);

			Network net;
			net.fromImage(vertices, edges, &measurements, &points);

			if (flags.pixelCount() != net.edges.size())
				throw ITLException("Invalid number of flags. There must be one flag for each edge.");

			vector<coord_t> edgeIndices;

			for (coord_t n = 0; n < flags.pixelCount(); n++)
			{
				if (flags(n) != 0)
					edgeIndices.push_back(n);
			}

			net.removeEdges(edgeIndices, disconnectStraightThrough, removeIsolated, true);

			net.toImage(vertices, edges, &measurements, &points);
		}
	};

	class PruneSkeletonCommand : public Command
	{
	protected:
		friend class CommandList;

		PruneSkeletonCommand() : Command("pruneskeleton", "Prunes a traced skeleton. Removes all edges that end in a node with no other branches connected to it (degree = 1) and whose length is less than specified value.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "vertices", "Image where vertex coordinates are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::InOut, "edges", "Image where vertex indices corresponding to each edge are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "edge measurements", "Image where properties of each edge are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<int32_t> >(ParameterDirection::InOut, "edge points", "Image that stores some points on each edge. See `tracelineskeleton` command."),
				CommandArgument<double>(ParameterDirection::In, "maximum length", "Edges shorter than this are pruned if they are not connected to other edges in both ends."),
				CommandArgument<bool>(ParameterDirection::In, "disconnect straight-through nodes", "If set to true, all straight-through nodes (nodes with degree = 2) are removed after pruning. This operation might change the network even if no edges are pruned."),
				CommandArgument<bool>(ParameterDirection::In, "remove isolated nodes", "If set to true, all isolated nodes (nodes with degree = 0) are removed after pruning. This operation might change the network even if no edges are pruned.")
			},
			tracedSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<float32_t>& vertices = *pop<Image<float32_t>* >(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>* >(args);
			Image<float32_t>& measurements = *pop<Image<float32_t>* >(args);
			Image<int32_t>& points = *pop<Image<int32_t>*>(args);
			double maxLength = pop<double>(args);
			bool disconnectStraightThrough = pop<bool>(args);
			bool removeIsolated = pop<bool>(args);

			Network net;
			net.fromImage(vertices, edges, &measurements, &points);
			net.prune((float32_t)maxLength, disconnectStraightThrough, removeIsolated, true);
			net.toImage(vertices, edges, &measurements, &points);
		}
	};


	class GetPointsAndLinesCommand : public Command
	{
	protected:
		friend class CommandList;

		GetPointsAndLinesCommand() : Command("getpointsandlines", "Converts a traced skeleton to points-and-lines format, where all points on the skeleton are stored in one array and all branches are stored as lists of indices of points that form the branch.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "vertices", "Image where vertex coordinates are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::In, "edges", "Image where vertex indices corresponding to each edge are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "edge measurements", "Image where properties of each edge are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<int32_t> >(ParameterDirection::In, "edge points", "Image that stores some points on each edge. See `tracelineskeleton` command."),
				
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "points", "Image that stores coordinates of all points in the network. These include points on edges and centroids of intersection regions. Each row of image stores (x, y, z) coordinates of one point."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::Out, "lines", "Image that stores a list of point indices for each edge. This image is stored in compressed format: [count of edges][count of points in 1st edge][index of point 1][index of point 2]...[count of points in 2nd edge][index of point 1][index of point 2]..."),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "If smooth lines are requested instead of jagged lines (consisting of pixel locations), specify positive value here. The value is standard deviation of a Gaussian kernel used to smooth the lines in an anchored convolution smoothing. The end points of each line are not changed by smoothing. Values in range [0.5, 1.5] often give suitable amount of smoothing. The smoothing algorithm is described in Suhadolnik - An anchored discrete convolution algorithm for measuring length in digital images.", 0.0),
				CommandArgument<double>(ParameterDirection::In, "max displacement", "Maximum displacement of points in anchored convolution smoothing of the lines.", 0.5),
			},
			"writevtk, " + tracedSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<float32_t>& vertices = *pop<Image<float32_t>* >(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>* >(args);
			Image<float32_t>& measurements = *pop<Image<float32_t>* >(args);
			Image<int32_t>& edgePoints = *pop<Image<int32_t>*>(args);
			
			Image<float32_t>& points = *pop<Image<float32_t>*>(args);
			Image<uint64_t>& lines = *pop<Image<uint64_t>*>(args);

			double sigma = pop<double>(args);
			double maxDisplacement = pop<double>(args);


			Network net;
			net.fromImage(vertices, edges, &measurements, &edgePoints);

			vector<Vec3f> pointsv;
			vector<vector<size_t>> linesv;
			getPointsAndLines(net, pointsv, linesv, sigma, maxDisplacement);

			// Convert lists to images
			points.ensureSize(3, pointsv.size());
			for(size_t n = 0; n < pointsv.size(); n++)
			{
				const auto& v = pointsv[n];
				points(0, n) = v.x;
				points(1, n) = v.y;
				points(2, n) = v.z;
			}

			// Count total number of entries in lines lists
			size_t total = 0;
			for (size_t n = 0; n < linesv.size(); n++)
				total += 1 + linesv[n].size(); // count + items

			lines.ensureSize(total + 1);

			coord_t i = 0;
			lines(i++) = linesv.size();
			for (size_t n = 0; n < linesv.size(); n++)
			{
				const auto& list = linesv[n];
				lines(i++) = list.size();
				for (size_t m = 0; m < list.size(); m++)
				{
					lines(i++) = list[m];
				}
			}
		}
	};


	class WriteVtkCommandBase : public Command
	{
	protected:

		using Command::Command;

		void writeVtk(const Image<float32_t>& points, const Image<uint64_t>& lines, const string& filename, const string& pointDataNames = "", const Image<float32_t>* pointData = nullptr, const string& lineDataNames = "", Image<float32_t>* lineData = nullptr) const
		{
			if (lines.pixelCount() < 1)
				throw ITLException("Empty image passed as lines.");

			// Convert images to lists
			vector<Vec3f> pointsv;
			pointsv.reserve(points.height());
			for (coord_t n = 0; n < points.height(); n++)
			{
				Vec3f v(points(0, n), points(1, n), points(2, n));
				pointsv.push_back(v);
			}

			coord_t i = 0;
			coord_t edgeCount = (coord_t)lines(i++);
			vector<vector<size_t>> linesv;
			linesv.reserve(edgeCount);

			if (lines.pixelCount() < edgeCount + 1)
				throw ITLException("Lines image contains less pixels than there are elements in the lines.");

			for (coord_t n = 0; n < edgeCount; n++)
			{
				coord_t pointCount = (coord_t)lines(i++);
				vector<size_t> line;
				line.reserve(pointCount);

				if (lines.pixelCount() < i + pointCount)
					throw ITLException("Lines image is inconsistent, it contains too few pixels.");

				for (coord_t m = 0; m < pointCount; m++)
				{
					line.push_back(lines(i++));
				}

				linesv.push_back(line);
			}

			vector<std::tuple<string, vector<float32_t>>> pointDataArrays;
			if (pointDataNames != "" && pointData != nullptr)
			{
				vector<string> headers = split(pointDataNames, true, ',', true);

				if (headers.size() != pointData->width())
					throw ITLException("Count of headers does not match column count of point data image.");

				for (size_t n = 0; n < headers.size(); n++)
				{
					vector<float32_t> datav;
					datav.reserve(pointData->height());
					for (coord_t m = 0; m < pointData->height(); m++)
						datav.push_back((*pointData)(n, m));

					pointDataArrays.push_back(make_tuple(headers[n], datav));
				}
			}

			vector<std::tuple<string, vector<float32_t>>> lineDataArrays;
			if (lineDataNames != "" && lineData != nullptr)
			{
				vector<string> headers = split(lineDataNames, true, ',', true);

				if (headers.size() != lineData->width())
					throw ITLException("Count of headers does not match column count of line data image.");

				for (size_t n = 0; n < headers.size(); n++)
				{
					vector<float32_t> datav;
					datav.reserve(lineData->height());
					for (coord_t m = 0; m < lineData->height(); m++)
						datav.push_back((*lineData)(n, m));

					lineDataArrays.push_back(make_tuple(headers[n], datav));
				}
			}

			vtk::writed(pointsv, linesv, filename, &pointDataArrays, &lineDataArrays);

			//vector<float32_t> pointDatav;
			//if (pointDataNames != "" && pointData != nullptr)
			//{
			//	if (pointData->pixelCount() != pointsv.size())
			//		throw ITLException("Point data image does not contain one pixel for each point.");

			//	pointDatav.reserve(points.height());
			//	for (coord_t n = 0; n < pointData->pixelCount(); n++)
			//		pointDatav.push_back((*pointData)(n));
			//}

			//vector<float32_t> lineDatav;
			//if (lineDataNames != "" && lineData != nullptr)
			//{
			//	if (lineData->pixelCount() != linesv.size())
			//		throw ITLException("Line data image does not contain one pixel for each line.");

			//	lineDatav.reserve(lines.height());
			//	for (coord_t n = 0; n < lineData->pixelCount(); n++)
			//		lineDatav.push_back((*lineData)(n));
			//}

			//vtk::writed(pointsv, linesv, filename, pointDataName, &pointDatav, lineDataName, &lineDatav);
		}
	};

	class WriteVtkCommand : public WriteVtkCommandBase
	{
	protected:
		friend class CommandList;

		WriteVtkCommand() : WriteVtkCommandBase("writevtk", "Writes network in points-and-lines format to a .vtk file.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "points", "Image that stores coordinates of all points in the network. See `getpointsandlines` command."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::In, "lines", "Image that stores a list of point indices for each edge. See `getpointsandlines` command."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name of file to write. Suffix .vtk will be added to the file name."),
				CommandArgument<string>(ParameterDirection::In, "point data name", "Comma-separated list of names of data fields for each point.", ""),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "point data", "Image that has one row for each point. The data in columns will be saved in the .vtk file using names given in argument 'point data name'."),
				CommandArgument<string>(ParameterDirection::In, "line data name", "Comma-separated list of names of data fields for each line.", ""),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "cell data", "Image that has one row for each line. The data in columns will be saved in the .vtk file using names given in argument 'line data name'.")
			},
			"getpointsandlines")
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<float32_t>& points = *pop<Image<float32_t>*>(args);
			Image<uint64_t>& lines = *pop<Image<uint64_t>*>(args);
			string filename = pop<string>(args);

			string pointDataName = pop<string>(args);
			Image<float32_t>* pointData = pop<Image<float32_t>*>(args);

			string lineDataName = pop<string>(args);
			Image<float32_t>* lineData = pop<Image<float32_t>*>(args);

			writeVtk(points, lines, filename, pointDataName, pointData, lineDataName, lineData);
		}
	};

	class WriteVtk3Command : public WriteVtkCommandBase
	{
	protected:
		friend class CommandList;

		WriteVtk3Command() : WriteVtkCommandBase("writevtk", "Writes network in points-and-lines format to a .vtk file.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "points", "Image that stores coordinates of all points in the network. See `getpointsandlines` command."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::In, "lines", "Image that stores a list of point indices for each edge. See `getpointsandlines` command."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name of file to write. Suffix .vtk will be added to the file name."),
				CommandArgument<string>(ParameterDirection::In, "point data name", "Comma-separated list of names of data fields for each point.", ""),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "point data", "Image that has one row for each point. The data in columns will be saved in the .vtk file using names given in argument 'point data name'."),
			})
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<float32_t>& points = *pop<Image<float32_t>*>(args);
			Image<uint64_t>& lines = *pop<Image<uint64_t>*>(args);
			string filename = pop<string>(args);

			string pointDataName = pop<string>(args);
			Image<float32_t>* pointData = pop<Image<float32_t>*>(args);

			writeVtk(points, lines, filename, pointDataName, pointData);
		}
	};

	class WriteVtk2Command : public WriteVtkCommandBase
	{
	protected:
		friend class CommandList;

		WriteVtk2Command() : WriteVtkCommandBase("writevtk", "Writes network in points-and-lines format to a .vtk file.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "points", "Image that stores coordinates of all points in the network. See `getpointsandlines` command."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::In, "lines", "Image that stores a list of point indices for each edge. See `getpointsandlines` command."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name of file to write. Suffix .vtk will be added to the file name."),
			})
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<float32_t>& points = *pop<Image<float32_t>*>(args);
			Image<uint64_t>& lines = *pop<Image<uint64_t>*>(args);
			string filename = pop<string>(args);

			writeVtk(points, lines, filename);
		}
	};




	template<typename pixel_t> class FillSkeletonCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		FillSkeletonCommand() : Command("fillskeleton", "Fills a traced skeleton with one of values measured from the skeleton. E.g. fills each branch by its length.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::InOut, "skeleton", "Image containing the skeleton. Point belonging to the skeleton must be marked with non-zero value. The pixel data type should be able to store all values of the fill quantity, or otherwise clipping will occur. (Note that traceskeleton command clears the input image.)"),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "vertices", "Image where vertex coordinates are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<uint64_t> >(ParameterDirection::In, "edges", "Image where vertex indices corresponding to each edge are stored. See `tracelineskeleton` command."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "edge measurements", "Image that stores properties of each edge. See `tracelineskeleton` command"),
				CommandArgument<Image<int32_t> >(ParameterDirection::In, "edge points", "Image that stores some points on each edge. See `tracelineskeleton` command"),
				CommandArgument<size_t>(ParameterDirection::In, "fill index", "Defines column in the measurement image where fill values are grabbed from. Zero corresponds to first column in edge measurements, etc."),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "If processing a block of original image, origin of the block in coordinates of the full image. This parameter is used internally in distributed processing and should usually be set to its default value.", Distributor::BLOCK_ORIGIN_ARG_TYPE(0, 0, 0)),
			},
			tracedSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& vertices = *pop<Image<float32_t>* >(args);
			Image<uint64_t>& edges = *pop<Image<uint64_t>* >(args);
			Image<float32_t>& measurements = *pop<Image<float32_t>* >(args);
			Image<int32_t>& points = *pop<Image<int32_t>*>(args);
			size_t fillIndex = pop<size_t>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE origin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

			Network net;
			net.fromImage(vertices, edges, &measurements, &points);
			net.vertices -= Vec3f(origin);
			for (size_t n = 0; n < net.edges.size(); n++)
			{
				net.edges[n].properties.edgePoints -= Vec3sc(origin);
			}

			fillSkeleton(in, net, fillIndex);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<uint64_t>& edges = *std::get<DistributedImage<uint64_t>* >(args[2]);
			DistributedImage<float32_t>& measurements = *std::get<DistributedImage<float32_t>* >(args[3]);
			
			measurements.checkSize(Vec3c(3, edges.height(), 1));

			// Fill
			distributor.distribute(this, args);

			// Grow regions towards branch color until no changes occur.
			CommandList::get<GrowLabelsCommand<pixel_t> >().runDistributed(distributor, {&in, (double)std::numeric_limits<pixel_t>::max(), (double)0, Connectivity::AllNeighbours});

			return vector<string>();
		}

		virtual Vec3c getMargin(const vector<ParamVariant>& args) const override
		{
			// TODO: I'm not sure if this overlap is enough in all cases?
			return Vec3c(5, 5, 5);
		}

		virtual size_t getRefIndex(const vector<ParamVariant>& args) const override
		{
			// The input/output image is always the reference
			return 0;
		}

		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			// Load images 1-4 (network definition) always
			if (argIndex == 1)
			{
				DistributedImage<float32_t>& img = *std::get<DistributedImage<float32_t>*>(args[1]);
				readStart = Vec3c(0, 0, 0);
				readSize = img.dimensions();
			}
			else if (argIndex == 2)
			{
				DistributedImage<uint64_t>& img = *std::get<DistributedImage<uint64_t>*>(args[2]);
				readStart = Vec3c(0, 0, 0);
				readSize = img.dimensions();
			}
			else if (argIndex == 3)
			{
				DistributedImage<float32_t>& img = *std::get<DistributedImage<float32_t>*>(args[3]);
				readStart = Vec3c(0, 0, 0);
				readSize = img.dimensions();
			}
			else if (argIndex == 4)
			{
				DistributedImage<int32_t>& img = *std::get<DistributedImage<int32_t>*>(args[4]);
				readStart = Vec3c(0, 0, 0);
				readSize = img.dimensions();
			}
		}
	};


}
