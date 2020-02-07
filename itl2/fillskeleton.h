#pragma once

#include "image.h"
#include "network.h"
#include "traceskeleton.h"
#include "pointprocess.h"


namespace itl2
{
	const int EDGE_POINT_COUNT = 0;
	const int EDGE_LENGTH = 1;
	const int EDGE_AREA = 2;
	const int EDGE_X = 3;
	const int EDGE_Y = 4;
	const int EDGE_Z = 5;
	const int EDGE_CUSTOM = 6;

	/**
	Fills each branch of a skeleton with a measured value.
	Does not fill edges that start or end in the edge of the image.
	@param skeleton Image containing the skeleton as non-zero pixels.
	@param network Network structure returned from skeleton tracing.
	*/
	template<typename pixel_t> void fillSkeleton(Image<pixel_t>& skeleton, const Network& network, size_t propertyIndex, bool indicateProgress = true)
	{
		// Find intersection regions
		//internals::classifyForTracing<pixel_t>(skeleton, (pixel_t)internals::BRANCHING);
		internals::classifyForTracing<pixel_t>(skeleton);

		// Replace color internals::CURVE by numeric_limits<pixel_t>::max()
		replace<pixel_t>(skeleton, Vec2<pixel_t>(internals::CURVE, std::numeric_limits<pixel_t>::max()));

		// Replace branch points by background
		replace<pixel_t>(skeleton, Vec2<pixel_t>(internals::BRANCHING, 0));

		// Flood fill from each seed point.
		// There is one flood fill per edge, so the filled regions never overlap and we can do a multithreaded fill.
		bool hasInvalid = false;
		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)network.edges.size(); n++)
		{
			const Edge& edge = network.edges[n];
			pixel_t fillColor = edge.properties.get<pixel_t>(propertyIndex);
			for (const auto& startPoint : edge.properties.edgePoints)
			{
				if (skeleton.isInImage(startPoint) && !skeleton.isOnEdge(startPoint))
				{
					// In distributed processing we have zeros in the edges of each block, in
					// addition to the edges of the whole image.
					// Thus this test is requried not to fill background if startPoint happens to
					// be located in the block edge.
					if (skeleton(startPoint) != 0)
						floodfill(skeleton, Vec3c(startPoint), fillColor, fillColor, Connectivity::AllNeighbours);
				}
			}

			if (edge.properties.edgePoints.size() <= 0)
				hasInvalid = true;

			showThreadProgress(counter, network.edges.size(), indicateProgress);
		}
		
		// Inform user if there was a problem.
		if (hasInvalid)
			std::cout << "Warning: Some edges do not contain location of any point on the edge, and remain having value " << std::numeric_limits<pixel_t>::max() << std::endl;
	}

	namespace tests
	{
		void fillSkeleton();
	}
}
