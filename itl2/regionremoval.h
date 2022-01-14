#pragma once

#include "image.h"
#include "floodfill.h"
#include "particleanalysis.h"

namespace itl2
{

	/**
	Removes all nonzero regions smaller than the given volume limit.
	@param img Image to process.
	@param volumeLimit Nonzero regions smaller than this value are removed.
	@param preserveEdges Set to true to skip processing of regions that touch image edge.
	@param connectivity Connectivity of pixels.
	*/
	template<typename pixel_t> void regionRemoval(Image<pixel_t>& img, size_t volumeLimit, bool preserveEdges = false, Connectivity connectivity = Connectivity::NearestNeighbours, bool multiThreaded = true)
	{
		// Analyze particles
		Results results;
		AnalyzerSet<Vec3sc, pixel_t> analyzers;
		analyzers.push_back(new analyzers::Coordinates3D<Vec3sc, pixel_t>());
		analyzers.push_back(new analyzers::Volume<Vec3sc, pixel_t>());
		if (preserveEdges)
			analyzers.push_back(new analyzers::IsOnEdge<Vec3sc, pixel_t>(img.dimensions()));

		std::cout << "Searching for particles..." << std::endl;
		if(multiThreaded)
			analyzeParticles<pixel_t>(img, analyzers, results, connectivity, volumeLimit);
		else
			analyzeParticlesSingleThreaded<pixel_t>(img, analyzers, results, connectivity, volumeLimit);
		
		if (preserveEdges)
		{
			// Remove all particles that touch edges.
			size_t n = 0;
			while (n < results.size())
			{
				if (results[n][results[n].size() - 1] != 0)
				{
					results.erase(results.begin() + n);
					n--;
				}
				n++;
			}
		}

		// Fill those with too small volume
		std::cout << "Filling small particles..." << std::endl;
		//fillParticles<pixel_t>(img, results, 0, 3, (double)volumeLimit, connectivity, false);
		fillParticles<pixel_t>(img, results, 0, connectivity);

		threshold(img, (pixel_t)0);
	}


	namespace tests
	{
		void regionRemoval();
	}
}
