#pragma once

#include "datatypes.h"
#include "image.h"
#include "neighbourhood.h"
#include "progress.h"

namespace itl2
{
	namespace internals
	{
		/**
		Retrieve or calculate tables for Danielsson's algorithm, up to maximum squared radius R2max.
		*/
		void getDanielssonTables(std::vector<coord_t>& table1, std::vector<coord_t>& table2, std::vector<coord_t>& table3, coord_t R2max);
	}

	/**
	Calculates centers of locally maximal disks using Danielsson algorithm.
	The output image can be interpreted as medial axis or distance ridge. Drawing a sphere on
	each nonzero pixel, with radius and color equal to pixel value, results in thickness map of the structure
	(assuming that larger colors replace smaller ones).
	See e.g.
	Yaorong Ge and J. Michael Fitzpatrick - On the Generation of Skeletons from Discrete Euclidean Distance Maps
	@param dmap2 Squared Euclidean distance map of the input geometry.
	@param out Distance ridge. Squared distance values of pixels that correspond to centers of locally maximal disks are set to this image. Other values are not changed. Usually this image should be empty before calling this function.
	*/
	template<typename pixel_t> void centersOfLocallyMaximalSpheres(const Image<pixel_t>& dmap2, Image<pixel_t>& out)
	{
		out.mustNotBe(dmap2);
		out.ensureSize(dmap2);


		// Find maximum squared distance value
		//cout << "Finding maximum radius..." << endl;
		coord_t maxr2 = (coord_t)max(dmap2);

		// Fill lookup table for three types of neighbours.
		//cout << "Calculating Danielsson tables (cached)..." << endl;
		std::vector<coord_t> table1;
		std::vector<coord_t> table2;
		std::vector<coord_t> table3;

		internals::getDanielssonTables(table1, table2, table3, maxr2);

		//cout << "Finding centers of locally maximal disks..." << endl;
		ProgressIndicator progress(dmap2.depth());
		#pragma omp parallel if(dmap2.pixelCount() >= PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			Image<coord_t> nb(3, 3, 3);

			// Process all points in the image
			#pragma omp for schedule(dynamic)
			for (coord_t z = 0; z < dmap2.depth(); z++)
			{
				for (coord_t y = 0; y < dmap2.height(); y++)
				{
					for (coord_t x = 0; x < dmap2.width(); x++)
					{
						coord_t c = (coord_t)dmap2(x, y, z);
						if (c != 0)
						{
							getNeighbourhood<pixel_t, coord_t>(dmap2, Vec3c(x, y, z), Vec3c(1, 1, 1), nb, BoundaryCondition::Zero);

							// Check all neighbours
							if (!(// 6-neighbours, one coordinate changes by one pixel.
								table1[nb(0, 1, 1)] >= c ||
								table1[nb(2, 1, 1)] >= c ||
								table1[nb(1, 0, 1)] >= c ||
								table1[nb(1, 2, 1)] >= c ||
								table1[nb(1, 1, 0)] >= c ||
								table1[nb(1, 1, 2)] >= c ||
								// 18-neighbours but not 6-neighbours, two coordinates change by one pixel.
								table2[nb(0, 0, 1)] >= c ||
								table2[nb(0, 2, 1)] >= c ||
								table2[nb(2, 0, 1)] >= c ||
								table2[nb(2, 2, 1)] >= c ||
								table2[nb(1, 0, 0)] >= c ||
								table2[nb(1, 0, 2)] >= c ||
								table2[nb(1, 2, 0)] >= c ||
								table2[nb(1, 2, 2)] >= c ||
								table2[nb(0, 1, 0)] >= c ||
								table2[nb(0, 1, 2)] >= c ||
								table2[nb(2, 1, 0)] >= c ||
								table2[nb(2, 1, 2)] >= c ||
								// Corners, three coordinates change by one pixel.
								table3[nb(0, 0, 0)] >= c ||
								table3[nb(0, 2, 0)] >= c ||
								table3[nb(2, 0, 0)] >= c ||
								table3[nb(2, 2, 0)] >= c ||
								table3[nb(0, 0, 2)] >= c ||
								table3[nb(0, 2, 2)] >= c ||
								table3[nb(2, 0, 2)] >= c ||
								table3[nb(2, 2, 2)] >= c
								))
							{
								// This is center of locally maximal sphere
								out(x, y, z) = dmap2(x, y, z);
							}
							else
							{
								// Not a center of locally maximal sphere
								out(x, y, z) = 0;
							}
						}
					}
				}

				progress.step();
			}
		}
	}
	

	namespace tests
	{
		void danielssonTableSpeedTest();
		void fullDanielssonTables();
		void danielsson();
	}
}
