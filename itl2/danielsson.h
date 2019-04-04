#pragma once

#include "datatypes.h"
#include "image.h"

namespace itl2
{
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
	void centersOfLocallyMaximalDisks(const Image<int32_t>& dmap2, Image<int32_t>& out);
	

	namespace tests
	{
		void danielssonTableSpeedTest();
		void fullDanielssonTables();
		void danielsson();
	}
}
