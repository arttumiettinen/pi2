
#include <map>
#include <mutex>
#include <shared_mutex>

#include "dmap.h"
#include "io/raw.h"
#include "pointprocess.h"
#include "projections.h"
#include "neighbourhood.h"
#include "testutils.h"

using namespace std;

namespace itl2
{
	

	namespace tests
	{
		void dmap(const string& binFile, const string& dmapOutFile, const string& dmapGTFile)
		{
			Image<uint8_t> head;
			raw::read(head, binFile);

			//raw::writed(head, "./dmap/head_bin");

			Image<float32_t> dmap;

			distanceTransform(head, dmap);

			raw::writed(dmap, dmapOutFile);


			Image<float32_t> gt;
			raw::read(gt, dmapGTFile);

			checkDifference(dmap, gt, "Difference between ground truth and distance map.");
		}

		void dmap1()
		{
			dmap("./t1-head_bin_256x256x129.raw", "./dmap/head_dmap", "./t1-head_bin_dmap_256x256x129.raw");
			dmap("./test_piece_bin_256x256x256.raw", "./dmap/test_piece_result", "./test_piece_dmap_GT_256x256x256.raw");
			dmap("./test_piece_bin_512x512x512.raw", "./dmap/test_piece_result", "./test_piece_dmap_GT_512x512x512.raw");
		}
		
	}

}