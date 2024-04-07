
#include <map>

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
		template<typename dmap_t> void dmap(const string& binFile, const string& dmapOutFile, const string& dmapGTFile, double tolerance)
		{
			Image<uint8_t> head;
			raw::read(head, binFile);

			//raw::writed(head, "./dmap/head_bin");

			Image<dmap_t> dmap;

			distanceTransform(head, dmap);

			raw::writed(dmap, dmapOutFile);


			Image<float32_t> gt;
			raw::read(gt, dmapGTFile);

			checkDifference(dmap, gt, "Difference between ground truth and distance map. (" + toString(imageDataType<dmap_t>()) + ")", tolerance);
		}

		void dmap1()
		{
			dmap<float32_t>("../test_input_data/t1-head_bin_256x256x129.raw", "./dmap/head_dmap", "../test_input_data/t1-head_bin_dmap_256x256x129.raw", 1e-6);
			dmap<float32_t>("../test_input_data/test_piece_bin_256x256x256.raw", "./dmap/test_piece_result", "../test_input_data/test_piece_dmap_GT_256x256x256.raw", 1e-6);
			dmap<float32_t>("../test_input_data/test_piece_bin_512x512x512.raw", "./dmap/test_piece_result", "../test_input_data/test_piece_dmap_GT_512x512x512.raw", 1e-6);

			dmap<int32_t>("../test_input_data/t1-head_bin_256x256x129.raw", "./dmap/head_dmap", "../test_input_data/t1-head_bin_dmap_256x256x129.raw", 0.5);
			dmap<int32_t>("../test_input_data/test_piece_bin_256x256x256.raw", "./dmap/test_piece_result", "../test_input_data/test_piece_dmap_GT_256x256x256.raw", 0.5);
			dmap<int32_t>("../test_input_data/test_piece_bin_512x512x512.raw", "./dmap/test_piece_result", "../test_input_data/test_piece_dmap_GT_512x512x512.raw", 0.5);

			dmap<uint32_t>("../test_input_data/t1-head_bin_256x256x129.raw", "./dmap/head_dmap", "../test_input_data/t1-head_bin_dmap_256x256x129.raw", 0.5);
			dmap<uint32_t>("../test_input_data/test_piece_bin_256x256x256.raw", "./dmap/test_piece_result", "../test_input_data/test_piece_dmap_GT_256x256x256.raw", 0.5);
			dmap<uint32_t>("../test_input_data/test_piece_bin_512x512x512.raw", "./dmap/test_piece_result", "../test_input_data/test_piece_dmap_GT_512x512x512.raw", 0.5);
		}
		
	}

}