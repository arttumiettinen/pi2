
#include <map>
#include <mutex>
#include <shared_mutex>

#include "dmap.h"
#include "io/raw.h"
#include "pointprocess.h"
#include "projections.h"
#include "neighbourhood.h"

using namespace std;

namespace itl2
{
	
	namespace tests
	{
		void dmap(const string& binFile, const string& dmapOutFile, const string& dmapGTFile, const string& dmapErrFile)
		{
			Image<uint8_t> head;
			raw::readd(head, binFile);

			//raw::writed(head, "./dmap/head_bin");

			Image<float32_t> dmap;

			distanceTransform(head, dmap);

			raw::writed(dmap, dmapOutFile);


			Image<float32_t> gt;
			raw::readd(gt, dmapGTFile);

			subtract(dmap, gt);
			abs(dmap);
			double err = sum(dmap);
			cout << "sum|dmap - gt| = " << err << endl;

			raw::writed(dmap, dmapErrFile);

			testAssert(NumberUtils<double>::equals(err, 0), "Difference between ground truth and distance map.");
		}

		void dmap1()
		{
			dmap("./t1-head_bin_256x256x129.raw", "./dmap/head_dmap", "./t1-head_bin_dmap_256x256x129.raw", "./dmap/head_dmap_err");
			dmap("./test_piece_bin_256x256x256.raw", "./dmap/test_piece_result", "./test_piece_dmap_GT_256x256x256.raw", "./dmap/test_piece_err");
		}

	}

}