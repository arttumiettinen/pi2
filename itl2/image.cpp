
#include "image.h"

#include "io/raw.h"
#include "dmap.h"
#include "pointprocess.h"
#include "projections.h"

#include <iostream>
using namespace std;

namespace itl2
{
	namespace tests
	{
		void image()
		{
			Image<uint8_t> img(10, 10);

			testAssert(img.dimensionality() == 2, "Dimensionality");
			testAssert(img.pixelCount() == 10 * 10, "Pixel count");

		}

		void buffers()
		{
			// Generate distance map using memory mapped file
			{
				Image<uint8_t> head(256, 256, 129);
				raw::read(head, "./input_data/t1-head_bin_256x256x129.raw");

				Image<float32_t> dmap("./buffers/head_dmap_256x256x129.raw", false, 256, 256, 129);

				distanceTransform(head, dmap);

				Image<float32_t> gt;
				raw::read(gt, "./input_data/t1-head_bin_dmap_256x256x129.raw");

				subtract(gt, dmap);
				abs(gt);
				double err = sum(gt);
				cout << "sum|dmap - gt| = " << err << endl;
			}


			// Check that the memory mapping works for reading, too.
			{
				Image<float32_t> dmap("./buffers/head_dmap_256x256x129.raw", true, 256, 256, 129);

				Image<float32_t> gt;
				raw::read(gt, "./input_data/t1-head_bin_dmap_256x256x129.raw");

				subtract(gt, dmap);
				abs(gt);
				double err = sum(gt);
				cout << "sum|dmap - gt| = " << err << endl;
			}
		}
	}
}