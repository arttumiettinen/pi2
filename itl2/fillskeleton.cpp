
#include "fillskeleton.h"

#include "io/raw.h"
#include "traceskeleton.h"
#include "testutils.h"

namespace itl2
{

	namespace tests
	{
		void fillSkeleton()
		{
			Image<uint8_t> skele, skele2, orig;
			raw::read(skele, "./input_data/real_skele_200x200x200.raw");
			Network net;
			traceLineSkeleton(skele, false, 1.0, 1.0, net);

			Image<float32_t> verts;
			Image<size_t> edges;
			Image<float32_t> meas;
			Image<int32_t> points;
			net.toImage(verts, edges, &meas, &points);

			raw::read(skele, "./input_data/real_skele_200x200x200.raw");
			itl2::fillSkeleton(skele, net, EDGE_LENGTH);
			raw::writed(skele, "./fillskeleton/filled_length");

			Network net2;
			net2.fromImage(verts, edges, &meas, &points);

			raw::read(skele2, "./input_data/real_skele_200x200x200.raw");
			itl2::fillSkeleton(skele2, net2, EDGE_LENGTH);
			raw::writed(skele2, "./fillskeleton/filled_length_comp");

			checkDifference(skele, skele2, "Filled skeleton differs before and after conversion of the network to image.");
		}
	}
}