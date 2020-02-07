
#include "regionremoval.h"
#include "io/raw.h"
#include "pointprocess.h"

namespace itl2
{
	namespace tests
	{
		void regionRemoval()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16(256, 256, 129);
			//raw::read(head16, "./input_data/t1-head_noisy_256x256x129.raw");
			raw::read(head16, "./input_data/t1-head_256x256x129.raw");

			//threshold(head16, 128);

			//regionRemoval(head16, 75);
			regionRemoval(head16, 1000);

			raw::writed(head16, "./regionremoval/result");

		}
	}
}
