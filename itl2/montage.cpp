
#include "montage.h"
#include "io/raw.h"

namespace itl2
{
	namespace tests
	{
		void montage()
		{
			Image<uint16_t> head;
			raw::read(head, "../test_input_data/t1-head");

			Image<uint16_t> montage;

			itl2::montage(head, montage, 6, 5, 0.75, 0, (size_t)head.depth(), 0, 2, (uint16_t)500);

			raw::writed(montage, "./montage/montage");
		}
	}
}