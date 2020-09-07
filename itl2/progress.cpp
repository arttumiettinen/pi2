
#include "progress.h"

namespace itl2
{
	namespace tests
	{
		void progress()
		{
			{
				ProgressIndicator prog(5);
				for (size_t n = 0; n < 5; n++)
					prog.step();
			}

			{
				ProgressIndicator prog(450 * 450 * 450);
				for (size_t n = 0; n < 450 * 450 * 450; n++)
					prog.step();
			}

			float MAX = 2235450000.0f;
			{
				ProgressIndicator prog(MAX);
				for (float n = 0; n < 3 * MAX / 4; n += 173)
					prog.step();
			}
			std::cout << "this message interrupts the progress indicator..." << std::endl;
		}
	}
}