
#include "progress.h"

namespace itl2
{
	namespace tests
	{
		void progress()
		{
			float MAX = 2235450000.0f;
			{
				ProgressIndicator prog(MAX);
				for (float n = 0; n < 3 * MAX / 4; n += 173)
					prog.Show(n);
			}
			std::cout << "this message interrupts the progress indicator..." << std::endl;
		}
	}
}