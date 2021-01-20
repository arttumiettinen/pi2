
#include "eval.h"
#include "testutils.h"

namespace itl2
{
	namespace tests
	{
		void eval()
		{
			Image<uint8_t> out(100, 200, 300);
			Image<float32_t> param1(100, 200, 300);
			Image<uint64_t> param2(100, 200, 300);

			float32_t x0 = 7.45f;
			uint64_t x1 = 76;
			uint8_t answer = pixelRound<uint8_t>((double)x0 + 2 * (double)x1);

			Image<uint8_t> ans(100, 200, 300);
			setValue(ans, answer);
			setValue(param1, x0);
			setValue(param2, x1);

			itl2::eval("x0 + 2 * x1", out, std::vector<ImageBase*>{&param1, & param2});
			checkDifference(out, ans, "eval result");
		}
	}
}