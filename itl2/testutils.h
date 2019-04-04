#pragma once

#include "image.h"
#include "test.h"
#include "conversions.h"
#include "pointprocess.h"
#include "math/numberutils.h"
#include "projections.h"

#include <string>

using std::string;

namespace itl2
{
	template<typename pixel1_t, typename pixel2_t> void checkDifference(const Image<pixel1_t>& a, const Image<pixel2_t>& b, const string& errorMessage)
	{
		//a.checkSize(b);
		if (!testAssert(a.sizeEquals(b), errorMessage))
			return;

		Image<float32_t> af(a.dimensions());
		convert(a, af);

		Image<float32_t> bf(b.dimensions());
		convert(b, bf);

		subtract(af, bf);
		abs(af);
		double diff = sum(af);

		//cout << "|a-b| = " << diff << endl;

		testAssert(NumberUtils<double>::equals(diff, 0.0), errorMessage);
	}
}