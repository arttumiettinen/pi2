#pragma once

#include "image.h"
#include "test.h"
#include "conversions.h"
#include "pointprocess.h"
#include "math/numberutils.h"
#include "projections.h"

#include <string>
#include <iostream>


namespace itl2
{
	template<typename T1, typename T2> void typeAssert()
	{
		std::string t1 = typeid(T1).name();
		std::string t2 = typeid(T2).name();
		std::cout << t1 << ", " << t2 << std::endl;
		testAssert(std::is_same_v<T1, T2>, std::string("Got ") + t1 + ", should be " + t2);
	}

	template<typename pixel1_t, typename pixel2_t> void checkDifference(const Image<pixel1_t>& a, const Image<pixel2_t>& b, const std::string& errorMessage, double tolerance = NumberUtils<double>::tolerance())
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
		double maxDiff = max(af);

		//cout << "|a-b| = " << diff << endl;

		testAssert(NumberUtils<double>::equals(maxDiff, 0.0, tolerance), errorMessage + " Max difference = " + toString(maxDiff));
	}
}
