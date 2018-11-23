
#include <cmath>
#include "pointprocess.h"
#include "projections.h"

namespace itl2
{
	namespace tests
	{
		void pointProcess()
		{
			Image<uint8_t> img(10, 10);

			setValue(img, 5);
			testAssert(sum(img) == 5 * 10 * 10, "sum after set");

			add(img, 1);
			testAssert(sum(img) == 6 * 10 * 10, "sum after add");

			subtract(img, 2);
			testAssert(sum(img) == 4 * 10 * 10, "sum after sub");

			multiply(img, 2);
			testAssert(sum(img) == 8 * 10 * 10, "sum after multiply");

			divide(img, 4);
			testAssert(sum(img) == 2 * 10 * 10, "sum after divide");

			setValue(img, 1);
			add(img, img);
			testAssert(sum(img) == 2 * 10 * 10, "sum after add img");

			multiply(img, img);
			testAssert(sum(img) == 4 * 10 * 10, "sum after multiply img");

			subtract(img, img);
			testAssert(sum(img) == 0 * 10 * 10, "sum after sub img");

			setValue(img, 7);
			divide(img, img);
			testAssert(sum(img) == 1 * 10 * 10, "sum after div img");

			setValue(img, 20);
			testAssert(sum(img) == 20 * 10 * 10, "sum before negate");
			negate(img);
			testAssert(sum(img) == 0 * 10 * 10, "sum after negate");

			setValue(img, 12);
			threshold(img, 7);
			testAssert(sum(img) == 1 * 10 * 10, "sum after threshold");

			threshold(img, 7);
			testAssert(sum(img) == 0 * 10 * 10, "sum after threshold 2");
		}

		void pointProcessComplex()
		{
			Image<complex32_t> img(10, 10);

			setValue(img, 5.0f + 3if);
			testAssert(img(1, 1) == 5.0f + 3.0if, "set");

			negate(img);
			testAssert(img(1, 1) == -5.0f - 3.0if, "negate");

			conjugate(img);
			testAssert(img(1, 1) == -5.0f + 3.0if, "conjugate");

			real(img);
			testAssert(img(1, 1) == -5.0f + 0.0if, "real");

			setValue(img, 5.0f + 3.0if);
			imag(img);
			testAssert(img(1, 1) == 3.0f, "imag");

			setValue(img, 5.0f + 3.0if);
			abs(img);
			testAssert(img(1, 1) == sqrt(5.0f*5.0f+3.0f*3.0f), "abs");

			setValue(img, 5.0f + 3.0if);
			normSquared(img);
			cout << img(1, 1).real() << " = " << (5.0f*5.0f + 3.0f*3.0f) << endl;
			//testAssert(NumberUtils<float32_t>::equals(img(1, 1).real(), 5.0f*5.0f + 3.0f*3.0f), "norm");
			testAssert(::abs(img(1, 1).real() - (5.0f*5.0f + 3.0f*3.0f)) < 0.001, "norm");


			setValue(img, 5.0f + 0.0if);
			add(img, 0.0f + 5.0if);
			testAssert(img(1, 1) == (5.0f + 5.0if), "add constant");

			add(img, img);
			testAssert(img(1, 1) == (10.0f + 10.0if), "add image");

			divide(img, 2.0f + 0.0if);
			testAssert(img(1, 1) == (5.0f + 5.0if), "add image");

		}
	}
}
