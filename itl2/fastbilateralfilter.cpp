
#include "fastbilateralfilter.h"
#include "io/raw.h"
#include "transform.h"

#include <iostream>
using namespace std;

namespace itl2
{
	namespace tests
	{
		//void fastBilateralGaussPolynomial()
		//{
		//	Image<uint16_t> head16(256, 256, 129);
		//	raw::read(head16, "../test_input_data/t1-head_noisy_256x256x129.raw");

		//	// This is needed if test below is enabled.
		//	//Image<uint16_t> filtered;
		//	//Vec3c dims;
		//	//ImageDataType dt;
		//	//const string slowName = "./fast_bilateral/slow";
		//	//if (!raw::getInfo(slowName, dims, dt))
		//	//{
		//	//	bilateralFilter(head16, filtered, 2, 25);
		//	//	raw::writed(filtered, slowName);
		//	//}
		//	//else
		//	//{
		//	//	raw::read(filtered, slowName);
		//	//}

		//	bilateralFilterGaussPolynomial(head16, 2, 25, 250, -1);
		//	raw::writed(head16, "./fast_bilateral/fast_gauss_polynomial");

		//	// This test will never succeed as the approximation is not good enough.
		//	//testAssert(equals(head16, filtered), "Difference between slow and fast bilateral filtering.");
		//}

		void fastBilateralSampling()
		{
			Image<uint16_t> head16Full(256, 256, 129);
			raw::read(head16Full, "../test_input_data/t1-head_noisy_256x256x129.raw");

			Image<uint16_t> head16(256, 256, 10);
			crop(head16Full, head16, Vec3c(0, 0, 60));

			float32_t spatialSigma = 3;
			float32_t rangeSigma = 150;

			//{
			//	Image<uint16_t> filtered;
			//	Timer timer;
			//	timer.start();
			//	bilateralFilter(head16, filtered, spatialSigma, rangeSigma);
			//	timer.stop();
			//	cout << "Slow: " << timer.getSeconds() << " s." << endl;

			//	raw::writed(filtered, "./fast_bilateral/slow");
			//}

			{
				Image<uint16_t> filtered;
				Timer timer;
				timer.start();
				//bilateralFilterSampling(head16, filtered, 2, 25);
				bilateralFilterSampling(head16, filtered, spatialSigma, rangeSigma);
				timer.stop();
				cout << "Sampling: " << timer.getSeconds() << " s." << endl;

				raw::writed(filtered, "./fast_bilateral/fast_sampling");
			}
		}
	}
}