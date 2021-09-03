
#include "filters.h"
#include "io/raw.h"
#include "conversions.h"
#include "projections.h"
#include "generation.h"

namespace itl2
{
	namespace tests
	{
		void gaussFilters()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16(256, 256, 129);
			raw::read(head16, "./input_data/t1-head_noisy_256x256x129.raw");
			divide(head16, 4);
			raw::writed(head16, "./filters/gauss_orig");

			Image<uint8_t> head8(head16.dimensions());
			convert(head16, head8);

			Image<float32_t> head32(head16.dimensions());
			convert(head16, head32);


			Image<uint8_t> headg8;
			Image<uint16_t> headg16;
			Image<float32_t> headg32;

			
			gaussFilter(head8, headg8, 2);
			raw::writed(headg8, "./filters/gauss8");
			
			//gaussDerivative(head8, headg8, 2, 0, 1);
			//raw::writed(headg8, "./filters/gauss8_dx");
			//gaussDerivative(head8, headg8, 2, 1, 1);
			//raw::writed(headg8, "./filters/gauss8_dy");
			//gaussDerivative(head8, headg8, 2, 2, 1);
			//raw::writed(headg8, "./filters/gauss8_dz");
			
			
			gaussFilter(head16, headg16, 2);
			raw::writed(headg16, "./filters/gauss16");

			//gaussDerivative(head16, headg16, 2, 0, 1);
			//raw::writed(headg16, "./filters/gauss16_dx");
			//gaussDerivative(head16, headg16, 2, 1, 1);
			//raw::writed(headg16, "./filters/gauss16_dy");
			//gaussDerivative(head16, headg16, 2, 2, 1);
			//raw::writed(headg16, "./filters/gauss16_dz");
			

			gaussFilter(head32, headg32, 2);
			raw::writed(headg32, "./filters/gauss32");
			gaussDerivative(head32, headg32, 2, 0, 1);
			raw::writed(headg32, "./filters/gauss32_dx");
			gaussDerivative(head32, headg32, 2, 1, 1);
			raw::writed(headg32, "./filters/gauss32_dy");
			gaussDerivative(head32, headg32, 2, 2, 1);
			raw::writed(headg32, "./filters/gauss32_dz");
			
			gaussFilter(head32, headg32, 2, false);
			raw::writed(headg32, "./filters/gauss32_no_opt");
		}

		void separableOptimization()
		{
			Image<uint16_t> head(256, 256, 129);
			raw::read(head, "./input_data/t1-head_noisy_256x256x129.raw");

			Image<float32_t> head32(head.dimensions());
			convert(head, head32);

			Image<uint16_t> filtered;
			Image<float32_t> filtered2;


			// This is rect mean filtering without optimization
			meanFilter(head, filtered, 3, NeighbourhoodType::Rectangular);
			raw::writed(filtered, "./filters/mean_rect_filtered");

			// This should use separable optimization (Rectangular neighbourhood and float image)
			meanFilter(head32, filtered2, 3, NeighbourhoodType::Rectangular);
			raw::writed(filtered2, "./filters/mean_rect_filtered_separable");

			// The results should agree up to rounding error
			Image<float32_t> ad;
			convert(filtered, ad);
			subtract(ad, filtered2);
			abs(ad);
			testAssert(max(ad) <= 0.5, "separable optimization (mean filter)");
		}

		void stddevuint16()
		{
			// NOTE: This test shows numerical inaccuracy in the right half of the 1000 pixel wide image.
			// There, the true result is 2.5, but for pixel values near 1000, the calculations might result in
			// 2.49, 2.51 etc. that will be rounded to 2 and 3 --> stripes in the integer image.
			using DT = uint16_t;
			Image<DT> img(1000, 1000);
			ramp(img, 0);
			draw(img, Sphere(Vec3f(50, 50, 0), 25.0f), (DT)10000);
			raw::writed(img, "./filters/stddev_orig");	
			Image<DT> filtered(img.dimensions());
			stddevFilter(img, filtered, 5);
			raw::writed(filtered, "./filters/stddev_filtered");
		}

		void filters()
		{
			// NOTE: No asserts!

			Image<uint16_t> head(256, 256, 129);
			raw::read(head, "./input_data/t1-head_noisy_256x256x129.raw");

			Image<float32_t> head32(head.dimensions());
			convert(head, head32);

			Image<uint16_t> filtered;
			Image<float32_t> filtered2;

			medianFilter(head, filtered, 2);
			raw::writed(filtered, "./filters/median_filtered");

			meanFilter(head, filtered, 2);
			raw::writed(filtered, "./filters/mean_filtered");

			varianceFilter(head, filtered, 2);
			raw::writed(filtered, "./filters/variance_filtered");

			stddevFilter(head, filtered, 2);
			raw::writed(filtered, "./filters/stddev_filtered");

			minFilter(head, filtered, 2);
			raw::writed(filtered, "./filters/min_filtered");

			// This should use separable optimization
			minFilter(head, filtered, 2, NeighbourhoodType::Rectangular);
			raw::writed(filtered, "./filters/min_rect_filtered");

			maxFilter(head, filtered, 2);
			raw::writed(filtered, "./filters/max_filtered");

			vaweFilter(head, filtered, 2, 25);
			raw::writed(filtered, "./filters/vawe_filtered");

			//bilateralFilter(head, filtered, 2, 25);
			//raw::writed(filtered, "./filters/bilateral_filtered");

			//minFilter(head32, filtered, 2);
			//min(filtered, head32, 2);

			Image<float32_t> tr(100, 100);
			raw::read(tr, "./input_data/test_rect_100x100x1.raw");

			meanFilter(tr, filtered2, Vec3c(10, 10, 0));
			raw::writed(filtered2, "./filters/test_rect_rect_filtered_xy");
		}

		void bilateral()
		{
			// NOTE: No asserts!

			Image<uint16_t> head(256, 256, 129);
			Image<uint16_t> filtered;
			raw::read(head, "./input_data/t1-head_noisy_256x256x129.raw");
			bilateralFilter(head, filtered, 2, 25);
			raw::writed(filtered, "./filters/bilateral_filtered2");
		}
	}
}