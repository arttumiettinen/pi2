
#include "projections.h"
#include "pointprocess.h"
#include "io/raw.h"

namespace itl2
{

	namespace tests
	{
		void projections()
		{

			Image<uint8_t> img(30, 30, 30);

			setValue(img, 2);
			testAssert(sum(img) == 2 * 30 * 30 * 30, "sum");
			testAssert(mean(img) == 2, "mean");

			img(10, 10, 10) = 1;
			testAssert(min(img) == 1, "min");

			img(20, 20, 20) = 10;
			testAssert(max(img) == 10, "max");


		}

		/*
		Checks test result by comparing to ground truth image, used in projections2 test method.
		*/
		void checkResult(const Image<float32_t>& result, const string& file, const string& opname)
		{
			Vec3c dimensions;
			ImageDataType dt;

			raw::getInfo(file, dimensions, dt);

			Image<float32_t> gt(dimensions);

			gt.checkSize(result);

			raw::read(gt, file);

			subtract(gt, result);
			abs(gt);
			testAssert(max(gt) < 1e-3, opname);
		}

		void projections2()
		{

			Image<uint16_t> head(256, 256, 129);
			raw::read(head, "t1-head_256x256x129.raw");

			Image<float32_t> maxy;
			max(head, 1, maxy, false);
			raw::writed(maxy, "./projections/maxy");

			Image<float32_t> zprojection;
			Image<float32_t> yprojection;
			Image<float32_t> xprojection;
			sum(head, 2, zprojection);
			sum(head, 1, yprojection);
			sum(head, 0, xprojection);
			raw::writed(zprojection, "./projections/zproj_sum");
			raw::writed(yprojection, "./projections/yproj_sum");
			raw::writed(xprojection, "./projections/xproj_sum");

			checkResult(zprojection, "./zproj_sum_true_256x256x1.raw", "z projection");
			checkResult(yprojection, "./yproj_sum_true_256x129x1.raw", "y projection");
			checkResult(xprojection, "./xproj_sum_true_129x256x1.raw", "x projection");

			double count;
			double sum = maskedsum(zprojection, (float32_t)0, count);
			double sum2 = maskedsum(zprojection, numeric_limits<float32_t>::signaling_NaN(), count);

			Image<float32_t> minproj, maxproj, meanproj;
			min(head, 2, minproj);
			max(head, 2, maxproj);
			mean(head, 2, meanproj);

			raw::writed(minproj, "./projections/zproj_min");
			raw::writed(maxproj, "./projections/zproj_max");
			raw::writed(meanproj, "./projections/zproj_mean");

			checkResult(minproj, "./zproj_min_true_256x256x1.raw", "min projection");
			checkResult(maxproj, "./zproj_max_true_256x256x1.raw", "max projection");
			checkResult(meanproj, "./zproj_mean_true_256x256x1.raw", "mean projection");

		}
	}
}