
#include <iostream>

#include "histogram.h"
#include "math/numberutils.h"
#include "io/raw.h"
#include "test.h"
#include "projections.h"
#include "testutils.h"

using namespace std;


namespace itl2
{

	namespace tests
	{
		void histogramIntermediateType()
		{
			typeAssert<histogram_intermediate_type<uint8_t, uint8_t>::type, uint64_t>();
			typeAssert<histogram_intermediate_type<uint8_t, float32_t>::type, double>();
			typeAssert<histogram_intermediate_type<float32_t, float32_t>::type, double>();
			typeAssert<histogram_intermediate_type<float32_t, uint8_t>::type, double>();
			typeAssert<histogram_intermediate_type<uint64_t, int32_t>::type, uint64_t>();
		}

		void histogram()
		{
			Image<uint16_t> head(256, 256, 129);
			raw::read(head, "../test_input_data/t1-head_256x256x129.raw");

			Image<int> hist(100);

			histogram(head, hist, Vec2d(0, 1000));

			// Read ground truth file

			ifstream file("../test_input_data/t1-head_histogram_0_1000_100.txt");
			string str;
			int n = 0;
			getline(file, str); // Reads header line
			while (getline(file, str))
			{
				testAssert(n <= hist.pixelCount(), "Not enough bins.");

				stringstream ss;
				ss << str;

				string part1, part2;
				getline(ss, part1, '\t');
				getline(ss, part2);

				double val = fromString<double>(part2);

				//cout << part1 << "\t" << part2 << "\t" << hist(n) << endl;

				testAssert(NumberUtils<double>::equals(val, hist(n)), string("bin ") + toString(n) + " value: true = " + toString(val) + ", measured = " + toString(hist(n)));

				n++;
			}
			testAssert(n == hist.pixelCount(), "Bin count");
		}

		void histogram2d()
		{
			Image<uint16_t> head(256, 256, 129);
			raw::read(head, "../test_input_data/t1-head_256x256x129.raw");

			Image<uint8_t> headbin(256, 256, 129+10);
			raw::read(headbin, "../test_input_data/t1-head_bin_256x256x129.raw");

			Image<float32_t> hist(100, 5);

			//itl2::histogram({ head, headbin }, hist, { Vec2d(0, 1000), Vec2d(0, 2) });
			//itl2::histogram(tuple<Image<uint16_t>&, Image<uint8_t>&> { head, headbin });
			//itl2::histogram(make_tuple(&head, &headbin), hist, { Vec2d(0, 1000), Vec2d(0, 2) });
			//ImageAndRange r2 = { head, Vec2d(0, 1000) };
			//itl2::histogram(make_tuple(ImageAndRange( head, Vec2d(0, 1000), 100 ), ImageAndRange( headbin, Vec2d(0, 2), 5 )), hist);
			
			itl2::multiHistogram(hist, 0, ImageAndRange(head, Vec2d(0, 1000), 100), ImageAndRange(headbin, Vec2d(0, 2), 5));
			raw::writed(hist, "./histogram/bivariate");

			Image<float32_t> hist2(hist.dimensions());
			Image<int8_t> weight(head.dimensions());
			setValue(weight, 1);
			
			itl2::multiHistogram(hist2, 0, &weight, ImageAndRange(head, Vec2d(0, 1000), 100), ImageAndRange(headbin, Vec2d(0, 2), 5));
			raw::writed(hist2, "./histogram/bivariate2");

			testAssert(itl2::equals(hist, hist2), "Weighted and non weighted bivariate histogram.");
		}
	}

}