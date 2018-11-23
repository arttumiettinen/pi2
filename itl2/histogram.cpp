
#include <iostream>

#include "histogram.h"
#include "math/numberutils.h"
#include "io/raw.h"
#include "test.h"

using namespace std;
using namespace math;

namespace itl2
{

	namespace tests
	{
		void histogram()
		{
			Image<uint16_t> head(256, 256, 129);
			raw::read(head, "./t1-head_256x256x129.raw");

			Image<int> hist(100);

			histogram(head, hist, Vec2d(0, 1000));

			// Read ground truth file

			ifstream file("./t1-head_histogram_0_1000_100.txt");
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
	}

}