
#include "minhash.h"

using namespace std;

namespace itl2
{

	namespace tests
	{
		void hashPow()
		{
			using hash_t = uint8_t;
			for (size_t n = 0; n < sizeof(hash_t); n++)
			{
				hash_t res = internals::hashPow<hash_t>(n);
				double tru = std::pow(2, n);
				//cout << "2^" << n << " = " << res << " = " << tru << endl;
				testAssert(res == tru, "hashPow");
			}
		}

		void nbHash()
		{
			using hash_t = uint8_t;
			const size_t hashSize = sizeof(hash_t) * 8;
			Image<uint8_t> nb(hashSize);

			// Empty neighbourhood should result in 0 hash
			for (size_t m = 0; m < hashSize; m++)
				nb(m) = 0;

			testAssert(internals::nbhash<uint8_t, hash_t>(nb) == 0, "empty nb");

			// Full-1 neighbourhood should result in maximum value of hash type
			for (size_t m = 0; m < hashSize; m++)
				nb(m) = 1;

			testAssert(internals::nbhash<uint8_t, hash_t>(nb) == numeric_limits<hash_t>::max(), "full nb");

			// Single 1 in the neighbourhood should result in the result of single call to hashPow
			for (size_t n = 0; n < hashSize; n++)
			{
				for (size_t m = 0; m < hashSize; m++)
					nb(m) = 0;

				nb(n) = 1;

				hash_t res = internals::nbhash<uint8_t, hash_t>(nb);
				hash_t tru = internals::hashPow<hash_t>(n);

				testAssert(res == tru, "single 1 hash");
			}

		}

		Image<uint8_t>& put(uint8_t nb[], Image<uint8_t>& img)
		{
			for (coord_t n = 0; n < img.pixelCount(); n++)
				img(n) = nb[n];

			return img;
		}

		void minHash()
		{
			uint8_t nb1[] = { 1, 1,
				1, 0,
				1, 0,
				0, 0 };
			//uint8_t nb1i[] = { 0, 0,
			//	0, 1,
			//	0, 1,
			//	1, 1 };
			uint8_t nb2[] = { 1, 1, // mirror of nb1
				0, 1,
				0, 1,
				0, 0 };
			uint8_t nb3[] = { 1, 0,
				0, 0,
				1, 1,
				1, 0 };
			uint8_t nb4[] = { 0, 1,
				0, 0,
				1, 1,
				0, 1 };


			uint8_t nb5[] = { 1, 0,
				1, 1,
				0, 0,
				1, 0 };
			uint8_t nb6[] = { 0, 1,
				1, 1,
				0, 0,
				0, 1 };
			uint8_t nb7[] = { 0, 0,
				1, 0,
				1, 0,
				1, 1 };
			uint8_t nb8[] = { 0, 0,
				0, 1,
				0, 1,
				1, 1 };

			uint8_t nb9[] = { 0, 0,
				1, 1,
				0, 1,
				1, 1 };

			uint8_t nb10[] = { 0, 0,
				0, 0,
				0, 1,
				0, 1 };

			uint8_t nb11[] = { 1, 1,
				0, 0,
				0, 0,
				0, 0 };

			uint8_t nb12[] = { 0, 0,
				0, 0,
				0, 0,
				1, 1 };

			using hash_t = uint8_t;

			Image<uint8_t> tmp(2, 2, 2);
			hash_t h1 = itl2::minHash<2, uint8_t, hash_t>(put(nb1, tmp));
			//hash_t h1i = itl2::minHash<2, uint8_t, hash_t>(put(nb1i, tmp));
			hash_t h2 = itl2::minHash<2, uint8_t, hash_t>(put(nb2, tmp));
			hash_t h3 = itl2::minHash<2, uint8_t, hash_t>(put(nb3, tmp));
			hash_t h4 = itl2::minHash<2, uint8_t, hash_t>(put(nb4, tmp));
			hash_t h5 = itl2::minHash<2, uint8_t, hash_t>(put(nb5, tmp));
			hash_t h6 = itl2::minHash<2, uint8_t, hash_t>(put(nb6, tmp));
			hash_t h7 = itl2::minHash<2, uint8_t, hash_t>(put(nb7, tmp));
			hash_t h8 = itl2::minHash<2, uint8_t, hash_t>(put(nb8, tmp));
			hash_t h9 = itl2::minHash<2, uint8_t, hash_t>(put(nb9, tmp));
			hash_t h10 = itl2::minHash<2, uint8_t, hash_t>(put(nb10, tmp));
			hash_t h11 = itl2::minHash<2, uint8_t, hash_t>(put(nb11, tmp));
			hash_t h12 = itl2::minHash<2, uint8_t, hash_t>(put(nb12, tmp));

			cout << "nb, hash (should be the same in all the cases)" << endl;
			cout << 1 << ", " << (uint64_t)h1 << endl;
			//cout << "1i" << ", " << (uint64_t)h1i << endl;
			cout << 2 << ", " << (uint64_t)h2 << endl;
			cout << 3 << ", " << (uint64_t)h3 << endl;
			cout << 4 << ", " << (uint64_t)h4 << endl;
			cout << 5 << ", " << (uint64_t)h5 << endl;
			cout << 6 << ", " << (uint64_t)h6 << endl;
			cout << 7 << ", " << (uint64_t)h7 << endl;
			cout << 8 << ", " << (uint64_t)h8 << endl;

			cout << "nb, hash (should be different from above)" << endl;
			cout << 9 << ", " << (uint64_t)h9 << endl;

			cout << "nb, hash (should be the same but different from all above)" << endl;
			cout << 10 << ", " << (uint64_t)h10 << endl;
			cout << 11 << ", " << (uint64_t)h11 << endl;
			cout << 12 << ", " << (uint64_t)h12 << endl;

			//testAssert(h1 == h1i, "h1i");
			testAssert(h1 == h2, "h2");
			testAssert(h1 == h3, "h3");
			testAssert(h1 == h4, "h4");
			testAssert(h1 == h5, "h5");
			testAssert(h1 == h6, "h6");
			testAssert(h1 == h7, "h7");
			testAssert(h1 == h8, "h8");

			testAssert(h1 != h9, "h9");

			testAssert(h9 != h10, "h10 1");
			testAssert(h1 != h10, "h10 2");

			testAssert(h10 == h11, "h11");
			testAssert(h10 == h12, "h12");
		}
	}
}