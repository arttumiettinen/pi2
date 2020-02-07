#pragma once

#include "image.h"
#include "surfaceskeleton.h"
#include "math/vec3.h"

#include <deque>

namespace itl2
{

	namespace internals
	{
		template<typename pixel_t> bool isBorderPoint(int direction, const Image<pixel_t>& img, coord_t x, coord_t y, coord_t z)
		{
			switch (direction)
			{
			case 1:
				if(y > 0)
					return img(x, y - 1, z) == (pixel_t)0;
				return false;
			case 2:
				if(y < img.height() - 1)
					return img(x, y + 1, z) == (pixel_t)0;
				return false;
			case 3:
				if(x > 0)
					return img(x - 1, y, z) == (pixel_t)0;
				return false;
			case 4:
				if(x < img.width() - 1)
					return img(x + 1, y, z) == (pixel_t)0;
				return false;
			case 5:
				if(z > 0)
					return img(x, y, z - 1) == (pixel_t)0;
				return false;
			case 6:
				if(z < img.depth() - 1)
					return img(x, y, z + 1) == (pixel_t)0;
				return false;
			default:
				throw std::runtime_error("Impossible direction value.");
			}
		}

		/**
		Calculates Euler characteristic of the given neighbourhood.
		*/
		inline int eulerChar(const uint8_t Np[26])
		{
			int zeroCells[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
			int oneCells[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
			int twoCells[] = { 0, 0, 0, 0, 0, 0 };

			//zeroCells[0] = 0;
			//zeroCells[1] = 0;
			//zeroCells[2] = 0;
			//zeroCells[3] = 0;
			//zeroCells[4] = 0;
			//zeroCells[5] = 0;
			//zeroCells[6] = 0;
			//zeroCells[7] = 0;

			//oneCells[0] = 0;
			//oneCells[1] = 0;
			//oneCells[2] = 0;
			//oneCells[3] = 0;
			//oneCells[4] = 0;
			//oneCells[5] = 0;
			//oneCells[6] = 0;
			//oneCells[7] = 0;
			//oneCells[8] = 0;
			//oneCells[9] = 0;
			//oneCells[10] = 0;
			//oneCells[11] = 0;

			//twoCells[0] = 0;
			//twoCells[1] = 0;
			//twoCells[2] = 0;
			//twoCells[3] = 0;
			//twoCells[4] = 0;
			//twoCells[5] = 0;

			int zeroCount = 0;
			int oneCount = 0;
			int twoCount = 0;

			if (Np[4] == 1)
			{
				twoCells[4] = 1;

				oneCells[8] = 1;
				oneCells[9] = 1;
				oneCells[10] = 1;
				oneCells[11] = 1;

				zeroCells[4] = 1;
				zeroCells[5] = 1;
				zeroCells[6] = 1;
				zeroCells[7] = 1;
			}

			if (Np[10] == 1)
			{
				twoCells[0] = 1;

				oneCells[3] = 1;
				oneCells[4] = 1;
				oneCells[7] = 1;
				oneCells[8] = 1;

				zeroCells[0] = 1;
				zeroCells[3] = 1;
				zeroCells[4] = 1;
				zeroCells[5] = 1;
			}

			if (Np[12] == 1)
			{
				twoCells[3] = 1;

				oneCells[2] = 1;
				oneCells[6] = 1;
				oneCells[7] = 1;
				oneCells[11] = 1;

				zeroCells[2] = 1;
				zeroCells[3] = 1;
				zeroCells[4] = 1;
				zeroCells[7] = 1;
			}

			if (Np[13] == 1)
			{
				twoCells[1] = 1;

				oneCells[0] = 1;
				oneCells[5] = 1;
				oneCells[9] = 1;
				oneCells[4] = 1;

				zeroCells[0] = 1;
				zeroCells[1] = 1;
				zeroCells[6] = 1;
				zeroCells[5] = 1;
			}

			if (Np[15] == 1)
			{
				twoCells[2] = 1;

				oneCells[5] = 1;
				oneCells[1] = 1;
				oneCells[6] = 1;
				oneCells[10] = 1;

				zeroCells[1] = 1;
				zeroCells[2] = 1;
				zeroCells[7] = 1;
				zeroCells[6] = 1;
			}

			if (Np[21] == 1)
			{
				twoCells[5] = 1;

				oneCells[0] = 1;
				oneCells[1] = 1;
				oneCells[2] = 1;
				oneCells[3] = 1;

				zeroCells[0] = 1;
				zeroCells[1] = 1;
				zeroCells[2] = 1;
				zeroCells[3] = 1;
			}

			if (Np[1] == 1)
			{
				oneCells[8] = 1;

				zeroCells[4] = 1;
				zeroCells[5] = 1;
			}

			if (Np[3] == 1)
			{
				oneCells[11] = 1;

				zeroCells[4] = 1;
				zeroCells[7] = 1;
			}

			if (Np[5] == 1)
			{
				oneCells[9] = 1;

				zeroCells[5] = 1;
				zeroCells[6] = 1;
			}

			if (Np[7] == 1)
			{
				oneCells[10] = 1;

				zeroCells[6] = 1;
				zeroCells[7] = 1;
			}

			if (Np[11] == 1)
			{
				oneCells[4] = 1;

				zeroCells[0] = 1;
				zeroCells[5] = 1;
			}

			if (Np[16] == 1)
			{
				oneCells[5] = 1;

				zeroCells[1] = 1;
				zeroCells[6] = 1;
			}

			if (Np[14] == 1)
			{
				oneCells[6] = 1;

				zeroCells[2] = 1;
				zeroCells[7] = 1;
			}

			if (Np[9] == 1)
			{
				oneCells[7] = 1;

				zeroCells[3] = 1;
				zeroCells[4] = 1;
			}

			if (Np[18] == 1)
			{
				oneCells[3] = 1;

				zeroCells[3] = 1;
				zeroCells[0] = 1;
			}

			if (Np[22] == 1)
			{
				oneCells[0] = 1;

				zeroCells[0] = 1;
				zeroCells[1] = 1;
			}

			if (Np[24] == 1)
			{
				oneCells[1] = 1;

				zeroCells[1] = 1;
				zeroCells[2] = 1;
			}

			if (Np[20] == 1)
			{
				oneCells[2] = 1;

				zeroCells[2] = 1;
				zeroCells[3] = 1;
			}

			if (Np[0] == 1)
				zeroCells[4] = 1;
			if (Np[2] == 1)
				zeroCells[5] = 1;
			if (Np[6] == 1)
				zeroCells[7] = 1;
			if (Np[8] == 1)
				zeroCells[6] = 1;
			if (Np[17] == 1)
				zeroCells[3] = 1;
			if (Np[19] == 1)
				zeroCells[0] = 1;
			if (Np[23] == 1)
				zeroCells[2] = 1;
			if (Np[25] == 1)
				zeroCells[1] = 1;

			for (int i = 0; i < 8; i++)
			{
				if (zeroCells[i] == 1)
					zeroCount++;
			}

			for (int i = 0; i < 12; i++)
			{
				if (oneCells[i] == 1)
					oneCount++;
			}

			for (int i = 0; i < 6; i++)
			{
				if (twoCells[i] == 1)
					twoCount++;
			}

			// euler number = zeroCount - oneCount + twoCount;
			return zeroCount - oneCount + twoCount;
		}

		inline bool isolated0cell(const uint8_t Np[26])
		{
			static const int zeroPoints[] = { 0, 2, 6, 8, 17, 19, 23, 25 };
			static const int zeroCorners[8][6] = {
				{ 1,  3,  4,  9,  10, 12 },
				{ 1,  4,  5,  10, 11, 13 },
				{ 3,  4,  7,  12, 14, 15 },
				{ 4,  5,  7,  13, 15, 16 },
				{ 9,  10, 12, 18, 20, 21 },
				{ 10, 11, 13, 18, 21, 22 },
				{ 12, 14, 15, 20, 21, 24 },
				{ 13, 15, 16, 21, 22, 24 }
				};

			for (int i = 0; i < 8; i++)
			{
				if (Np[zeroPoints[i]] == 1)
				{
					int counter = 0;
					for (int j = 0; j < 6; j++)
					{
						if (Np[zeroCorners[i][j]] == 1)
							counter++;
					}

					if (counter == 0)
						return true;
				}
			}

			return false;
		}

		inline bool isolated1cell(const uint8_t Np[26])
		{
			static const int onePoints[] = { 1, 3, 5, 7, 9, 11, 14, 16, 18, 20, 22, 24 };
			static const int oneCorners[12][10] = {
				{ 0,  2,  3,  4,  5,  9,  10, 11, 12, 13 },
				{ 0,  1,  4,  6,  7,  9,  10, 12, 14, 15 },
				{ 1,  2,  4,  7,  8,  10, 11, 13, 15, 16 },
				{ 3,  4,  5,  6,  8,  12, 13, 14, 15, 16 },

				{ 0,  1,  3,  4,  10, 12, 17, 18, 20, 21 },
				{ 1,  2,  4,  5,  10, 13, 18, 19, 21, 22 },
				{ 3,  4,  6,  7,  12, 15, 20, 21, 23, 24 },
				{ 4,  5,  7,  8,  13, 15, 21, 22, 24, 25 },

				{ 9,  10, 11, 12, 13, 17, 19, 20, 21, 22 },
				{ 9,  10, 12, 14, 15, 17, 18, 21, 23, 24 },
				{ 10, 11, 13, 15, 16, 18, 19, 21, 24, 25 },
				{ 12, 13, 14, 15, 16, 20, 21, 22, 23, 25 }
				};

			for (int i = 0; i < 12; i++)
			{
				if (Np[onePoints[i]] == 1)
				{
					int counter = 0;
					for (int j = 0; j < 10; j++)
					{
						if (Np[oneCorners[i][j]] == 1)
							counter++;
					}

					if (counter == 0)
						return true;
				}
			}

			return false;
		}

		inline bool isolatedInv2cell(const uint8_t Np[26])
		{
			static const int twoCorners[6][16] = {
				{ 0, 1,  2,  3,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16 },
				{ 0, 1,  2,  3,  4,  5,  9,  11, 12, 13, 17, 18, 19, 20, 21, 22 },
				{ 0, 1,  3,  4,  6,  7,  9,  10, 14, 15, 17, 18, 20, 21, 23, 24 },
				{ 1, 2,  4,  5,  7,  8,  10, 11, 15, 16, 18, 19, 21, 22, 24, 25 },
				{ 3, 4,  5,  6,  7,  8,  12, 13, 14, 16, 20, 21, 22, 23, 24, 25 },
				{ 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25 }
				};

			static const int twoPoints[] = { 4, 10, 12, 13, 15, 21 };


			for (int i = 0; i < 6; i++)
			{
				if (Np[twoPoints[i]] == 0)
				{
					for (int j = 0; j < 16; j++)
					{
						if (Np[twoCorners[i][j]] == 0)
							return false;
					}
				}
			}

			return true;
		}


		template<typename pixel_t> bool isSimplePointLine(const Image<pixel_t>& nb)
		{
			// Initialize neighbor list. The list will not contain the center pixel.
			uint8_t Np[26];
			for (coord_t i = 0; i < 13; i++)  // i =  0..12 -> cube[0..12]
				Np[i] = nb(i) != (pixel_t)0 ? 1 : 0;
			for (coord_t i = 14; i < 27; i++) // i = 14..26 -> cube[13..25]
				Np[i - 1] = nb(i) != (pixel_t)0 ? 1 : 0;

			// Calculate count of foreground 6-connected neighbours
			static const size_t N6[] = { 4, 10, 12, 13, 15, 21 };
			int N6Sum = 6;
			for (int i = 0; i < 6; i++)
			{
				N6Sum = N6Sum - Np[N6[i]];
			}

			if (N6Sum == 1)
				return true;

			if (N6Sum == 2)
			{
				if (Np[4] == 0 && Np[21] == 0)
					return false;
				if (Np[12] == 0 && Np[13] == 0)
					return false;
				if (Np[15] == 0 && Np[10] == 0)
					return false;

				return true;
			}

			if (N6Sum > 2)
			{
				if (eulerChar(Np) != 1)
				{
					return false;
				}
				else
				{
					if (isolated0cell(Np))
						return false;
					if (isolated1cell(Np))
						return false;
					if (isolatedInv2cell(Np))
						return false;
				}
			}

			return true;
		}

		/*
		This function uses similar getNeighbourhood pixel ordering than ImageJ implementation,
		but it should not be necessary to use exactly this ordering...
		template<typename pixel_t> void getNeighbourhoodTest(Image<pixel_t>& img, coord_t x, coord_t y, coord_t z, Image<pixel_t>& result)
		{
			if (x <= 1 || y <= 1 || z <= 1 ||
				x >= img.width() - 1 || y >= img.height() - 1 || z >= img.depth() - 1)
			{
				setValue(result, 0);
				return;
			}

			result(0) = img(x - 1,	y - 1, z + 1) != 0 ? 1 : 0;
			result(1) = img(x,		y - 1, z + 1) != 0 ? 1 : 0;
			result(2) = img(x + 1,	y - 1, z + 1) != 0 ? 1 : 0;
			result(3) = img(x - 1,	y - 1, z) != 0 ? 1 : 0;
			result(4) = img(x,		y - 1, z) != 0 ? 1 : 0;
			result(5) = img(x + 1,	y - 1, z) != 0 ? 1 : 0;
			result(6) = img(x - 1,	y - 1, z - 1) != 0 ? 1 : 0;
			result(7) = img(x,		y - 1, z - 1) != 0 ? 1 : 0;
			result(8) = img(x + 1,	y - 1, z - 1) != 0 ? 1 : 0;

			result(9) = img(x - 1,	y,		z + 1) != 0 ? 1 : 0;
			result(10) = img(x,		y,		z + 1) != 0 ? 1 : 0;
			result(11) = img(x + 1,	y,		z + 1) != 0 ? 1 : 0;
			result(12) = img(x - 1,	y,		z) != 0 ? 1 : 0;
			result(13) = img(x,		y,		z) != 0 ? 1 : 0;
			result(14) = img(x + 1,	y,		z) != 0 ? 1 : 0;
			result(15) = img(x - 1,	y,		z - 1) != 0 ? 1 : 0;
			result(16) = img(x,		y,		z - 1) != 0 ? 1 : 0;
			result(17) = img(x + 1,	y,		z - 1) != 0 ? 1 : 0;

			result(18) = img(x - 1,	y + 1,	z + 1) != 0 ? 1 : 0;
			result(19) = img(x,		y + 1,	z + 1) != 0 ? 1 : 0;
			result(20) = img(x + 1, y + 1,	z + 1) != 0 ? 1 : 0;
			result(21) = img(x - 1, y + 1,	z) != 0 ? 1 : 0;
			result(22) = img(x,		y + 1,	z) != 0 ? 1 : 0;
			result(23) = img(x + 1, y + 1,	z) != 0 ? 1 : 0;
			result(24) = img(x - 1, y + 1,	z - 1) != 0 ? 1 : 0;
			result(25) = img(x,		y + 1,	z - 1) != 0 ? 1 : 0;
			result(26) = img(x + 1, y + 1,	z - 1) != 0 ? 1 : 0;
		}
		*/
	}

	/*
	Performs one thinning iteration and returns count of pixels removed from the image.
	Repeating this thinning process until the image does not change results in line skeleton.
	@param img Image that should be thinned.
	*/
	template<typename pixel_t> size_t lineThin(Image<pixel_t>& img)
	{
		size_t changed = 0;
		std::deque<Vec3c> points;
		size_t counter = 0;

		for(int direction = 1; direction <= 6; direction++)
		{
			points.clear();

			// Collect initial list of points
			Image<pixel_t> nbPrivate(3, 3, 3);
			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						pixel_t c = img(x, y, z);
						if (c != (pixel_t)0)
						{
							if (internals::isBorderPoint(direction, img, x, y, z))
							{
								getNeighbourhood(img, Vec3c(x, y, z), Vec3c(1, 1, 1), nbPrivate, BoundaryCondition::Zero);

								// Note: This uses the same isEndPoint than hybridSkeleton!
								if (!internals::isEndPoint(nbPrivate) &&
									internals::isSimplePointLine(nbPrivate))
								{
									points.push_back(Vec3c(x, y, z));
								}
							}
						}
					}
				}

				showThreadProgress(counter, 6 * img.depth());
			}

			// This is required to make the result exactly the same than in Fiji (non-threaded version)
			// Otherwise, the point removal order may change the skeleton points.
			std::sort(points.begin(), points.end(), vecComparer<coord_t>);

			// Process all points
			size_t counter = 0;
			Image<pixel_t> nb(3, 3, 3);
			while (!points.empty())
			{
				Vec3c p;
				if (counter % 2 == 0)
				{
					p = points.front();
					points.pop_front();
				}
				else
				{
					p = points.back();
					points.pop_back();
				}

				getNeighbourhood(img, p, Vec3c(1, 1, 1), nb, BoundaryCondition::Zero);
				if (!internals::isEndPoint(nb) &&
					internals::isSimplePointLine(nb))
				{
					img(p) = 0;
					changed++;
				}

				counter++;
			}
		}

		return changed;
	}

	/*
	Calculates skeleton of image by thinning it until no pixels can be removed.
	Background is assumed to have value 0 and all nonzero pixels are assumed to be foreground.
	Uses algorithm published in
	Klette - Six iteration thinning with simple voxel definition.
	The implementation is based on one by Tuomas Turpeinen.
	The resulting skeleton consists of lines only so it is referred to as "line skeleton".
	Note that if line skeleton is required, it might be better idea to fill holes in the structure
	and use hybrid skeleton as that algorithm seems to produce cleaner looking skeletons.
	@param img Image that should be skeletonized.
	*/
	template<typename pixel_t> void lineSkeleton(Image<pixel_t>& img, size_t maxIterations = std::numeric_limits<size_t>::max())
	{

		size_t it = 0;
		size_t changes;
		do
		{
			changes = lineThin(img);

			std::cout << changes << " pixels removed." << std::endl;

			it++;
			//raw::writed(img, string("./skeleton/line_skeleton_sequence") + toString(it));
			if (it >= maxIterations)
				break;

		} while (changes > 0);

	}

	namespace tests
	{
		void lineSkeleton();
	}

}
