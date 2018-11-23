#pragma once

#include "image.h"
#include "neighbourhood.h"
#include "utilities.h"

#include <algorithm>


namespace itl2
{

	namespace internals
	{
		template<typename pixel_t> bool N(const Image<pixel_t>& img, coord_t x, coord_t y, coord_t z)
		{
			if (y > 0)
				return img(x, y - 1, z) == (pixel_t)0;
			return false;
		}

		template<typename pixel_t> bool S(const Image<pixel_t>& img, coord_t x, coord_t y, coord_t z)
		{
			if (y < img.height() - 1)
				return img(x, y + 1, z) == (pixel_t)0;
			return false;
		}

		template<typename pixel_t> bool E(const Image<pixel_t>& img, coord_t x, coord_t y, coord_t z)
		{
			if(x < img.width() - 1)
				return img(x + 1, y, z) == (pixel_t)0;
			return false;
		}

		template<typename pixel_t> bool W(const Image<pixel_t>& img, coord_t x, coord_t y, coord_t z)
		{
			if(x > 0)
				return img(x - 1, y, z) == (pixel_t)0;
			return false;
		}

		template<typename pixel_t> bool U(const Image<pixel_t>& img, coord_t x, coord_t y, coord_t z)
		{
			if(z < img.depth() - 1)
				return img(x, y, z + 1) == (pixel_t)0;
			return false;
		}

		template<typename pixel_t> bool B(const Image<pixel_t>& img, coord_t x, coord_t y, coord_t z)
		{
			if(z > 0)
				return img(x, y, z - 1) == (pixel_t)0;
			return false;
		}

		/*
		Checks if a point is at the end of an arc.
		*/
		template<typename pixel_t> bool isEndPoint(const Image<pixel_t>& nb)
		{
			coord_t numberOfNeighbors = -1;   // -1 (not 0) because the center pixel will be counted as well
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (nb(n) != (pixel_t)0)
					numberOfNeighbors++;
			}


			return numberOfNeighbors == 1;
		}

		/*
		Lookup table used in function isEulerInvariant.
		*/
		const int eulerLUT[] = {
			0, 1,
			0, -1,
			0, -1,
			0, 1,
			0, -3,
			0, -1,
			0, -1,
			0, 1,
			0, -1,
			0, 1,
			0, 1,
			0, -1,
			0, 3,
			0, 1,
			0, 1,
			0, -1,
			0, -3,
			0, -1,
			0, 3,
			0, 1,
			0, 1,
			0, -1,
			0, 3,
			0, 1,
			0, -1,
			0, 1,

			0, 1,
			0, -1,
			0, 3,
			0, 1,
			0, 1,
			0, -1,
			0, -3,
			0, 3,
			0, -1,
			0, 1,
			0, 1,
			0, 3,
			0, -1,
			0, 1,
			0, -1,
			0, 1,
			0, 1,
			0, -1,
			0, 3,
			0, 1,
			0, 1,
			0, -1,
			0, 1,
			0, 3,
			0, 3,
			0, 1,

			0, 5,
			0, 3,
			0, 3,
			0, 1,
			0, -1,
			0, 1,
			0, 1,
			0, -1,
			0, 3,
			0, 1,
			0, 1,
			0, -1,
			0, -7,
			0, -1,
			0, -1,
			0, 1,
			0, -3,
			0, -1,
			0, -1,
			0, 1,
			0, -1,
			0, 1,
			0, 1,
			0, -1,
			0, 3,
			0, 1,

			0, 1,
			0, -1,
			0, -3,
			0, -1,
			0, 3,
			0, 1,
			0, 1,
			0, -1,
			0, 3,
			0, 1,
			0, -1,
			0, 1,
			0, 1,
			0, -1,
			0, 3,
			0, 1,
			0, 1,
			0, -1,
			0, -3,
			0, 3,
			0, -1,
			0, 1,
			0, 1,
			0, 3,
			0, -1,
			0, 1,

			0, -1,
			0, 1,
			0, 1,
			0, -1,
			0, 3,
			0, 1,
			0, 1,
			0, -1,
			0, 1,
			0, 3,
			0, 3,
			0, 1,
			0, 5,
			0, 3,
			0, 3,
			0, 1,
			0, -1,
			0, 1,
			0, 1,
			0, -1,
			0, 3,
			0, 1,
			0, 1,
			0, -1 };

		/*
		Checks if a point is Euler invariant.
		*/
		template<typename pixel_t> bool isEulerInvariant(const Image<pixel_t>& nb)
		{
			// calculate Euler characteristic for each octant and sum up
			constexpr pixel_t ZERO = (pixel_t)0;
			int eulerChar = 0;
			unsigned char n;

			// Octant SWU
			n = 1;
			if (nb(24) != ZERO)
				n |= 128;
			if (nb(25) != ZERO)
				n |= 64;
			if (nb(15) != ZERO)
				n |= 32;
			if (nb(16) != ZERO)
				n |= 16;
			if (nb(21) != ZERO)
				n |= 8;
			if (nb(22) != ZERO)
				n |= 4;
			if (nb(12) != ZERO)
				n |= 2;
			eulerChar += eulerLUT[n];

			// Octant SEU
			n = 1;
			if (nb(26) != ZERO)
				n |= 128;
			if (nb(23) != ZERO)
				n |= 64;
			if (nb(17) != ZERO)
				n |= 32;
			if (nb(14) != ZERO)
				n |= 16;
			if (nb(25) != ZERO)
				n |= 8;
			if (nb(22) != ZERO)
				n |= 4;
			if (nb(16) != ZERO)
				n |= 2;
			eulerChar += eulerLUT[n];

			// Octant NWU
			n = 1;
			if (nb(18) != ZERO)
				n |= 128;
			if (nb(21) != ZERO)
				n |= 64;
			if (nb(9) != ZERO)
				n |= 32;
			if (nb(12) != ZERO)
				n |= 16;
			if (nb(19) != ZERO)
				n |= 8;
			if (nb(22) != ZERO)
				n |= 4;
			if (nb(10) != ZERO)
				n |= 2;
			eulerChar += eulerLUT[n];

			// Octant NEU
			n = 1;
			if (nb(20) != ZERO)
				n |= 128;
			if (nb(23) != ZERO)
				n |= 64;
			if (nb(19) != ZERO)
				n |= 32;
			if (nb(22) != ZERO)
				n |= 16;
			if (nb(11) != ZERO)
				n |= 8;
			if (nb(14) != ZERO)
				n |= 4;
			if (nb(10) != ZERO)
				n |= 2;
			eulerChar += eulerLUT[n];

			// Octant SWB
			n = 1;
			if (nb(6) != ZERO)
				n |= 128;
			if (nb(15) != ZERO)
				n |= 64;
			if (nb(7) != ZERO)
				n |= 32;
			if (nb(16) != ZERO)
				n |= 16;
			if (nb(3) != ZERO)
				n |= 8;
			if (nb(12) != ZERO)
				n |= 4;
			if (nb(4) != ZERO)
				n |= 2;
			eulerChar += eulerLUT[n];

			// Octant SEB
			n = 1;
			if (nb(8) != ZERO)
				n |= 128;
			if (nb(7) != ZERO)
				n |= 64;
			if (nb(17) != ZERO)
				n |= 32;
			if (nb(16) != ZERO)
				n |= 16;
			if (nb(5) != ZERO)
				n |= 8;
			if (nb(4) != ZERO)
				n |= 4;
			if (nb(14) != ZERO)
				n |= 2;
			eulerChar += eulerLUT[n];

			// Octant NWB
			n = 1;
			if (nb(0) != ZERO)
				n |= 128;
			if (nb(9) != ZERO)
				n |= 64;
			if (nb(3) != ZERO)
				n |= 32;
			if (nb(12) != ZERO)
				n |= 16;
			if (nb(1) != ZERO)
				n |= 8;
			if (nb(10) != ZERO)
				n |= 4;
			if (nb(4) != ZERO)
				n |= 2;
			eulerChar += eulerLUT[n];

			// Octant NEB
			n = 1;
			if (nb(2) != ZERO)
				n |= 128;
			if (nb(1) != ZERO)
				n |= 64;
			if (nb(11) != ZERO)
				n |= 32;
			if (nb(10) != ZERO)
				n |= 16;
			if (nb(5) != ZERO)
				n |= 8;
			if (nb(4) != ZERO)
				n |= 4;
			if (nb(14) != ZERO)
				n |= 2;
			eulerChar += eulerLUT[n];

			return eulerChar == 0;
		}

        /*
		Calculates the number of connected components in the 3D neighbourhood after the center pixel would have been removed.
		*/
		inline void octreeLabeling(int octant, int label, uint8_t *cube)
		{
			if (octant == 1)
			{
				// Set points in this octant to current label
				// and perform recursive labeling of adjacent octants
				if (cube[0] == 1)
					cube[0] = label;
				if (cube[1] == 1)
				{
					cube[1] = label;
					octreeLabeling(2, label, cube);
				}
				if (cube[3] == 1)
				{
					cube[3] = label;
					octreeLabeling(3, label, cube);
				}
				if (cube[4] == 1)
				{
					cube[4] = label;
					octreeLabeling(2, label, cube);
					octreeLabeling(3, label, cube);
					octreeLabeling(4, label, cube);
				}
				if (cube[9] == 1)
				{
					cube[9] = label;
					octreeLabeling(5, label, cube);
				}
				if (cube[10] == 1)
				{
					cube[10] = label;
					octreeLabeling(2, label, cube);
					octreeLabeling(5, label, cube);
					octreeLabeling(6, label, cube);
				}
				if (cube[12] == 1)
				{
					cube[12] = label;
					octreeLabeling(3, label, cube);
					octreeLabeling(5, label, cube);
					octreeLabeling(7, label, cube);
				}
			}

			if (octant == 2)
			{
				if (cube[1] == 1)
				{
					cube[1] = label;
					octreeLabeling(1, label, cube);
				}
				if (cube[4] == 1)
				{
					cube[4] = label;
					octreeLabeling(1, label, cube);
					octreeLabeling(3, label, cube);
					octreeLabeling(4, label, cube);
				}
				if (cube[10] == 1)
				{
					cube[10] = label;
					octreeLabeling(1, label, cube);
					octreeLabeling(5, label, cube);
					octreeLabeling(6, label, cube);
				}
				if (cube[2] == 1)
					cube[2] = label;
				if (cube[5] == 1)
				{
					cube[5] = label;
					octreeLabeling(4, label, cube);
				}
				if (cube[11] == 1)
				{
					cube[11] = label;
					octreeLabeling(6, label, cube);
				}
				if (cube[13] == 1)
				{
					cube[13] = label;
					octreeLabeling(4, label, cube);
					octreeLabeling(6, label, cube);
					octreeLabeling(8, label, cube);
				}
			}

			if (octant == 3)
			{
				if (cube[3] == 1)
				{
					cube[3] = label;
					octreeLabeling(1, label, cube);
				}
				if (cube[4] == 1)
				{
					cube[4] = label;
					octreeLabeling(1, label, cube);
					octreeLabeling(2, label, cube);
					octreeLabeling(4, label, cube);
				}
				if (cube[12] == 1)
				{
					cube[12] = label;
					octreeLabeling(1, label, cube);
					octreeLabeling(5, label, cube);
					octreeLabeling(7, label, cube);
				}
				if (cube[6] == 1)
					cube[6] = label;
				if (cube[7] == 1)
				{
					cube[7] = label;
					octreeLabeling(4, label, cube);
				}
				if (cube[14] == 1)
				{
					cube[14] = label;
					octreeLabeling(7, label, cube);
				}
				if (cube[15] == 1)
				{
					cube[15] = label;
					octreeLabeling(4, label, cube);
					octreeLabeling(7, label, cube);
					octreeLabeling(8, label, cube);
				}
			}

			if (octant == 4)
			{
				if (cube[4] == 1)
				{
					cube[4] = label;
					octreeLabeling(1, label, cube);
					octreeLabeling(2, label, cube);
					octreeLabeling(3, label, cube);
				}
				if (cube[5] == 1)
				{
					cube[5] = label;
					octreeLabeling(2, label, cube);
				}
				if (cube[13] == 1)
				{
					cube[13] = label;
					octreeLabeling(2, label, cube);
					octreeLabeling(6, label, cube);
					octreeLabeling(8, label, cube);
				}
				if (cube[7] == 1)
				{
					cube[7] = label;
					octreeLabeling(3, label, cube);
				}
				if (cube[15] == 1)
				{
					cube[15] = label;
					octreeLabeling(3, label, cube);
					octreeLabeling(7, label, cube);
					octreeLabeling(8, label, cube);
				}
				if (cube[8] == 1)
					cube[8] = label;
				if (cube[16] == 1)
				{
					cube[16] = label;
					octreeLabeling(8, label, cube);
				}
			}

			if (octant == 5)
			{
				if (cube[9] == 1)
				{
					cube[9] = label;
					octreeLabeling(1, label, cube);
				}
				if (cube[10] == 1)
				{
					cube[10] = label;
					octreeLabeling(1, label, cube);
					octreeLabeling(2, label, cube);
					octreeLabeling(6, label, cube);
				}
				if (cube[12] == 1)
				{
					cube[12] = label;
					octreeLabeling(1, label, cube);
					octreeLabeling(3, label, cube);
					octreeLabeling(7, label, cube);
				}
				if (cube[17] == 1)
					cube[17] = label;
				if (cube[18] == 1)
				{
					cube[18] = label;
					octreeLabeling(6, label, cube);
				}
				if (cube[20] == 1)
				{
					cube[20] = label;
					octreeLabeling(7, label, cube);
				}
				if (cube[21] == 1)
				{
					cube[21] = label;
					octreeLabeling(6, label, cube);
					octreeLabeling(7, label, cube);
					octreeLabeling(8, label, cube);
				}
			}

			if (octant == 6)
			{
				if (cube[10] == 1)
				{
					cube[10] = label;
					octreeLabeling(1, label, cube);
					octreeLabeling(2, label, cube);
					octreeLabeling(5, label, cube);
				}
				if (cube[11] == 1)
				{
					cube[11] = label;
					octreeLabeling(2, label, cube);
				}
				if (cube[13] == 1)
				{
					cube[13] = label;
					octreeLabeling(2, label, cube);
					octreeLabeling(4, label, cube);
					octreeLabeling(8, label, cube);
				}
				if (cube[18] == 1)
				{
					cube[18] = label;
					octreeLabeling(5, label, cube);
				}
				if (cube[21] == 1)
				{
					cube[21] = label;
					octreeLabeling(5, label, cube);
					octreeLabeling(7, label, cube);
					octreeLabeling(8, label, cube);
				}
				if (cube[19] == 1)
					cube[19] = label;
				if (cube[22] == 1)
				{
					cube[22] = label;
					octreeLabeling(8, label, cube);
				}
			}

			if (octant == 7)
			{
				if (cube[12] == 1)
				{
					cube[12] = label;
					octreeLabeling(1, label, cube);
					octreeLabeling(3, label, cube);
					octreeLabeling(5, label, cube);
				}
				if (cube[14] == 1)
				{
					cube[14] = label;
					octreeLabeling(3, label, cube);
				}
				if (cube[15] == 1)
				{
					cube[15] = label;
					octreeLabeling(3, label, cube);
					octreeLabeling(4, label, cube);
					octreeLabeling(8, label, cube);
				}
				if (cube[20] == 1)
				{
					cube[20] = label;
					octreeLabeling(5, label, cube);
				}
				if (cube[21] == 1)
				{
					cube[21] = label;
					octreeLabeling(5, label, cube);
					octreeLabeling(6, label, cube);
					octreeLabeling(8, label, cube);
				}
				if (cube[23] == 1)
					cube[23] = label;
				if (cube[24] == 1)
				{
					cube[24] = label;
					octreeLabeling(8, label, cube);
				}
			}

			if (octant == 8)
			{
				if (cube[13] == 1)
				{
					cube[13] = label;
					octreeLabeling(2, label, cube);
					octreeLabeling(4, label, cube);
					octreeLabeling(6, label, cube);
				}
				if (cube[15] == 1)
				{
					cube[15] = label;
					octreeLabeling(3, label, cube);
					octreeLabeling(4, label, cube);
					octreeLabeling(7, label, cube);
				}
				if (cube[16] == 1)
				{
					cube[16] = label;
					octreeLabeling(4, label, cube);
				}
				if (cube[21] == 1)
				{
					cube[21] = label;
					octreeLabeling(5, label, cube);
					octreeLabeling(6, label, cube);
					octreeLabeling(7, label, cube);
				}
				if (cube[22] == 1)
				{
					cube[22] = label;
					octreeLabeling(6, label, cube);
				}
				if (cube[24] == 1)
				{
					cube[24] = label;
					octreeLabeling(7, label, cube);
				}
				if (cube[25] == 1)
					cube[25] = label;
			}
		}

		template<typename pixel_t> bool isSimplePointHybrid(const Image<pixel_t>& nb)
		{
			// Initialize neighbor list for labeling. The list will not contain center pixel.
			// Value 1 denotes unlabeled foreground voxels, value 0 denotes background.
			uint8_t cube[26];
			for (coord_t i = 0; i < 13; i++)  // i =  0..12 -> cube[0..12]
				cube[i] = nb(i) != (pixel_t)0 ? 1 : 0;
			for (coord_t i = 14; i < 27; i++) // i = 14..26 -> cube[13..25]
				cube[i - 1] = nb(i) != (pixel_t)0 ? 1 : 0;
			
			int label = 2;
			for (coord_t i = 0; i < 26; i++)
			{
				if (cube[i] == 1)
				{
					// Start recursion with any octant that contains the point i
					switch (i)
					{
					case 0:
					case 1:
					case 3:
					case 4:
					case 9:
					case 10:
					case 12:
						octreeLabeling(1, label, cube);
						break;
					case 2:
					case 5:
					case 11:
					case 13:
						octreeLabeling(2, label, cube);
						break;
					case 6:
					case 7:
					case 14:
					case 15:
						octreeLabeling(3, label, cube);
						break;
					case 8:
					case 16:
						octreeLabeling(4, label, cube);
						break;
					case 17:
					case 18:
					case 20:
					case 21:
						octreeLabeling(5, label, cube);
						break;
					case 19:
					case 22:
						octreeLabeling(6, label, cube);
						break;
					case 23:
					case 24:
						octreeLabeling(7, label, cube);
						break;
					case 25:
						octreeLabeling(8, label, cube);
						break;
					}
					label++;
					if (label - 2 >= 2)
					{
						return false;
					}
				}
			}

			// label-2 is the number of connected compontents
			return true;
		}

		inline bool pointSorter(const math::Vec3c& a, const math::Vec3c& b)
		{
			if (a.z != b.z)
				return a.z < b.z;
			if (a.y != b.y)
				return a.y < b.y;
			return a.x < b.x;
		}
	}

	/*
	Performs one thinning iteration and returns count of pixels removed from the image.
	Repeating this thinning process until the image does not change results in hybrid skeleton.
	@param img Image that should be thinned.
	*/
	template<typename pixel_t> size_t hybridThin(Image<pixel_t>& img)
	{
		size_t changed = 0;
		size_t counter = 0;
		vector<math::Vec3c> pointsToRemove;
		pointsToRemove.reserve(10000);

		for (int currentBorder = 1; currentBorder <= 6; currentBorder++)
		{
			pointsToRemove.clear();

			#pragma omp parallel
			{
				Image<pixel_t> nbPrivate(3, 3, 3);

				#pragma omp for
				for (coord_t z = 0; z < img.depth(); z++)
				{
					for (coord_t y = 0; y < img.height(); y++)
					{
						for (coord_t x = 0; x < img.width(); x++)
						{
							pixel_t c = img(x, y, z);
							if (c != (pixel_t)0)
							{

								// check 6-neighbors if point is a border point of type currentBorder
								bool isBorderPoint = false;
								switch (currentBorder)
								{
								case 1: isBorderPoint = internals::N(img, x, y, z); break;
								case 2: isBorderPoint = internals::S(img, x, y, z); break;
								case 3: isBorderPoint = internals::E(img, x, y, z); break;
								case 4: isBorderPoint = internals::W(img, x, y, z); break;
								case 5: isBorderPoint = internals::U(img, x, y, z); break;
								case 6: isBorderPoint = internals::B(img, x, y, z); break;
								}

								if (isBorderPoint)
								{
									getNeighbourhood(img, math::Vec3c(x, y, z), math::Vec3c(1, 1, 1), nbPrivate, Zero);

									if (!internals::isEndPoint(nbPrivate) &&
										internals::isEulerInvariant(nbPrivate) &&
										internals::isSimplePointHybrid(nbPrivate))
									{
										#pragma omp critical(thin_insert)
											pointsToRemove.push_back(math::Vec3c(x, y, z));
									}
								}
							}
						}

					}

					showThreadProgress(counter, 6 * img.depth());
				}
			}

			// This is required to make the result exactly the same than in Fiji (non-threaded version)
			// Otherwise, the point removal order may change the skeleton points.
			std::sort(pointsToRemove.begin(), pointsToRemove.end(), internals::pointSorter);

			// Re-check while removing points so that connectivity is preserved.
			Image<pixel_t> nb(3, 3, 3);
			for (size_t n = 0; n < pointsToRemove.size(); n++)
			{
				math::Vec3c p = pointsToRemove[n];

				getNeighbourhood(img, p, math::Vec3c(1, 1, 1), nb, Zero);
				//nb(1, 1, 1) = 0; // This is redundant as isSimplePointHybrid does not consider the center pixel anyway.

				if (internals::isSimplePointHybrid(nb))
				{
					img(p) = 0;
					changed++;
				}
			}
		}

		return changed;
	}

	/**
	Calculates skeleton of image by thinning it until no pixels can be removed.
	Background is assumed to have value 0 and all nonzero pixels are assumed to be foreground.
	Uses algorithm published in
	Lee et al. Building skeleton models via 3-D medial surface/axis thinning algorithms. Computer Vision, Graphics, and Image Processing, 56(6):462–478, 1994.
	The implementation is based on Fiji Skeletonize3D_ plugin (GPL) and ITK implementation,
	but this version tries to save memory and is multithreaded.
	The resulting skeleton consists of lines and plates so it is referred to as "hybrid skeleton".
	@param img Image that should be skeletonized.
	*/
	template<typename pixel_t> void hybridSkeleton(Image<pixel_t>& img)
	{

		size_t changes;
		do
		{
			changes = hybridThin(img);

			cout << changes << " pixels removed." << endl;

		} while (changes > 0);

	}

	namespace tests
	{
		void hybridSkeleton();
		void cavities();
	}
}
