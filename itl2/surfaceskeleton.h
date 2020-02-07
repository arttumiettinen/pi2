#pragma once

#include "image.h"
#include "neighbourhood.h"
#include "utilities.h"
#include "minhash.h"
#include "math/vec3.h"

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

		/**
		Counts foreground pixels in the given neighbourhood
		*/
		template<typename pixel_t> coord_t foregroundCount(const Image<pixel_t>& nb)
		{
			coord_t count = 0;
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (nb(n) != (pixel_t)0)
					count++;
			}
			return count;
		}

		/*
		Checks if a point is at the end of an arc.
		*/
		template<typename pixel_t> bool isEndPoint(const Image<pixel_t>& nb)
		{
			coord_t numberOfNeighbors = foregroundCount(nb) - 1;   // -1 because the center pixel will be counted as well
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

		/**
		Counts number of connected components in the neighbourhood, setting center pixel to background, and returns true if there are zero or one connected components.
		*/
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

		/**
		Gets 2x2x2 'cube' from 'nb'.
		The first corner of the cube is placed at 'start'.
		Does not do bounds checking.
		*/
		template<typename pixel_t> void getCube(const Image<pixel_t>& nb, Image<pixel_t>& cube, const Vec3c& start)
		{
			cube(0, 0, 0) = nb(start.x, start.y, start.z);
			cube(1, 0, 0) = nb(start.x + 1, start.y, start.z);
			cube(0, 1, 0) = nb(start.x, start.y + 1, start.z);
			cube(1, 1, 0) = nb(start.x + 1, start.y + 1, start.z);

			cube(0, 0, 1) = nb(start.x, start.y, start.z + 1);
			cube(1, 0, 1) = nb(start.x + 1, start.y, start.z + 1);
			cube(0, 1, 1) = nb(start.x, start.y + 1, start.z + 1);
			cube(1, 1, 1) = nb(start.x + 1, start.y + 1, start.z + 1);
		}


		/**
		List that contains minhashes of all 2x2x2 neighbourhoods that represent come kind of surface.
		*/
		inline std::vector<uint8_t> surfaceCubeHashes;

		/**
		Checks if a 2x2x2 pixel cube represents a surface.
		*/
		template<typename pixel_t> bool isSurfaceCube(const Image<pixel_t>& N2)
		{
			coord_t fgCount = foregroundCount(N2);

			// All neighbourhoods with less than 3 foreground points are considered planar.
			if (fgCount < 3)
				return true;

			// All neighbourhoods with more than 4 points are considered non-planar
			if (fgCount > 4)
				return false;

			// Here fgCount == 3 || fgCount == 4

			uint8_t hash = minHash<2, pixel_t, uint8_t>(N2);
			return std::binary_search(surfaceCubeHashes.begin(), surfaceCubeHashes.end(), hash);
		}

		/**
		Flag that indicates if a 2x2x2 binary neighbourhood represents a surface.
		Index this array with nbhash(neighbourhood)
		*/
		inline std::array<bool, 256> isSurfaceFlags;

		/**
		Fills surfaceCubeHashes list and isSurfaceFlags array if they have not been filled yet.
		*/
		inline void createAllowedHashList()
		{
			//vector<uint8_t> surfaceCubeHashes;
			if (surfaceCubeHashes.size() <= 0)
			{
				// 4-point plane
				surfaceCubeHashes.push_back(minHash<2, uint8_t, uint8_t>(Image<uint8_t>(2, 2, 2,
					{ 1, 1,
					1, 1,
					0, 0,
					0, 0 })));

				// 4-point plane, one point displaced
				surfaceCubeHashes.push_back(minHash<2, uint8_t, uint8_t>(Image<uint8_t>(2, 2, 2,
					{ 1, 1,
					1, 0,
					0, 0,
					0, 1 })));

				// 4-point diagonal plane
				surfaceCubeHashes.push_back(minHash<2, uint8_t, uint8_t>(Image<uint8_t>(2, 2, 2,
					{ 1, 0,
					0, 1,
					1, 0,
					0, 1 })));

				// 3-point plane edge
				surfaceCubeHashes.push_back(minHash<2, uint8_t, uint8_t>(Image<uint8_t>(2, 2, 2,
					{ 1, 1,
					1, 0,
					0, 0,
					0, 0 })));

				// 3-point plane, one point displaced
				surfaceCubeHashes.push_back(minHash<2, uint8_t, uint8_t>(Image<uint8_t>(2, 2, 2,
					{ 0, 1,
					1, 0,
					1, 0,
					0, 0 })));

				// 3-point plane edge, diagonal
				surfaceCubeHashes.push_back(minHash<2, uint8_t, uint8_t>(Image<uint8_t>(2, 2, 2,
					{ 1, 0,
					0, 1,
					1, 0,
					0, 0 })));

				std::sort(surfaceCubeHashes.begin(), surfaceCubeHashes.end());

				// Test every possible 2x2x2 neighbourhood and store a flag indicating whether it is a surface block or not.
				Image<uint8_t> N(2, 2, 2);
				for (size_t n = 0; n < 256; n++)
				{
					N(0) = ((n >> 0) & 0x1) != 0 ? 1 : 0;
					N(1) = ((n >> 1) & 0x1) != 0 ? 1 : 0;
					N(2) = ((n >> 2) & 0x1) != 0 ? 1 : 0;
					N(3) = ((n >> 3) & 0x1) != 0 ? 1 : 0;
					N(4) = ((n >> 4) & 0x1) != 0 ? 1 : 0;
					N(5) = ((n >> 5) & 0x1) != 0 ? 1 : 0;
					N(6) = ((n >> 6) & 0x1) != 0 ? 1 : 0;
					N(7) = ((n >> 7) & 0x1) != 0 ? 1 : 0;

					size_t ind = nbhash<uint8_t, size_t>(N);
					if (ind != n)
						throw std::logic_error("Invalid bitwise operations.");

					isSurfaceFlags[n] = isSurfaceCube(N);
				}
			}	
		}

		/**
		Checks if 2x2x2 pixel cube at start represents surface.
		Uses isSurfaceFlags array for fast processing.
		*/
		template<typename pixel_t> bool isSurfaceCubeCached(const Image<pixel_t>& nb, const Vec3c& start)
		{
			Image<pixel_t> N2(2, 2, 2);
			getCube(nb, N2, start);

			//return isSurfaceCube(N2);

			size_t ind = nbhash<uint8_t, size_t>(N2);
			return isSurfaceFlags[ind];
		}

		/**
		Test if the point in the center of the given 3x3x3 neighbourhood is a surface point.
		*/
		template<typename pixel_t> bool isSurfacePoint(const Image<pixel_t>& nb)
		{
			// Divide the 3x3x3 neighbourhood to 8 2x2x2 neighbourhoods, and see if all of them look like surface
			return isSurfaceCubeCached(nb, Vec3c(0, 0, 0)) &&
				isSurfaceCubeCached(nb, Vec3c(1, 0, 0)) &&
				isSurfaceCubeCached(nb, Vec3c(0, 1, 0)) &&
				isSurfaceCubeCached(nb, Vec3c(1, 1, 0)) &&
				isSurfaceCubeCached(nb, Vec3c(0, 0, 1)) &&
				isSurfaceCubeCached(nb, Vec3c(1, 0, 1)) &&
				isSurfaceCubeCached(nb, Vec3c(0, 1, 1)) &&
				isSurfaceCubeCached(nb, Vec3c(1, 1, 1));
		}
	}


	/*
	Performs one thinning iteration and returns count of pixels removed from the image.
	Repeating this thinning process until the image does not change results in skeleton.
	@param retainSurfaces If true, surfaces are not thinned to lines.
	@param img Image that should be thinned.
	*/
	template<typename pixel_t> size_t thin(Image<pixel_t>& img, bool retainSurfaces)
	{
		internals::createAllowedHashList();

		size_t changed = 0;
		size_t counter = 0;
		std::vector<Vec3c> pointsToRemove;
		pointsToRemove.reserve(10000);

		for (int currentBorder = 1; currentBorder <= 6; currentBorder++)
		{
			pointsToRemove.clear();

			#pragma omp parallel
			{
				Image<uint8_t> nbPrivate(3, 3, 3);

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

								// check 6-neighbors to see if the point is a border point of type currentBorder
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
									getNeighbourhood(img, Vec3c(x, y, z), Vec3c(1, 1, 1), nbPrivate, BoundaryCondition::Zero);

									if (!internals::isEndPoint(nbPrivate) &&
										internals::isEulerInvariant(nbPrivate) &&
										internals::isSimplePointHybrid(nbPrivate))
									{
										#pragma omp critical(surfacethin_insert)
											pointsToRemove.push_back(Vec3c(x, y, z));
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
			std::sort(pointsToRemove.begin(), pointsToRemove.end(), vecComparer<coord_t>);

			// Re-check while removing points so that connectivity is preserved.
			Image<uint8_t> nb(3, 3, 3);
			for (size_t n = 0; n < pointsToRemove.size(); n++)
			{
				Vec3c p = pointsToRemove[n];

				getNeighbourhood(img, p, Vec3c(1, 1, 1), nb, BoundaryCondition::Zero);
				//nb(1, 1, 1) = 0; // This must not be done. Function isSimplePointHybrid does not consider the center pixel, but it is needed for isSurfacePoint.

				if (internals::isSimplePointHybrid(nb) &&
					((retainSurfaces && !internals::isSurfacePoint(nb)) || !retainSurfaces)
					)
				{
					img(p) = 0;
					changed++;
				}
			}
		}

		return changed;
	}

	/*
	Performs one thinning iteration and returns count of pixels removed from the image.
	Repeating this thinning process until the image does not change results in surface skeleton.
	@param img Image that should be thinned.
	*/
	template<typename pixel_t> size_t surfaceThin(Image<pixel_t>& img)
	{
		return thin(img, true);
	}

	/**
	Calculates skeleton of image by thinning it until no pixels can be removed.
	Background is assumed to have value 0 and all nonzero pixels are assumed to be foreground.
	Uses algorithm published in
	Lee et al. Building skeleton models via 3-D medial surface/axis thinning algorithms. Computer Vision, Graphics, and Image Processing, 56(6):462–478, 1994.
	The implementation is based on Fiji Skeletonize3D_ plugin (GPL) and ITK implementation, but this version tries to save memory and is multithreaded.
	This version also uses additional conditions to avoid skeletonizing that are not present in Lee's paper. See internals::isSurfacePoint
	and internals::createAllowedHashList for details.
	The resulting skeleton consists of planes (= surfaces) (not lines) so it is referred to as "surface skeleton".
	@param retainSurfaces If true, surfaces are not thinned to lines. If false, as many surfaces as possible are turned into lines, and only those surfaces will be left that surround a cavity.
	@param img Image that should be skeletonized.
	*/
	template<typename pixel_t> void surfaceSkeleton(Image<pixel_t>& img, bool retainSurfaces = true, size_t maxIterations = std::numeric_limits<size_t>::max())
	{

		size_t it = 0;
		size_t changes;
		do
		{
			changes = thin(img, retainSurfaces);

			std::cout << std::endl << changes << " pixels removed." << std::endl;

			it++;
			//raw::writed(img, string("./skeleton/surface_skeleton_sequence") + toString(it));
			if (it >= maxIterations)
				break;

		} while (changes > 0);

	}

	namespace tests
	{
		void surfaceSkeleton();
	}
}
