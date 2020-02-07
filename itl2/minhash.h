#pragma once

#include "image.h"
#include "pointprocess.h"

#include <array>

namespace itl2
{
	namespace internals
	{
	
	    template<class T> struct dependent_false : std::false_type {};
	
		/**
		Template definition for neighbourhood rotating functions.
		This function rotates clockwise around z-axis such that x becomes y.
		*/
		template<int size, typename pixel_t> inline void rotateright3(const Image<pixel_t>& original, Image<pixel_t>& rotated)
		{
			if constexpr (size == 2)
			{
				// Original
				//0 1 | 4 5
				//2 3 | 6 7

				// Rotated
				//2 0 | 6 4
				//3 1 | 7 5

				rotated(0) = original(2);
				rotated(1) = original(0);
				rotated(2) = original(3);
				rotated(3) = original(1);

				rotated(4) = original(6);
				rotated(5) = original(4);
				rotated(6) = original(7);
				rotated(7) = original(5);
			}
			else if constexpr (size == 3)
			{
				// Original
				//0 1 2 | 9  10 11 | 18 19 20
				//3 4 5 | 12 13 14 | 21 22 23
				//6 7 8 | 15 16 17 | 24 25 26

				// Rotated
				//6 3 0 | 15 12 9  | 24 21 18
				//7 4 1 | 16 13 10 | 25 22 19
				//8 5 2 | 17 14 11 | 26 23 20

				rotated(0) = original(6);
				rotated(1) = original(3);
				rotated(2) = original(0);
				rotated(3) = original(7);
				rotated(4) = original(4);
				rotated(5) = original(1);
				rotated(6) = original(8);
				rotated(7) = original(5);
				rotated(8) = original(2);
				rotated(9) = original(15);
				rotated(10) = original(12);
				rotated(11) = original(9);
				rotated(12) = original(16);
				rotated(13) = original(13);
				rotated(14) = original(10);
				rotated(15) = original(17);
				rotated(16) = original(14);
				rotated(17) = original(11);
				rotated(18) = original(24);
				rotated(19) = original(21);
				rotated(20) = original(18);
				rotated(21) = original(25);
				rotated(22) = original(22);
				rotated(23) = original(19);
				rotated(24) = original(26);
				rotated(25) = original(23);
				rotated(26) = original(20);
			}
			else if constexpr (size == 4)
			{
				// Original
				//0   1  2  3 | 16 17 18 19 | 32 33 34 35 | 48 49 50 51
				//4   5  6  7 | 20 21 22 23 | 36 37 38 39 | 52 53 54 55
				//8   9 10 11 | 24 25 26 27 | 40 41 42 43 | 56 57 58 59
				//12 13 14 15 | 28 29 30 31 | 44 45 46 47 | 60 61 62 63

				// Rotated:
				//12  8  4  0 | --
				//13  9  5  1 | --
				//14 10  6  2 | --
				//15 11  7  3 | --

				for (size_t i = 0; i < 64; i += 16)
				{
					rotated(i + 0) = original(i + 12);
					rotated(i + 1) = original(i + 8);
					rotated(i + 2) = original(i + 4);
					rotated(i + 3) = original(i + 0);
					rotated(i + 4) = original(i + 13);
					rotated(i + 5) = original(i + 9);
					rotated(i + 6) = original(i + 5);
					rotated(i + 7) = original(i + 1);
					rotated(i + 8) = original(i + 14);
					rotated(i + 9) = original(i + 10);
					rotated(i + 10) = original(i + 6);
					rotated(i + 11) = original(i + 2);
					rotated(i + 12) = original(i + 15);
					rotated(i + 13) = original(i + 11);
					rotated(i + 14) = original(i + 7);
					rotated(i + 15) = original(i + 3);
				}
			}
			else
			{
				static_assert(dependent_false<pixel_t>::value, "Not configured");
			}
		}

		template<int size, typename pixel_t> inline void reslicetop3(const Image<pixel_t>& original, Image<pixel_t>& rotated)
		{
			if constexpr (size == 2)
			{
				// Original
				//0 1 | 4 5
				//2 3 | 6 7

				// Rotated
				//4 5 | 6 7
				//0 1 | 2 3

				rotated(0) = original(4);
				rotated(1) = original(5);
				rotated(2) = original(0);
				rotated(3) = original(1);

				rotated(4) = original(6);
				rotated(5) = original(7);
				rotated(6) = original(2);
				rotated(7) = original(3);
			}
			else if constexpr (size == 3)
			{
				// Original
				//0 1 2 | 9  10 11 | 18 19 20
				//3 4 5 | 12 13 14 | 21 22 23
				//6 7 8 | 15 16 17 | 24 25 26

				// Rotated
				//18 19 20 | 21 22 23 | 24 25 26
				// 9 10 11 | 12 13 14 | 15 16 17
				// 0  1  2 |  3  4  5 |  6  7  8

				rotated(0) = original(18);
				rotated(1) = original(19);
				rotated(2) = original(20);
				rotated(3) = original(9);
				rotated(4) = original(10);
				rotated(5) = original(11);
				rotated(6) = original(0);
				rotated(7) = original(1);
				rotated(8) = original(2);
				rotated(9) = original(21);
				rotated(10) = original(22);
				rotated(11) = original(23);
				rotated(12) = original(12);
				rotated(13) = original(13);
				rotated(14) = original(14);
				rotated(15) = original(3);
				rotated(16) = original(4);
				rotated(17) = original(5);
				rotated(18) = original(24);
				rotated(19) = original(25);
				rotated(20) = original(26);
				rotated(21) = original(15);
				rotated(22) = original(16);
				rotated(23) = original(17);
				rotated(24) = original(6);
				rotated(25) = original(7);
				rotated(26) = original(8);
			}
			else if constexpr (size == 4)
			{
				// Original
				//0   1  2  3 | 16 17 18 19 | 32 33 34 35 | 48 49 50 51
				//4   5  6  7 | 20 21 22 23 | 36 37 38 39 | 52 53 54 55
				//8   9 10 11 | 24 25 26 27 | 40 41 42 43 | 56 57 58 59
				//12 13 14 15 | 28 29 30 31 | 44 45 46 47 | 60 61 62 63

				// Rotated
				//48 49 50 51 | 52 53 54 55 | 56 57 58 59 | 60 61 62 63
				//32 33 34 35 | 36 37 38 39 | 40 41 42 43 | 44 45 46 47
				//16 17 18 19 | 20 21 22 23 | 24 25 26 27 | 28 29 30 31
				// 0  1  2  3 |  4  5  6  7 |  8  9 10 11 | 12 13 14 15


				rotated(0) = original(48);
				rotated(1) = original(49);
				rotated(2) = original(50);
				rotated(3) = original(51);
				rotated(4) = original(32);
				rotated(5) = original(33);
				rotated(6) = original(34);
				rotated(7) = original(35);
				rotated(8) = original(16);
				rotated(9) = original(17);
				rotated(10) = original(18);
				rotated(11) = original(19);
				rotated(12) = original(0);
				rotated(13) = original(1);
				rotated(14) = original(2);
				rotated(15) = original(3);
				rotated(16) = original(52);
				rotated(17) = original(53);
				rotated(18) = original(54);
				rotated(19) = original(55);
				rotated(20) = original(36);
				rotated(21) = original(37);
				rotated(22) = original(38);
				rotated(23) = original(39);
				rotated(24) = original(20);
				rotated(25) = original(21);
				rotated(26) = original(22);
				rotated(27) = original(23);
				rotated(28) = original(4);
				rotated(29) = original(5);
				rotated(30) = original(6);
				rotated(31) = original(7);

				// Rotated
				//48 49 50 51 | 52 53 54 55 | 56 57 58 59 | 60 61 62 63
				//32 33 34 35 | 36 37 38 39 | 40 41 42 43 | 44 45 46 47
				//16 17 18 19 | 20 21 22 23 | 24 25 26 27 | 28 29 30 31
				// 0  1  2  3 |  4  5  6  7 |  8  9 10 11 | 12 13 14 15
				rotated(32) = original(56);
				rotated(33) = original(57);
				rotated(34) = original(58);
				rotated(35) = original(59);
				rotated(36) = original(40);
				rotated(37) = original(41);
				rotated(38) = original(42);
				rotated(39) = original(43);
				rotated(40) = original(24);
				rotated(41) = original(25);
				rotated(42) = original(26);
				rotated(43) = original(27);
				rotated(44) = original(8);
				rotated(45) = original(9);
				rotated(46) = original(10);
				rotated(47) = original(11);
				rotated(48) = original(60);
				rotated(49) = original(61);
				rotated(50) = original(62);
				rotated(51) = original(63);
				rotated(52) = original(44);
				rotated(53) = original(45);
				rotated(54) = original(46);
				rotated(55) = original(47);
				rotated(56) = original(28);
				rotated(57) = original(29);
				rotated(58) = original(30);
				rotated(59) = original(31);
				rotated(60) = original(12);
				rotated(61) = original(13);
				rotated(62) = original(14);
				rotated(63) = original(15);
			}
			else
			{
				static_assert(dependent_false<pixel_t>::value, "Not configured");
			}
		}

		template<int size, typename pixel_t> inline void resliceright3(const Image<pixel_t>& original, Image<pixel_t>& rotated)
		{
			if constexpr (size == 2)
			{
				// Original
				//0 1 | 4 5
				//2 3 | 6 7

				// Rotated
				//1 5 | 0 4
				//3 7 | 2 6

				rotated(0) = original(1);
				rotated(1) = original(5);
				rotated(2) = original(3);
				rotated(3) = original(7);

				rotated(4) = original(0);
				rotated(5) = original(4);
				rotated(6) = original(2);
				rotated(7) = original(6);
			}
			else if constexpr (size == 3)
			{
				// Original
				//0 1 2 | 9  10 11 | 18 19 20
				//3 4 5 | 12 13 14 | 21 22 23
				//6 7 8 | 15 16 17 | 24 25 26

				// Rotated
				//2 11 20 |  1 10 19 |  0  9 18
				//5 14 23 |  4 13 22 |  3 12 21
				//8 17 26 |  7 16 25 |  6 15 24

				rotated(0) = original(2);
				rotated(1) = original(11);
				rotated(2) = original(20);
				rotated(3) = original(5);
				rotated(4) = original(14);
				rotated(5) = original(23);
				rotated(6) = original(8);
				rotated(7) = original(17);
				rotated(8) = original(26);
				rotated(9) = original(1);
				rotated(10) = original(10);
				rotated(11) = original(19);
				rotated(12) = original(4);
				rotated(13) = original(13);
				rotated(14) = original(22);
				rotated(15) = original(7);
				rotated(16) = original(16);
				rotated(17) = original(25);
				rotated(18) = original(0);
				rotated(19) = original(9);
				rotated(20) = original(18);
				rotated(21) = original(3);
				rotated(22) = original(12);
				rotated(23) = original(21);
				rotated(24) = original(6);
				rotated(25) = original(15);
				rotated(26) = original(24);
			}
			else if constexpr (size == 4)
			{
				// Original
				//0   1  2  3 | 16 17 18 19 | 32 33 34 35 | 48 49 50 51
				//4   5  6  7 | 20 21 22 23 | 36 37 38 39 | 52 53 54 55
				//8   9 10 11 | 24 25 26 27 | 40 41 42 43 | 56 57 58 59
				//12 13 14 15 | 28 29 30 31 | 44 45 46 47 | 60 61 62 63

				// Rotated
				//3  19 35 51 | 2  18 34 50 | 1  17 33 49 | 0  16 32 48
				//7  23 39 55 | 6  22 38 54 | 5  21 37 53 | 4  20 36 52
				//11 27 43 59 | 10 26 42 58 | 9  25 41 57 | 8  24 40 56
				//15 31 47 63 | 14 30 46 62 | 13 29 45 61 | 12 28 44 60

				rotated(0) = original(3);
				rotated(1) = original(19);
				rotated(2) = original(35);
				rotated(3) = original(51);
				rotated(4) = original(7);
				rotated(5) = original(23);
				rotated(6) = original(39);
				rotated(7) = original(55);
				rotated(8) = original(11);
				rotated(9) = original(27);
				rotated(10) = original(43);
				rotated(11) = original(59);
				rotated(12) = original(15);
				rotated(13) = original(31);
				rotated(14) = original(47);
				rotated(15) = original(63);
				rotated(16) = original(2);
				rotated(17) = original(18);
				rotated(18) = original(34);
				rotated(19) = original(50);
				rotated(20) = original(6);
				rotated(21) = original(22);
				rotated(22) = original(38);
				rotated(23) = original(54);
				rotated(24) = original(10);
				rotated(25) = original(26);
				rotated(26) = original(42);
				rotated(27) = original(58);
				rotated(28) = original(14);
				rotated(29) = original(30);
				rotated(30) = original(46);
				rotated(31) = original(62);

				// Rotated
				//3  19 35 51 | 2  18 34 50 | 1  17 33 49 | 0  16 32 48
				//7  23 39 55 | 6  22 38 54 | 5  21 37 53 | 4  20 36 52
				//11 27 43 59 | 10 26 42 58 | 9  25 41 57 | 8  24 40 56
				//15 31 47 63 | 14 30 46 62 | 13 29 45 61 | 12 28 44 60
				rotated(32) = original(1);
				rotated(33) = original(17);
				rotated(34) = original(33);
				rotated(35) = original(49);
				rotated(36) = original(5);
				rotated(37) = original(21);
				rotated(38) = original(37);
				rotated(39) = original(53);
				rotated(40) = original(9);
				rotated(41) = original(25);
				rotated(42) = original(41);
				rotated(43) = original(57);
				rotated(44) = original(13);
				rotated(45) = original(29);
				rotated(46) = original(45);
				rotated(47) = original(61);
				rotated(48) = original(0);
				rotated(49) = original(16);
				rotated(50) = original(32);
				rotated(51) = original(48);
				rotated(52) = original(4);
				rotated(53) = original(20);
				rotated(54) = original(36);
				rotated(55) = original(52);
				rotated(56) = original(8);
				rotated(57) = original(24);
				rotated(58) = original(40);
				rotated(59) = original(56);
				rotated(60) = original(12);
				rotated(61) = original(28);
				rotated(62) = original(44);
				rotated(63) = original(60);
			}
			else
			{
				static_assert(dependent_false<pixel_t>::value, "Not configured");
			}
		}

		template<int size, typename pixel_t> inline void mirror3(const Image<pixel_t>& original, Image<pixel_t>& rotated)
		{
			if constexpr (size == 2)
			{
				// Original
				//0 1 | 4 5
				//2 3 | 6 7

				// Rotated
				//1 0 | 5 4
				//3 2 | 7 6

				rotated(0) = original(1);
				rotated(1) = original(0);
				rotated(2) = original(3);
				rotated(3) = original(2);

				rotated(4) = original(5);
				rotated(5) = original(4);
				rotated(6) = original(7);
				rotated(7) = original(6);
			}
			else if constexpr (size == 3)
			{
				// Original
				//0 1 2 | 9  10 11 | 18 19 20
				//3 4 5 | 12 13 14 | 21 22 23
				//6 7 8 | 15 16 17 | 24 25 26

				// Rotated
				//2 1 0 | 11 10 9  | 20 19 18
				//5 4 3 | 14 13 12 | 23 22 21
				//8 7 6 | 17 16 15 | 26 25 24

				rotated(0) = original(2);
				rotated(1) = original(1);
				rotated(2) = original(0);
				rotated(3) = original(5);
				rotated(4) = original(4);
				rotated(5) = original(3);
				rotated(6) = original(8);
				rotated(7) = original(7);
				rotated(8) = original(6);
				rotated(9) = original(11);
				rotated(10) = original(10);
				rotated(11) = original(9);
				rotated(12) = original(14);
				rotated(13) = original(13);
				rotated(14) = original(12);
				rotated(15) = original(17);
				rotated(16) = original(16);
				rotated(17) = original(15);
				rotated(18) = original(20);
				rotated(19) = original(19);
				rotated(20) = original(18);
				rotated(21) = original(23);
				rotated(22) = original(22);
				rotated(23) = original(21);
				rotated(24) = original(26);
				rotated(25) = original(25);
				rotated(26) = original(24);
			}
			else if constexpr (size == 4)
			{
				// Original
				//0   1  2  3 | 16 17 18 19 | 32 33 34 35 | 48 49 50 51
				//4   5  6  7 | 20 21 22 23 | 36 37 38 39 | 52 53 54 55
				//8   9 10 11 | 24 25 26 27 | 40 41 42 43 | 56 57 58 59
				//12 13 14 15 | 28 29 30 31 | 44 45 46 47 | 60 61 62 63

				//Rotated
				//3   2  1  0 | --
				//7   6  5  4 | --
				//11 10  9  8 | --
				//15 14 13 12 | --

				for (size_t i = 0; i < 64; i += 16)
				{
					rotated(i + 0) = original(i + 3);
					rotated(i + 1) = original(i + 2);
					rotated(i + 2) = original(i + 1);
					rotated(i + 3) = original(i + 0);
					rotated(i + 4) = original(i + 7);
					rotated(i + 5) = original(i + 6);
					rotated(i + 6) = original(i + 5);
					rotated(i + 7) = original(i + 4);
					rotated(i + 8) = original(i + 11);
					rotated(i + 9) = original(i + 10);
					rotated(i + 10) = original(i + 9);
					rotated(i + 11) = original(i + 8);
					rotated(i + 12) = original(i + 15);
					rotated(i + 13) = original(i + 14);
					rotated(i + 14) = original(i + 13);
					rotated(i + 15) = original(i + 12);
				}
			}
			else
			{
				static_assert(dependent_false<pixel_t>::value, "Not configured");
			}
		}


		/*
		Calculates 2^n for n in [0, sizeof(hash_t)].
		Uses simple shift operation for the calculation.
		*/
		template<typename hash_t> hash_t hashPow(size_t n)
		{
			hash_t result = 1;
			return result << n;
		}

		/**
		Calculate hash of neighbourhood.
		Assumes that the neighbourhood array contains only values 0 and 1.
		*/
		template<typename pixel_t, typename hash_t> hash_t nbhash(const Image<pixel_t>& bytes)
		{
			hash_t hash = 0;
			for (coord_t n = 0; n < bytes.pixelCount(); n++)
			{
				//hash += (bytes(n) != 0 ? 1 : 0) * hashPow(n); // Slowest, safest
				hash += (bytes(n) & 0x1) * hashPow<hash_t>(n); // Faster, should be safe
				//hash += bytes(n) * hashPow(n); // Fastest
			}
			return hash;
		}


		/**
		Calculate hashes for the four rotations of the neighbourhood.
		@param temp1 The original neighbourhood, will be overwritten with 180 deg rotated one.
		@param temp2 Temporary space. Will contain 270 deg rotated version of original at output.
		@param h List of hashes
		@param h0 Index in the hash list where the first hash should be placed.
		*/
		template<size_t size, typename pixel_t, typename hash_t> void performRotations(Image<pixel_t>& temp1, Image<pixel_t>& temp2, std::array<hash_t, 48>& h, size_t h0)
		{
			// Calculate hashes
			h[h0] = nbhash<pixel_t, hash_t>(temp1);	// original

			rotateright3<size>(temp1, temp2);
			h[h0 + 1] = nbhash<pixel_t, hash_t>(temp2); // 90 deg rotated

			rotateright3<size>(temp2, temp1);
			h[h0 + 2] = nbhash<pixel_t, hash_t>(temp1); // 180 deg rotated

			rotateright3<size>(temp1, temp2);
			h[h0 + 3] = nbhash<pixel_t, hash_t>(temp2); // 270 deg rotated
		}

		/**
		Calculates 24 hashes of different orientations of the neighbourhood.
		@param neighbourhood Original neighbourhood.
		@param temp1, temp2 Temporary arrays, the same size as the neighbourhood.
		@param h Hash array.
		@param h0 Index where the first hash in the hash array h should be placed.
		*/
		template<size_t size, typename pixel_t, typename hash_t> void performReslices(const Image<pixel_t>& neighbourhood, Image<pixel_t>& temp1, Image<pixel_t>& temp2, std::array<hash_t, 48>& h, size_t h0)
		{
			const size_t count = size * size * size;

			//original
			//memcpy((void*)&temp1[0], (void*)&neighbourhood[0], count * sizeof(uint8_t));
			setValue(temp1, neighbourhood);
			performRotations<size, pixel_t, hash_t>(temp1, temp2, h, h0 + 0);

			//original+reslice
			//memcpy((void*)&temp1[0], (void*)&neighbourhood[0], count * sizeof(uint8_t));
			setValue(temp1, neighbourhood);
			reslicetop3<size, pixel_t>(temp1, temp2);
			performRotations<size, pixel_t, hash_t>(temp2, temp1, h, h0 + 4);

			//original+reslice
			//memcpy((void*)&temp1[0], (void*)&neighbourhood[0], count * sizeof(uint8_t));
			setValue(temp1, neighbourhood);
			reslicetop3<size, pixel_t>(temp1, temp2);
			reslicetop3<size, pixel_t>(temp2, temp1);
			performRotations<size, pixel_t, hash_t>(temp1, temp2, h, h0 + 8);

			//original+reslice+reslice+reslice
			//memcpy((void*)&temp1[0], (void*)&neighbourhood[0], count * sizeof(uint8_t));
			setValue(temp1, neighbourhood);
			reslicetop3<size, pixel_t>(temp1, temp2);
			reslicetop3<size, pixel_t>(temp2, temp1);
			reslicetop3<size, pixel_t>(temp1, temp2);
			performRotations<size, pixel_t, hash_t>(temp2, temp1, h, h0 + 12);

			//original+reslice right
			//memcpy((void*)&temp1[0], (void*)&neighbourhood[0], count * sizeof(uint8_t));
			setValue(temp1, neighbourhood);
			resliceright3<size, pixel_t>(temp1, temp2);
			performRotations<size, pixel_t, hash_t>(temp2, temp1, h, h0 + 16);

			//original+reslice right+reslice right+reslice right
			//memcpy((void*)&temp1[0], (void*)&neighbourhood[0], count * sizeof(uint8_t));
			setValue(temp1, neighbourhood);
			resliceright3<size, pixel_t>(temp1, temp2);
			resliceright3<size, pixel_t>(temp2, temp1);
			resliceright3<size, pixel_t>(temp1, temp2);
			performRotations<size, pixel_t, hash_t>(temp2, temp1, h, h0 + 20);
		}

		/*
		Calculates hashes of the neighbourhood and all possible (total of 48) geometric transformations of it.
		*/
		template<size_t size, typename pixel_t, typename hash_t> void hashGeometricTransformations3(Image<pixel_t>& neighbourhood, Image<pixel_t>& temp1, Image<pixel_t>& temp2, std::array<hash_t, 48>& h, size_t firstIndex)
		{
			performReslices<size, pixel_t, hash_t>(neighbourhood, temp1, temp2, h, firstIndex + 0);

			//memcpy((void*)&temp1[0], (void*)&neighbourhood[0], (size * size * size) * sizeof(uint8_t));
			setValue(temp1, neighbourhood);
			mirror3<size, pixel_t>(temp1, neighbourhood);

			performReslices<size, pixel_t, hash_t>(neighbourhood, temp1, temp2, h, firstIndex + 24);
		}
	}


	/**
	Classify size x size x size binary neighbourhood using minhash algorithm.
	Calculates all possible rotations and mirrorings of a binary neighbourhood, hashes each of the transformed neighbourhoods
	and returns minimum of the hashes.
	As a result, all geometrically similar neighbourhoods return the same value.
	'neighbourhood' must contain size*size*size pixels, zero marking background and nonzero marking foreground.
	The 'neighbourhood' image will be modified during the call to this function.
	*/
	template<size_t size, typename pixel_t, typename hash_t> hash_t minHash(Image<pixel_t>& neighbourhood)
	{
		static_assert(sizeof(hash_t) * 8 >= size * size * size, "Hash type is not wide enough to hold all possible hash values. Please use hash_t with larger sizeof(hash_t).");

		if (neighbourhood.width() != size || neighbourhood.height() != size || neighbourhood.depth() != size)
			throw ITLException("Neighbourhood size does not match what was specified.");

		Image<pixel_t> temp1(size, size, size);
		Image<pixel_t> temp2(size, size, size);

		std::array<hash_t, 48> h;

		internals::hashGeometricTransformations3<size, pixel_t, hash_t>(neighbourhood, temp1, temp2, h, 0);
		// These lines can be used to include inverted neighbourhoods in the hash.
		//internals::invert(neighbourhood, count);
		//internals::hashGeometricTransformations3<size, pixel_t, hash_t>(neighbourhood, temp1, temp2, h, 48);

		return *std::min_element(h.begin(), h.end());
	}
	
	template<size_t size, typename pixel_t, typename hash_t> hash_t minHash(const Image<pixel_t>& neighbourhood)
	{
	    Image<pixel_t> tmp;
	    setValue(tmp, neighbourhood);
	    return minHash<size, pixel_t, hash_t>(tmp);
	}


	namespace tests
	{

		/*
		Tests hashPow function.
		*/
		void hashPow();

		/*
		Tests nbhash function.
		*/
		void nbHash();

		/*
		Test method for classify3 method.
		*/
		void minHash();
	}
}
