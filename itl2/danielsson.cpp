
#include "danielsson.h"
#include "math/vec3.h"
#include "pointprocess.h"
#include "io/raw.h"
#include "conversions.h"
#include "projections.h"
#include "dmap.h"
#include "testutils.h"
#include "math/mathutils.h"
#include "math/vectoroperations.h"
#include "io/vectorio.h"

#include "thickmap.h"

#include <vector>
#include <algorithm>

using namespace std;

namespace itl2
{
// Define this to remove disk caching of lookup tables from Danielsson algorithm implementation
//#define NO_DANIELSSON_CACHE

	namespace internals
	{
		vector<coord_t> squareTable;

		/*
		Init table of largest integers whose square is less than table index.
		*/
		void initSquareTable(size_t maxSquare)
		{
			if (squareTable.size() < maxSquare)
			{
				squareTable.reserve(maxSquare);
				while (squareTable.size() <= maxSquare)
					squareTable.push_back(largestIntWhoseSquareIsLessThan(squareTable.size()));
			}
		}

		/*
		Get largest integer whose square is less than given value. To use this function, the cache must be first initialized with initSquareTable.
		*/
		coord_t largestIntWhoseSquareIsLessThanCached(coord_t square)
		{
			return squareTable[square];
		}


		/*
		Helper for getMaxSphereRadius.
		Tests if a sphere of squared radius Rdot2 centered at (cx, cy, cz) fits inside a sphere of radius R2 centered at origin.
		*/
		bool testFit(coord_t cx, coord_t cy, coord_t cz, coord_t R2, coord_t size, coord_t Rdot2)
		{
			for (coord_t z = 0; z < size; z++)
			{
				coord_t dz = z - cz;
				for (int32_t y = 0; y < size; y++)
				{
					// Squared x-coordinate of the surface of a sphere centered at (0, cy, cz) and having squared radiut Rdot2.
					coord_t dy = y - cy;
					coord_t test = Rdot2 - dy * dy - dz * dz;

					if (test >= 0)
					{
						test = largestIntWhoseSquareIsLessThanCached(test) + cx;


						// Squared x-coordinate of the surface of a sphere centered at origin and having squared radius R2.
						coord_t Rx = R2 - y * y - z * z;

						if (Rx >= 0)
							Rx = largestIntWhoseSquareIsLessThanCached(Rx);
						else
							Rx = -1;

						if (test > Rx)
						{
							// The sphere does not fit
							return false;
						}
					}
				}
			}

			// Here sphere with radius of Rdot fits inside the central sphere.
			return true;
		}



		coord_t search(coord_t R2Index, const vector<coord_t>& radii2, coord_t R2, coord_t cx, coord_t cy, coord_t cz, coord_t size)
		{
			// Binary search for last index for which testFit(... R2, size, radii2[index]) gives true.
			// This handles also -1 elements in the radii2 array.

			coord_t first = R2Index;
			while (first > 0 && (radii2[first] < 0 || largestIntWhoseSquareIsLessThanCached(radii2[first]) >= largestIntWhoseSquareIsLessThanCached(R2) - 2))
				first--;

			coord_t last = R2Index;

			while (last - first > 1)
			{
				coord_t mid = (last + first) / 2;

				// Handle -1 elements in the array
				while (mid < last && radii2[mid] < 0)
					mid++;

				if (mid >= last)
				{
					// last is the first non -1 element after (and including) mid.
					last = (last + first) / 2;
				}
				else
				{

					coord_t Rdot2 = radii2[mid];

					if (testFit(cx, cy, cz, R2, size, Rdot2))
					{
						// Fits. Set start = mid
						first = mid;
					}
					else
					{
						// Does not fit. Set end = mid
						last = mid;
					}
				}
			}

			return first;
		}



		/*
		Calculates squared radius of largest sphere that is centered at (cx, cy, cz) and fits inside sphere of radius sqrt(r2) centered at (0, 0, 0).
		Works for positive (cx, cy, cz) only.
		indices lookup table gives index for squared radius.
		radii lookup table gives radius for index.
		NOTE: This is optimized version that does not use separate mask array, and does only 1/8 of processing of the unoptimized version.
		*/
		coord_t getMaxSphereRadius(coord_t cx, coord_t cy, coord_t cz, const vector<coord_t>& radii2, coord_t R2Index)
		{
			coord_t R2 = radii2[R2Index];
			if (R2 < 0)
				throw ITLException("Determining Danielsson tables for squared radius that is not a sum of three squares.");

			
			coord_t Rint = itl2::ceil(sqrt(R2));
			coord_t size = Rint + 1;


			// Find the largest r that still fits inside the central sphere
			// by testing mask of sphere centered at (cx, cy, cz).
			// This uses linear search. It's a little bit slower than the search(..) function below.
			//coord_t startInd = R2Index;
			//coord_t trial = 0;
			//for (coord_t RdotInd = startInd - 1; RdotInd >= 0; RdotInd--)
			//{
			//	coord_t Rdot2 = radii2[RdotInd];

			//	if (Rdot2 >= 0)
			//	{
			//		if (testFit(cx, cy, cz, R2, size, Rdot2))
			//		{
			//			trial = Rdot2;
			//			break;
			//		}
			//	}
			//}

			coord_t trial = search(R2Index, radii2, R2, cx, cy, cz, size);

			return trial;
		}


		

		///*
		//Calculates squared radius of largest sphere that is centered at (cx, cy, cz) and fits inside sphere of radius sqrt(r2) centered at (0, 0, 0).
		//Works for positive (cx, cy, cz) only.
		//indices lookup table gives index for squared radius.
		//radii lookup table gives radius for index.
		//NOTE: This is optimized version and does only 1/8 of processing of the unoptimized version.
		//*/
		//int32_t getMaxSphereRadius(int32_t cx, int32_t cy, int32_t cz, const vector<int32_t>& radii2, coord_t R2Index)
		//{
		//	int32_t R2 = radii2[R2Index];
		//	if (R2 < 0)
		//		throw ITLException("Determining Danielsson tables for squared radius that is not a sum of three squares.");

		//	// Create mask of sphere centered at origin, having radius R
		//	int32_t Rint = (int32_t)ceil(sqrt(R2));

		//	coord_t size = Rint + 1;
		//	Image<int32_t> xMask(size, size, 1, -1);
		//	for (int32_t dz = 0; dz < xMask.height(); dz++)
		//	{
		//		for (int32_t dy = 0; dy < xMask.width(); dy++)
		//		{
		//			int32_t dx2 = R2 - dy * dy - dz * dz;
		//			if (dx2 >= 0)
		//			{
		//				xMask(dy, dz) = largestIntWhoseSquareIsLessThan(dx2);
		//			}
		//		}
		//	}

		//	//raw::writed(xMask, "./danielsson/xMask");


		//	// Find the largest r that still fits inside the central sphere
		//	// by testing mask of sphere centered at (cx, cy, cz).
		//	coord_t startInd = R2Index;
		//	for (coord_t RdotInd = startInd - 1; RdotInd >= 0; RdotInd--)
		//	{
		//		int32_t Rdot2 = radii2[RdotInd];

		//		if (Rdot2 >= 0)
		//		{
		//			for (int32_t z = 0; z < xMask.height(); z++)
		//			{
		//				int32_t dz = z - cz;
		//				for (int32_t y = 0; y < xMask.width(); y++)
		//				{
		//					int32_t dy = y - cy;
		//					int32_t test = Rdot2 - dy * dy - dz * dz;

		//					if (test >= 0)
		//					{
		//						test = largestIntWhoseSquareIsLessThan(test) + cx;

		//						int32_t comp = xMask(y, z);

		//						if (test > comp)
		//						{
		//							// The sphere does not fit, we can skip processing of both inner loops.
		//							goto doesNotFit;
		//						}
		//					}
		//				}
		//			}

		//			// Here sphere with radius of Rdot fits inside the central sphere.
		//			return Rdot2;
		//		}

		//	doesNotFit:;
		//	}

		//	return 0;
		//}


		///*
		//Calculates squared radius of largest sphere that is centered at (cx, cy, cz) and fits inside sphere of radius sqrt(r2) centered at (0, 0, 0).
		//Works for positive (cx, cy, cz) only.
		//indices lookup table gives index for squared radius.
		//radii lookup table gives radius for index.
		//NOTE: This is an unoptimized but working version.
		//*/
		//int32_t getMaxSphereRadius(int32_t R2, coord_t cx, coord_t cy, coord_t cz, const vector<coord_t>& indices, const vector<int32_t>& radii2)
		//{
		//	// Create mask of sphere centered at origin, having radius R
		//	int32_t Rint = (int32_t)ceil(sqrt(R2));

		//	coord_t size = 2 * Rint + 1;
		//	Image<int32_t> xMask(size, size, 1, -1);
		//	for (int32_t z = 0; z < xMask.height(); z++)
		//	{
		//		int32_t dz = z - Rint;
		//		for (int32_t y = 0; y < xMask.width(); y++)
		//		{
		//			int32_t dy = y - Rint;
		//			int32_t dx2 = R2 - dy * dy - dz * dz;
		//			if (dx2 >= 0)
		//			{
		//				xMask(y, z) = largestIntWhoseSquareIsLessThan(dx2);
		//			}
		//		}
		//	}

		//	// Find the largest r that still fits inside the central sphere
		//	// by testing mask of sphere centered at (cx, cy, cz).
		//	coord_t startInd = indices[R2];
		//	for (coord_t RdotInd = startInd - 1; RdotInd >= 0; RdotInd--)
		//	{
		//		int32_t Rdot2 = radii2[RdotInd];

		//		for (int32_t z = 0; z < xMask.height(); z++)
		//		{
		//			int32_t dz = z - Rint - cz;
		//			for (int32_t y = 0; y < xMask.width(); y++)
		//			{
		//				int32_t dy = y - Rint - cy;
		//				int32_t test = Rdot2 - dy * dy - dz * dz;

		//				if (test >= 0)
		//				{
		//					test = largestIntWhoseSquareIsLessThan(test) + cx;

		//					int32_t comp = xMask(y, z);

		//					if (test > comp)
		//					{
		//						// The sphere does not fit, we can skip processing of both inner loops.
		//						goto doesNotFit;
		//					}
		//				}
		//			}
		//		}

		//		// Here sphere with radius of Rdot fits inside the central sphere.
		//		return Rdot2;

		//	doesNotFit:;
		//	}

		//	return 0;
		//}
	}



#if defined(NO_DANIELSSON_CACHE)

	namespace internals
	{
		/*
		Get a value from 'table' lookup table that corresponds to radius 'dist'.
		*/
		int32_t get(const Image<int32_t>& table, int32_t dist2, const vector<coord_t>& indices)
		{
			coord_t ind = indices[dist2];
			return table(ind);
		}
	}

	/**
	Calculates centers of locally maximal disks using Danielsson algorithm.
	NOTE: This version works but does not use cached Danielsson lookup tables.
	See e.g.
	Yaorong Ge and J. Michael Fitzpatrick - On the Generation of Skeletons from Discrete Euclidean Distance Maps
	@param dmap2 Squared Euclidean distance map of the input geometry.
	@param out Squared distance values of pixels that correspond to centers of locally maximal disks are set to this image. Other values are not changed. Normally pass empty image.
	*/
	void centersOfLocallyMaximalSpheres(const Image<int32_t>& dmap2, Image<int32_t>& out)
	{
		out.mustNotBe(dmap2);
		out.ensureSize(dmap2);


		// Find maximum squared distance value
		cout << "Finding maximum radius..." << endl;
		int32_t maxr2 = max(dmap2);

		cout << "Maximal radius = " << sqrt(maxr2) << endl;

		// Create lookup table that converts squared radius to index in lookup table,
		// i.e., indices[squared radius] = index in other lookup tables.
		// Smallest squared radius gets index 0, second smallest gets index 1, etc.
		// The index is used to access lookup tables.
		cout << "Finding unique squared radius values and building lookup tables..." << endl;
		vector<coord_t> indices(maxr2 + 1, -1);

		// First mark which indices appear in the image
#pragma omp parallel for if(dmap2.pixelCount() >= PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t n = 0; n < dmap2.pixelCount(); n++)
		{
			int32_t r2 = dmap2(n);
			indices[r2] = 0;
		}

		// Then create index for each occurence
		coord_t radiusCount = 0;
		for (coord_t n = 0; n < (coord_t)indices.size(); n++)
		{
			if (indices[n] == 0)
			{
				indices[n] = radiusCount;
				radiusCount++;
			}
		}

		// Create list of diameter values in the image
		vector<int32_t> radii2(radiusCount);
#pragma omp parallel for if(indices.size() >= PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t n = 0; n < (coord_t)indices.size(); n++)
		{
			coord_t index = indices[n];
			if (index >= 0)
			{
				radii2[index] = (int32_t)n;
			}
		}

		cout << "Found " << radiusCount << " unique radii." << endl;

		//raw::writed(radii2, "radii2");
		//raw::writed(indices, "indices");

		// Fill lookup table for three types of neighbours.
		cout << "Calculating Danielsson tables..." << endl;
		Image<int32_t> table1(radiusCount);
		Image<int32_t> table2(radiusCount);
		Image<int32_t> table3(radiusCount);

		size_t totalProcessed = 0;
#pragma omp parallel for if(3*radiusCount >= PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t n = 0; n < radiusCount; n++)
		{
			int32_t r2 = radii2[n];
			//table1(n) = internals::getMaxSphereRadius(r2, 1, 0, 0, indices, radii2);
			//table2(n) = internals::getMaxSphereRadius(r2, 1, 1, 0, indices, radii2);
			//table3(n) = internals::getMaxSphereRadius(r2, 1, 1, 1, indices, radii2);
			coord_t index = indices[r2];
			table1(n) = internals::getMaxSphereRadius(1, 0, 0, radii2, index);
			table2(n) = internals::getMaxSphereRadius(1, 1, 0, radii2, index);
			table3(n) = internals::getMaxSphereRadius(1, 1, 1, radii2, index);

			showThreadProgress(totalProcessed, radiusCount);
		}

		cout << "Processing image..." << endl;
		totalProcessed = 0;
		//#pragma omp parallel if(dmap.pixelCount() >= PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			Image<int32_t> nb(3, 3, 3);

			// Process all points in the image
			//#pragma omp for
			for (coord_t z = 0; z < dmap2.depth(); z++)
			{
				for (coord_t y = 0; y < dmap2.height(); y++)
				{
					for (coord_t x = 0; x < dmap2.width(); x++)
					{
						int32_t c = dmap2(x, y, z);
						if (c != 0)
						{
							getNeighbourhood(dmap2, Vec3c(x, y, z), Vec3c(1, 1, 1), nb, BoundaryCondition::Zero);

							// Check all neighbours
							if (!(// 6-neighbours, one coordinate changes by one pixel.
								internals::get(table1, nb(0, 1, 1), indices) >= c ||
								internals::get(table1, nb(2, 1, 1), indices) >= c ||
								internals::get(table1, nb(1, 0, 1), indices) >= c ||
								internals::get(table1, nb(1, 2, 1), indices) >= c ||
								internals::get(table1, nb(1, 1, 0), indices) >= c ||
								internals::get(table1, nb(1, 1, 2), indices) >= c ||
								// 18-neighbours but not 6-neighbours, two coordinates change by one pixel.
								internals::get(table2, nb(0, 0, 1), indices) >= c ||
								internals::get(table2, nb(0, 2, 1), indices) >= c ||
								internals::get(table2, nb(2, 0, 1), indices) >= c ||
								internals::get(table2, nb(2, 2, 1), indices) >= c ||
								internals::get(table2, nb(1, 0, 0), indices) >= c ||
								internals::get(table2, nb(1, 0, 2), indices) >= c ||
								internals::get(table2, nb(1, 2, 0), indices) >= c ||
								internals::get(table2, nb(1, 2, 2), indices) >= c ||
								internals::get(table2, nb(0, 1, 0), indices) >= c ||
								internals::get(table2, nb(0, 1, 2), indices) >= c ||
								internals::get(table2, nb(2, 1, 0), indices) >= c ||
								internals::get(table2, nb(2, 1, 2), indices) >= c ||
								// Corners, three coordinates change by one pixel.
								internals::get(table3, nb(0, 0, 0), indices) >= c ||
								internals::get(table3, nb(0, 2, 0), indices) >= c ||
								internals::get(table3, nb(2, 0, 0), indices) >= c ||
								internals::get(table3, nb(2, 2, 0), indices) >= c ||
								internals::get(table3, nb(0, 0, 2), indices) >= c ||
								internals::get(table3, nb(0, 2, 2), indices) >= c ||
								internals::get(table3, nb(2, 0, 2), indices) >= c ||
								internals::get(table3, nb(2, 2, 2), indices) >= c
								))
							{
								// This is center of locally maximal disk
								out(x, y, z) = c;
							}
						}
					}
				}

				showThreadProgress(totalProcessed, dmap2.depth());
			}
		}
	}
#else


	namespace internals
	{

		/**
		Expands Danielsson lookup tables so that they cover at least squared radius R2max.
		*/
		void expandDanielssonTables(vector<coord_t>& table1, vector<coord_t>& table2, vector<coord_t>& table3, coord_t R2max)
		{
			// Test if the tables are already complete.
			if ((coord_t)table1.size() >= R2max)
				return;

			// Determine which elements of the tables need to be determined
			coord_t Rmax = largestIntWhoseSquareIsLessThan(R2max) + 1;
			constexpr coord_t invalidValue = -1;
			constexpr coord_t toBeDeterminedValue = numeric_limits<coord_t>::max();
			{
				ProgressIndicator progress(Rmax);
				for (coord_t z = 0; z < Rmax; z++)
				{
					for (coord_t y = 0; y < Rmax; y++)
					{
						for (coord_t x = 0; x < Rmax; x++)
						{
							coord_t R2 = x * x + y * y + z * z;
							if (R2 <= R2max)
							{
								// radius is less than maximum, add it to the tables.
								while ((coord_t)table1.size() <= R2)
								{
									table1.push_back(invalidValue);
									table2.push_back(invalidValue);
									table3.push_back(invalidValue);
								}

								if (table1[R2] == invalidValue)
								{
									table1[R2] = toBeDeterminedValue;
									table2[R2] = toBeDeterminedValue;
									table3[R2] = toBeDeterminedValue;
								}
							}
						}
					}

					progress.step();
				}
			}

			// Create a list of possible radius^2 values.
			vector<coord_t> radii2;
			for (coord_t r2 = 0; r2 < (coord_t)table1.size(); r2++)
			{
				if (table1[r2] != invalidValue) // Value is possible if it is not marked as invalid in the tables.
					radii2.push_back(r2);
				else
					radii2.push_back(-1);
			}



			// Calculate count of items
			//cout << "Total size of tables = 3*" << table1.size() << endl;

			//size_t totalCount = 0;
			//size_t toBeDetermined = 0;
			//for (size_t n = 0; n < table1.size(); n++)
			//{
			//	if (table1[n] != invalidValue)
			//		totalCount++;
			//	if (table1[n] == toBeDeterminedValue)
			//		toBeDetermined++;
			//}

			//cout << "Count of distinct R^2 values = " << totalCount << " when R < " << Rmax << endl;
			//cout << "Count of values whose table entries are to be determined = " << toBeDetermined << endl;
			
			initSquareTable(max(radii2));

			// Calculate the missing entries
			{
				coord_t tableSize = (coord_t)table1.size();
				ProgressIndicator progress(tableSize);
#pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
				for (coord_t r2 = 0; r2 < tableSize; r2++)
				{
					if (table1[r2] == toBeDeterminedValue)
					{
						table1[r2] = internals::getMaxSphereRadius(1, 0, 0, radii2, r2);
						table2[r2] = internals::getMaxSphereRadius(1, 1, 0, radii2, r2);
						table3[r2] = internals::getMaxSphereRadius(1, 1, 1, radii2, r2);
					}
					progress.step();
				}
			}
		}

		/**
		Reads lookup table from disk.
		*/
		void readTable(vector<coord_t>& table, const string& filename)
		{
			try
			{
				readListFile(filename, table);
			}
			catch (const ITLException&)
			{
				// The file does not exist - just return without loading the table
				return;
			}
			//std::ifstream in(filename, ios_base::in | ios_base::binary);
			//if (!in.is_open())
			//	return;

			//while (in.good())
			//{
			//	coord_t value;
			//	in.read((char*)&value, sizeof(coord_t));
			//	if (in.eof())
			//		break;
			//	table.push_back(value);
			//}
		}

		/**
		Writes lookup table to disk.
		*/
		void writeTable(vector<coord_t>& table, const string& filename)
		{
			writeListFile(filename, table);
			//createFoldersFor(filename);
			//std::ofstream out(filename, ios_base::out | ios_base::trunc | ios_base::binary);
			//if (!out.is_open())
			//	return;

			//for (size_t n = 0; n < table.size(); n++)
			//	out.write((char*)&table[n], sizeof(coord_t));
		}

		/**
		Creates Danielsson lookup tables that cover at least squared radius R2max.
		Uses old tables if they are found.
		*/
		void getDanielssonTables(vector<coord_t>& table1, vector<coord_t>& table2, vector<coord_t>& table3, coord_t R2max)
		{
			readTable(table1, "danielsson_table_1.dat");
			readTable(table2, "danielsson_table_2.dat");
			readTable(table3, "danielsson_table_3.dat");

			if ((coord_t)table1.size() < R2max)
			{
				expandDanielssonTables(table1, table2, table3, R2max);
				writeTable(table1, "danielsson_table_1.dat");
				writeTable(table2, "danielsson_table_2.dat");
				writeTable(table3, "danielsson_table_3.dat");
			}
		}
	}





	namespace tests
	{
		void danielssonTableSpeedTest()
		{
			coord_t R2max = 150 * 150;
			//coord_t R2max = 300 * 300;

			// Get tables from disk. This assumes they have been saved before.
			vector<coord_t> table1GT, table2GT, table3GT;
			internals::getDanielssonTables(table1GT, table2GT, table3GT, R2max);

			// Measure how long it takes to calculate the tables
			vector<coord_t> table1, table2, table3;
			
			Timer timer;
			timer.start();
			internals::expandDanielssonTables(table1, table2, table3, R2max);
			timer.stop();
			cout << "Calculation took " << timer.getTime() << " ms" << endl;

			// Initial version: 15 s for R2max = 150^2
			// Binary search: 12 s
			// With int square cache + linear search: 2.1 s
			// With int square cache + binary search: 1.8 s

			// Compare tables
			testAssert(table1 == table1GT, "table 1");
			testAssert(table2 == table2GT, "table 2");
			testAssert(table3 == table3GT, "table 3");
		}


		void fullDanielssonTables()
		{

			vector<coord_t> testTable;
			testTable.push_back(1);
			testTable.push_back(2);
			testTable.push_back(3);
			testTable.push_back(4);
			testTable.push_back(8);
			testTable.push_back('\n');
			testTable.push_back(100);
			testTable.push_back('\r');
			testTable.push_back(200);
			testTable.push_back(numeric_limits<int32_t>::max());
			testTable.push_back(-200);
			internals::writeTable(testTable, "./danielsson/testtable.dat");

			vector<coord_t> readTable;
			internals::readTable(readTable, "./danielsson/testtable.dat");

			testAssert(testTable == readTable, "danielsson IO");

			vector<coord_t> table1PartialFill, table2PartialFill, table3PartialFill;
			internals::expandDanielssonTables(table1PartialFill, table2PartialFill, table3PartialFill, 10 * 10);
			internals::expandDanielssonTables(table1PartialFill, table2PartialFill, table3PartialFill, 10 * 10);
			internals::expandDanielssonTables(table1PartialFill, table2PartialFill, table3PartialFill, 11 * 11);
			internals::expandDanielssonTables(table1PartialFill, table2PartialFill, table3PartialFill, 50 * 50);

			vector<coord_t> table1, table2, table3;
			internals::expandDanielssonTables(table1, table2, table3, 50 * 50);

			testAssert(table1.size() >= 50 * 50, "Danielsson table size");

			testAssert(table1 == table1PartialFill, "Danielsson tables differ between full and partial fill");
			testAssert(table2 == table2PartialFill, "Danielsson tables differ between full and partial fill");
			testAssert(table3 == table3PartialFill, "Danielsson tables differ between full and partial fill");


			vector<coord_t> table1Cached, table2Cached, table3Cached;
			internals::getDanielssonTables(table1Cached, table2Cached, table3Cached, 10 * 10);
			table1Cached.clear();
			table2Cached.clear();
			table3Cached.clear();
			internals::getDanielssonTables(table1Cached, table2Cached, table3Cached, 10 * 10);
			table1Cached.clear();
			table2Cached.clear();
			table3Cached.clear();
			internals::getDanielssonTables(table1Cached, table2Cached, table3Cached, 11 * 11);
			table1Cached.clear();
			table2Cached.clear();
			table3Cached.clear();
			internals::getDanielssonTables(table1Cached, table2Cached, table3Cached, 50 * 50);

			// If cached tables are larger than expected, erase the extra items.
			if (table1Cached.size() > table1.size())
				table1Cached.erase(table1Cached.begin() + table1.size(), table1Cached.end());
			if (table2Cached.size() > table2.size())
				table2Cached.erase(table2Cached.begin() + table2.size(), table2Cached.end());
			if (table3Cached.size() > table3.size())
				table3Cached.erase(table3Cached.begin() + table3.size(), table3Cached.end());

			testAssert(table1 == table1Cached, "Danielsson tables differ between normal and cached fills");
			testAssert(table2 == table2Cached, "Danielsson tables differ between normal and cached fills");
			testAssert(table3 == table3Cached, "Danielsson tables differ between normal and cached fills");

		}
	}

#endif


	namespace tests
	{

		template<typename dmap_t> void singleDanielssonTest()
		{
			Image<float32_t> geom;
			//raw::read(geom, "../test_input_data/simple_structures_128x128x128.raw");
			raw::read(geom, "../test_input_data/simple_structures_256x256x256.raw");

			double nonzeroCount = 0;
			for (coord_t n = 0; n < geom.pixelCount(); n++)
				if (geom(n) > 0)
					nonzeroCount++;

			Image<dmap_t> dmap2;
			distanceTransform2(geom, dmap2);
			raw::writed(dmap2, "./danielsson/dmap2");



			Image<dmap_t> ridgeApprox2;
			centersOfLocallyMaximalSpheres(dmap2, ridgeApprox2);

			raw::writed(ridgeApprox2, "./danielsson/ridge2");

			// Calculate thickness map from ridge
			cout << "From ridge:" << endl;
			Image<dmap_t> thicknessFromRidge;
			
			Vec3d sum;
			itl2::dimredsuper::thickmap2(ridgeApprox2, thicknessFromRidge, nullptr, &sum);
			raw::writed(thicknessFromRidge, "./danielsson/local_thickness_from_ridge");
			cout << "Number of ri items per fg pixel after dimension 1: " << sum[0] / nonzeroCount << endl;
			cout << "Number of ri items per fg pixel after dimension 2: " << sum[1] / nonzeroCount << endl;
			cout << "Number of ri items per fg pixel after dimension 3: " << sum[2] / nonzeroCount << endl;





			// Calculate thickness map from distance map
			cout << "From distance map:" << endl;
			Image<dmap_t> thicknessFromDmap;
			itl2::dimredsuper::thickmap2(dmap2, thicknessFromDmap, nullptr, &sum);
			raw::writed(thicknessFromDmap, "./danielsson/local_thickness_from_dmap");

			cout << "Number of ri items per fg pixel after dimension 1: " << sum[0] / nonzeroCount << endl;
			cout << "Number of ri items per fg pixel after dimension 2: " << sum[1] / nonzeroCount << endl;
			cout << "Number of ri items per fg pixel after dimension 3: " << sum[2] / nonzeroCount << endl;


			cout << "Sanity check with simple algorithm..." << endl;
			Image<dmap_t> thicknessSimple;
			setValue(thicknessSimple, dmap2);
			itl2::optimized::thickmap2(thicknessSimple);
			raw::writed(thicknessSimple, "./danielsson/local_thickness_simple_dmap");


			// Check that thickness maps are the same
			checkDifference(thicknessFromDmap, thicknessFromRidge, "difference between local thickness calculated from distance map and from distance ridge");
			checkDifference(thicknessSimple, thicknessFromRidge, "difference between local thickness calculated using simple algorithm and dimred");


			// Check that both thickness maps have nonzero pixels at the same locations than in the original
			threshold(thicknessFromRidge, 0);
			checkDifference(thicknessFromRidge, geom, "difference between original geometry and nonzero points of thickness map from distance ridge.");

			threshold(thicknessFromDmap, 0);
			checkDifference(thicknessFromDmap, geom, "difference between original geometry and nonzero points of thickness map from distance map.");

		}

		void danielsson()
		{
			singleDanielssonTest<int32_t>();
			singleDanielssonTest<uint32_t>();
			singleDanielssonTest<float32_t>();
		}
	}
}