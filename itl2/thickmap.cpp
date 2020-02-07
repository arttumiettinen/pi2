
#include "thickmap.h"

#include "io/raw.h"

#include "sphere.h"
#include "generation.h"
#include "dmap.h"
#include "testutils.h"
#include "danielsson.h"
#include "transform.h"

#include <map>

using namespace std;



namespace itl2
{

	namespace internals
	{
		/**
		Tests if discretized circle whose squared radius is r1Square fit completely into discretized circle whose squared radius is r2Square.
		*/
		bool doesDiscretizedCircle1FitInto2(int32_t r1Square, int32_t r2Square)
		{
			if (r1Square == r2Square)
				return true;

			int32_t r1 = largestIntWhoseSquareIsLessThan(r1Square);
			int32_t r2 = largestIntWhoseSquareIsLessThan(r2Square);

			// If the r1 > r2, the projection of the circle 1 to x-axis completely overlaps projection of circle 2.
			// Circle 1 cannot fit into circle 2.
			if (r1 > r2)
				return false;

			// If r1 < r2, the circle 1 surely fits into circle 2.
			if (r1 < r2)
				return true;

			// Projection of both circles to x-axis is the same. Now test if the circles are the same in each
			// y-directional pixel column.
			for (int32_t x = 0; x <= r1; x++)
			{
				int32_t y1Square = r1Square - x * x;
				int32_t y2Square = r2Square - x * x;
				if (largestIntWhoseSquareIsLessThan(y1Square) > largestIntWhoseSquareIsLessThan(y2Square))
					return false;
			}

			return true;
		}

		/**
		Lookup table for doesDiscretizedCircle1FitInto2Cached.
		Element circleLookup[r2] stores the maximal squared radius of a circle that fits into a circle of squared radius r2.
		*/
		vector<int32_t> circleLookup;

		/**
		Builds lookup table that doesDiscretizedCircle1FitInto2Cached uses.
		@param maxrSquare Maximum squared radius found in the image.
		*/
		void buildCircleLookup(int32_t maxrSquare)
		{
			int32_t firstr2 = (int32_t)circleLookup.size();
			circleLookup.resize(maxrSquare + 1);

			for (int32_t r2 = firstr2; r2 < maxrSquare + 1; r2++)
			{
				// Find the maximal squared radius of a circle that fits into a circle of squared radius r2.
				circleLookup[r2] = r2;
				int32_t rdot2 = r2;
				while (true)
				{
					rdot2++;
					if (!doesDiscretizedCircle1FitInto2(rdot2, r2))
						break;
					circleLookup[r2] = rdot2;
				}
			}
		}
	}

















	namespace internals
	{
		/**
		Draws maximal sphere corresponding to center point and squared radius.
		*/
		inline size_t drawMax2(Image<int32_t>& image, const Vec3c& center, int32_t r2)
		{
			Vec3c minPos = round(Vec3d(center) - sqrt(r2) * Vec3d(1, 1, 1));
			Vec3c maxPos = round(Vec3d(center) + sqrt(r2) * Vec3d(1, 1, 1));

			clamp(minPos, Vec3c(0, 0, 0), image.dimensions() - Vec3c(1, 1, 1));
			clamp(maxPos, Vec3c(0, 0, 0), image.dimensions() - Vec3c(1, 1, 1));

			size_t filledCount = 0;
			//#pragma omp parallel for if(!omp_in_parallel() && AABox(minPos, maxPos).volume() > PARALLELIZATION_THRESHOLD) reduction(+:filledCount)
			for (coord_t z = minPos.z; z <= maxPos.z; z++)
			{
				for (coord_t y = minPos.y; y <= maxPos.y; y++)
				{
					for (coord_t x = minPos.x; x <= maxPos.x; x++)
					{
						Vec3c d = center - Vec3c(x, y, z);

						if (d.x * d.x + d.y * d.y + d.z * d.z < r2)
						{
							int32_t& p = image(x, y, z);
							p = std::max(p, r2);
							filledCount++;
						}
					}
				}
			}

			return filledCount;
		}
	}


	namespace standard
	{
		/**
		Calculates squared local radius from squared distance map
		Uses standard Hildebrand & Ruegsegger algorithm.
		Plots maximal spheres corresponding to squared distance map, larger distance values replacing smaller ones.
		*/
		void thickmap2(const Image<int32_t>& dmap2, Image<int32_t>& result, bool showProgressInfo)
		{
			result.mustNotBe(dmap2);
			result.ensureSize(dmap2);

			size_t counter = 0;
			for (coord_t z = 0; z < result.depth(); z++)
			{
				for (coord_t y = 0; y < result.height(); y++)
				{
					for (coord_t x = 0; x < result.width(); x++)
					{
						int32_t r2 = dmap2(x, y, z);
						internals::drawMax2(result, Vec3c(x, y, z), r2);
					}
				}
				showThreadProgress(counter, result.depth(), showProgressInfo);
			}
		}
	}







	namespace internals
	{
#pragma pack(push, 1)
		/*
		Stores 16-bit block index and 48-bit start index in 64-bits.
		*/
		struct IndexItem
		{
			uint64_t data;

			uint16_t getBlockIndex() const
			{
				return (uint16_t)(data & 0x000000000000ffff);
			}

			uint64_t getStartIndex() const
			{
				return (data & 0xffffffffffff0000) >> 16;
			}

			void setBlockAndStart(uint16_t blockIndex, uint64_t startIndex)
			{
				data = (uint64_t)blockIndex + ((startIndex & 0xffffffffffff) << 16);
			}

		};
#pragma pack(pop)


		string createDatFileName(const string& indexFile, uint16_t blockIndex)
		{
			return indexFile + "_block" + toString(blockIndex) + ".dat";
		}

		/**
		Writes ri image block to file.
		*/
		void writeRiBlock(Image<internals::RiStorageSet>& ri, const string& indexFilePrefix, uint16_t blockIndex, const Vec3c& filePosition, const Vec3c& fileDimensions)
		{
			Image<IndexItem> index(ri.dimensions());

			string indexFile = concatDimensions(indexFilePrefix, fileDimensions);

			createFoldersFor(indexFile);

			string blockFileName = createDatFileName(indexFilePrefix, blockIndex);

			// Create ofstream with large buffer size.
			std::vector<char> buf;
			buf.resize(1024 * 1024 * 10);
			ofstream data;
			data.rdbuf()->pubsetbuf(&buf.front(), buf.size());
			data.open(blockFileName.c_str(), ios_base::out | ios_base::binary);

			if (!data)
				throw ITLException(string("Unable to create output file ") + blockFileName);

			uint64_t startIndex = 0;
			for (coord_t z = 0; z < ri.depth(); z++)
			{
				for (coord_t y = 0; y < ri.height(); y++)
				{
					for (coord_t x = 0; x < ri.width(); x++)
					{
						index(x, y, z).setBlockAndStart(blockIndex, startIndex);

						const internals::RiStorageSet& s = ri(x, y, z);

						// Write size
						uint16_t val = s.size();
						data.write((char*)&val, sizeof(uint16_t));

						// Write items
						for (uint16_t m = 0; m < s.size(); m++)
						{
							val = s[m].srcX;
							data.write((char*)&val, sizeof(uint16_t));
							
							val = s[m].srcY;
							data.write((char*)&val, sizeof(uint16_t));
						}

						startIndex += 2 * s.size() + 1;
					}
				}
			}

			raw::writeBlock(index, indexFile, filePosition, fileDimensions);
		}

		/**
		Reads ri image block from file.
		*/
		void readRiBlock(Image<internals::RiStorageSet>& ri, std::string indexFilePrefix, const Vec3c& start)
		{
			// index(x, y, z) gives the start index to read from the data file dat.
			// dat[index(x, y, z)] gives the count of elements.
			Image<IndexItem> index(ri.dimensions());
			raw::readBlock(index, indexFilePrefix, start);

			map<uint16_t, unique_ptr<DiskMappedBuffer<int16_t> > > datFiles;
			for (coord_t z = 0; z < ri.depth(); z++)
			{
				for (coord_t y = 0; y < ri.height(); y++)
				{
					for (coord_t x = 0; x < ri.width(); x++)
					{
						const IndexItem& startItem = index(x, y, z);
						uint16_t blockIndex = startItem.getBlockIndex();
						uint64_t startIndex = startItem.getStartIndex();

						// Map dat file if it is not open.
						auto it = datFiles.find(blockIndex);
						if (it == datFiles.end())
						{
							string datFileName = createDatFileName(indexFilePrefix, blockIndex);
							it = datFiles.insert(it, make_pair(blockIndex, make_unique<DiskMappedBuffer<int16_t> >(0, datFileName, true)));
						}

						const int16_t* dat = it->second->getBufferPointer();

						uint16_t count = reinterpret_cast<const uint16_t*>(dat)[startIndex];

						vector<RiStorageItem> vals;
						vals.reserve(count);
						for (uint16_t i = 0; i < count; i++)
						{
							RiStorageItem item;
							item.srcX = dat[(startIndex + 1) + 2 * i];
							item.srcY = dat[(startIndex + 1) + 2 * i + 1];
							vals.push_back(item);
						}
						ri(x, y, z).setValues(vals);
					}
				}

			}
		}
	}





	namespace tests
	{

		/**
		Tests IndexItem class.
		*/
		void indexItem()
		{
			internals::IndexItem item;

			item.setBlockAndStart(1, 7);

			testAssert(item.getBlockIndex() == 1, "block index");
			testAssert(item.getStartIndex() == 7, "start index");


			item.setBlockAndStart(12, 56373);

			testAssert(item.getBlockIndex() == 12, "block index");
			testAssert(item.getStartIndex() == 56373, "start index");


			item.setBlockAndStart(858, 478787);

			testAssert(item.getBlockIndex() == 858, "block index");
			testAssert(item.getStartIndex() == 478787, "start index");
			
			
			// There was some interesting bug in either IndexItem or GCC:
			// block and start index could be set in stand alone declared index item,
			// but block index could not be set if the index item was in a container.
			// This code was used to find the bug.
			Image<internals::IndexItem> index(10, 1, 1);
	        uint16_t blockIndex = 7;
	        uint64_t startIndex = 0;
	        coord_t y = 0;
	        coord_t z = 0;
	        internals::IndexItem standalone;
	        vector<internals::IndexItem> v;
	        for (coord_t x = 0; x < index.width(); x++)
	        {
	            v.push_back(internals::IndexItem());
	        
		        index(x, y, z).setBlockAndStart(blockIndex, startIndex);
		        
	            v[x].setBlockAndStart(blockIndex, startIndex);
		        
		        standalone.setBlockAndStart(blockIndex, startIndex);

                testAssert(standalone.data == index(x, y, z).data && standalone.data == v[x].data, "data difference 1");
                testAssert(v[x].getStartIndex() == startIndex && v[x].getBlockIndex() == blockIndex, "vector get* difference");
				testAssert(index(x, y, z).getStartIndex() == startIndex && index(x, y, z).getBlockIndex() == blockIndex, "vector get* difference");

				//cout << standalone.data << " == " << index(x, y, z).data << " == " << v[x].data << endl;
				//cout << "vector:\t\tstartIndex = " << startIndex << " == " << v[x].getStartIndex() << ", blockIndex = " << blockIndex << " == " << v[x].getBlockIndex() << endl;
                //cout << "standalone:\tstartIndex = " << startIndex << " == " << standalone.getStartIndex() << ", blockIndex = " << blockIndex << " == " << standalone.getBlockIndex() << endl;
                //cout << "image:\t\tstartIndex = " << startIndex << " == " << index(x, y, z).getStartIndex() << ", blockIndex = " << blockIndex << " == " << index(x, y, z).getBlockIndex() << endl;
                
                
		        startIndex++;
            }
		}

		/**
		Tests writing and reading of ri image.
		*/
		void readWriteRi()
		{
			Image<internals::RiStorageSet> ri(30, 20, 10);
			//Image<internals::RiStorageSet> ri(10, 1, 1);

            vector<uint16_t> checkCounts;

			// Create ri image
			for (coord_t z = 0; z < ri.depth(); z++)
			{
				for (coord_t y = 0; y < ri.height(); y++)
				{
					for (coord_t x = 0; x < ri.width(); x++)
					{
						uint16_t count = (uint16_t)randc(15);
						
						checkCounts.push_back(count);
						
						vector<internals::RiStorageItem> vals;
						for (int n = 0; n < count; n++)
						{
							internals::RiStorageItem it;
							it.srcX = (int16_t)(10);
							it.srcY = (int16_t)(20);
							vals.push_back(it);
						}
						
						ri(x, y, z).setValues(vals);
					}
				}
			}

			// Write
			string file = "./ri/ri";
			Vec3c writeBlockSize(10, 10, 10);
			//Vec3c writeBlockSize(5, 1, 1);
			//Vec3c writeBlockSize(ri.dimensions());
			int n = 0;
			for (coord_t z = 0; z < ri.depth(); z += writeBlockSize.z)
			{
				for (coord_t y = 0; y < ri.height(); y += writeBlockSize.y)
				{
					for (coord_t x = 0; x < ri.width(); x += writeBlockSize.x)
					{
						Image<internals::RiStorageSet> riBlock(writeBlockSize);
						Vec3c pos(x, y, z);
						crop(ri, riBlock, pos);
						internals::writeRiBlock(riBlock, file, n, pos, ri.dimensions());
						n++;
					}
				}
			}


			// Read
			Image<internals::RiStorageSet> riRead(ri.dimensions());
			Vec3c readBlockSize(5, 7, 8);
			//Vec3c readBlockSize(ri.dimensions());
			for (coord_t z = 0; z < ri.depth(); z += readBlockSize.z)
			{
				for (coord_t y = 0; y < ri.height(); y += readBlockSize.y)
				{
					for (coord_t x = 0; x < ri.width(); x += readBlockSize.x)
					{
						Image<internals::RiStorageSet> riBlock(readBlockSize);
						Vec3c pos(x, y, z);
						internals::readRiBlock(riBlock, file, pos);
						copyValues(riRead, riBlock, pos);
					}
				}
			}

			testAssert(equals(ri, riRead), "ri read from disk does not equal the one saved to disk.");
			
			
			n = 0;
			for (coord_t z = 0; z < ri.depth(); z++)
			{
				for (coord_t y = 0; y < ri.height(); y++)
				{
					for (coord_t x = 0; x < ri.width(); x++)
					{
    					const internals::RiStorageSet& r = riRead(x, y, z);
    					
    					testAssert(r.size() == checkCounts[n], string("Count difference at ") + toString(n) + ". Should be " + toString(checkCounts[n]) + ", is " + toString(r.size()));
    					n++;
    					
    					for(coord_t i = 0; i < r.size(); i++)
    					{
        					testAssert(r[(uint16_t)i].srcX == 10, string("Wrong srcX at ") + toString(Vec3c(x, y, z)) + " at pos " + toString(i) + ": " + toString(r[(uint16_t)i].srcX));
        					testAssert(r[(uint16_t)i].srcY == 20, string("Wrong srcY at ") + toString(Vec3c(x, y, z)) + " at pos " + toString(i) + ": " + toString(r[(uint16_t)i].srcY));
    					}
    					
					}
				}
			}
			
			
			
		}


		void readInParts(Image<internals::RiStorageSet>& ri, const string& indexPrefix)
		{
			Vec3c readBlockSize(ri.dimensions());
			readBlockSize[0] = 6;
			for (coord_t z = 0; z < ri.depth(); z += readBlockSize.z)
			{
				for (coord_t y = 0; y < ri.height(); y += readBlockSize.y)
				{
					for (coord_t x = 0; x < ri.width(); x += readBlockSize.x)
					{
						Image<internals::RiStorageSet> riBlock(readBlockSize);
						Vec3c pos(x, y, z);
						internals::readRiBlock(riBlock, indexPrefix, pos);
						copyValues(ri, riBlock, pos);
					}
				}
			}
		}



		/**
		Tests equality of doesDiscretizedCircle1FitInto2 and doesDiscretizedCircle1FitInto2Cached.
		*/
		void discretizedCircles()
		{
			int32_t MAX = 100 * 100;

			Image<uint8_t> img1(210, 210);
			Image<uint8_t> img2(210, 210);

			//int32_t r2_1 = 390;
			//int32_t r2_2 = 393;
			//internals::draw2(img1, Vec3c(50, 50, 0), r2_1, (uint8_t)128);
			//internals::draw2(img2, Vec3c(50, 50, 0), r2_2, (uint8_t)200);

			//raw::writed(img1, "./discretized_circles/1");
			//raw::writed(img2, "./discretized_circles/2");

			//cout << "Fits = " << dimred::doesDiscretizedCircle1FitInto2(r2_2, r2_1) << endl;

			internals::buildCircleLookup(MAX);
			for (int32_t r2_1 = 0; r2_1 < MAX; r2_1++)
			{
				int32_t r = largestIntWhoseSquareIsLessThan(r2_1);
				for (int32_t r2_2 = r2_1; r2_2 <= (r + 1)*(r + 1); r2_2++)
				{
					setValue(img1, 0);
					setValue(img2, 0);
					draw(img1, Sphere2(Vec3sc(50, 50, 0), r2_1), (uint8_t)255);
					draw(img2, Sphere2(Vec3sc(50, 50, 0), r2_2), (uint8_t)255);
					bool imagesEqual = equals(img1, img2);
					bool fits = internals::doesDiscretizedCircle1FitInto2(r2_2, r2_1);
					bool fitsCached = internals::doesDiscretizedCircle1FitInto2Cached(r2_2, r2_1);
					if (!testAssert(fits == imagesEqual, "discretized circle equality"))
					{
						raw::writed(img1, "./discretized_circles/1");
						raw::writed(img2, "./discretized_circles/2");
						cout << "Fits = " << fits << endl;
					}

					testAssert(fits == fitsCached, "difference between cached and non-cached result");
				}

				showProgress(r2_1, MAX);
			}
		}


		/**
		Tests that dimred version of local thickness calculation gives the same result than simple version.
		*/
		void testThickmapEquality(const Image<uint8_t>& geom)
		{
			Timer timer;

			Image<int32_t> dmap2;
			timer.start();
			distanceTransform2(geom, dmap2);
			timer.stop();
			cout << "Distance map took " << timer.getTime() << " ms." << endl;
			

			// Calculate distance ridge
			Image<int32_t> ridge2;
			timer.start();
			centersOfLocallyMaximalSpheres(dmap2, ridge2);
			timer.stop();
			cout << "Distance ridge took " << timer.getTime() << " ms." << endl;

			raw::writed(ridge2, "./localthickness/ridge2");

			

			// Calculate local thickness using simple sphere plotting algorithm from distance map
			Image<int32_t> thicknessSimpleDmap;
			timer.start();
			itl2::standard::thickmap2(dmap2, thicknessSimpleDmap);
			timer.stop();
			cout << "Simple algorithm using dmap took " << timer.getTime() << " ms." << endl;
			raw::writed(thicknessSimpleDmap, "./localthickness/thickness2_simple_dmap");


			// Calculate local thickness using simple sphere plotting algorithm from distance ridge
			Image<int32_t> thicknessSimpleRidge;
			timer.start();
			itl2::standard::thickmap2(ridge2, thicknessSimpleRidge);
			timer.stop();
			cout << "Simple algorithm using ridge took " << timer.getTime() << " ms." << endl;
			raw::writed(thicknessSimpleRidge, "./localthickness/thickness2_simple_ridge");



			// Calculate local thickness using optimized sphere plotting algorithm from distance map
			Image<int32_t> thicknessOptDmap;
			timer.start();
			setValue(thicknessOptDmap, dmap2);
			itl2::optimized::thickmap2(thicknessOptDmap);
			timer.stop();
			cout << "Optimized algorithm using dmap took " << timer.getTime() << " ms." << endl;
			raw::writed(thicknessOptDmap, "./localthickness/thickness2_optimized_dmap");


			// Calculate local thickness using optimized sphere plotting algorithm from distance ridge
			Image<int32_t> thicknessOptRidge;
			timer.start();
			setValue(thicknessOptRidge, ridge2);
			itl2::optimized::thickmap2(thicknessOptRidge);
			timer.stop();
			cout << "Optimized algorithm using ridge took " << timer.getTime() << " ms." << endl;
			raw::writed(thicknessOptRidge, "./localthickness/thickness2_optimized_ridge");



			// Calculate local thickness using dimensionality reduction algorithm from distance map
			Image<int32_t> thicknessDimRedSuperDmap;
			timer.start();
			itl2::dimredsuper::thickmap2(dmap2, thicknessDimRedSuperDmap);
			timer.stop();
			cout << "DimRedSuper algorithm using dmap took " << timer.getTime() << " ms." << endl;
			raw::writed(thicknessDimRedSuperDmap, "./localthickness/thickness2_dimredsuper_dmap");

			// Calculate local thickness using dimensionality reduction algorithm from distance ridge
			Image<int32_t> thicknessDimRedSuperRidge;
			timer.start();
			itl2::dimredsuper::thickmap2(ridge2, thicknessDimRedSuperRidge);
			timer.stop();
			cout << "DimRedSuper algorithm using ridge took " << timer.getTime() << " ms." << endl;
			raw::writed(thicknessDimRedSuperRidge, "./localthickness/thickness2_dimredsuper_ridge");



			// Check difference between results of various algorithms
			checkDifference(thicknessSimpleDmap, thicknessSimpleRidge, "distance ridge does not contain all required points.");
			checkDifference(thicknessSimpleDmap, thicknessOptDmap, "simple local thickness algorithm and optimized algorithm (dmap) do not give the same results.");
			checkDifference(thicknessSimpleDmap, thicknessOptRidge, "simple local thickness algorithm and optimized algorithm (dmap) do not give the same results.");
			checkDifference(thicknessSimpleDmap, thicknessDimRedSuperDmap, "simple local thickness algorithm and dimensionality reduction super algorithm (dmap) do not give the same results.");
			checkDifference(thicknessSimpleDmap, thicknessDimRedSuperRidge, "simple local thickness algorithm and dimensionality reduction super algorithm (ridge) do not give the same results.");
			

			// DEBUG: output difference automatically
			//subtract(thicknessDimRedBlocksRidge, thicknessSimpleDmap);
			//raw::writed(thicknessDimRedBlocksRidge, "./localthickness/difference");


			// Check that thickness maps have nonzero pixels at the same locations than in the original
			// It's sufficient to test one of the dmaps as we have checked that they are equal.
			threshold(thicknessSimpleDmap, 0);
			checkDifference(thicknessSimpleDmap, geom, "difference between original geometry and nonzero points of thickness map.");
		}



		void thickmapsEquality()
		{
			//throwOnFailedAssertion(true);

			//for (size_t n = 0; n < 100; n+=2)
			size_t n = 0;
			{
				cout << "Test " << n << endl;

				Image<uint8_t> geom(200, 200, 200);
				generateSimpleGeometry(geom, (unsigned int)n);

				raw::writed(geom, "./localthickness/geom");
				testThickmapEquality(geom);


				//cout << "Test " << (n + 1) << endl;

				//linearMap(geom, Vec4d(0, 1, 1, 0));
				//raw::writed(geom, "./localthickness/geom");
				//testThickmapEquality(geom);

			}
		}





		/**
		Calculates block size for block-wise processing such that there will be 8 blocks (in 3D).
		*/
		Vec3c getBlockSize(const Vec3c& originalDimensions)
		{
			coord_t subdivisions = 3;
			// take ceiling using equation: 1 + ((x - 1) / y)
			return Vec3c(1, 1, 1) + (originalDimensions - Vec3c(1, 1, 1)) / subdivisions;
			//return originalDimensions;
		}
		
		struct ThickmapStatistics
		{
			/**
			Algorithm title.
			*/
			string title;

			/**
			Random seed used to generate the image.
			*/
			int seed;

			/**
			Linear image size.
			*/
			int imageSize;

			/**
			Timings.
			*/
			double dmapTime;
			double ridgeTime;
			double thickmapTime;
			double finalizeTime;

			/**
			Volume fraction of foreground pixels.
			*/
			double volumeFraction;

			/**
			Volume fraction of non-zero pixels in distance ridge.
			*/
			double ridgeVolumeFraction;

			/**
			Mean r, r^2 and r^3 in the ridge, ZERO POINTS NOT ACCOUNTED FOR!
			*/
			double ridgeMeanR;
			double ridgeMeanR2;
			double ridgeMeanR3;

			/**
			Total memory required as fraction of input image size.
			*/
			double memReq;

			/**
			Total temporary disk space required as a fraction of input image size.
			*/
			double diskReq;

			/**
			Gets titles of fields printed by operator <<.
			*/
			static string titles()
			{
				return "algorithm, image size [pix], seed, object volume fraction [1], ridge volume fraction [1], <r>' [pix], <r^2>' [pix^2], <r^3>' [pix^3], dmap time [ms], ridge time [ms], thickmap time [ms], finalize time [ms], mem req [fraction of dmap size], disk req [fraction of dmap size]";
			}

			friend ostream& operator<<(ostream& out, const ThickmapStatistics& s)
			{
				out << s.title << ", "
					<< s.imageSize << ", "
					<< s.seed << ", "

					<< s.volumeFraction << ", "
					<< s.ridgeVolumeFraction << ", "

					<< s.ridgeMeanR << ", "
					<< s.ridgeMeanR2 << ", "
					<< s.ridgeMeanR3 << ", "

					<< s.dmapTime << ", "
					<< s.ridgeTime << ", "
					<< s.thickmapTime << ", "
					<< s.finalizeTime << ", "

					<< s.memReq << ", "
					<< s.diskReq;
				return out;
			}
		};

		/**
		Converts algorithm index to name.
		*/
		string algorithmName(size_t algorithm)
		{
			if (algorithm == 0)
			{
				return "standard (from dmap)";
			}
			else if (algorithm == 1)
			{
				return "standard (from ridge)";
			}
			else if (algorithm == 2)
			{
				return "optimized (from dmap)";
			}
			else if (algorithm == 3)
			{
				return "optimized (from ridge)";
			}
			else if (algorithm == 12)
			{
				return "dimred super (from dmap)";
			}
			else if (algorithm == 13)
			{
				return "dimred super (from ridge)";
			}
			else
			{
				throw ITLException("Invalid thickness map algorithm index.");
			}
		}

		/**
		Replaces spaces by underscore.
		*/
		void replaceSpaces(string& s)
		{
			std::replace(s.begin(), s.end(), ' ', '_');
		}


		/**
		Runs one thickness map test (given distance map) and returns statistics about the processing.
		*/
		ThickmapStatistics testThickmap(const Image<int32_t>& dmap2, size_t algorithm, const string& identifier, bool doRounding = false)
		{
			Timer timer;
			ThickmapStatistics result;

			result.title = algorithmName(algorithm);

			// Ridge or not?
			Image<int32_t> ridge2(dmap2.dimensions());
			if (algorithm == 1 || algorithm == 3 || algorithm == 13)
			{
				timer.start();
				centersOfLocallyMaximalSpheres(dmap2, ridge2);
				timer.stop();
				result.ridgeTime = timer.getTime();
			}
			else
			{
				setValue(ridge2, dmap2);
				result.ridgeTime = 0;
			}

			// Round?
			if (doRounding)
			{
				roundDistanceRidge2(ridge2);
			}


			Image<int32_t> thickness2(ridge2.dimensions());
			setValue(thickness2, ridge2);

			// Calculate statistics of original distance values.
			{
				Image<float32_t> tmp;

				convert(ridge2, tmp);
				squareRoot(tmp);
				result.ridgeMeanR = maskedmean(tmp, 0.0f);

				pow(tmp, 2);
				result.ridgeMeanR2 = maskedmean(tmp, 0.0f);

				convert(ridge2, tmp);
				squareRoot(tmp);
				pow(tmp, 3);
				result.ridgeMeanR3 = maskedmean(tmp, 0.0f);

				convert(ridge2, tmp);
				threshold(tmp, 0);
				double nonZeroCount = sum(tmp);
				result.ridgeVolumeFraction = nonZeroCount / (double)ridge2.pixelCount();


				//double meanDistance = sum(dmap);
				//double maxDistance = max(dmap);
				//threshold(dmap, 0);
				//double nonZeroCount = sum(dmap);
				//meanDistance /= nonZeroCount;
			}

			
			Vec3c blockSize = getBlockSize(ridge2.dimensions());

			// Squared local thickness.
			double inputMem = (double)ridge2.pixelCount() * sizeof(int32_t);
			double extraMem = 0;
			double diskReq = 0;
			timer.start();
			if (algorithm == 0 || algorithm == 1)
			{
				itl2::standard::thickmap2(ridge2, thickness2);
				extraMem = (double)(thickness2.pixelCount() * sizeof(int32_t)); // We need separate output image.
			}
			else if (algorithm == 2 || algorithm == 3)
			{
				itl2::optimized::thickmap2(thickness2, &extraMem);
			}
			else if (algorithm == 12 || algorithm == 13)
			{
				itl2::dimredsuper::thickmap2<int32_t>(ridge2, thickness2, &extraMem);
			}
			else
			{
				throw ITLException("Invalid thickness map algorithm index.");
			}
			
			timer.stop();
			result.thickmapTime = timer.getTime();

			// Calculate memory requirement as a fraction of input size
			// (input block size + extra data size) / (total image size)
			result.memReq = (inputMem + extraMem) / (ridge2.pixelCount() * sizeof(int32_t));
			result.diskReq = diskReq / (ridge2.pixelCount() * sizeof(int32_t));

			// Take square root and multiply by 2
			Image<float32_t> thickness;
			timer.start();
			finalizeThickmap(thickness2, thickness);
			timer.stop();
			result.finalizeTime = timer.getTime();


			

			raw::writed(thickness2, "./localthickness/tmap2_" + identifier);

			return result;
		}

		/**
		Runs one thickness map test, given geometry.
		*/
		ThickmapStatistics testThickmap(const Image<uint8_t>& geom, size_t algorithm, string identifier, bool doRounding = false)
		{
			double objectVolume = sum(geom);
			double volumeFraction = objectVolume / geom.pixelCount();

			Timer timer;

			Image<int32_t> dmap2;
			timer.start();
			distanceTransform2(geom, dmap2);
			timer.stop();
			double dmapTime = timer.getTime();
			raw::writed(dmap2, "./localthickness/dmap2_" + identifier);

			ThickmapStatistics result = testThickmap(dmap2, algorithm, identifier, doRounding);

			result.dmapTime = dmapTime;
			result.imageSize = (int)geom.width();
			result.volumeFraction = volumeFraction;

			return result;
		}

		string makeIdentifier(size_t algorithm, size_t imageSize, bool invert)
		{
			string identifier = algorithmName(algorithm) + "_" + toString(imageSize);
			replaceSpaces(identifier);

			if(invert)
				identifier += "_inv";
			else
				identifier += "_normal";

			return identifier;
		}

		/**
		Generates test geometry by creating a large image and then downscaling it to given size.
		This allows testing with the same geometry but multiple resolutions.
		*/
		void generateDownscaledThickmapGeometry(Image<uint8_t>& target, unsigned int seed, double targetVolumeFraction, size_t maxSize, size_t imageSize)
		{
			// Generate full-scale image first.
			Image<uint8_t> geom;

			// If it has been generated already, just read it.
			string fullFile = "./localthickness/geom_source_" + toString(maxSize) + "_" + toString(targetVolumeFraction) + "_" + toString(seed);
			Vec3c tempDims;
			ImageDataType tempDt;
			string reason;
			if (raw::getInfo(fullFile, tempDims, tempDt, reason))
			{
				cout << "Reading full-scale geometry from " << fullFile << endl;
				raw::read(geom, fullFile);
			}
			else
			{
				cout << "Generating full-scale geometry..." << endl;
				geom.ensureSize(maxSize, maxSize, maxSize);
				generateGeometry(geom, seed, targetVolumeFraction);
				raw::writed(geom, fullFile);
				cout << "Done" << endl;
			}

			target.ensureSize(imageSize, imageSize, imageSize);
			scale(geom, target, false, InterpolationMode::Nearest);
		}

		/**
		Runs one thickness map test, given seed to generate geometry.
		*/
		ThickmapStatistics thickmapSpeed(size_t algorithm, size_t imageSize, size_t maxSize, unsigned int seed, bool invert, double targetVolumeFraction, bool doRounding)
		{
			Image<uint8_t> geom;
			generateDownscaledThickmapGeometry(geom, seed, targetVolumeFraction, maxSize, imageSize);

			int userSeed = (int)seed;

			string identifier = makeIdentifier(algorithm, imageSize, invert);

			if (invert)
			{
				linearMap(geom, Vec4d(0, 1, 1, 0));
				userSeed *= -1;
			}

			ThickmapStatistics result = testThickmap(geom, algorithm, identifier, doRounding);
			result.seed = seed;

			return result;
		}

		/**
		Checks that two algorithms gave the same thickness map and throws exception if they are not the same.
		*/
		void checkResults(size_t algo1, size_t algo2, size_t imageSize, bool invert)
		{
			string id1 = makeIdentifier(algo1, imageSize, invert);
			string id2 = makeIdentifier(algo1, imageSize, invert);

			string file1 = "./localthickness/tmap2_" + id1;
			string file2 = "./localthickness/tmap2_" + id2;

			Image<int32_t> img1, img2;
			raw::read(img1, file1);
			raw::read(img2, file2);

			if (differs(img1, img2))
				throw ITLException("Difference between thickness maps " + file1 + " and " + file2);
		}

		void print(const ThickmapStatistics& s, ofstream& out)
		{
			cout << s << endl;
			out << s << endl;
			out.flush();
		}

		/**
		Tests scaling behaviour of the thickness map implementations with computer-generated images.
		Checks also that all the tested algorithms give the same result.
		*/
		void thickmapScaling()
		{
			vector<double> volumeFractions = { 0.5, 0.4, 0.3, 0.2, 0.1 }; // No need to test > 0.5 as we invert the images anyway.
			vector<size_t> imageSizes = { 50, 100, 200, 300, 400, 600, 800, 1000, 1500 };
			//vector<size_t> imageSizes = { 50, 500 };
			//vector<size_t> imageSizes = { 100 };
			//vector<unsigned int> seeds = {11, 21, 31, 41, 51, 61, 71, 81, 91, 101};
			vector<unsigned int> seeds = {11, 21, 31, 41};

			size_t maxSize = max(imageSizes);

			bool doRounding = true;

			string outname;
			if (doRounding)
				outname = "./localthickness/thickmap_results_with_rounding.txt";
			else
				outname = "./localthickness/thickmap_results.txt";

			ofstream out(outname);

			cout << ThickmapStatistics::titles() << endl;
			out << ThickmapStatistics::titles() << endl;

			for (size_t imageSize : imageSizes)
			{
			
				for (double targetVolumeFraction : volumeFractions)
				{
				
					for (unsigned int seed : seeds)
					{
					
						//print(thickmapSpeed(1, imageSize, seed, false, targetVolumeFraction), out);  // Standard
						//print(thickmapSpeed(1, imageSize, seed, true, targetVolumeFraction), out);   // Standard
						//print(thickmapSpeed(3, imageSize, maxSize, seed, false, targetVolumeFraction, doRounding), out);  // Optimized
						//print(thickmapSpeed(3, imageSize, maxSize, seed, true, targetVolumeFraction, doRounding), out);   // Optimized
						print(thickmapSpeed(13, imageSize, maxSize, seed, false, targetVolumeFraction, doRounding), out); // Dimred super
						print(thickmapSpeed(13, imageSize, maxSize, seed, true, targetVolumeFraction, doRounding), out);  // Dimred super

						//checkResults(5, 1, imageSize, false);
						//checkResults(5, 1, imageSize, true);
						//checkResults(3, 13, imageSize, false);
						//checkResults(3, 13, imageSize, true);
					}
				}
			}

		}


		void dimred(bool do2D)
		{
			Image<uint8_t> geom;
			
			if(do2D)
				geom.ensureSize(100, 100, 1);
			else
				geom.ensureSize(100, 100, 100);

			//draw(geom, Sphere(Vec3f(50, 45, 50), 10.0f), (uint8_t)1);
			//draw(geom, Sphere(Vec3f(58, 45, 50), 8.0f), (uint8_t)1);
			generateGeometry(geom, 234, 0.6);

			raw::writed(geom, "./dimred/geom");


			Image<int32_t> thicknessDimRedInt32;
			Image<uint32_t> thicknessDimRedUInt32;
			Image<int32_t> thicknessSimpleInt32;
			Image<uint32_t> thicknessSimpleUInt32;

			{
				Image<int32_t> dmap2;
				itl2::distanceTransform2(geom, dmap2);

				raw::writed(dmap2, "./dimred/dmap2_int32");

				// Calculate local thickness using dimensionality reduction algorithm
				Vec3d counts;
				itl2::dimredsuper::thickmap2(dmap2, thicknessDimRedInt32, nullptr, &counts);
				raw::writed(thicknessDimRedInt32, "./dimred/thickness2_dimredsuper_int32");

				cout << "rilist sizes: " << counts << endl;
			}


			{
				Image<uint32_t> dmap2;
				itl2::distanceTransform2(geom, dmap2);

				raw::writed(dmap2, "./dimred/dmap2_uint32");

				// Calculate local thickness using dimensionality reduction algorithm
				Vec3d counts;
				itl2::dimredsuper::thickmap2(dmap2, thicknessDimRedUInt32, nullptr, &counts);
				raw::writed(thicknessDimRedUInt32, "./dimred/thickness2_dimredsuper_uint32");

				cout << "rilist sizes: " << counts << endl;
			}

			{
				itl2::distanceTransform2(geom, thicknessSimpleInt32);
				itl2::optimized::thickmap2(thicknessSimpleInt32);
				raw::writed(thicknessSimpleInt32, "./dimred/thickness2_optimized_int32");
			}

			{
				itl2::distanceTransform2(geom, thicknessSimpleUInt32);
				itl2::optimized::thickmap2(thicknessSimpleUInt32);
				raw::writed(thicknessSimpleUInt32, "./dimred/thickness2_optimized_uint32");
			}


			checkDifference(thicknessDimRedInt32, thicknessDimRedUInt32, "int and uint dimred");
			checkDifference(thicknessSimpleInt32, thicknessSimpleUInt32, "int and uint optimized");
			checkDifference(thicknessSimpleInt32, thicknessDimRedInt32, "optimized and dimred");
		}

		void dimred2D()
		{
			dimred(true);
		}

		void dimred3D()
		{
			dimred(false);
		}


		void testDataForFiji()
		{
			// This generates comparison data for Fiji implementation
			Image<uint8_t> geometry;
			raw::read(geometry, "./input_data/t1-head_bin_256x256x129.raw");
			raw::writed(geometry, "./tmap_processing_phases/geometry");

			Image<float32_t> dmap2;
			itl2::distanceTransform2(geometry, dmap2);
			raw::writed(dmap2, "./tmap_processing_phases/dmap2");

			Image<float32_t> ridge2;
			itl2::centersOfLocallyMaximalSpheres(dmap2, ridge2);
			raw::writed(ridge2, "./tmap_processing_phases/ridge2");

			Image<float32_t> r2;
			itl2::dimredsuper::thickmap2(dmap2, r2);
			raw::writed(r2, "./tmap_processing_phases/rmap2");

			Image<float32_t> tmap;
			itl2::finalizeThickmap(r2, tmap);
			raw::writed(tmap, "./tmap_processing_phases/tmap");
		}


		void thickmapRounding()
		{
			Image<uint8_t> geometry;

			//int algo = 3; // optimized
			int algo = 13; // dimred super

			{
				raw::read(geometry, "./input_data/t1-head_bin_256x256x129.raw");
				linearMap(geometry, Vec4d(0, 1, 1, 0));
				
				auto stats = testThickmap(geometry, algo, algorithmName(algo) + "_no_rounding", false);
				cout << "No rounding: " << stats.thickmapTime << " ms" << endl;
			}

			{
				raw::read(geometry, "./input_data/t1-head_bin_256x256x129.raw");
				linearMap(geometry, Vec4d(0, 1, 1, 0));

				auto stats = testThickmap(geometry, algo, algorithmName(algo) + "_with_rounding", true);
				cout << "With rounding: " << stats.thickmapTime << " ms" << endl;
			}

		}
	}

}
