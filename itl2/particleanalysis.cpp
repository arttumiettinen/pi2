
#include "particleanalysis.h"
#include "generation.h"
#include "pointprocess.h"


using namespace std;

namespace itl2
{
	namespace tests
	{

		bool resultsComparer(const vector<double>& v1, const vector<double>& v2)
		{
			for (size_t n = 0; n < v1.size(); n++)
			{
				if (v1[n] != v2[n])
					return v1[n] < v2[n];
			}
			return false;
		}

		void checkResults(Results& resultsSingle, Results& resultsMulti)
		{
			size_t N = std::min(resultsSingle.size(), resultsMulti.size());
			bool different = false;
			for (size_t n = 0; n < N; n++)
			{
				for (size_t m = 0; m < resultsSingle[n].size(); m++)
				{
					if (!NumberUtils<double>::equals(resultsSingle[n][m], resultsMulti[n][m], 0.001))
					{
						cout << "First difference at " << n << endl;
						cout << resultsSingle.headers() << endl;
						cout << resultsSingle.str(n) << endl;
						cout << resultsMulti.str(n) << endl;
						different = true;
						break;
					}
				}

				if (different)
					break;
			}

			if (!different && resultsSingle.size() != resultsMulti.size())
				different = true;

			if (different)
			{
				//raw::writed(img, "particleanalysis/geometry");
				writeText("./particleanalysis/results_single.txt", resultsSingle.str());
				writeText("./particleanalysis/results_multi.txt", resultsMulti.str());
				throw ITLException("Difference found.");

			}

			testAssert(!different, "multi- and single-threaded particle analysis");
		}

		Results checkThreading(const Image<uint8_t>& img, Connectivity conn, coord_t volumeLimit)
		{

			Image<uint8_t> img1;
			setValue(img1, img);

			Image<uint8_t> img2;
			setValue(img2, img);

			Image<uint8_t> img3;
			setValue(img3, img);


			auto analyzers = allAnalyzers(img);
			Results results_multi, results_single, results_alternate;

			Timer timer;


			timer.start();
			{
				vector<Results> allResults;
				vector<vector<vector<Vec3sc> > > allBlockLargeEdgePoints;
				vector<vector<vector<Vec3sc> > > allBlockIncompleteParticles;
				vector<vector<coord_t> > allBlockEdgeZ;

				coord_t blockSize = 100;
				for (coord_t minZ = 0; minZ < img3.depth(); minZ += blockSize)
				{
					coord_t maxZ = minZ + blockSize - 1;
					if (maxZ >= img3.depth())
						maxZ = img3.depth() - 1;

					Image<uint8_t> block(img3, minZ, maxZ);

					Results blockResults;
					vector<vector<Vec3sc> > blockLargeEdgePoints;
					vector<vector<Vec3sc> > blockIncompleteParticles;
					vector<coord_t> blockEdgeZ;
					size_t counter = 0;
					uint8_t fillColor = internals::SpecialColors<uint8_t>::fillColor();
					uint8_t largeColor = internals::SpecialColors<uint8_t>::largeColor();
					internals::analyzeParticlesSingleBlock(block, analyzers, blockResults, &blockIncompleteParticles, &blockLargeEdgePoints, conn, volumeLimit, fillColor, largeColor, counter, img3.depth() * img3.height(), Vec3sc(0, 0, (int32_t)minZ));

					allResults.push_back(blockResults);

					allBlockLargeEdgePoints.push_back(blockLargeEdgePoints);
					allBlockIncompleteParticles.push_back(blockIncompleteParticles);

					blockEdgeZ.push_back(minZ);
					blockEdgeZ.push_back(maxZ);

					allBlockEdgeZ.push_back(blockEdgeZ);
				}


				Results finalResults = allResults[0];
				finalResults.headers() = analyzers.headers();
				vector<vector<Vec3sc> > finalLargeEdgePoints = allBlockLargeEdgePoints[0];
				vector<vector<Vec3sc> > finalIncompleteParticles = allBlockIncompleteParticles[0];
				vector<coord_t> finalEdgeZ = allBlockEdgeZ[0];

				if (allResults.size() > 1)
				{
					for (size_t n = 1; n < allResults.size(); n++)
					{
						finalResults.insert(finalResults.end(), allResults[n].begin(), allResults[n].end());
						finalLargeEdgePoints.insert(finalLargeEdgePoints.end(), allBlockLargeEdgePoints[n].begin(), allBlockLargeEdgePoints[n].end());
						finalIncompleteParticles.insert(finalIncompleteParticles.end(), allBlockIncompleteParticles[n].begin(), allBlockIncompleteParticles[n].end());
						finalEdgeZ.insert(finalEdgeZ.end(), allBlockEdgeZ[n].begin(), allBlockEdgeZ[n].end());

						internals::combineParticleAnalysisResults(analyzers, finalResults, finalLargeEdgePoints, finalIncompleteParticles, volumeLimit, conn, finalEdgeZ, n < allResults.size() - 1);
					}
				}
				else
				{
					internals::combineParticleAnalysisResults(analyzers, finalResults, finalLargeEdgePoints, finalIncompleteParticles, volumeLimit, conn, finalEdgeZ, false);
				}

				results_alternate = finalResults;

			}
			timer.stop();
			cout << "Alternate block combine method test (single-threaded) took " << timer.getTime() << " ms." << endl;

			timer.start();
			itl2::analyzeParticles(img1, analyzers, results_multi, conn, volumeLimit);
			timer.stop();
			cout << "Multithreaded particle analysis took " << timer.getTime() << " ms." << endl;


			timer.start();
			itl2::analyzeParticlesSingleThreaded(img2, analyzers, results_single, conn, volumeLimit);
			timer.stop();
			cout << "Singlethreaded particle analysis took " << timer.getTime() << " ms." << endl;

			// Remove x, y, and z from the results as they might not be the same for single- and multithreaded versions.
			// (They are just arbitrary location inside the particle)
			results_multi.removeColumn("x");
			results_multi.removeColumn("y");
			results_multi.removeColumn("z");
			results_single.removeColumn("x");
			results_single.removeColumn("y");
			results_single.removeColumn("z");
			results_alternate.removeColumn("x");
			results_alternate.removeColumn("y");
			results_alternate.removeColumn("z");

			// Artificial difference to test that code below works.
			//results_single[3][3] = 7;

			// Sort the results so that threading does not affect the comparisons
			std::sort(results_multi.begin(), results_multi.end(), resultsComparer);
			std::sort(results_single.begin(), results_single.end(), resultsComparer);
			std::sort(results_alternate.begin(), results_alternate.end(), resultsComparer);

			//cout << results_single.str() << endl;
			//cout << results_multi.str() << endl;

			checkResults(results_single, results_multi);
			checkResults(results_single, results_alternate);

			return results_multi;
		}

		void checkThreading(const Image<uint8_t>& img)
		{
			cout << "------- AllNeighbours, no volume limit -------" << endl;
			checkThreading(img, Connectivity::AllNeighbours, 0);
			cout << "------- NearestNeighbours, no volume limit -------" << endl;
			checkThreading(img, Connectivity::NearestNeighbours, 0);

			cout << "------- AllNeighbours, with volume limit -------" << endl;
			checkThreading(img, Connectivity::AllNeighbours, 100);
			cout << "------- NearestNeighbours, with volume limit -------" << endl;
			checkThreading(img, Connectivity::NearestNeighbours, 100);
		}

		void analyzeParticlesThreading()
		{
			Image<uint8_t> img(150, 150, 150);

			AABox box1 = AABox(Vec3c(30, 30, 30), Vec3c(50, 50, 50));
			AABox box2 = AABox(Vec3c(60, 60, 60), Vec3c(150, 150, 150));
			draw(img, box1, (uint8_t)1);
			draw(img, box2, (uint8_t)1);

			checkThreading(img);
		}

		void analyzeParticlesThreadingBig()
		{
			for (size_t trial = 0; trial < 50; trial++)
			{
				cout << "Trial " << trial << endl;
				cout << "*****************************************" << endl;

				Image<uint8_t> img(600, 600, 600);

				for (coord_t n = 0; n < 1000*8; n++)
				{
					Vec3d pos(frand((double)img.width()), frand((double)img.height()), frand((double)img.depth()));
					double r = frand(1, 10);

					draw(img, Sphere(pos, r), (uint8_t)1);
				}

				raw::writed(img, "./particleanalysis/geometry");

				checkThreading(img);
			}
		}

		void analyzeParticlesSanity2()
		{
			{
				Image<uint8_t> img;
				raw::read(img, "./input_data/complicated_particles_1");

				// General check
				checkThreading(img);

				// Detailed check
				Results results;
				itl2::analyzeParticles(img, "volume", results, Connectivity::AllNeighbours);

				sort(results.begin(), results.end(), resultsComparer);

				testAssert(results.size() == 3, "particle count");
				testAssert(results.get("volume", 0) == 70, "particle 1");
				testAssert(results.get("volume", 1) == 2685, "particle 2");
				testAssert(results.get("volume", 2) == 3221, "particle 3");
			}

			{
				Image<uint8_t> img;
				raw::read(img, "./input_data/complicated_particles_1");

				// Detailed check
				Results results;
				itl2::analyzeParticles(img, "volume", results, Connectivity::NearestNeighbours);

				sort(results.begin(), results.end(), resultsComparer);

				testAssert(results.size() == 4, "particle count");
				testAssert(results.get("volume", 0) == 70, "particle 1");
				testAssert(results.get("volume", 1) == 212, "particle 2");
				testAssert(results.get("volume", 2) == 2685, "particle 3");
				testAssert(results.get("volume", 3) == 3009, "particle 4");
			}

			// Draw particles to the image
			//for (size_t n = 0; n < results.size(); n++)
			//{
			//	coord_t x = results.get("x", n);
			//	coord_t y = results.get("y", n);
			//	coord_t z = results.get("z", n);

			//	itl2::floodfill(img, Vec3c(x, y, z), (uint8_t)(n + 1), (uint8_t)(n + 1), Connectivity::AllNeighbours);
			//}

			//raw::writed(img, "./particleanalysis/labels");
		}

		void analyzeParticlesVolumeLimit()
		{
			{
				Image<uint8_t> img;
				raw::read(img, "./input_data/complicated_particles_1");

				// General check
				checkThreading(img, Connectivity::AllNeighbours, 3000);

				// Detailed check
				Results results;
				itl2::analyzeParticles(img, "volume", results, Connectivity::AllNeighbours, 3000);

				sort(results.begin(), results.end(), resultsComparer);

				testAssert(results.size() == 2, "particle count");
				testAssert(results.get("volume", 0) == 70, "particle 1");
				testAssert(results.get("volume", 1) == 2685, "particle 2");
				//testAssert(results.get("volume", 2) == 3221, "particle 3"); This particle must be missing as its size is above volume limit.
			}

			{
				Image<uint8_t> img;
				raw::read(img, "./input_data/complicated_particles_1");

				// General check
				checkThreading(img, Connectivity::NearestNeighbours, 3000);

				// Detailed check
				Results results;
				itl2::analyzeParticles(img, "volume", results, Connectivity::NearestNeighbours, 3000);

				// Detailed check
				sort(results.begin(), results.end(), resultsComparer);

				testAssert(results.size() == 3, "particle count");
				testAssert(results.get("volume", 0) == 70, "particle 1");
				testAssert(results.get("volume", 1) == 212, "particle 2");
				testAssert(results.get("volume", 2) == 2685, "particle 3");
				//testAssert(results.get("volume", 3) == 3009, "particle 3"); This particle must be missing as its size is above volume limit.
			}
		}

		void analyzeParticlesSanity()
		{
			Image<uint8_t> img(150, 150, 150);

			AABox box1 = AABox(Vec3c(30, 30, 30), Vec3c(50, 50, 50));
			AABox box2 = AABox(Vec3c(60, 60, 60), Vec3c(150, 150, 150));
			draw(img, box1, (uint8_t)1);
			draw(img, box2, (uint8_t)1);

			auto analyzers = allAnalyzers(img);
			Results results;
			itl2::analyzeParticles(img, analyzers, results);

			cout << results << endl;


			if (testAssert(results.size() == 2, "count of particles found"))
			{

				testAssert(results.get("volume", 0) == box1.volume(), "Volume 1");
				testAssert(results.get("volume", 1) == box2.volume(), "Volume 2");

				testAssert(results.get("is on edge", 0) == 0.0, "edge 1");
				testAssert(results.get("is on edge", 1) == 1.0, "edge 2");

				testAssert(results.get("minx", 0) == box1.minc.x, "minx 1");
				testAssert(results.get("miny", 0) == box1.minc.y, "miny 1");
				testAssert(results.get("minz", 0) == box1.minc.z, "minz 1");

				testAssert(results.get("minx", 1) == box2.minc.x, "minx 2");
				testAssert(results.get("miny", 1) == box2.minc.y, "miny 2");
				testAssert(results.get("minz", 1) == box2.minc.z, "minz 2");

				testAssert(results.get("maxx", 0) == box1.maxc.x - 1, "maxx 1");
				testAssert(results.get("maxy", 0) == box1.maxc.y - 1, "maxy 1");
				testAssert(results.get("maxz", 0) == box1.maxc.z - 1, "maxz 1");

				testAssert(results.get("maxx", 1) == box2.maxc.x - 1, "maxx 2");
				testAssert(results.get("maxy", 1) == box2.maxc.y - 1, "maxy 2");
				testAssert(results.get("maxz", 1) == box2.maxc.z - 1, "maxz 2");
			}
		}

	}
}