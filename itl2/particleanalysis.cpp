
#include "particleanalysis.h"
#include "generation.h"
#include "pointprocess.h"

using namespace math;

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

		Results checkThreading(const Image<uint8_t>& img, Connectivity conn, coord_t volumeLimit)
		{

			Image<uint8_t> img1;
			setValue(img1, img);

			Image<uint8_t> img2;
			setValue(img2, img);

			auto analyzers = allAnalyzers(img);
			Results results1, results2;
			
			Timer timer;
			timer.start();
			itl2::analyzeParticles(img1, analyzers, results1, conn, volumeLimit);
			timer.stop();
			cout << "Multithreaded particle analysis took " << timer.getTime() << " ms." << endl;

			timer.start();
			itl2::analyzeParticlesSingleThreaded(img2, analyzers, results2, conn, volumeLimit);
			timer.stop();
			cout << "Singlethreaded particle analysis took " << timer.getTime() << " ms." << endl;

			// Remove x, y, and z from the results as they might not be the same for single- and multithreaded versions.
			// (They are just arbitrary location inside the particle)
			results1.removeColumn("x");
			results1.removeColumn("y");
			results1.removeColumn("z");
			results2.removeColumn("x");
			results2.removeColumn("y");
			results2.removeColumn("z");

			// Sort the results so that threading does not affect the comparisons
			std::sort(results1.begin(), results1.end(), resultsComparer);
			std::sort(results2.begin(), results2.end(), resultsComparer);

			//cout << results1.size() << endl;
			//cout << results2.size() << endl;

			size_t N = math::min(results1.size(), results2.size());
			for (size_t n = 0; n < N; n++)
			{
				if (results1.str(n) != results2.str(n))
				{
					cout << "First difference at " << n << endl;
					cout << results1.headers() << endl;
					cout << results1.str(n) << endl;
					cout << results2.str(n) << endl;
					break;
				}
			}
			
			bool different = results1.str() != results2.str();
			if (different)
			{
				raw::writed(img, "particleanalysis/geometry");
				writeText("./particleanalysis/results_multi.txt", results1.str());
				writeText("./particleanalysis/results_single.txt", results2.str());
				//cout << results1 << endl;
				//cout << "----" << endl;
				//cout << results2 << endl;
			}

			testAssert(!different, "multi- and single-threaded particle analysis");

			return results1;
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

			Box box1 = Box(Vec3c(30, 30, 30), Vec3c(50, 50, 50));
			Box box2 = Box(Vec3c(60, 60, 60), Vec3c(150, 150, 150));
			draw(img, box1, (uint8_t)1);
			draw(img, box2, (uint8_t)1);

			checkThreading(img);
		}

		void analyzeParticlesThreadingBig()
		{
			for (size_t trial = 0; trial < 50; trial++)
			{
				cout << "Trial " << trial << endl;
				cout << "-----------------------------------------" << endl;

				Image<uint8_t> img(600, 600, 600);

				for (coord_t n = 0; n < 1000*8; n++)
				{
					Vec3d pos(frand((double)img.width()), frand((double)img.height()), frand((double)img.depth()));
					double r = frand(1, 10);

					draw(img, Sphere(pos, r), (uint8_t)1);
				}

				//if (trial == 2)
				//{
					checkThreading(img);
				//	return;
				//}
			}
		}

		void analyzeParticlesSanity2()
		{
			{
				Image<uint8_t> img;
				raw::read(img, "complicated_particles_1");

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
				raw::read(img, "complicated_particles_1");

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
				raw::read(img, "complicated_particles_1");

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
				raw::read(img, "complicated_particles_1");

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

			Box box1 = Box(Vec3c(30, 30, 30), Vec3c(50, 50, 50));
			Box box2 = Box(Vec3c(60, 60, 60), Vec3c(150, 150, 150));
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