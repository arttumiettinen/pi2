
#include "fastmaxminfilters.h"

#include "filters.h"
#include "conversions.h"
#include "stringutils.h"
#include "pointprocess.h"
#include "projections.h"
#include "io/raw.h"

#include "testutils.h"

#include <map>
using namespace std;


namespace itl2
{

	namespace internals
	{
		/**
		Adds all possible vectors of the form (a, b, c) to list dirs, where each of a, b, and c is in range [-stepR, stepR] excluding 0.
		*/
		void addDirs(vector<Vec3c>& dirs, coord_t stepR)
		{
			for (coord_t dz = stepR; dz >= -stepR; dz--)
			{
				for (coord_t dy = stepR; dy >= -stepR; dy--)
				{
					for (coord_t dx = stepR; dx >= -stepR; dx--)
					{
						Vec3c dir(dx, dy, dz);
						if (dir != Vec3c(0, 0, 0))
						{
							if (std::find(dirs.begin(), dirs.end(), dir) == dirs.end() &&
								std::find(dirs.begin(), dirs.end(), -dir) == dirs.end())
							{
								// The dir is new
								// Check if its unit length version is also new
								bool add = true;

								Vec3d newUnit(dir);
								newUnit.normalize();
								for (Vec3c oldDir : dirs)
								{
									Vec3d oldUnit(oldDir);
									oldUnit.normalize();
									if (newUnit.equals(oldUnit) || newUnit.equals(-oldUnit))
									{
										add = false;
										break;
									}
								}

								if (add)
									dirs.push_back(dir);
							}
						}
					}
				}
			}
		}

		/**
		Helps in sorting direction vectors according to their length.
		*/
		bool dirComparer(const Vec3c& a, const Vec3c& b)
		{
			return NumberUtils<float>::lessThan(a.normSquared(), b.normSquared());
		}

		/**
		Creates standard DecomposedSphere with empty rls array.
		*/
		DecomposedSphere createDirections()
		{
			DecomposedSphere result;

			// Find all possible directions
			addDirs(result.dirs, 1);
			addDirs(result.dirs, 2);

			sort(result.dirs.begin(), result.dirs.end(), dirComparer);

			// Find count of directions with the same length.
			for (size_t n = 1; n < result.dirs.size(); n++)
			{
				if (!NumberUtils<float>::equals(result.dirs[n].norm(), result.dirs[n - 1].norm()))
				{
					result.Ns.push_back(n);
				}
			}
			result.Ns.push_back(result.dirs.size());

			// Create empty rls array.
			result.rls.insert(result.rls.begin(), result.Ns.size(), 0);

			return result;
		}

		/**
		Finds optimal structuring element composed of multiple periodic lines that corresponding to a sphere as well as possible.
		*/
		DecomposedSphere optimizeStructuringElement(coord_t r)
		{
			cout << "Determining approximate structuring element decomposition for sphere of radius " << r << endl;

			// Generate ground truth image
			// TODO: This is actually not required, we may use the analytical expression of sphere instead of image.
			coord_t radj = itl2::round(1.1 * r + 2);
			coord_t size = 2 * radj + 1;
			Vec3c c(radj, radj, radj);
			Image<uint8_t> gt(size, size, size);

			for (coord_t z = 0; z < size; z++)
			{
				for (coord_t y = 0; y < size; y++)
				{

					for (coord_t x = 0; x < size; x++)
					{
						if ((Vec3c(x, y, z) - c).norm() < r)
						{
							gt(x, y, z) = 1;
						}
					}
				}
			}

			//raw::writed(gt, string("./tests/ground_truth_sphere_r") + toString(r));


			DecomposedSphere result = createDirections();

			// Unknown variables are rl (line radius) for each set of dirs.
			//vector<coord_t> rls(Ns.size(), 0);
			Image<uint8_t> img(size, size, size);

			// Strategy:
			// Set rls[0] = 1
			// Start from rls[4], increase until approximate sphere is too large. Set to largest value resulting in not too large sphere.
			// Continue similarly with rls[3], then rls[2], etc.

			result.rls[0] = 1;
			coord_t topInd = 4;
			for (coord_t topInd = 4; topInd >= 0; topInd--)
			{
				cout << "Optimizing element group " << topInd << "\r" << flush;
				result.rls[topInd] = 1;

				double lastError = numeric_limits<double>::max();
				while (true)
				{

					// Generate test image
					setValue(img, 0);
					img(c) = 1;

					maxFilterSphereApprox(img, result, BoundaryCondition::Zero);

					if (topInd >= 1)
					{
						// Multiplicative rls adjustment, test for too big sphere

						bool tooBig = false;
						for (coord_t z = 0; z < size; z++)
						{
							for (coord_t y = 0; y < size; y++)
							{

								for (coord_t x = 0; x < size; x++)
								{
									if ((Vec3c(x, y, z) - c).norm() >= r && img(x, y, z) != 0)
									{
										tooBig = true;
									}
								}
							}
						}

						if (tooBig)
						{
							// Undo last rls change
							if (result.rls[topInd] > 1)
								result.rls[topInd] /= 2;
							else
								result.rls[topInd] = 0;

							break;
						}

						// Increase rls[topInd]
						result.rls[topInd] *= 2;
					}
					else
					{
						// Additive rls adjustment for fine tuning, test for total error in the sphere

						// Calculate absolute difference to original
						double error = absDifference(img, gt);

						// Break when error starts to increase
						if (error > lastError)
						{
							// Undo rls change
							result.rls[topInd]--;
							break;
						}
						lastError = error;

						result.rls[topInd]++;
					}
				}

			}
			cout << "                           \r" << flush;

			return result;

			//// Create image of the optimization result
			//setValue(img, 0);
			//img(c) = 1;
			//maxFilterSphereApprox(img, result);

			//raw::writed(img, string("./tests/optimized_sphere_r") + toString(r));

			//double totalError = absDifference(img, gt);
			//double fracError = totalError / sum(gt);

			//cout << "Best rls: ";
			//for (size_t n = 0; n < result.rls.size(); n++)
			//{
			//	cout << result.rls[n];
			//	if (n < result.rls.size() - 1)
			//		cout << ", ";
			//}
			//cout << " (error = " << totalError << " pixels = " << (fracError * 100) << " %)" << endl;

		}

		/**
		Structuring element decomposition cache.
		Used to store decompositions that have been optimized as the optimization process may be slow for large structuring elements.
		*/
		map<coord_t, DecomposedSphere> sphereCache;

		/**
		Reads structuring element decompositions from file into cache map.
		*/
		void readCache()
		{
			string data = readText("spherical_structuring_element_cache.txt", false);

			vector<string> lines = split(data, false);
			sphereCache.clear();
			for (const string& line : lines)
			{
				vector<string> items = split(line, false, ',');
				coord_t r = itl2::fromString<coord_t>(items[0]);

				DecomposedSphere elem = createDirections();
				if (items.size() == elem.rls.size() + 1)
				{
					for (size_t n = 0; n < items.size() - 1; n++)
						elem.rls[n] = itl2::fromString<coord_t>(items[n + 1]);
				}

				sphereCache[r] = elem;
			}
		}

		/**
		Writes structuring element decompositions from cache map to a file.
		*/
		void writeCache()
		{
			stringstream s;
			for (const auto& elem : sphereCache)
			{
				s << elem.first << ", ";
				for (size_t n = 0; n < elem.second.rls.size(); n++)
				{
					s << elem.second.rls[n];
					if (n < elem.second.rls.size() - 1)
						s << ", ";
				}
				s << endl;
			}

			writeText("spherical_structuring_element_cache.txt", s.str());
		}

		/**
		Finds optimal decomposition of sphere of radius r from cache, or if the cache does not contain it, calculates the decomposition.
		*/
		DecomposedSphere optimizeStructuringElementCached(coord_t r)
		{
			readCache();

			auto iter = sphereCache.find(r);
			if (iter != sphereCache.end())
				return iter->second;

			auto elem = optimizeStructuringElement(r);
			sphereCache[r] = elem;
			writeCache();

			return elem;
		}
	}




	namespace tests
	{

		template<typename pixel_t> void testOneLineMax(Image<pixel_t>& img, coord_t r, BoundaryCondition bc)
		{
			raw::writed(img, "tests/linemax_original");

			Image<pixel_t> imgGT;
			setValue(imgGT, img);


			lineMax(img, r, Vec3c(1, 0, 0), bc);
			raw::writed(img, "tests/linemax");

			maxFilter(imgGT, Vec3c(r, 0, 0), bc);
			raw::writed(imgGT, "tests/linemaxGT");

			checkDifference(img, imgGT, string("linemax r = ") + toString(r) + string(", size = ") + toString(img.dimensions()));
		}

		template<typename pixel_t> void testOneLineMin(Image<pixel_t>& img, coord_t r, BoundaryCondition bc)
		{
			raw::writed(img, "tests/linemin_original");

			Image<pixel_t> imgGT;
			setValue(imgGT, img);


			lineMin(img, r, Vec3c(1, 0, 0), bc);
			raw::writed(img, "tests/linemin");

			minFilter(imgGT, Vec3c(r, 0, 0), bc);
			raw::writed(imgGT, "tests/lineminGT");

			checkDifference(img, imgGT, string("linemin r = ") + toString(r) + string(", size = ") + toString(img.dimensions()));
		}

		void testSetLineMax(BoundaryCondition bc)
		{
			for (coord_t r = 1; r < 10; r++)
				//coord_t r = 2;
			{
				for (coord_t w = 1; w < 20; w++)
					//coord_t w = 7;
				{
					Image<float32_t> img(w, 1, 1);
					for (coord_t x = 0; x < img.width(); x++)
						img(x, 0, 0) = -(float32_t)x;
					testOneLineMax(img, r, bc);

					for (coord_t x = 0; x < img.width(); x++)
						img(x, 0, 0) = -(float32_t)(img.width() - x);
					testOneLineMax(img, r, bc);
				}
			}
		}

		void lineMax()
		{
			cout << "bc = Zero" << endl;
			testSetLineMax(BoundaryCondition::Zero);
			cout << "bc = Nearest" << endl;
			testSetLineMax(BoundaryCondition::Nearest);
		}

		void testSetLineMin(BoundaryCondition bc)
		{
			for (coord_t r = 1; r < 10; r++)
				//coord_t r = 2;
			{
				for (coord_t w = 1; w < 20; w++)
					//coord_t w = 7;
				{
					Image<float32_t> img(w, 1, 1);
					for (coord_t x = 0; x < img.width(); x++)
						img(x, 0, 0) = (float32_t)x;
					testOneLineMin(img, r, bc);

					for (coord_t x = 0; x < img.width(); x++)
						img(x, 0, 0) = (float32_t)(img.width() - x);
					testOneLineMin(img, r, bc);
				}
			}
		}

		void lineMin()
		{
			cout << "bc = Zero" << endl;
			testSetLineMin(BoundaryCondition::Zero);
			cout << "bc = Nearest" << endl;
			testSetLineMin(BoundaryCondition::Nearest);
		}

		void sphereMaxSpeed()
		{
			Image<uint16_t> head;
			raw::read(head, "./input_data/t1-head_256x256x129.raw");

			Image<uint16_t> result(head.dimensions());

			coord_t r = 5;

			Timer timer;

			timer.start();
			maxFilterSphereApprox(head, r, BoundaryCondition::Nearest);
			timer.stop();
			cout << "Decomposed approximate max takes " << timer.getSeconds() << " s" << endl;
			raw::writed(head, "filters/max_sphere_decomp");

			raw::read(head, "./input_data/t1-head_256x256x129.raw");
			timer.start();
			//maxFilter(head, result, r, NeighbourhoodType::Ellipsoidal, BoundaryCondition::Nearest);
			filter<uint16_t, uint16_t, internals::maxOp<uint16_t> >(head, result, Vec3c(r, r, r), NeighbourhoodType::Ellipsoidal, BoundaryCondition::Nearest);
			timer.stop();
			cout << "Normal max takes " << timer.getSeconds() << " s" << endl;
			raw::writed(result, "filters/max_sphere_exact");

		}
	}
}