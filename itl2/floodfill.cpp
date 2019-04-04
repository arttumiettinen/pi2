
#include "floodfill.h"
#include "io/raw.h"
#include "pointprocess.h"
#include "transform.h"
#include "generation.h"

#include "testutils.h"

#include <algorithm>

using namespace math;
using namespace std;

namespace itl2
{
	namespace tests
	{
		/**
		Flood fill beginning from the given seed points.
		@param origColor Original color that we are filling. (the color of the region where the fill is allowed to proceed)
		@param fillColor Fill color. The filled pixels will be colored with this color.
		@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered.
		@return Count of pixels filled if the fill was terminated naturally; -1 times count of pixels filled if the fill was terminated by reaching fillLimit in filled pixel count or by encountering pixel with stopColor value.
		*/
		template<typename pixel_t> coord_t slowFloodfill(Image<pixel_t>& image, const vector<Vec3sc>& seeds, pixel_t origColor, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, std::vector<math::Vec3sc>* pFilledPoints = 0, size_t fillLimit = numeric_limits<size_t>::max(), std::set<pixel_t>* pNeighbouringColors = 0)
		{
			std::queue<Vec3sc> points;
			for(const Vec3sc& p : seeds)
				points.push(p);

			size_t lastPrinted = 0;
			coord_t filledPoints = 0;
			while (!points.empty())
			{
				const math::Vec3sc p = points.front();

				pixel_t pixelValue = image(p);
				if (pixelValue == origColor)
				{
					filledPoints++;
					image(p) = fillColor;
					if (pFilledPoints)
						pFilledPoints->push_back(p);

					// Add items to queue
					if (connectivity == Connectivity::NearestNeighbours)
					{
						for (size_t n = 0; n < 3; n++)
						{
							if (p[n] > 0)
							{
								Vec3sc np = p;
								np[n]--;
								points.push(np);
							}

							if (p[n] < image.dimension(n) - 1)
							{
								Vec3sc np = p;
								np[n]++;
								points.push(np);
							}
						}
					}
					else
					{
						// All neighbours
						for (int32_t dx = -1; dx <= 1; dx++)
						{
							for (int32_t dy = -1; dy <= 1; dy++)
							{
								for (int32_t dz = -1; dz <= 1; dz++)
								{
									Vec3sc np(p.x + dx, p.y + dy, p.z + dz);
									if (np.x >= 0 && np.y >= 0 && np.z >= 0 &&
										np.x < image.dimension(0) &&
										np.y < image.dimension(1) &&
										np.z < image.dimension(2)
										)
										points.push(np);
								}
							}
						}
					}
				}
				else if (fillColor != stopColor && pixelValue == stopColor)
				{
					// Stop color has been encountered.
					return -filledPoints;
				}
				else
				{
					if (pNeighbouringColors != 0)
						pNeighbouringColors->insert(pixelValue);
				}

				points.pop();

				if ((size_t)filledPoints > fillLimit)
					return -filledPoints; // Fill volume limit has been encountered.

				// Progress report for large fills
				size_t s = points.size();
				if (s > 0 && s % 50000 == 0 && lastPrinted != s)
				{
					lastPrinted = s;
					cout << s << " seed points...\r" << flush;
				}
			}

			return filledPoints;
		}

		template<typename pixel_t> coord_t slowFloodfill(Image<pixel_t>& image, const math::Vec3c& start, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, std::vector<math::Vec3sc>* pFilledPoints = 0, size_t fillLimit = numeric_limits<size_t>::max(), std::set<pixel_t>* pNeighbouringColors = 0)
		{
			pixel_t origColor = image(start);
			if (origColor == fillColor)
				return true;

			vector<Vec3sc> seeds;
			seeds.push_back(Vec3sc(start));

			return slowFloodfill(image, seeds, origColor, fillColor, stopColor, connectivity, pFilledPoints, fillLimit, pNeighbouringColors);
		}

/*
		bool vecComp(const Vec3c&a, const Vec3c& b)
		{
			if (a.z < b.z)
			{
				return true;
			}
			else if (a.z == b.z)
			{
				if (a.y < b.y)
				{
					return true;
				}
				else if (a.y == b.y)
				{
					if (a.x < b.x)
						return true;
				}
			}

			return false;
		}
*/
		void singleTest(Image<uint8_t>& image, const Vec3c& start, Connectivity conn)
		{
			raw::writed(image, "./floodfill/test_image");

			Image<uint8_t> filled1, filled2;
			setValue(filled1, image);
			setValue(filled2, image);

			vector<Vec3sc> filledPoints1, filledPoints2;
			set<uint8_t> neighbours1, neighbours2;

			Timer timer;

			timer.start();
			size_t count1;
			itl2::floodfill(filled1, start, (uint8_t)128, (uint8_t)128, conn, &count1, &filledPoints1, numeric_limits<size_t>::max(), &neighbours1);
			timer.stop();
			cout << "Fast version: " << timer.getTime() << " ms" << endl;

			timer.start();
			coord_t count2 = itl2::tests::slowFloodfill(filled2, start, (uint8_t)128, (uint8_t)128, conn, &filledPoints2, numeric_limits<size_t>::max(), &neighbours2);
			timer.stop();
			cout << "Slow version: " << timer.getTime() << " ms" << endl;
			

			testAssert(count1 == count2, "filled point count");

			sort(filledPoints1.begin(), filledPoints1.end(), math::vecComparer<int32_t>);
			sort(filledPoints2.begin(), filledPoints2.end(), math::vecComparer<int32_t>);
			testAssert(filledPoints1 == filledPoints2, "filled points");

			// The 'new' version does not return fillColor in neighbouring points list!
			//testAssert(neighbours1 == neighbours2, "neighbouring points");

			checkDifference(filled1, filled2, "filled images");

			raw::writed(filled1, "./floodfill/filled");
			raw::writed(filled2, "./floodfill/filled_true");
		}

		void floodfillSanityChecks()
		{
			{
				// Generate grid of alternating 0- and 255-pixels.
				Image<uint8_t> image(200, 200, 200);

				for (coord_t dim = 0; dim < 3; dim++)
				{
					for (coord_t x = 0; x < image.dimension(dim); x += 2)
					{
						Vec3c p(0, 0, 0);
						p[dim] = x;
						image(p) = (uint8_t)x;
					}
				}

				for (coord_t z = 1; z < image.depth(); z++)
				{
					for (coord_t y = 1; y < image.height(); y++)
					{
						for (coord_t x = 1; x < image.width(); x++)
						{
							bool isGap = y >= 8 && y <= 12;
							bool prevXWhite = image(x - 1, y, z) != 0;
							bool prevYWhite = image(x, y - 1, z) != 0;
							bool prevZWhite = image(x, y, z - 1) != 0;

							if (!isGap && !prevXWhite && !prevYWhite && !prevZWhite)
								image(x, y, z) = (uint8_t)x;
						}
					}
				}

				singleTest(image, Vec3c(10, 10, 10), Connectivity::NearestNeighbours);
				singleTest(image, Vec3c(10, 10, 10), Connectivity::AllNeighbours);
			}

			{
				// Two cubes
				Image<uint8_t> image(20, 20, 20);
				for (coord_t z = 0; z < 10; z++)
				{
					for (coord_t y = 0; y < 10; y++)
					{
						for (coord_t x = 0; x < 10; x++)
						{
							image(x, y, z) = 255;
							image(image.width() - 1 - x, image.height() - 1 - y, image.depth() - 1 - z) = 255;
						}
					}
				}

				singleTest(image, Vec3c(14, 14, 14), Connectivity::NearestNeighbours);
				singleTest(image, Vec3c(14 / 2, 14 / 2, 14 / 2), Connectivity::AllNeighbours);
			}

		}

		


		void floodfill()
		{
			// NOTE: No asserts!

			Image<uint8_t> head;
			raw::read(head, "t1-head_bin_256x256x129.raw");

			itl2::floodfill(head, Vec3c(110, 110, 25), (uint8_t)128, (uint8_t)128);

			raw::writed(head, "./floodfill/filled");
		}


		void floodfillLeaks()
		{
			{
				Image<uint8_t> img(3, 2);
				draw(img, Box(Vec3c(), Vec3c(1, 1, 1)), (uint8_t)255);
				draw(img, Box(Vec3c(2, 1, 0), img.dimensions()), (uint8_t)255);

				singleTest(img, Vec3c(0, 0, 0), Connectivity::AllNeighbours);
				singleTest(img, Vec3c(0, 0, 0), Connectivity::NearestNeighbours);
			}

			{
				Image<uint8_t> img(3, 2);
				draw(img, Box(Vec3c(2, 0, 0), Vec3c(3, 1, 1)), (uint8_t)255);
				draw(img, Box(Vec3c(0, 1, 0), Vec3c(1, 2, 1)), (uint8_t)255);
				

				singleTest(img, Vec3c(0, 0, 0), Connectivity::AllNeighbours);
				singleTest(img, Vec3c(0, 0, 0), Connectivity::NearestNeighbours);
			}

			{
				Image<uint8_t> full;
				raw::read(full, "complicated_particles_1");

				Image<uint8_t> img(full.width(), full.height(), 10);
				itl2::crop(full, img, Vec3c(0, 0, 10));

				singleTest(img, Vec3c(7, 0, 0), Connectivity::AllNeighbours);
				singleTest(img, Vec3c(7, 0, 0), Connectivity::NearestNeighbours);
			}
			
			{
				Image<uint8_t> full;
				raw::read(full, "complicated_particles_1");

				Image<uint8_t> img(full, 10, 19); // view of full image

				singleTest(img, Vec3c(7, 0, 0), Connectivity::AllNeighbours);
				singleTest(img, Vec3c(7, 0, 0), Connectivity::NearestNeighbours);
			}
		}
	}
}
