
#include "floodfill.h"
#include "io/raw.h"
#include "io/itlpng.h"
#include "pointprocess.h"
#include "transform.h"
#include "generation.h"
#include "dmap.h"
#include "iteration.h"

#include "testutils.h"

#include <algorithm>


using namespace std;

namespace itl2
{
	

	namespace internals
	{
		template<typename weight_t> class SimpleSeed
		{
		private:
			Vec3c pos;

			weight_t myWeight;
			size_t birthday;

		public:

			/**
			Constructor
			@param p The point.
			@param w Weight of the point. Used to prioritize points_shared with larger weight before points_shared with smaller weight.
			@param birthday The filling round number. Used to prioritize older points_shared before newer points_shared so that points_shared near seeds are filled first.
			*/
			SimpleSeed(const Vec3c& p, weight_t w, size_t birthday) :
				pos(p),
				myWeight(w),
				birthday(birthday)
			{
			}

			/**
			Gets position.
			*/
			const Vec3c& position() const
			{
				return pos;
			}

			/**
			Compares weights.
			*/
			bool operator < (const SimpleSeed& right) const
			{
				// Large weights first, small birthdays first
				if (!NumberUtils<weight_t>::equals(myWeight, right.myWeight))
					return NumberUtils<weight_t>::lessThan(myWeight, right.myWeight);
				else
					return birthday > right.birthday;
			}

			const weight_t weight() const
			{
				return myWeight;
			}

		};
	}


	/**
	This function is used as a baseline implementation for region grow (speed) testing.
	 
	This function is inspired by OpenCV implementation of watershed, found in file \opencv\modules\imgproc\src\segmentation.cpp.
	This was chosen as it is quite fast according to Kornilov - An Overview of Watershed Algorithm Implementations in Open Source Libraries.
	OpenCV implementation is licensed with the following text:

								   License Agreement
					For Open Source Computer Vision Library

	 Copyright (C) 2000, Intel Corporation, all rights reserved.
	 Copyright (C) 2013, OpenCV Foundation, all rights reserved.
	 Third party copyrights are property of their respective owners.

	 Redistribution and use in source and binary forms, with or without modification,
	 are permitted provided that the following conditions are met:

	   * Redistribution's of source code must retain the above copyright notice,
		 this list of conditions and the following disclaimer.

	   * Redistribution's in binary form must reproduce the above copyright notice,
		 this list of conditions and the following disclaimer in the documentation
		 and/or other materials provided with the distribution.

	   * The name of the copyright holders may not be used to endorse or promote products
		 derived from this software without specific prior written permission.

	 This software is provided by the copyright holders and contributors "as is" and
	 any express or implied warranties, including, but not limited to, the implied
	 warranties of merchantability and fitness for a particular purpose are disclaimed.
	 In no event shall the Intel Corporation or contributors be liable for any direct,
	 indirect, incidental, special, exemplary, or consequential damages
	 (including, but not limited to, procurement of substitute goods or services;
	 loss of use, data, or profits; or business interruption) however caused
	 and on any theory of liability, whether in contract, strict liability,
	 or tort (including negligence or otherwise) arising in any way out of
	 the use of this software, even if advised of the possibility of such damage.


	*/
	template<typename label_t, typename weight_t> void growCVTypeAlgorithm(Image<label_t>& labels, const Image<weight_t>& weights, label_t inQueue = numeric_limits<label_t>::max())
	{
		weights.checkSize(labels);
		weights.mustNotBe(labels);

		priority_queue<internals::SimpleSeed<weight_t> > points;
		//queue<internals::SimpleSeed<weight_t> > pits;
		

		// Add all neighbours of seed points_shared to the priority queue
		for (coord_t z = 0; z < labels.depth(); z++)
		{
			for (coord_t y = 0; y < labels.height(); y++)
			{
				for (coord_t x = 0; x < labels.width(); x++)
				{
					Vec3 p(x, y, z);

					if (labels(p) == inQueue)
						throw ITLException(string("Value ") + toString(inQueue) + " must not be used in the labels image, it is needed as a temporary value.");

					if (labels(p) == 0) // Only background points_shared can be neighbours of seeds
					{
						bool hasSeedNeighbour = false;

						for (size_t n = 0; n < labels.dimensionality(); n++)
						{
							for (int delta = -1; delta <= 1; delta += 2)
							{
								Vec3c np = p;
								np[n] += delta;

								if (np[n] >= 0 && np[n] < labels.dimension(n))
								{
									label_t lbl = labels(np);
									weight_t w = weights(np);
									if (w > 0 && lbl != 0 && lbl != inQueue)
									{
										hasSeedNeighbour = true;
									}
								}
							}
						}

						if (hasSeedNeighbour)
						{
							points.push(internals::SimpleSeed<weight_t>(p, weights(p), 0));
							labels(p) = inQueue;
						}
					}
				}
			}
		}


		long round = 0;
		//size_t maxPriorityQueueDepth = points_shared.size();
		//size_t maxQueueDepth = points_shared.size();
		// Grow from the point p to all directions if they are not filled yet.
		while (!points.empty() /*|| !pits.empty()*/)
		{
			//maxQueueDepth = std::max(maxQueueDepth, pits.size());
			//maxPriorityQueueDepth = std::max(maxPriorityQueueDepth, points_shared.size());

			round++;

			weight_t myw;
			Vec3c p;
			//if (!pits.empty())
			//{
			//	const auto& obj = pits.front();
			//	p = obj.position();
			//	myw = obj.weight();
			//	pits.pop();
			//}
			//else
			//{
				const auto& obj = points.top();
				p = obj.position();
				myw = obj.weight();
				points.pop();
			//}

			// Find label for this pixel from the surrounding pixels.
			// Add adjacent unlabeled pixels to the queue.
			label_t currPixelLabel = 0;
			weight_t currPixelLabelWeight = 0;
			for (size_t n = 0; n < labels.dimensionality(); n++)
			{
				for (int delta = -1; delta <= 1; delta += 2)
				{
					Vec3c np = p;
					np[n] += delta;

					if (np[n] >= 0 && np[n] < labels.dimension(n))
					{
						label_t& lbl = labels(np);
						
						if (lbl == 0)
						{
							// This is a non-labeled neighbour. It goes to the queue.
							weight_t w = weights(np);
							if (w > 0)
							{
								//if(w > myw)
								//	pits.push(internals::SimpleSeed<weight_t>(np, myw, round));
								//else
									points.push(internals::SimpleSeed<weight_t>(np, w, round));

								lbl = inQueue;
							}
						}
						else
						{
							// This is a labeled neighbour, find if the label of this point should be the label of the neighbour.
							if (lbl != inQueue)
							{
								weight_t w = weights(np);
								if (currPixelLabelWeight == 0 || w > currPixelLabelWeight) // If there are multiple possibilities, choose the one with the highest weight.
								{
									currPixelLabel = lbl;
									currPixelLabelWeight = w;
								}
							}
						}
					}
				}
			}

			// Fill current pixel
			labels(p) = currPixelLabel;

			//if (round % 10 == 0)
			//	png::writed(labels, "./grow_comparison/animation/round_" + toString(round), 0);
		}

		//cout << "Max pits queue depth = " << maxQueueDepth << endl;
		//cout << "Max priority queue depth = " << maxPriorityQueueDepth << endl;
	}




	namespace tests
	{
		

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
			itl2::floodfillSingleThreaded(filled1, start, (uint8_t)128, (uint8_t)128, conn, &count1, &filledPoints1, numeric_limits<size_t>::max(), &neighbours1);
			timer.stop();
			cout << "Fast version: " << timer.getTime() << " ms" << endl;

			timer.start();
			coord_t count2 = itl2::slowFloodfill(filled2, start, (uint8_t)128, (uint8_t)128, conn, &filledPoints2, numeric_limits<size_t>::max(), &neighbours2);
			timer.stop();
			cout << "Slow version: " << timer.getTime() << " ms" << endl;
			

			testAssert(count1 == count2, "filled point count");

			sort(filledPoints1.begin(), filledPoints1.end(), vecComparer<int32_t>);
			sort(filledPoints2.begin(), filledPoints2.end(), vecComparer<int32_t>);
			testAssert(filledPoints1 == filledPoints2, "filled points");

			// The 'new' version does not return fillColor in neighbouring points_shared list!
			//testAssert(neighbours1 == neighbours2, "neighbouring points_shared");

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
			raw::read(head, "input_data/t1-head_bin_256x256x129.raw");

			Timer t;
			t.start();
			itl2::floodfillSingleThreaded(head, Vec3c(110, 110, 25), (uint8_t)128, (uint8_t)128);
			t.stop();
			cout << "Flood fill (single-threaded) took " << t.getTime() << " ms" << endl;

			raw::writed(head, "./floodfill/filled");
		}

		void floodfillThreading()
		{
			string infile = "input_data/t1-head_bin_256x256x129.raw";
			Vec3c startPoint(110, 110, 25);

			//string infile = "input_data/test_piece_bin_512x512x512.raw";
			//Vec3c startPoint(150, 150, 0);

			//{
			//	Image<uint8_t> orig(750, 750, 750);
			//	raw::writed(orig, "./floodfill/geometry");
			//}
			//string infile = "./floodfill/geometry";
			//Vec3c startPoint(500, 500, 500);
			
			// Multi-threaded
			Image<uint8_t> head;
			raw::read(head, infile);

			size_t filledPointCountMT;
			vector<Vec3sc> filledPointsMT;
			set<uint8_t> nbColorsMT;

			Timer t;
			t.start();
			// We use very small minimum block size to test multithreading
			floodfillBlocks(head, startPoint, (uint8_t)128, (uint8_t)128, Connectivity::AllNeighbours, &filledPointCountMT, &filledPointsMT, 0, &nbColorsMT, true, 20);
			t.stop();
			cout << "Flood fill (multi-threaded) took " << t.getTime() << " ms" << endl;
			raw::writed(head, "./floodfill/filled_multithreaded");


			// Single-threaded
			Image<uint8_t> headTrue;
			raw::read(headTrue, infile);

			size_t filledPointCountST;
			vector<Vec3sc> filledPointsST;
			set<uint8_t> nbColorsST;

			t.start();
			itl2::floodfillSingleThreaded(headTrue, startPoint, (uint8_t)128, (uint8_t)128, Connectivity::AllNeighbours, &filledPointCountST, &filledPointsST, 0, &nbColorsST);
			t.stop();
			cout << "Flood fill (single-threaded) took " << t.getTime() << " ms" << endl;
			raw::writed(headTrue, "./floodfill/filled_singlethreaded");


			checkDifference(head, headTrue, "flood fill vs multithreaded flood fill");
			testAssert(filledPointCountMT == filledPointCountST, "filled point count");
			std::sort(filledPointsMT.begin(), filledPointsMT.end(), vecComparer<int32_t>);
			std::sort(filledPointsST.begin(), filledPointsST.end(), vecComparer<int32_t>);
			testAssert(filledPointsMT == filledPointsST, "filled points");
			testAssert(nbColorsMT == nbColorsST, "neighbouring colors");
		}


		void floodfillLeaks()
		{
			{
				Image<uint8_t> img(3, 2);
				draw(img, AABoxc::fromMinMax(Vec3c(), Vec3c(1, 1, 1)), (uint8_t)255);
				draw(img, AABoxc::fromMinMax(Vec3c(2, 1, 0), img.dimensions()), (uint8_t)255);

				singleTest(img, Vec3c(0, 0, 0), Connectivity::AllNeighbours);
				singleTest(img, Vec3c(0, 0, 0), Connectivity::NearestNeighbours);
			}

			{
				Image<uint8_t> img(3, 2);
				draw(img, AABoxc::fromMinMax(Vec3c(2, 0, 0), Vec3c(3, 1, 1)), (uint8_t)255);
				draw(img, AABoxc::fromMinMax(Vec3c(0, 1, 0), Vec3c(1, 2, 1)), (uint8_t)255);
				

				singleTest(img, Vec3c(0, 0, 0), Connectivity::AllNeighbours);
				singleTest(img, Vec3c(0, 0, 0), Connectivity::NearestNeighbours);
			}

			{
				Image<uint8_t> full;
				raw::read(full, "../test_input_data/complicated_particles_1");

				Image<uint8_t> img(full.width(), full.height(), 10);
				itl2::crop(full, img, Vec3c(0, 0, 10));

				singleTest(img, Vec3c(7, 0, 0), Connectivity::AllNeighbours);
				singleTest(img, Vec3c(7, 0, 0), Connectivity::NearestNeighbours);
			}
			
			{
				Image<uint8_t> full;
				raw::read(full, "../test_input_data/complicated_particles_1");

				Image<uint8_t> img(full, 10, 19); // view of full image

				singleTest(img, Vec3c(7, 0, 0), Connectivity::AllNeighbours);
				singleTest(img, Vec3c(7, 0, 0), Connectivity::NearestNeighbours);
			}
		}

		void growPriority()
		{
			Image<uint16_t> weights;
			raw::read(weights, "../test_input_data/t1-head_256x256x129.raw");

			Image<uint8_t> labels(weights.dimensions());

			labels(110, 90, 63) = 100;
			labels(182, 165, 63) = 200;

			raw::writed(labels, "./grow_priority/labels");

			itl2::grow(labels, weights);

			raw::writed(labels, "./grow_priority/grown");
		}

		void growAll()
		{
			Image<uint8_t> img;
			//raw::read(img, "../test_input_data/t1-head_bin_256x256x129.raw");
			img.ensureSize(1000, 1000, 1000);
			setValue(img, 255);

			img(110, 90, 63) = 100;
			img(182, 165, 63) = 200;

			raw::writed(img, "./grow_all/before_grow");

			itl2::growAll(img, (uint8_t)255, (uint8_t)0);

			raw::writed(img, "./grow_all/after_grow");
		}

		void growComparison()
		{
			//hheap();
			//return;

			
			//coord_t m = 2;
			coord_t m = 4;
			coord_t d = 50;
			
			float32_t mf = (float32_t)m;
			float32_t df = (float32_t)d;
			//Vec3c delta(30, 30, 0);
			Vec3c delta(0, 0, 0);

			Image<uint8_t> geometry(m*100, m*100, m*2*d+1);
			draw(geometry, Sphere<float32_t>(mf*Vec3f(30, 30, df) - Vec3f(delta), mf*11.0f), (uint8_t)255);
			draw(geometry, Sphere<float32_t>(mf*Vec3f(50, 30, df) - Vec3f(delta), mf*11.0f), (uint8_t)255);
			draw(geometry, Sphere<float32_t>(mf*Vec3f(60, 70, df) - Vec3f(delta), mf*30.0f), (uint8_t)255);
			draw(geometry, Sphere<float32_t>(mf*Vec3f(15, 60, df) - Vec3f(delta), mf*20.0f), (uint8_t)255);

			Image<uint8_t> labels(geometry.dimensions());
			labels(m * Vec3c(25, 30, d) - delta) = 80;
			labels(m * Vec3c(52, 30, d) - delta) = 120;
			labels(m * Vec3c(60, 50, d) - delta) = 200;

			labels(m * Vec3c(51, 80, d) - delta) = 250;

			Image<float32_t> weights;
			distanceTransform(geometry, weights);

			raw::writed(geometry, "./grow_comparison/geometry");
			

			/*
			Image<uint16_t> weights;
			raw::read(weights, "C:/mytemp/big_vessel_segmentation_tests/bin4_piece1.2_100x100x100.raw");
			Image<uint8_t> labels;
			raw::read(labels, "C:/mytemp/big_vessel_segmentation_tests/bin4_piece1.2_seeds_100x100x100.raw");
			divide(labels, 2);
			*/
			
			raw::writed(weights, "./grow_comparison/weights");
			raw::writed(labels, "./grow_comparison/labels");

			Timer t;

			Image<uint8_t> labelsMeyer;
			setValue(labelsMeyer, labels);
			t.start();
			growOld(labelsMeyer, weights);
			t.stop();
			cout << "Original version took " << t.getSeconds() << " s." << endl;
			raw::writed(labelsMeyer, "./grow_comparison/watershed_orig");

			Image<uint8_t> labelsOpt;
			setValue(labelsOpt, labels);
			t.start();
			growCVTypeAlgorithm(labelsOpt, weights);
			t.stop();
			cout << "CV type algorithm took " << t.getSeconds() << " s." << endl;
			raw::writed(labelsOpt, "./grow_comparison/watershed_cv");

			Image<uint8_t> labelsHH;
			setValue(labelsHH, labels);
			t.start();
			grow(labelsHH, weights);
			t.stop();
			cout << "CV type algorithm + custom heap took " << t.getSeconds() << " s." << endl;
			raw::writed(labelsHH, "./grow_comparison/watershed_hh");


			checkDifference(labelsMeyer, labelsOpt, "meyer -> CV type");
			checkDifference(labelsMeyer, labelsHH, "meyer -> CV type + custom heap");
			checkDifference(labelsOpt, labelsHH, "CV type -> CV type + custom heap");
		}

		/**
		This is a test method that was used to check performance of various grow function implementations
		and varieties. The data was generated using the generate_particles method from pi2py2_examples,
		and then inverted to create a computationally more challenging problem.
		The results are here:
		Algorithm,						Runtime,	Max RAM usage
		CV-type + bucket map parallel,	34 s,		911 MB
		CV-type + bucket map + pits		39 s,		913 MB
		CV-type + bucket map,			39 s,		913 MB
		CV-type + gheap,				74 s,		1.1 GB
		CV-type + hheap,				59 s,		5.8 GB
		original grow,					311 s,		1.1 GB
		*/
		void growCustomHeap()
		{
			Image<uint16_t> labels;
			Image<float32_t> weights;

			raw::read(labels, "./pi2py2/inv_particles_labeled_500x500x500.raw");
			raw::read(weights, "./pi2py2/inv_particles_priority_500x500x500.raw");

			Timer t;
			t.start();
			grow(labels, weights);
			//growSingleThreaded(labels, weights);
			//growOld(labels, weights);
			t.stop();
			cout << "CV type algorithm + custom heap took " << t.getSeconds() << " s." << endl;
			raw::writed(labels, "./pi2py2/watershed_custom_heap");
		}


	}
}
