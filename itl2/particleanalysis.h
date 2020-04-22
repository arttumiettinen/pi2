#pragma once

#include <iostream>
#include <algorithm>

#include "image.h"
#include "floodfill.h"
#include "utilities.h"

#include "math/vec2.h"
#include "math/vec3.h"

#include "indexforest.h"
#include "aabox.h"
#include "misc.h"

#include "resultstable.h"
#include "analyzers.h"

#include "math/vectoroperations.h"

#include "generation.h"


namespace itl2
{

	///**
	//Number labels such that the numbering begins from one and is consecutive.
	//*/
	//template<typename T> void renumberLabels(Image<T>& image)
	//{
	//	T maxLabel = 0;
	//	map<T, T> oldtonew;	// Maps color of current image to its color in result.
	//	for(size_t n = 0; n < image.pixelCount(); n++)
	//	{
	//		T oldval = **p;
	//		T newval = oldtonew[oldval];
	//		if(newval == 0)
	//		{
	//			// Take unused color (unused in labels and unused in new labels)
	//			maxLabel++;
	//			newval = maxLabel;
	//			oldtonew[oldval] = newval;
	//		}
	//		**p = newval;
	//	}
	//}

	///**
	//Re-colors particles such that particles in the same disjoint set are colored with the same color.
	//@param image The image containing the particles.
	//@param particleRelations Forest describing the relations of the particles.
	//*/
	//template<typename T> void recolorParticles(Image<T>& image, DisjointSetForest<T>& particleRelations, T backgroundColor = 0)
	//{
	//	CursorPointer<T, LinearCursor<T> > p(image, image.createLinearCursor());
	//	for(; p->hasNext(); p->proceed())
	//	{
	//		T oldColor = **p;
	//		T newColor;

	//		if(oldColor != backgroundColor)
	//			newColor = particleRelations.find_set(oldColor);
	//		else
	//			newColor = 0;

	//		if(oldColor != newColor)
	//			**p = newColor;
	//	}
	//}

	/**
	Label all particles with distinct colors beginning from one.
	It is assumed that background pixels are set to zero.
	Uses simple flood fill algorithm.
	If the pixel data type does not support large enough values to label all particles, an exception is thrown.
	@param image Image containing the particles.
	@param particleColor Color of particles. This color will be skipped in labeling, and regions having some other color than this are not labeled. Pass zero to label particles of any color.
	@param firstLabelValue Label value for the first particle encountered.
	@returns The largest label value used in the image.
	*/
	template<typename pixel_t> pixel_t labelParticles(Image<pixel_t>& image, pixel_t particleColor = 0, pixel_t firstLabelValue = 1, Connectivity connectivity = Connectivity::NearestNeighbours, bool showProgress = true)
	{
		size_t counter = 0;

		// Handle "color all particles" case first
		if (particleColor == 0)
		{
			particleColor = std::numeric_limits<pixel_t>::max();

			#pragma omp parallel for if (image.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t n = 0; n < image.pixelCount(); n++)
			{
				if (image(n) != 0)
					image(n) = particleColor;
			}
		}

		pixel_t particleNumber = firstLabelValue;
		if(particleNumber == particleColor)
			particleNumber++;

		for (coord_t z = 0; z < image.depth(); z++)
		{
			for (coord_t y = 0; y < image.height(); y++)
			{
				for (coord_t x = 0; x < image.width(); x++)
				{
					pixel_t pixel = image(x, y, z);

					if (pixel == particleColor)
					{
						if (particleNumber >= std::numeric_limits<pixel_t>::max() - 2)
							throw ITLException("Unable to label because there are not enough distinct pixel values available. Consider converting input image to higher bitdepth (e.g. uint8 to uint16).");

						floodfill(image, Vec3c(x, y, z), particleNumber, particleNumber, connectivity);

						particleNumber++;
						if (particleNumber == particleColor)
							particleNumber++;
					}

					showThreadProgress(counter, image.pixelCount(), showProgress);
				}
			}
		}

		return particleNumber - 1;
	}

    

	namespace internals
	{

		template<typename T> class SpecialColors
		{
		public:
			static T fillColor()
			{
				return std::numeric_limits<T>::max();
			}

			static T largeColor()
			{
				return std::numeric_limits<T>::max() - 1;
			}
		};

		template<> class SpecialColors<uint8_t>
		{
		public:
			static uint8_t fillColor()
			{
				return 128;
			}

			static uint8_t largeColor()
			{
				return 230;
			}
		};

		template<> class SpecialColors<int8_t>
		{
		public:
			static int8_t fillColor()
			{
				return 64;
			}

			static int8_t largeColor()
			{
				return 100;
			}
		};

		template<> class SpecialColors<float32_t>
		{
		public:
			static float32_t fillColor()
			{
				return 2.0f;
			}

			static float32_t largeColor()
			{
				return 3.0f;
			}
		};


		/**
		Tests if Z coordinate of any particle point is 0 or maxZ.
		*/
		inline bool isOnZEdge(const std::vector<Vec3sc>& particle, coord_t maxZ)
		{
			for (size_t n = 0; n < particle.size(); n++)
			{
				if (particle[n].z == 0 || particle[n].z == maxZ)
					return true;
			}

			return false;
		}

		/**
		Perform particle analysis for single image block.
		Has no restrictions on particle count. (uses flood fill algorithm)
		Needs image that can store at least four distinct colours - original particle color, temporary color, big particle color and background color.
		Particles colored with big particle color are skipped. Additionally, particles that are larger than volumeLimit will be colored with big
		particle color. Other particles will be colored with temporary color.
		@param image The image containing the particles.
		@param analyzers List of analyzers to apply to the particles.
		@param results List for result matrix.
		@param largeEdgePoints Pointer to array that will be filled with coordinates of image edge points belonging to large particle.
		@param volumeLimit Only include particles smaller than this value in the results.
		@param fillColor The analyzed particles will be colored with this color.
		@param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
		@param backgroundColor Color of background.
		*/
		template<typename pixel_t> void analyzeParticlesSingleBlock(Image<pixel_t>& image, AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results, std::vector<std::vector<Vec3sc> >* pIncompleteParticles, std::vector<std::vector<Vec3sc> >* pLargeEdgePoints, Connectivity connectivity, size_t volumeLimit, pixel_t fillColor, pixel_t largeColor, size_t& counter, size_t counterMax, const Vec3sc& coordinateShift)
		{
			std::vector<Vec3sc> particlePoints;
			particlePoints.reserve(1000);

			for (coord_t z = 0; z < image.depth(); z++)
			{
				for (coord_t y = 0; y < image.height(); y++)
				{
					for (coord_t x = 0; x < image.width(); x++)
					{
						pixel_t pixel = image(x, y, z);

						if (pixel != fillColor && pixel != 0 && pixel != largeColor)
						{
							bool isOk;
							if (volumeLimit <= 0)
							{
								// No volume limit, use simple fill
								isOk = floodfill(image, Vec3c(x, y, z), fillColor, fillColor, connectivity, 0, &particlePoints);
							}
							else
							{
								// Volume limit (we may have large color region neighbouring the current particle)
								// Capture neighbourColors so that large particles are not collected as 1-pixel regions but larger ones
								// (good for combineParticleAnalysisResults).
								std::set<pixel_t> neighbourColors;
								isOk = floodfill(image, Vec3c(x, y, z), fillColor, fillColor, connectivity, 0, &particlePoints, volumeLimit, &neighbourColors);
								if (isOk) // Check if the region touches region with large color
									isOk = neighbourColors.count(largeColor) <= 0;
							}

							if (isOk)
							{
								// Particle was filled and it was not too large.

								// Test if particle is on image edge. Result is only needed if incomplete particles are to be dealt separately.
								// The test must be done before shifting!
								bool isIncomplete = pIncompleteParticles && isOnZEdge(particlePoints, image.depth() - 1);
									//analyzers::IsOnEdge<Vec3c, pixel_t>::isOnEdge(particlePoints, image.dimensions());

								// Apply shift to particlePoints
								if (coordinateShift != Vec3sc(0, 0, 0))
								{
									for (size_t n = 0; n < particlePoints.size(); n++)
										particlePoints[n] += coordinateShift;
								}

								if (isIncomplete)
								{
									// The particle touches image edge, store it as incomplete point set.
									pIncompleteParticles->push_back(particlePoints);
								}
								else
								{
									// The particle does not touch image edge, store only analysis results.
									std::vector<double> resultLine;
									analyzers.analyze(particlePoints, resultLine);
									results.push_back(resultLine);
								}
							}
							else
							{
								// The fill was ended by size limit or the filled region touches large region
								// --> the filled region is part of large particle.

								if (particlePoints.size() <= 0)
									particlePoints.push_back(Vec3sc((int32_t)x, (int32_t)y, (int32_t)z));

								// Mark the filled points as belonging to a large particle.
								for (size_t n = 0; n < particlePoints.size(); n++)
									image(particlePoints[n]) = largeColor;

								// Store image edge points belonging to large particle
								if (pLargeEdgePoints)
								{
									std::vector<Vec3sc> edgePoints;

									// Find edge points and apply shift
									for (size_t n = 0; n < particlePoints.size(); n++)
									{
										if (particlePoints[n].z <= 0 || particlePoints[n].z >= image.depth() - 1)
											edgePoints.push_back(particlePoints[n] + coordinateShift);
									}

									if(edgePoints.size() > 0)
										pLargeEdgePoints->push_back(edgePoints);

									//for (size_t n = 0; n < particlePoints.size(); n++)
									//{
									//	if (particlePoints[n].z <= 0 || particlePoints[n].z >= image.depth() - 1)
									//		pLargeEdgePoints->push_back(particlePoints[n] + coordinateShift);
									//}
								}
							}

							particlePoints.clear();
						}


					}

					showThreadProgress(counter, counterMax);
				}
			}

		}


		/**
		Perform particle analysis in blocks, one thread per block.
		Has no restrictions on particle count. (uses flood fill algorithm)
		Needs image that can store at least four distinct colours - original particle color, temporary color, big particle color and background color.
		Particles colored with big particle color are skipped. Additionally, particles that are larger than volumeLimit will be colored with big
		particle color. Other particles will be colored with temporary color.
		@param image The image containing the particles.
		@param analyzers List of analyzers to apply to the particles.
		@param results List for result matrix.
		@param pLargeEdgePoints Pointer to array that will be filled with coordinates of image edge points belonging to large particle.
		@param volumeLimit Only include particles smaller than this value in the results.
		@param fillColor The analyzed particles will be colored with this color.
		@param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
		*/
		template<typename pixel_t> void analyzeParticlesBlocks(Image<pixel_t>& image, AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results, std::vector<std::vector<Vec3sc> >& largeEdgePoints, std::vector<std::vector<Vec3sc> >& incompleteParticles, Connectivity connectivity, size_t volumeLimit, pixel_t fillColor, pixel_t largeColor, const Vec3sc& origin, std::vector<coord_t>& blockEdgeZ)
		{
			size_t counter = 0;

			//int count = 4;
			//int idx[] = { 0, 1, 2, 3 };
			//int idxs[] = { 0, 2, 1, 3 };
			//int idxs[] = { 0, 3, 1, 2 };
			//int idxs[] = { 1, 0, 2, 3 };
			//int idxs[] = { 1, 2, 0, 3 };
			//int idxs[] = { 1, 2, 3, 0 };
			//for(int i = 0; i < 4; i++)
			//{
			//	int idx = idxs[i];
			//int count = 4;
			//for(int idx = 0; idx < count; idx++)
			//{
			#pragma omp parallel
			{
				int idx = omp_get_thread_num();
				int count = omp_get_num_threads();

				if (idx == 0)
				{
					#pragma omp critical
					{
						std::cout << "Analyzing in " << count << " blocks." << std::endl;
					}
				}

				// Calculate amount of slices single thread should process.
				// The value is capped so that we don't divide the work unnecessarily too much if the image is small.
				coord_t size = image.depth() / count;
				if (size < 100)
					size = 100;

				// If there's nothing to do for all the threads, the excess threads will just skip processing.
				coord_t minZ = idx * size;
				if (minZ < 0)
					minZ = 0;
				if (minZ < image.depth())
				{

					coord_t maxZ = minZ + size - 1;
					if (maxZ >= image.depth())
						maxZ = image.depth() - 1;

					if (idx == count - 1)
					{
						// The last thread processes possible "rounding error" slices
						maxZ = image.depth() - 1;
					}

//#pragma omp critical
//					{
//						std::cout << "Thread " << idx << " analyzes z-range [" << minZ << ", " << maxZ << "[" << std::endl;
//					}

					// Get view of part of the image
					Image<pixel_t> block(image, minZ, maxZ);

					// Analyze the block. All coordinates are shifted in analyzeParticlesSingleBlock function to global original image coordinates.
					Results blockResults;
					std::vector<std::vector<Vec3sc> > blockLargeEdgePoints;
					std::vector<std::vector<Vec3sc> > blockIncompleteParticles;
					internals::analyzeParticlesSingleBlock(block, analyzers, blockResults, &blockIncompleteParticles, &blockLargeEdgePoints, connectivity, volumeLimit, fillColor, largeColor, counter, image.depth() * image.height(), origin + Vec3sc(0, 0, (int32_t)minZ));

					//raw::writed(block, "./particleanalysis/block");
					//raw::writed(image, "./particleanalysis/image");

					#pragma omp critical(analyzeParticles)
					{
						results.insert(results.end(), blockResults.begin(), blockResults.end());
						largeEdgePoints.insert(largeEdgePoints.end(), blockLargeEdgePoints.begin(), blockLargeEdgePoints.end());
						incompleteParticles.insert(incompleteParticles.end(), blockIncompleteParticles.begin(), blockIncompleteParticles.end());
						blockEdgeZ.push_back(minZ + origin.z);
						blockEdgeZ.push_back(maxZ + origin.z);
					}
				}
			}
		}

		/**
		Tests if the two point lists have neighbouring elements.
		*/
		inline bool isNeighbouring(const std::vector<Vec3sc>& points1, const std::vector<Vec3sc>& points2, Connectivity connectivity)
		{
			switch (connectivity)
			{
			case Connectivity::NearestNeighbours:
				for (size_t n = 0; n < points1.size(); n++)
				{
					for (size_t m = 0; m < points2.size(); m++)
					{
						Vec3sc ad = (points1[n] - points2[m]).abs();
						if (ad.max() <= 1 && ad.sum() <= 1) // Max difference in any coordinate direction is 1 and only one coordinate direction may have absolute value 1 -> sum of coordinates is 1 or 0.
							return true;
					}
				}
				break;
			case Connectivity::AllNeighbours:
				for (size_t n = 0; n < points1.size(); n++)
				{
					for (size_t m = 0; m < points2.size(); m++)
					{
						if ((points1[n] - points2[m]).abs().max() <= 1)
							return true;
					}
				}
				break;
			default:
				throw ITLException("Unsupported connectivity.");
			}
			
			return false;
		}

		/**
		Finds neighbouring incomplete particle regions and combines their indices in the forest.
		At output the forest will contain all indices of incomplete particles list, and indices of neighbouring particles belong to the same root.
		@param edgePoints Points of incomplete particles that are at calculation block edges.
		*/
		inline void findNeighboringParticles(const std::vector<std::vector<Vec3sc> >& edgePoints, IndexForest& forest, Connectivity connectivity)
		{
			// Add all indices to the forest
			std::cout << "Build forest..." << std::endl;
			forest.initialize(edgePoints.size());

			// Calculate bounding boxes for incomplete vertex sets
			std::cout << "Determine bounding boxes..." << std::endl;
			std::vector<AABox<int32_t> > boundingBoxes;
			boundingBoxes.reserve(edgePoints.size());
			Vec3sc m = Vec3sc(std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max());
			Vec3sc M = Vec3sc(std::numeric_limits<int32_t>::lowest(), std::numeric_limits<int32_t>::lowest(), std::numeric_limits<int32_t>::lowest());
			for (size_t n = 0; n < edgePoints.size(); n++)
			{
				AABox<int32_t> b = AABox<int32_t>::boundingBox(edgePoints[n]);
				if (edgePoints[n].size() > 0)
				{
					b.inflate(1);

					m = min(m, b.minc);
					M = max(M, b.maxc);
				}
				boundingBoxes.push_back(b);
			}

			Image<std::vector<size_t> > grid(100, 100, 100);

			std::cout << "Divide boxes to a grid..." << std::endl;
			size_t counter = 0;
			for (size_t n = 0; n < boundingBoxes.size(); n++)
			{
				if (edgePoints[n].size() > 0)
				{
					Vec3c mcell(itl2::floor((double)(boundingBoxes[n].minc.x - m.x) / (double)(M.x - m.x) * (grid.width() - 1)),
								itl2::floor((double)(boundingBoxes[n].minc.y - m.y) / (double)(M.y - m.y) * (grid.height() - 1)),
								itl2::floor((double)(boundingBoxes[n].minc.z - m.z) / (double)(M.z - m.z) * (grid.depth() - 1)));

					Vec3c Mcell(itl2::ceil((double)(boundingBoxes[n].maxc.x - m.x) / (double)(M.x - m.x) * (grid.width() - 1)),
								itl2::ceil((double)(boundingBoxes[n].maxc.y - m.y) / (double)(M.y - m.y) * (grid.height() - 1)),
								itl2::ceil((double)(boundingBoxes[n].maxc.z - m.z) / (double)(M.z - m.z) * (grid.depth() - 1)));

					for (coord_t zc = mcell.z; zc <= Mcell.z; zc++)
					{
						for (coord_t yc = mcell.y; yc <= Mcell.y; yc++)
						{
							for (coord_t xc = mcell.x; xc <= Mcell.x; xc++)
							{
								grid(xc, yc, zc).push_back(n);
							}
						}
					}
				}

				showThreadProgress(counter, boundingBoxes.size());
			}

			std::cout << "Determine overlaps..." << std::endl;
			counter = 0;
			
			#pragma omp parallel if(!omp_in_parallel())
			{
				std::set<std::tuple<size_t, size_t> > nonOverlapping;

				#pragma omp for schedule(dynamic)
				for (coord_t i = 0; i < grid.pixelCount(); i++)
				{
					auto& list = grid(i);

					for (coord_t ni = 0; ni < (coord_t)list.size(); ni++)
					{
						size_t n = list[ni];
						const std::vector<Vec3sc>& v1 = edgePoints[n];
						const AABox<int32_t>& b1 = boundingBoxes[n];

						if (v1.size() > 0)
						{
							for (size_t mi = ni + 1; mi < list.size(); mi++)
							{
								size_t m = list[mi];
								const std::vector<Vec3sc>& v2 = edgePoints[m];
								const AABox<int32_t>& b2 = boundingBoxes[m];

								//if (std::min(n, m) == 308 && std::max(n, m) == 417)
								//{
								//	std::cout << "testing " << n << ", " << m << std::endl;
								//}

								if (v2.size() > 0)
								{
									// Only test if the particles have not been connected yet.
									// There is a race condition between union_sets and find_sets on this line, but
									// that may only result in the test below evaluating true even though it should be
									// false, and in that case we just do some extra work.
									if (forest.find_set(n) != forest.find_set(m))
									{

										if (b1.overlaps(b2)) // Do bounding boxes overlap?
										{
											auto key = std::make_tuple(std::min(n, m), std::max(n, m));

											// Only test points if they have not been determined to be non-overlapping.
											if(nonOverlapping.find(key) == nonOverlapping.end())
											{
												if (isNeighbouring(v1, v2, connectivity))
												{
													// v1 and v2 are the neighbouring, so the regions must be combined
													#pragma omp critical(forestInsert)
													{
														forest.union_sets(n, m);
													}
												}
												else
												{
													nonOverlapping.emplace(key);
												}
											}
										}
									}
								}
							}
						}
					}

					showThreadProgress(counter, grid.pixelCount());
				}
			}
		}

		inline void concatAndShrink(std::vector<Vec3sc>& target, std::vector<Vec3sc>& source)
		{
			target.insert(target.end(), source.begin(), source.end());
			source.clear();
			source.shrink_to_fit();
		}

		/**
		Combines incomplete particles and analyzes them.
		Can be called multiple times, adding more particles to the results, largeEdgePoints, and incompleteParticles arrays.
		The added particles must not be from already processed blocks.
		@param isSubBlock Set to false if this is the final call to this function. This ensures that also particles that touch z-block edges are included in the results table.
		*/
		template<typename pixel_t> void combineParticleAnalysisResults(const AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results, std::vector<std::vector<Vec3sc> >& largeEdgePoints, std::vector<std::vector<Vec3sc> >& incompleteParticles, size_t volumeLimit, Connectivity connectivity, std::vector<coord_t>& blockEdgeZ, bool isSubBlock)
		{
			/*
			At end of this function:
			largeEdgePoints = points of large particles that touch blockEdgeZ.front() or blockEdgeZ.back() (if blockEdgeZ.size() > 0)
			blockEdgeZ = (blockEdgeZ.front(), blockEdgeZ.back())
			incompleteParticles = those particles that touch blockEdgeZ.front() or blockEdgeZ.back()
			results = all other particles (combined ones)

			==> we don't need to load all data in distributed processing, it is enough to load the data sequentially and call this method for combination of previously
			combined points + newly read points.
			*/


			coord_t maxZ = max(blockEdgeZ);

			std::cout << "Processing " << incompleteParticles.size() << " particles on block edge boundaries..." << std::endl;

			// Create list of edge points for each incomplete particle
			std::cout << "Find edge points of particles..." << std::endl;
			std::vector<std::vector<Vec3sc> > edgePoints;
			size_t counter = 0;
			for (coord_t n = 0; n < (coord_t)incompleteParticles.size(); n++)
			{
				auto& points = incompleteParticles[n];
				std::vector<Vec3sc> currEdgePoints;
				for (size_t m = 0; m < points.size(); m++)
				{
					if (find(blockEdgeZ.begin(), blockEdgeZ.end(), points[m].z) != blockEdgeZ.end())
						currEdgePoints.push_back(points[m]);
				}

				// Order of insertions must be the same than order of particles in incompleteParticles list.
				edgePoints.push_back(currEdgePoints);

				showThreadProgress(counter, incompleteParticles.size());
			}

			// Add big particles to edge points
			// NOTE: After this the items {edgePoints[0], ..., edgePoints[incompleteParticles.size()-1]} correspond to the incomplete particles,
			// and the itemes {edgePoints[incompleteParticles.size()], ..., edgePoints[edgePoints.size()-1]} correspond to large particles.
			edgePoints.insert(edgePoints.end(), largeEdgePoints.begin(), largeEdgePoints.end());

			std::vector<bool> isBig(incompleteParticles.size(), false);

			// Find overlapping particles (based on edge points) and combine them
			{
				IndexForest forest;
				findNeighboringParticles(edgePoints, forest, connectivity);

				// Move all points to roots
				std::cout << "Combine particles..." << std::endl;
				size_t counter = 0;
				for (coord_t n = 0; n < (coord_t)incompleteParticles.size(); n++)
				{
					size_t base = forest.find_set(n);
					if (base != n)
					{
						if (base < incompleteParticles.size())
						{
							concatAndShrink(incompleteParticles[base], incompleteParticles[n]);
							concatAndShrink(edgePoints[base], edgePoints[n]);
						}
						else
						{
							// We are going to combine particle n to a big particle.
							// Mark particle n as big and clear its points (except those that touch image z edge).
							isBig[n] = true;
							

							edgePoints[n].clear();
							edgePoints[n].shrink_to_fit();

							incompleteParticles[n].erase(
								std::remove_if(incompleteParticles[n].begin(), incompleteParticles[n].end(),
									[=](Vec3sc& p)
									{
										return p.z != 0 && p.z != maxZ;
									}),
								incompleteParticles[n].end());
						}
					}

					showThreadProgress(counter, incompleteParticles.size());
				}

				if (volumeLimit > 0)
				{
					// Find which particles are big (they have been combined so they may be big themselves or touch
					// big particle not in incompleteParticles list)
					for (coord_t n = 0; n < (coord_t)incompleteParticles.size(); n++)
					{
						if (incompleteParticles[n].size() >= volumeLimit)
							isBig[n] = true;
					}

					// Particles that have been combined with any of the big particles are big.
					for (coord_t n = (coord_t)incompleteParticles.size(); n < (coord_t)edgePoints.size(); n++)
					{
						size_t base = forest.find_set(n);
						if (base < incompleteParticles.size())
							isBig[base] = true;
					}
				}
			}


			// Analyze combined particles
			std::cout << "Analyze combined particles..." << std::endl;
			counter = 0;

			std::vector<std::vector<Vec3sc> > remainingLargeEdgePoints;
			std::vector<std::vector<Vec3sc> > remainingIncompleteParticles;
			std::vector<coord_t> remainingBlockEdgeZ;

			remainingBlockEdgeZ.push_back(0);
			remainingBlockEdgeZ.push_back(maxZ);

			#pragma omp parallel for if(!omp_in_parallel() && incompleteParticles.size() > 1000) schedule(dynamic)
			for (coord_t n = 0; n < (coord_t)incompleteParticles.size(); n++)
			{

				auto& points = incompleteParticles[n];
				// NOTE: Empty point list is used to indicate that the particle has been combined with some other particle,
				// and thus must not be analyzed. (the other particle to which this one was combined contains all the points)
				if(points.size() > 0)
				{

					if(volumeLimit <= 0 || (volumeLimit > 0 && !isBig[n]))
					{
						if (isSubBlock && isOnZEdge(points, maxZ))
						{
							// The particle is incomplete.
							#pragma omp critical(remainingIncompleteParticlesInsert)
							{
								remainingIncompleteParticles.push_back(points);
							}
						}
						else
						{
							// The particle is complete and not big. Analyze it.
							// Note: the parts of particles originating from different calculation blocks
							// do not overlap so there are no duplicate points in the point list
							//Vec3c point = points[0];
							std::vector<double> resultLine;
							analyzers.analyze(points, resultLine);

							#pragma omp critical(resultsInsert)
							{
								results.push_back(resultLine);
							}
						}
					}
					else if (volumeLimit > 0 && isBig[n])
					{
						// The particle is big. Add its edge points to the new largeEdgePoints list.
						#pragma omp critical(remainingLargeEdgePointsInsert)
						{
							remainingLargeEdgePoints.push_back(points);
						}
					}

				}
				showThreadProgress(counter, incompleteParticles.size());
			}

			if (isSubBlock)
			{
				// Add large edge points from last call to the new large edge points list
				// This is necessary as otherwise edge points that are not combined with any particle will
				// be discarded.
				for (auto& points : largeEdgePoints)
				{
					std::vector<Vec3sc> filtered;

					for (Vec3sc& v : points)
					{
						if (v.z == 0 || v.z == maxZ)
							filtered.push_back(v);
					}

					if (filtered.size() > 0)
						remainingLargeEdgePoints.push_back(filtered);
				}

				largeEdgePoints = remainingLargeEdgePoints;
				incompleteParticles = remainingIncompleteParticles;
				blockEdgeZ = remainingBlockEdgeZ;
			}
			else
			{
				incompleteParticles.clear();
			}
		}
	}

	/**
	Makes sure that all particles have the same color.
	*/
	template<typename pixel_t> void prepareParticleAnalysis(Image<pixel_t>& image, pixel_t fillColor = internals::SpecialColors<pixel_t>::fillColor(), pixel_t largeColor = internals::SpecialColors<pixel_t>::largeColor())
	{
		// Find particle color
		pixel_t particleColor = 1;
		while (particleColor == fillColor || particleColor == largeColor)
			particleColor++;

		// Make sure that there are no extra colors in the image.
		#pragma omp parallel for if (image.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t n = 0; n < image.pixelCount(); n++)
		{
			pixel_t& pixel = image(n);
			if (pixel != fillColor && pixel != 0 && pixel != largeColor)
				pixel = particleColor;
		}
	}

	/**
	Measure properties of binary particles separated by background regions.
	Needs image that can store at least four distinct colours - original particle color, temporary color, large particle color and background color,
	and therefore is able to process many particles without need for labeling of all of them with individual values.
	Regions colored with large particle color are skipped.
	Additionally, particles that are larger than volumeLimit will be colored with large particle color.
	Other particles will be colored with fill color.
	@param image The image containing the particles.
	@param analyzers Analyzer objects to apply to the particles.
	@param results Results table.
	@param connectivity Connectivity of the particles.
	@param volumeLimit Only include particles smaller than this value in the results.
	@param fillColor The analyzed particles will be colored with this color.
	@param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
	*/
	template<typename pixel_t> void analyzeParticles(Image<pixel_t>& image, AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results, Connectivity connectivity = Connectivity::AllNeighbours, size_t volumeLimit = 0, pixel_t fillColor = internals::SpecialColors<pixel_t>::fillColor(), pixel_t largeColor = internals::SpecialColors<pixel_t>::largeColor())
	{
		if (fillColor == largeColor)
			throw ITLException("Fill color and large color must not be equal.");
		if (largeColor == 0)
			throw ITLException("Large color must not be zero as zero is background color.");

		std::vector<std::vector<Vec3sc> > incompleteParticles;
		std::vector<std::vector<Vec3sc> > largeEdgePoints;
		std::vector<coord_t> edgeZ;
		
		results.headers() = analyzers.headers();

		prepareParticleAnalysis(image, fillColor, largeColor);

		internals::analyzeParticlesBlocks(image, analyzers, results, largeEdgePoints, incompleteParticles, connectivity, volumeLimit, fillColor, largeColor, Vec3sc(), edgeZ);

		//std::cout << "Block edges are at" << std::endl;
		//sort(edgeZ.begin(), edgeZ.end());
		//for (coord_t z : edgeZ)
		//	std::cout << z << std::endl;

		//Image<uint16_t> labels(image.dimensions());
		//for (size_t n = 0; n < incompleteParticles.size(); n++)
		//{
		//	for (size_t m = 0; m < incompleteParticles[n].size(); m++)
		//	{
		//		labels(incompleteParticles[n][m]) = (uint16_t)(n);
		//	}
		//}
		//raw::writed(labels, "particleanalysis/labels");

		internals::combineParticleAnalysisResults(analyzers, results, largeEdgePoints, incompleteParticles, volumeLimit, connectivity, edgeZ, false);
	}

	/**
	Perform particle analysis.
	Needs image that can store at least four distinct colours - original particle color, temporary color, large particle color and background color.
	Regions colored with large particle color are skipped. Additionally, particles that are larger than volumeLimit will be colored with large
	particle color. Other particles will be colored with fill color.
	@param image The image containing the particles.
	@param analyzerNames List of names of analyzers to apply to the particles.
	@param results Results table.
	@param connectivity Connectivity of the particles.
	@param volumeLimit Only include particles smaller than this value in the results.
	@param fillColor The analyzed particles will be colored with this color.
	@param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
	*/
	template<typename pixel_t> void analyzeParticles(Image<pixel_t>& image, const string& analyzerNames, Results& results, Connectivity connectivity = Connectivity::AllNeighbours, size_t volumeLimit = 0, pixel_t fillColor = internals::SpecialColors<pixel_t>::fillColor(), pixel_t largeColor = internals::SpecialColors<pixel_t>::largeColor())
	{
		auto analyzers = createAnalyzers<pixel_t>(analyzerNames, image.dimensions());
		analyzeParticles(image, analyzers, results, connectivity, volumeLimit, fillColor, largeColor);
	}

	/**
	Perform particle analysis using single-threaded method.
	This one is faster for very small images or images containing only a very few particles.
	Needs image that can store at least four distinct colours - original particle color, temporary color, large particle color and background color.
	Regions colored with large particle color are skipped. Additionally, particles that are larger than volumeLimit will be colored with large
	particle color. Other particles will be colored with fill color.
	@param image The image containing the particles.
	@param analyzers Analyzer objects to apply to the particles.
	@param results Results table.
	@param connectivity Connectivity of the particles.
	@param volumeLimit Only include particles smaller than this value in the results.
	@param fillColor The analyzed particles will be colored with this color.
	@param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
	*/
	template<typename pixel_t> void analyzeParticlesSingleThreaded(Image<pixel_t>& image, AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results, Connectivity connectivity = Connectivity::AllNeighbours, size_t volumeLimit = 0, pixel_t fillColor = internals::SpecialColors<pixel_t>::fillColor(), pixel_t largeColor = internals::SpecialColors<pixel_t>::largeColor())
	{
		if (fillColor == largeColor)
			throw ITLException("Fill color and large color must not be equal.");
		if (largeColor == 0)
			throw ITLException("Large color must not be zero as zero is background color.");

		size_t counter;

		results.headers() = analyzers.headers();

		prepareParticleAnalysis(image, fillColor, largeColor);

		internals::analyzeParticlesSingleBlock(image, analyzers, results, nullptr, nullptr, connectivity, volumeLimit, fillColor, largeColor, counter, image.depth() * image.height(), Vec3sc());

		//raw::writed(image, "particleanalysis/labels");
	}

	/**
	Perform particle analysis using single-threaded method.
	This one is faster for very small images or images containing only a very few particles.
	Needs image that can store at least four distinct colours - original particle color, temporary color, large particle color and background color.
	Regions colored with large particle color are skipped. Additionally, particles that are larger than volumeLimit will be colored with large
	particle color. Other particles will be colored with fill color.
	@param image The image containing the particles.
	@param analyzerNames List of names of analyzers to apply to the particles.
	@param results Results table.
	@param connectivity Connectivity of the particles.
	@param volumeLimit Only include particles smaller than this value in the results.
	@param fillColor The analyzed particles will be colored with this color.
	@param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
	*/
	template<typename pixel_t> void analyzeParticlesSingleThreaded(Image<pixel_t>& image, const string& analyzerNames, Results& results, Connectivity connectivity = Connectivity::AllNeighbours, size_t volumeLimit = 0, pixel_t fillColor = internals::SpecialColors<pixel_t>::fillColor(), pixel_t largeColor = internals::SpecialColors<pixel_t>::largeColor())
	{
		auto analyzers = createAnalyzers<pixel_t>(analyzerNames, image.dimensions());
		analyzeParticlesSingleThreaded(image, analyzers, results, connectivity, volumeLimit, fillColor, largeColor);
	}

	/**
	Analyzes labeled regions of the input image.
	The regions do not need to be connected.
	Region having value zero is skipped.
	*/
	template<typename pixel_t> void analyzeLabels(const Image<pixel_t>& image, AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results)
	{
		// TODO: This can be done more efficiently by changing analyzers such that they don't require all the points at once.

		// Divide pixels into point sets based on their value
		std::map<pixel_t, std::vector<Vec3sc> > points;
		for (coord_t z = 0; z < image.depth(); z++)
		{
			for (coord_t y = 0; y < image.height(); y++)
			{
				for (coord_t x = 0; x < image.width(); x++)
				{
					pixel_t pixel = image(x, y, z);
					if (pixel != 0)
					{
						points[pixel].push_back(Vec3sc(Vec3c(x, y, z)));
					}
				}
			}
		}

		// Analyze each set
		results.headers() = analyzers.headers();
		for (auto& item : points)
		{
			std::vector<double> resultLine;
			analyzers.analyze(item.second, resultLine);
			results.push_back(resultLine);
		}
	}

	/**
	Analyzes labeled regions of the input image.
	The regions do not need to be connected.
	Region having value zero is skipped.
	*/
	template<typename pixel_t> void analyzeLabels(const Image<pixel_t>& image, const std::string& analyzerNames, Results& results)
	{
		auto analyzers = createAnalyzers<pixel_t>(analyzerNames, image.dimensions());
		analyzeLabels(image, analyzers, results);
	}

	/**
	Perform greedy coloring of regions.
	Colors each region in image such that its neighbours are colored with different colors, and uses as little colors as possible.
	Uses greedy algorithm so the count of colors used might not be minimal.
	Assumes background to have value 0.
    @param image The image containing the particles. Background is assumed to have value 0. The image must be able to store at least one value that it does not contain at input.
    */
    template<typename T> void greedyColoring(Image<T>& image, Connectivity connectivity = Connectivity::AllNeighbours)
    {
		
		// Find suitable fill color
		T fillColor = findUnusedValue(image);


		// Paint regions
        
		std::set<T> usedColors;
		std::set<T> neighbourColors;

		std::vector<Vec3sc> particlePoints;
        particlePoints.reserve(1000);

		for (coord_t z = 0; z < image.depth(); z++)
		{
			for (coord_t y = 0; y < image.height(); y++)
			{
				for (coord_t x = 0; x < image.width(); x++)
				{
					Vec3c p(x, y, z);
					T pixel = image(p);
					if (pixel != 0 && pixel != fillColor && usedColors.find(pixel) == usedColors.end())
					{
						floodfill(image, p, fillColor, fillColor, connectivity, nullptr, &particlePoints, 0, &neighbourColors);

						// Get color not in use by neighbours
						T newColor = 1;
						for (; newColor < std::numeric_limits<T>::max(); newColor++)
						{
							if (newColor != fillColor && neighbourColors.find(newColor) == neighbourColors.end())
								break;
						}

						usedColors.insert(newColor);

						// Fill the particle with that color
						for (size_t n = 0; n < particlePoints.size(); n++)
							image(particlePoints[n]) = newColor;

						particlePoints.clear();
						neighbourColors.clear();
					}

				}
			}

			showProgress(z, image.depth());
		}

    }


    /**
	Fill all particles that have some measurement value above or below some threshold.
	@param image The image.
	@param results Particle analysis results table. The table must contain results of 'coordinates' analyzer.
	@param fillColor Color to use when filling the particles.
	@param condi Index of column containing the value to test.
	@param condThreshold Threshold for the condition value.
	@param higher If set to true, fill if condition value for particle is larger than the threshold value; If set to false, fill if condition value is lower than threshold value.
	*/
    template<typename pixel_t> void fillParticles(Image<pixel_t>& image, const Results& results, pixel_t fillColor, int condi, double condThreshold, Connectivity conn, bool higher = false)
	{
		size_t xi = results.getColumnIndex("X");
		size_t yi = results.getColumnIndex("Y");
		size_t zi = results.getColumnIndex("Z");

		prepareParticleAnalysis(image);

		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)results.size(); n++)
		{
			double value = results[n][condi];
			bool fill = higher ? value > condThreshold : value <= condThreshold;

			if(fill)
			{
				Vec3c pos;
				pos.x = round(results[n][xi]);
				pos.y = round(results[n][yi]);
				pos.z = round(results[n][zi]);

				floodfill<pixel_t>(image, pos, fillColor, fillColor, conn);
			}

			showThreadProgress(counter, results.size());
		}
	}

	/**
	Fill all particles with their corresponding measurement value.
	@param image The image.
	@param results Particle analysis results table. The table must contain results of 'coordinates' analyzer.
	@param condi Index of column containing the fill value
	@param multiplier Multiplier for fill value.
	*/
    template<typename pixel_t> void fillParticles(Image<pixel_t>& image, const Results& results, int condi, double multiplier, Connectivity conn)
	{
		size_t xi = results.getColumnIndex("X");
		size_t yi = results.getColumnIndex("Y");
		size_t zi = results.getColumnIndex("Z");

		prepareParticleAnalysis(image);

		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)results.size(); n++)
		{
			pixel_t value = pixelRound<pixel_t>(results[n][condi] * multiplier);

			Vec3c pos;
			pos.x = round(results[n][xi]);
			pos.y = round(results[n][yi]);
			pos.z = round(results[n][zi]);

			floodfill<pixel_t>(image, pos, value, value, conn);

			showThreadProgress(counter, results.size());
		}
	}

	/**
	Fill all particles with given color.
	@param image The image.
	@param results Particle analysis results table. The table must contain results of 'coordinates' analyzer.
	@param fillColor Fill color.
	@param shift Shift that will be added to the particle coordinates before filling.
	*/
    template<typename pixel_t> void fillParticles(Image<pixel_t>& image, const Results& results, pixel_t fillColor, Connectivity conn, const Vec3c& shift = Vec3c())
	{
		size_t xi = results.getColumnIndex("X");
		size_t yi = results.getColumnIndex("Y");
		size_t zi = results.getColumnIndex("Z");

		prepareParticleAnalysis(image);

		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)results.size(); n++)
		{
			Vec3c pos;
			pos.x = round(results[n][xi]);
			pos.y = round(results[n][yi]);
			pos.z = round(results[n][zi]);

			floodfill<pixel_t>(image, pos + shift, fillColor, fillColor, conn);

			showThreadProgress(counter, results.size());
		}
	}

	/**
	Enumerates types of ellipsoids that can be used for visualization of particles.
	*/
	enum class EllipsoidType
	{
		/**
		Draw the principal ellipsoid as is without scaling.
		*/
		Principal,
		/**
		Scale the principal ellipsoid such that all particle points fit inside it.
		*/
		BoundingEllipsoid,
		/**
		Draw bounding sphere of the particle.
		*/
		BoundingSphere,
		/**
		Scale the principal ellipsoid such that its volume equals the volume of the particle
		*/
		CorrectVolume
	};

	template<>
	inline std::string toString(const EllipsoidType& x)
	{
		switch (x)
		{
		case EllipsoidType::Principal: return "Principal";
		case EllipsoidType::BoundingEllipsoid: return "BoundingEllipsoid";
		case EllipsoidType::BoundingSphere: return "BoundingSphere";
		case EllipsoidType::CorrectVolume: return "Volume";
		}
		throw ITLException("Invalid ellipsoid type.");
	}

	template<>
	inline EllipsoidType fromString(const string& dt)
	{
		string dt2 = dt;
		trim(dt2);
		toLower(dt2);
		if (dt2 == "principal")
			return EllipsoidType::Principal;
		if (dt2 == "boundingellipsoid" || dt2 == "bounding ellipsoid")
			return EllipsoidType::BoundingEllipsoid;
		if (dt2 == "boundingsphere" || dt2 == "bounding sphere")
			return EllipsoidType::BoundingSphere;
		if (dt2 == "correctvolume" || dt2 == "volume")
			return EllipsoidType::CorrectVolume;

		throw ITLException(string("Unsupported ellipsoid type: ") + dt);
	}

	/**
	Visualizes particles as ellipsoids.
	@param results Results table. If ellipsoid type is not BoundingSphere, the table must contain columns CX, CY, CZ, l1, l2, l3, phi1, theta1, phi2, theta2, and bounding scale (i.e. most of the results of the PCA analyzer). Otherwise, the table must contain 'bounding sphere X', 'bounding sphere Y', 'bounding sphere Z', and 'bounding sphere radius', i.e. the output of the boundingsphere analyzer.
	*/
	template<typename pixel_t> void drawEllipsoids(Image<pixel_t>& image, const Results& results, pixel_t color, EllipsoidType type)
	{
		size_t cxi = 0;
		size_t cyi = 0; 
		size_t czi = 0; 
		size_t l1i = 0; 
		size_t l2i = 0; 
		size_t l3i = 0; 
		size_t phi1i = 0; 
		size_t theta1i = 0; 
		size_t phi2i = 0;
		size_t theta2i = 0;
		if (type != EllipsoidType::BoundingSphere)
		{
			cxi = results.getColumnIndex("CX [pixel]");
			cyi = results.getColumnIndex("CY [pixel]");
			czi = results.getColumnIndex("CZ [pixel]");
			l1i = results.getColumnIndex("l1 [pixel]");
			l2i = results.getColumnIndex("l2 [pixel]");
			l3i = results.getColumnIndex("l3 [pixel]");
			phi1i = results.getColumnIndex("phi1 [rad]");
			theta1i = results.getColumnIndex("theta1 [rad]");
			phi2i = results.getColumnIndex("phi2 [rad]");
			theta2i = results.getColumnIndex("theta2 [rad]");
		}

		size_t boundingScalei = 0;
		if(type == EllipsoidType::BoundingEllipsoid)
			boundingScalei = results.getColumnIndex("bounding scale");

		size_t volumei = 0;
		if(type == EllipsoidType::CorrectVolume)
			volumei = results.getColumnIndex("Volume [pixel]");

		size_t bsxi = 0; 
		size_t bsyi = 0; 
		size_t bszi = 0; 
		size_t bsri = 0; 
		if (type == EllipsoidType::BoundingSphere)
		{
			bsxi = results.getColumnIndex("bounding sphere X [pixel]");
			bsyi = results.getColumnIndex("bounding sphere Y [pixel]");
			bszi = results.getColumnIndex("bounding sphere Z [pixel]");
			bsri = results.getColumnIndex("bounding sphere radius [pixel]");
		}

		for(size_t n = 0; n < results.size(); n++)
		{
			if (type != EllipsoidType::BoundingSphere)
			{
				Vec3d pos;
				pos.x = results[n][cxi];
				pos.y = results[n][cyi];
				pos.z = results[n][czi];

				double l1 = results[n][l1i];
				double l2 = results[n][l2i];
				double l3 = results[n][l3i];

				double phi1 = results[n][phi1i];
				double phi2 = results[n][phi2i];

				double theta1 = results[n][theta1i];
				double theta2 = results[n][theta2i];

				if (type == EllipsoidType::BoundingEllipsoid)
				{
					double boundingScale = results[n][boundingScalei];

					// Scale the ellipsoid so that it becomes bounding ellipsoid
					l1 *= boundingScale;
					l2 *= boundingScale;
					l3 *= boundingScale;
				}
				else if (type == EllipsoidType::CorrectVolume)
				{
					double particleVolume = results[n][volumei];
					double ellipsoidVolume = 4.0 / 3.0 * PI * l1 * l2 * l3;
					double scale = std::pow(particleVolume / ellipsoidVolume, 1.0 / 3.0);
					l1 *= scale;
					l2 *= scale;
					l3 *= scale;
				}

				// Only draw if the ellipsoid is sane!
				if (!std::isinf(l1) && !std::isnan(l1) && l1 != 0 &&
					!std::isinf(l2) && !std::isnan(l2) && l2 != 0 &&
					!std::isinf(l3) && !std::isnan(l3) && l3 != 0)
					draw(image, Ellipsoid(pos, Vec3d(l1, l2, l3), phi1, theta1, phi2, theta2), color);
			}
			else
			{
				// Draw bounding sphere
				Vec3d pos;
				pos.x = results[n][bsxi];
				pos.y = results[n][bsyi];
				pos.z = results[n][bszi];
				double r = results[n][bsri];

				draw(image, Sphere(pos, r), color);
			}
			showProgress(n, results.size());
		}
	}

	namespace tests
	{
		void analyzeParticlesThreading();
		void analyzeParticlesThreadingBig();
		void analyzeParticlesSanity();
		void analyzeParticlesSanity2();
		void analyzeParticlesVolumeLimit();
	}
}

