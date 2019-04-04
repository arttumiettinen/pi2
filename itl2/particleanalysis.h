#pragma once

#include <iostream>
#include <algorithm>

#include "image.h"
#include "floodfill.h"
#include "utilities.h"

#include "math/vec2.h"
#include "math/vec3.h"

#include "indexforest.h"
#include "box.h"

#include "resultstable.h"
#include "analyzers.h"


namespace itl2
{

	using math::Vec2c;
	using math::Vec3c;
	using math::Vec3sc;
	using math::pixelRound;

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
	@param image Image containing the particles.
	@param particleColor Color of particles. This color will be skipped in labeling. Pass zero to label particles of any color.
	@param firstLabelValue Value of the first label.
	@returns The largest label value used in the image.
	*/
	template<typename pixel_t> pixel_t labelParticles(Image<pixel_t>& image, pixel_t particleColor = 0, pixel_t firstLabelValue = 1, Connectivity connectivity = Connectivity::NearestNeighbours, bool showProgress = true)
	{
		size_t counter = 0;

		// Handle "color all particles" case first
		if (particleColor == 0)
		{
			particleColor = numeric_limits<pixel_t>::max();

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
						if (particleNumber >= numeric_limits<pixel_t>::max() - 2)
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
			static T fillColor();
			static T largeColor();
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

		template<> class SpecialColors<uint16_t>
		{
		public:
			static uint16_t fillColor()
			{
				return numeric_limits<uint16_t>::max();
			}

			static uint16_t largeColor()
			{
				return numeric_limits<uint16_t>::max() - 1;
			}
		};

		template<> class SpecialColors<uint32_t>
		{
		public:
			static uint32_t fillColor()
			{
				return numeric_limits<uint32_t>::max();
			}

			static uint32_t largeColor()
			{
				return numeric_limits<uint32_t>::max() - 1;
			}
		};

		template<> class SpecialColors<uint64_t>
		{
		public:
			static uint64_t fillColor()
			{
				return numeric_limits<uint64_t>::max();
			}

			static uint64_t largeColor()
			{
				return numeric_limits<uint64_t>::max() - 1;
			}
		};

		/**
		Tests if Z coordinate of any particle point is 0 or maxZ.
		*/
		inline bool isOnZEdge(const vector<Vec3sc>& particle, coord_t maxZ)
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
		template<typename pixel_t> void analyzeParticlesSingleBlock(Image<pixel_t>& image, AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results, vector<vector<Vec3sc> >* pIncompleteParticles, vector<vector<Vec3sc> >* pLargeEdgePoints, Connectivity connectivity, size_t volumeLimit, pixel_t fillColor, pixel_t largeColor, size_t& counter, size_t counterMax, const Vec3sc& coordinateShift)
		{
			vector<Vec3sc> particlePoints;
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
								set<pixel_t> neighbourColors;
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
									vector<double> resultLine;
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
									vector<Vec3sc> edgePoints;

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
		template<typename pixel_t> void analyzeParticlesBlocks(Image<pixel_t>& image, AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results, vector<vector<Vec3sc> >& largeEdgePoints, vector<vector<Vec3sc> >& incompleteParticles, Connectivity connectivity, size_t volumeLimit, pixel_t fillColor, pixel_t largeColor, const Vec3sc& origin, vector<coord_t>& blockEdgeZ)
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
						cout << "Analyzing in " << count << " blocks." << endl;
					}
				}

				// Calculate amount of slices single thread should process.
				// The value is capped so that we don't divide the work unnecessarily too much if the image is small.
				coord_t size = image.depth() / count;
				if (size < 50)
					size = 50;

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
//						cout << "Thread " << idx << " analyzes z-range [" << minZ << ", " << maxZ << "[" << endl;
//					}

					// Get view of part of the image
					Image<pixel_t> block(image, minZ, maxZ);

					// Analyze the block. All coordinates are shifted in analyzeParticlesSingleBlock function to global original image coordinates.
					Results blockResults;
					vector<vector<Vec3sc> > blockLargeEdgePoints;
					vector<vector<Vec3sc> > blockIncompleteParticles;
					internals::analyzeParticlesSingleBlock(block, analyzers, blockResults, &blockIncompleteParticles, &blockLargeEdgePoints, connectivity, volumeLimit, fillColor, largeColor, counter, image.depth() * image.height(), origin + math::Vec3sc(0, 0, (int32_t)minZ));

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
		inline bool isNeighbouring(const vector<Vec3sc>& points1, const vector<Vec3sc>& points2, Connectivity connectivity)
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
		inline void findNeighboringParticles(const vector<vector<Vec3sc> >& edgePoints, IndexForest& forest, Connectivity connectivity)
		{
			// Add all indices to the forest
			cout << "Build forest..." << endl;
			forest.initialize(edgePoints.size());

			// Calculate bounding boxes for incomplete vertex sets
			cout << "Determine bounding boxes..." << endl;
			vector<Box<int32_t> > boundingBoxes;
			boundingBoxes.reserve(edgePoints.size());
			Vec3sc m = Vec3sc(numeric_limits<int32_t>::max(), numeric_limits<int32_t>::max(), numeric_limits<int32_t>::max());
			Vec3sc M = Vec3sc(numeric_limits<int32_t>::lowest(), numeric_limits<int32_t>::lowest(), numeric_limits<int32_t>::lowest());
			for (size_t n = 0; n < edgePoints.size(); n++)
			{
				Box<int32_t> b = Box<int32_t>::boundingBox(edgePoints[n]);
				if (edgePoints[n].size() > 0)
				{
					b.inflate(1);

					m = min(m, b.minc);
					M = max(M, b.maxc);
				}
				boundingBoxes.push_back(b);
			}

			Image<vector<size_t> > grid(100, 100, 100);

			cout << "Divide boxes to a grid..." << endl;
			size_t counter = 0;
			for (size_t n = 0; n < boundingBoxes.size(); n++)
			{
				if (edgePoints[n].size() > 0)
				{
					Vec3c mcell((coord_t)floor((double)(boundingBoxes[n].minc.x - m.x) / (double)(M.x - m.x) * (grid.width() - 1)),
						(coord_t)floor((double)(boundingBoxes[n].minc.y - m.y) / (double)(M.y - m.y) * (grid.height() - 1)),
						(coord_t)floor((double)(boundingBoxes[n].minc.z - m.z) / (double)(M.z - m.z) * (grid.depth() - 1)));

					Vec3c Mcell((coord_t)ceil((double)(boundingBoxes[n].maxc.x - m.x) / (double)(M.x - m.x) * (grid.width() - 1)),
						(coord_t)ceil((double)(boundingBoxes[n].maxc.y - m.y) / (double)(M.y - m.y) * (grid.height() - 1)),
						(coord_t)ceil((double)(boundingBoxes[n].maxc.z - m.z) / (double)(M.z - m.z) * (grid.depth() - 1)));

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

			cout << "Determine overlaps..." << endl;
			counter = 0;
			
			#pragma omp parallel if(!omp_in_parallel())
			{
				set<tuple<size_t, size_t> > nonOverlapping;

				#pragma omp for schedule(dynamic)
				for (coord_t i = 0; i < grid.pixelCount(); i++)
				{
					auto& list = grid(i);

					for (coord_t ni = 0; ni < (coord_t)list.size(); ni++)
					{
						size_t n = list[ni];
						const vector<Vec3sc>& v1 = edgePoints[n];
						const Box<int32_t>& b1 = boundingBoxes[n];

						if (v1.size() > 0)
						{
							for (size_t mi = ni + 1; mi < list.size(); mi++)
							{
								size_t m = list[mi];
								const vector<Vec3sc>& v2 = edgePoints[m];
								const Box<int32_t>& b2 = boundingBoxes[m];

								//if (math::min(n, m) == 308 && math::max(n, m) == 417)
								//{
								//	cout << "testing " << n << ", " << m << endl;
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
											auto key = make_tuple(math::min(n, m), math::max(n, m));

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

		inline void concatAndShrink(vector<Vec3sc>& target, vector<Vec3sc>& source)
		{
			target.insert(target.end(), source.begin(), source.end());
			source.clear();
			source.shrink_to_fit();
		}

		/**
		Combines incomplete particles and analyzes them.
		*/
		template<typename pixel_t> void combineParticleAnalysisResults(const AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results, const vector<vector<Vec3sc> >& largeEdgePoints, vector<vector<Vec3sc> >& incompleteParticles, size_t volumeLimit, Connectivity connectivity, const vector<coord_t>& blockEdgeZ)
		{
			if (incompleteParticles.size() <= 0)
				return;

			cout << "Processing " << incompleteParticles.size() << " particles on block edge boundaries..." << endl;

			// Create list of edge points for each incomplete particle
			cout << "Find edge points of particles..." << endl;
			vector<vector<Vec3sc> > edgePoints;
			size_t counter = 0;
			for (coord_t n = 0; n < (coord_t)incompleteParticles.size(); n++)
			{
				auto& points = incompleteParticles[n];
				vector<Vec3sc> currEdgePoints;
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
			edgePoints.insert(edgePoints.end(), largeEdgePoints.begin(), largeEdgePoints.end());

			vector<bool> isBig(incompleteParticles.size(), false);

			// Find overlapping particles (based on edge points) and combine them
			{
				IndexForest forest;
				findNeighboringParticles(edgePoints, forest, connectivity);

				// Move all points to roots
				cout << "Combine particles..." << endl;
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
							// Mark particle n as big and clear its points.
							isBig[n] = true;
							incompleteParticles[n].clear();
							incompleteParticles[n].shrink_to_fit();
							edgePoints[n].clear();
							edgePoints[n].shrink_to_fit();
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


			// Analyze combined partiles
			cout << "Analyze combined particles..." << endl;
			counter = 0;
			#pragma omp parallel for if(!omp_in_parallel() && incompleteParticles.size() > 1000) schedule(dynamic)
			for (coord_t n = 0; n < (coord_t)incompleteParticles.size(); n++)
			{

				auto& points = incompleteParticles[n];
				//if (points.size() > 0)
				if(points.size() > 0)
				{
					//bool isLarge = points.size() >= volumeLimit;

					//if (!isLarge) // Test if the particle touches a position known to belong to large particle?
					//{
					//	
					//	isLarge = isNeighbouring(edgePoints[n], largeEdgePoints, connectivity);
					//}

					//if (!isLarge)
					//{

					if(volumeLimit <= 0 || (volumeLimit > 0 && !isBig[n]))
					{
						// Analyze the particle
						// Note: the parts of particles originating from different calculation blocks
						// do not overlap so there are no duplicate points in the point list
						//Vec3c point = points[0];
						vector<double> resultLine;
						analyzers.analyze(points, resultLine);

						#pragma omp critical(resultsInsert)
						{
							results.push_back(resultLine);
						}
					}

				}
				showThreadProgress(counter, incompleteParticles.size());
			}

			// Clear incomplete particles so that doing the same call again does not change anything.
			incompleteParticles.clear();
		}
	}

	/**
	Perform particle analysis.
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
	template<typename pixel_t> void analyzeParticles(Image<pixel_t>& image, AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results, Connectivity connectivity = Connectivity::AllNeighbours, size_t volumeLimit = 0, pixel_t fillColor = internals::SpecialColors<pixel_t>::fillColor(), pixel_t largeColor = internals::SpecialColors<pixel_t>::largeColor())
	{
		if (fillColor == largeColor)
			throw ITLException("Fill color and large color must not be equal.");
		if (largeColor == 0)
			throw ITLException("Large color must not be zero as zero is background color.");

		vector<vector<Vec3sc> > incompleteParticles;
		vector<vector<Vec3sc> > largeEdgePoints;
		vector<coord_t> edgeZ;
		
		results.headers() = analyzers.headers();

		internals::analyzeParticlesBlocks(image, analyzers, results, largeEdgePoints, incompleteParticles, connectivity, volumeLimit, fillColor, largeColor, Vec3sc(), edgeZ);

		//cout << "Block edges are at" << endl;
		//sort(edgeZ.begin(), edgeZ.end());
		//for (coord_t z : edgeZ)
		//	cout << z << endl;

		//Image<uint16_t> labels(image.dimensions());
		//for (size_t n = 0; n < incompleteParticles.size(); n++)
		//{
		//	for (size_t m = 0; m < incompleteParticles[n].size(); m++)
		//	{
		//		labels(incompleteParticles[n][m]) = (uint16_t)(n);
		//	}
		//}
		//raw::writed(labels, "particleanalysis/labels");

		internals::combineParticleAnalysisResults(analyzers, results, largeEdgePoints, incompleteParticles, volumeLimit, connectivity, edgeZ);
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
		auto analyzers = createAnalyzers<uint8_t>(analyzerNames, image.dimensions());
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

		internals::analyzeParticlesSingleBlock(image, analyzers, results, 0, 0, connectivity, volumeLimit, fillColor, largeColor, counter, image.depth() * image.height(), Vec3sc());

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
		auto analyzers = createAnalyzers<uint8_t>(analyzerNames, image.dimensions());
		analyzeParticlesSingleThreaded(image, analyzers, results, connectivity, volumeLimit, fillColor, largeColor);
	}

	///**
	//Perform greedy coloring of particles in 3D.
	//Colors each region in image such that its neighbours are colored with different colors.
	//Uses greedy algorithm.
	//Needs image that can store at least four distinct colours - original particle color, temporary color, big particle color and background color.
 //   Particles colored with big particle color are skipped. Additionally, particles that are larger than volumeLimit will be colored with big
 //   particle color. Other particles will be colored with colors such that neighbouring particles are colored with different colors.
 //   @param image The image containing the particles.
 //   @param volumeLimit Only include particles smaller than this value in the results.
 //   @param fillColor The analyzed particles will be colored with this color.
 //   @param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
 //   @param backgroundColor Color of background.
 //   */
 //   template<typename T> void colorMapGreedy3D(Image<T>& image, Connectivity connectivity = All, size_t volumeLimit = 0, T fillColor = internal::SpecialColors<T>::fillColor(), T largeColor = internal::SpecialColors<T>::largeColor(), T backgroundColor = 0)
 //   {
	//	if(image.dimensionality() != 3)
	//		throw ITLException("This method supports only 3-dimensional images.");

 //       if(fillColor == largeColor)
 //           throw ITLException("Fill color and large color must not be equal.");
 //       if(largeColor == backgroundColor)
 //           throw ITLException("Background color and large color must not be equal.");

 //       CursorPointer<T, LinearCursor<T> > cursor(image, image.createLinearCursor());

	//	set<T> usedColors;

	//	set<T> neighbourColors;

 //       vector<Vec3c> particlePoints;
 //       particlePoints.reserve(1000);

 //       size_t n = 0;
 //       for(; cursor->hasNext(); cursor->proceed())
 //       {
 //           T pixel = **cursor;
 //           if(pixel != fillColor && pixel != backgroundColor && pixel != largeColor && usedColors.find(pixel) == usedColors.end())
 //           {
 //               if(singlethreaded::floodfill3D(image, cursor->getPosition(), fillColor, largeColor, connectivity, &particlePoints, volumeLimit, &neighbourColors))
 //               {
 //                   // Analyze the points only if the volume of the particle is small enough.

	//				// Get color not in use by neighbours
	//				T newColor = 0;
	//				for(; newColor < numeric_limits<T>::max(); newColor++)
	//				{
	//					if(newColor != fillColor && newColor != backgroundColor && newColor != largeColor && neighbourColors.find(newColor) == neighbourColors.end())
	//						break;
	//				}

	//				usedColors.insert(newColor);

	//				// Fill the particle with that color
	//				for(size_t n = 0; n < particlePoints.size(); n++)
 //                       image.setPixel(particlePoints[n], newColor);

 //               }
 //               else
 //               {
 //                   // The fill was ended by size limit. Mark the filled points as belonging to a large particle.

 //                   for(size_t n = 0; n < particlePoints.size(); n++)
 //                       image.setPixel(particlePoints[n], largeColor);
 //               }

 //               particlePoints.clear();
	//			neighbourColors.clear();
 //           }

 //           n++;
	//		utilities::showProgress(n, image.pixelCount());
 //       }

 //   }


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

		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)results.size(); n++)
		{
			double value = results[n][condi];
			bool fill = higher ? value > condThreshold : value <= condThreshold;

			if(fill)
			{
				Vec3c pos;
				pos.x = math::round(results[n][xi]);
				pos.y = math::round(results[n][yi]);
				pos.z = math::round(results[n][zi]);

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

		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)results.size(); n++)
		{
			pixel_t value = pixelRound<pixel_t>(results[n][condi] * multiplier);

			Vec3c pos;
			pos.x = math::round(results[n][xi]);
			pos.y = math::round(results[n][yi]);
			pos.z = math::round(results[n][zi]);

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

		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)results.size(); n++)
		{
			Vec3c pos;
			pos.x = math::round(results[n][xi]);
			pos.y = math::round(results[n][yi]);
			pos.z = math::round(results[n][zi]);

			floodfill<pixel_t>(image, pos + shift, fillColor, fillColor, conn);

			showThreadProgress(counter, results.size());
		}
	}

	///**
	//Gets value of left side of ellipsoid equation at p for ellipsoid located at c that has
	//semi-axis lengths l1, l2 and l3	and whose semi-axis orientations are given by (phiN, thetaN), N=1..3.
	//Points inside the ellipsoid are characterized by return value <= 1.
	//*/
	//inline double getEllipsoidFunctionValue(const Vec3d& p,
	//							const Vec3d& c,
	//							double l1, double l2, double l3,
	//							double phi1, double theta1,
	//							double phi2, double theta2,
	//							double phi3, double theta3)
	//{
	//	// Make the ellipsoid centered to origin
	//	Vec3d pdot = p - c;

	//	// Rotate pdot such that ellipsoid axes are aligned with coordinate axes
	//	Vec3d u1 = tocartesian(1.0, phi1, theta1);
	//	Vec3d u2 = tocartesian(1.0, phi2, theta2);
	//	Vec3d u3 = tocartesian(1.0, phi3, theta3);

	//	Matrix3x3d R(u1.x, u2.x, u3.x,
	//				 u1.y, u2.y, u3.y,
	//				 u1.z, u2.z, u3.z);

	//	R.transpose();
	//	pdot = R * pdot;

	//	//Matrix3x3d Ri;
	//	//R.inverse(Ri);
	//	//pdot = Ri * pdot;

	//	// Use ellipsoid equation
	//	double f = (pdot.x * pdot.x) / (l1 * l1) + (pdot.y * pdot.y) / (l2 * l2) + (pdot.z * pdot.z) / (l3 * l3);

	//	return f;
	//}

	///**
	//Test if point p is in ellipsoid located at c that has semi-axis lengths l1, l2 and l3
	//and whose semi-axis orientations are given by (phiN, thetaN), N=1..3.
	//*/
	//inline bool isInEllipsoid(const Vec3d& p,
	//							const Vec3d& c,
	//							double l1, double l2, double l3,
	//							double phi1, double theta1,
	//							double phi2, double theta2,
	//							double phi3, double theta3)
	//{
	//	double f = getEllipsoidFunctionValue(p, c, l1, l2, l3, phi1, theta1, phi2, theta2, phi3, theta3);
	//	return f <= 1;
	//}

	///**
	//Test if point p is in ellipsoid located at c that has semi-axis lengths l1, l2 and l3
	//and that can be transformed to axis-aligned ellipsoid by rotating with rotation matrix Rinv
	//*/
	//inline bool isInEllipsoid(const Vec3d& p,
	//							const Vec3d& c,
	//							double l1, double l2, double l3,
	//							const Matrix3x3d& Rinv)
	//{
	//	// Make the ellipsoid centered to origin
	//	Vec3d pdot = p - c;

	//	// Rotate pdot such that ellipsoid axes are aligned with coordinate axes
	//	pdot = Rinv * pdot;

	//	// Use ellipsoid equation
	//	double f = (pdot.x * pdot.x) / (l1 * l1) + (pdot.y * pdot.y) / (l2 * l2) + (pdot.z * pdot.z) / (l3 * l3);

	//	return f <= 1;
	//}

	///**
	//Draws single ellipsoid.
	//@param pMask Positions where *pMask == 0 are not filled.
	//@return Count of filled pixels.
	//*/
	//template<typename T, typename Tmask> unsigned long long drawEllipsoid(Image<T>& image, const Vec3d& pos,
	//																		double l1, double l2, double l3,
	//																		double phi1, double theta1,
	//																		double phi2, double theta2,
	//																		double phi3, double theta3,
	//																		T color,
	//																		Image<Tmask>* pMask)
	//{
	//	double maxl = mathutils::max(l1, l2);
	//	maxl = mathutils::max(maxl, l3);

	//	Vec3c minPos = round(pos - Vec3d(maxl, maxl, maxl));
	//	Vec3c maxPos = round(pos + Vec3d(maxl, maxl, maxl));
	//	
	//	Vec3c dimensions((int)image.getDimension(0), (int)image.getDimension(1), (int)image.getDimension(2));
	//	clamp(minPos, Vec3c(0, 0, 0), dimensions);
	//	clamp(maxPos, Vec3c(0, 0, 0), dimensions);

	//	Vec3d u1 = tocartesian(1.0, phi1, theta1);
	//	Vec3d u2 = tocartesian(1.0, phi2, theta2);
	//	Vec3d u3 = tocartesian(1.0, phi3, theta3);

	//	Matrix3x3d R(u1.x, u2.x, u3.x,
	//				 u1.y, u2.y, u3.y,
	//				 u1.z, u2.z, u3.z);

	//	R.transpose();
	//	
	//	unsigned long long filledCount = 0;
	//	#pragma omp parallel for reduction(+:filledCount)
	//	for(coord_t z = minPos.z; z < maxPos.z; z++)
	//	{
	//		for(size_t y = minPos.y; y < maxPos.y; y++)
	//		{
	//			for(size_t x = minPos.x; x < maxPos.x; x++)
	//			{
	//				if(pMask == 0 || pMask->getPixel((coord_t)x, (coord_t)y, (coord_t)z) != 0)
	//				{
	//					if(isInEllipsoid(Vec3d((double)x, (double)y, (double)z), pos, l1, l2, l3, R))
	//					{
	//						image.setPixel((coord_t)x, (coord_t)y, (coord_t)z, color);
	//						filledCount++;
	//					}
	//				}
	//			}
	//		}
	//	}
	//	
	//	//line<T>(image, pos, pos + tocartesian(l1, phi1, theta1), 200);
	//	//line<T>(image, pos, pos + tocartesian(l2, phi2, theta2), 220);
	//	//line<T>(image, pos, pos + tocartesian(l3, phi3, theta3), 240);

	//	return filledCount;
	//}

	///**
	//Draws bounding PCA ellipsoids.
	//@param pMask Pointer to mask image. Only pixels where mask != 0 will be colored.
	//@param pResults Pointer to results array. If not null, count of filled pixels for each bounding ellipsoid is added to the array.
	//*/
	//template<typename T, typename Tmask> void drawEllipsoids(Image<T>& image, const vector<string>& headers, const vector<vector<double> >& results, T ellipsoidColor,
	//										Image<Tmask>* pMask = 0, vector<unsigned long long>* pResults = 0)
	//											/*int cxi = 5, int cyi = 6, int czi = 7,
	//											int l1i = 9, int l2i = 10, int l3i = 11,
	//											int phi1i = 12, int theta1i = 13,
	//											int phi2i = 14, int theta2i = 15,
	//											int phi3i = 16, int theta3i = 17)*/
	//{
	//	size_t cxi = getColumnIndex("CX [pixel]", headers);
	//	size_t cyi = getColumnIndex("CY [pixel]", headers);
	//	size_t czi = getColumnIndex("CZ [pixel]", headers);
	//	size_t l1i = getColumnIndex("l1 [pixel]", headers);
	//	size_t l2i = getColumnIndex("l2 [pixel]", headers);
	//	size_t l3i = getColumnIndex("l3 [pixel]", headers);
	//	size_t phi1i = getColumnIndex("phi1 [rad]", headers);
	//	size_t theta1i = getColumnIndex("theta1 [rad]", headers);
	//	size_t phi2i = getColumnIndex("phi2 [rad]", headers);
	//	size_t theta2i = getColumnIndex("theta2 [rad]", headers);
	//	size_t phi3i = getColumnIndex("phi3 [rad]", headers);
	//	size_t theta3i = getColumnIndex("theta3 [rad]", headers);
	//	size_t d1i = getColumnIndex("d1 [pixel]", headers);
	//	size_t d2i = getColumnIndex("d2 [pixel]", headers);
	//	size_t d3i = getColumnIndex("d3 [pixel]", headers);
	//	size_t boundingScalei = getColumnIndex("bounding scale", headers);

	//	Vec3d pos;
	//	for(size_t n = 0; n < results.size(); n++)
	//	{
	//		pos.x = results[n][cxi];
	//		pos.y = results[n][cyi];
	//		pos.z = results[n][czi];

	//		double l1 = results[n][l1i];
	//		double l2 = results[n][l2i];
	//		double l3 = results[n][l3i];

	//		double phi1 = results[n][phi1i];
	//		double phi2 = results[n][phi2i];
	//		double phi3 = results[n][phi3i];

	//		double theta1 = results[n][theta1i];
	//		double theta2 = results[n][theta2i];
	//		double theta3 = results[n][theta3i];

	//		double d1 = results[n][d1i];
	//		double d2 = results[n][d2i];
	//		double d3 = results[n][d3i];

	//		double boundingScale = results[n][boundingScalei];

	//		// Scale the ellipsoid so that it becomes bounding ellipsoid
	//		l1 = boundingScale * sqrt(l1);
	//		l2 = boundingScale * sqrt(l2);
	//		l3 = boundingScale * sqrt(l3);

	//		unsigned long long count = 0;

	//		// Only draw if the ellipsoid is sane!
	//		if(!isinf(l1) && !isnan(l1) && l1 != 0 &&
	//			!isinf(l2) && !isnan(l2) && l2 != 0 &&
	//			!isinf(l3) && !isnan(l3) && l3 != 0)
	//			count = drawEllipsoid<T, Tmask>(image, pos, l1, l2, l3, phi1, theta1, phi2, theta2, phi3, theta3, ellipsoidColor, pMask);

	//		if(pResults != 0)
	//			pResults->push_back(count);

	//		utilities::showProgress(n, results.size());
	//	}
	//}

	
	///**
	//Draws bounding spheres.
	//@param pMask Pointer to mask image. Only pixels where mask != 0 will be colored.
	//@param pResults Pointer to results array. If not null, contains count of filled pixels for each bounding sphere.
	//*/
	//template<typename T, typename Tmask> void drawBoundingSpheres(Image<T>& image, const vector<string>& headers, const vector<vector<double> >& results, T sphereColor,
	//												Image<Tmask>* pMask = 0, vector<unsigned long long>* pResults = 0)
	//{
	//	//int cxi = 25, int cyi = 26, int czi = 27, int radiusi = 28
	//	size_t cxi = getColumnIndex("bcenter X [pixel]", headers);
	//	size_t cyi = getColumnIndex("bcenter Y [pixel]", headers);
	//	size_t czi = getColumnIndex("bcenter Z [pixel]", headers);
	//	size_t radiusi = getColumnIndex("bradius [pixel]", headers);

	//	Vec3d pos;
	//	Sphere<double> sphere;
	//	for(size_t n = 0; n < results.size(); n++)
	//	{
	//		pos.x = results[n][cxi];
	//		pos.y = results[n][cyi];
	//		pos.z = results[n][czi];

	//		double radius = results[n][radiusi];
	//		
	//		sphere.center = pos;
	//		sphere.radius = radius;

	//		unsigned long long count = drawSphere<T, double, Tmask>(image, sphere, sphereColor, pMask);
	//		if(pResults != 0)
	//			pResults->push_back(count);

	//		utilities::showProgress(n, results.size());
	//	}
	//}

	namespace tests
	{
		void analyzeParticlesThreading();
		void analyzeParticlesThreadingBig();
		void analyzeParticlesSanity();
		void analyzeParticlesSanity2();
		void analyzeParticlesVolumeLimit();
	}
}

