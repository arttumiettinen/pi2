//
//#include "thickmap_experiments.h"
//
//#include "thickmap.h"
//
//#include "datatypes.h"
//#include "math/mathutils.h"
//#include "math/vec3.h"
//#include "image.h"
//#include "pointprocess.h"
//#include "generation.h"
//#include "dmap.h"
//#include "conversions.h"
//#include "testutils.h"
//#include "projections.h"
//#include "danielsson.h"
//#include "aabox.h"
//#include "transform.h"
//
//#include <vector>
//#include <algorithm>
//#include <queue>
//#include <functional>
//#include <set>
//
//#if defined(_WIN32)
//#include <execution>
//#endif
//
//namespace itl2
//{
//
//
//	namespace dimredsmartblocks
//	{
//		
//
//		/*
//		@param centers Image containing only the row to be processed before first call to this method.
//		@param ri Full ri image that will be updated.
//		@param start Start point of the row to be processed.
//		@param dim Dimension that we are processing.
//		@param step +1 or -1 to indicate the direction of the pass.
//		*/
//		inline void singlePass(const Image<internals::RiSet>& centers, Image<internals::RiSet>& ri, const Vec3c& rowStart, coord_t dim, coord_t step, vector<internals::ActiveSpheresItem>& activeSpheres)
//		{
//			// Stores the spheres that have been encountered and that have not been passed yet.
//			// Stores the center point, original radius, and ri.
//			//vector<internals::ActiveSpheresItem> activeSpheres;
//			vector<internals::ActiveSpheresItem> activeSpheresTmp;
//			vector<internals::ActiveSpheresItem> Ctmp;
//
//			internals::RiSet resultTmp;
//			internals::RiSet resultTmp2;
//
//			//activeSpheres.reserve(40);
//			activeSpheresTmp.reserve(40);
//			Ctmp.reserve(40);
//			resultTmp.reserve(40);
//			resultTmp2.reserve(40);
//
//			// Set start point to the start or to the end of the current row.
//			Vec3c p = rowStart;
//			if (step < 0)
//				p[dim] += (ri.dimension(dim) - 1);
//
//			coord_t dimensionality = ri.dimensionality();
//
//			for (coord_t i = 0; i < ri.dimension(dim); i++, p[dim] += step)
//			{
//				int32_t x = (int32_t)p[dim];
//
//				// If there is one or more sphere centers at the current location, add them to the set of active spheres.
//				const internals::RiSet& C = centers(x);
//
//				// The C list is sorted by R and so is activeSpheres list.
//				// Instead of finding place for each item, convert each item in C list to activeSpheres format (add x coordinate)
//				// and then merge the two sorted lists to construct the new activeSpheres list.
//				Ctmp.clear();
//				for (auto& item : C)
//					Ctmp.push_back({ item.R2, item.ri2, x });
//
//				activeSpheresTmp.clear();
//				activeSpheresTmp.resize(C.size() + activeSpheres.size());
//				merge(Ctmp.begin(), Ctmp.end(), activeSpheres.begin(), activeSpheres.end(), activeSpheresTmp.begin(), internals::activeSpheresSorter);
//				swap(activeSpheres, activeSpheresTmp);
//
//
//				// Iterate through all active spheres and calculate radius for next dimension.
//				resultTmp.clear();
//				auto it = activeSpheres.begin();
//				while (it != activeSpheres.end())
//				{
//					int32_t Rorig2 = it->R2;
//					int32_t R2 = it->ri2;
//					int32_t cx = it->xi;
//
//					int32_t dx = abs(x - cx);
//
//					// Calculate new ri^2
//					int32_t rn2 = R2 - dx * dx;
//					if (rn2 > 0)
//					{
//						// Insert ry2 into the list (TODO: Don't insert if the item is going to be removed anyway?)
//						resultTmp.push_back({ Rorig2, rn2 });
//
//						// Remove possible duplicate Rorig and replace it by (Rorig, max(duplicates))
//						if (resultTmp.size() > 1 && Rorig2 == resultTmp[resultTmp.size() - 2].R2)
//						{
//							if (rn2 > resultTmp[resultTmp.size() - 2].ri2)
//								resultTmp[resultTmp.size() - 2].ri2 = rn2;
//
//							resultTmp.pop_back();
//						}
//
//						it++;
//					}
//					else
//					{
//						// ry is non-positive, i.e. dx >= R
//						// This sphere is not active anymore, so remove it from the list of active spheres.
//						it = activeSpheres.erase(it);
//					}
//				}
//
//				// Rebuild ri list from resultTmp. Don't include those items that are hidden by other items.
//				if (resultTmp.size() > 0)
//				{
//					// Add ri from the previous pass to the ri list and save the result to resultTmp2.
//					// Note that both resultTmp2 and rilist are sorted so we can just merge them.
//					internals::RiSet& rilist = ri(p);
//					resultTmp2.clear();
//					resultTmp2.resize(rilist.size() + resultTmp.size());
//					merge(rilist.begin(), rilist.end(), resultTmp.begin(), resultTmp.end(), resultTmp2.begin(), internals::origRSorter);
//
//
//					// NOTE: This is the basic version, the memory saving version is below.
//					//// Linear time algorithm for finding relevant ri (those not hidden by other items).
//					//// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//					//// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//					//rilist.clear();
//					//rilist.push_back(resultTmp2[0]);
//					//for (size_t n = 1; n < resultTmp2.size(); n++)
//					//{
//					//	int32_t currri2 = rilist[rilist.size() - 1].ri2;
//					//	int32_t newri2 = resultTmp2[n].ri2;
//					//	if (newri2 > currri2)
//					//	{
//					//		rilist.push_back(resultTmp2[n]);
//					//	}
//					//}
//
//					// Linear time algorithm for finding relevant ri (those not hidden by other items).
//					// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//					// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//					rilist.clear();
//					rilist.push_back(resultTmp2[0]);
//					//if (dim < ri.dimensionality() - 1) // In the last dimension only the item with maximum R is needed, but we have separate function for that
//					//{
//					for (size_t n = 1; n < resultTmp2.size(); n++)
//					{
//						int32_t currri2 = rilist[rilist.size() - 1].ri2;
//						int32_t newri2 = resultTmp2[n].ri2;
//
//						if (newri2 > currri2) // This is the basic condition that works always (but may include unnecessary items in the rilist)
//						{
//							if (dim == dimensionality - 2)
//							{
//								// In the second last dimension only really visible spans are needed as there's no next dimension whose ri we would calculate based on the spans.
//								if (largestIntWhoseSquareIsLessThan(newri2) > largestIntWhoseSquareIsLessThan(currri2))
//									rilist.push_back(resultTmp2[n]);
//							}
//							else if (dim == dimensionality - 3)
//							{
//								// In the third last dimension we know that only those spans are required that produce visible circles in the output.
//								if (!internals::doesDiscretizedCircle1FitInto2Cached(newri2, currri2))
//									rilist.push_back(resultTmp2[n]);
//							}
//							else
//							{
//								// Here we could insert test if discretized spheres fit into each other etc. etc.
//								rilist.push_back(resultTmp2[n]);
//							}
//						}
//					}
//					//}
//
//					// Make sure we don't use extraneous memory by storing empty items in each pixel of the ri image.
//					rilist.shrink_to_fit();
//				}
//			}
//		}
//
//		/*
//		Same than singlePass but generates result image instead of ri image, and is thus a little bit faster.
//		@param centers Image containing only the row to be processed before first call to this method.
//		@param Result image. Must be set to zero before first call to this method.
//		@param start Start point of the row to be processed.
//		@param dim Dimension that we are processing.
//		@param step +1 or -1 to indicate the direction of the pass.
//		*/
//		inline void singlePassFinal(const Image<internals::RiSet>& centers, Image<int32_t>& result, const Vec3c& rowStart, size_t dim, coord_t step, vector<internals::ActiveSpheresItem>& initialActiveSpheres)
//		{
//			// Stores the spheres that have been encountered and that have not been passed yet.
//			// Stores the sphere with the largest R at the top of the priority queue.
//			priority_queue<internals::ActiveSpheresItem, vector<internals::ActiveSpheresItem>, std::function<decltype(internals::activeSpheresSorter2)> > activeSpheres(internals::activeSpheresSorter2);
//
//			for (const auto& item : initialActiveSpheres)
//				activeSpheres.push(item);
//
//
//			// Set start point to the start or end of the current row.
//			Vec3c p = rowStart;
//			if (step < 0)
//				p[dim] += (result.dimension(dim) - 1);
//
//			for (coord_t i = 0; i < result.dimension(dim); i++, p[dim] += step)
//			{
//				int32_t x = (int32_t)p[dim];
//
//				// If there is one or more sphere centers at the current location, add them to the set of active spheres.
//				const internals::RiSet& C = centers(x);
//				for (auto& item : C)
//				{
//					activeSpheres.push({ item.R2, item.ri2, x });
//				}
//
//				while (!activeSpheres.empty())
//				{
//					auto& item = activeSpheres.top();
//
//					int32_t Rorig2 = item.R2;
//					int32_t R2 = item.ri2;
//					int32_t cx = item.xi;
//					int32_t dx = abs(x - cx);
//
//					// Calculate new ri^2
//					int32_t rn2 = R2 - dx * dx;
//					if (rn2 > 0)
//					{
//						// Note that previous pass may have assigned larger value to the output.
//						if (Rorig2 > result(p))
//							result(p) = Rorig2;
//						break;
//					}
//					else
//					{
//						// ry is non-positive, i.e. dx >= R
//						// This sphere is not active anymore, so remove it from the list of active spheres.
//						activeSpheres.pop();
//					}
//
//				}
//
//			}
//
//			// Put stuff back to initialActiveSpheres list as it will be the initial active spheres list for possible adjacent block.
//			initialActiveSpheres.clear();
//			while (!activeSpheres.empty())
//			{
//				initialActiveSpheres.push_back(activeSpheres.top());
//				activeSpheres.pop();
//			}
//			reverse(initialActiveSpheres.begin(), initialActiveSpheres.end());
//			// TODO: How to avoid sorting here?
//			sort(initialActiveSpheres.begin(), initialActiveSpheres.end(), internals::activeSpheresSorter);
//			//if (!is_sorted(initialActiveSpheres.begin(), initialActiveSpheres.end(), internals::activeSpheresSorter))
//			//{
//			//	cout << "Not sorted" << endl;
//			//	sort(initialActiveSpheres.begin(), initialActiveSpheres.end(), internals::activeSpheresSorter);
//			//}
//		}
//
//		inline bool ridgeTestComparerri(const internals::RiItem& a, const internals::RiItem& b)
//		{
//			if (a.ri2 != b.ri2)
//				return a.ri2 < b.ri2;
//			return a.R2 < b.R2;
//		}
//
//		///*
//		//Makes forward and backward pass in single dimension.
//		//@param ri Image containing ri values from processing of last dimension.
//		//@param dim Dimension to process.
//		//@param dir 1 to make forward pass, -1 to make backward pass.
//		//@param result Result image, changed only if this is the final pass (processing the last dimension).
//		//*/
//		//inline void processDimension(Image<internals::RiSet>& ri, size_t dim, coord_t dir, Image<int32_t>& result, bool showProgressInfo)
//		//{
//		//	// Determine count of pixels to process
//		//	Vec3c reducedDimensions = ri.dimensions();
//		//	reducedDimensions[dim] = 1;
//		//	coord_t rowCount = reducedDimensions.x * reducedDimensions.y * reducedDimensions.z;
//
//
//
//		//	size_t counter = 0;
//		//	#pragma omp parallel if(!omp_in_parallel() && result.pixelCount() > PARALLELIZATION_THRESHOLD)
//		//	{
//		//		// Temporary buffer
//		//		Image<internals::RiSet> centers(ri.dimension(dim));
//
//		//		// Determine start points of pixel rows and process each row
//		//		#pragma omp for schedule(dynamic)
//		//		for (coord_t n = 0; n < rowCount; n++)
//		//		{
//		//			Vec3c start = indexToCoords(n, reducedDimensions);
//
//		//			// Make a copy of the current row as we update the row in the first pass but need
//		//			// the original data in the second pass.
//		//			Vec3c pos = start;
//		//			for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
//		//			{
//		//				centers(x) = ri(pos);
//		//				ri(pos).clear();
//		//			}
//
//		//			if (dim < ri.dimensionality() - 1)
//		//			{
//		//				singlePass(centers, ri, start, dim, dir);
//
//		//				//// Left to right pass
//		//				//singlePass(centers, ri, start, dim, 1);
//
//		//				//// Right to left pass
//		//				//singlePass(centers, ri, start, dim, -1);
//		//			}
//		//			else
//		//			{
//		//				singlePassFinal(centers, result, start, dim, dir);
//
//		//				//// TODO: Test if it is better/faster to combine singlePass and singlePassFinal as in the ridge version
//
//		//				//// Left to right pass
//		//				//singlePassFinal(centers, result, start, dim, 1);
//
//		//				//// Right to left pass
//		//				//singlePassFinal(centers, result, start, dim, -1);
//		//			}
//
//		//			showThreadProgress(counter, rowCount, showProgressInfo);
//		//		}
//		//	}
//		//}
//
//		Vec3c getReducedDimensions(Vec3c size, size_t dim)
//		{
//			size[dim] = 1;
//			return size;
//		}
//
//		/*
//		Makes one pass over image in specific dimension and direction.
//		@param ri Image containing ri values from processing of last dimension.
//		@param dim Dimension to process.
//		@param dir 1 to make forward pass, -1 to make backward pass.
//		@param result Result image, changed only if this is the final pass (processing the last dimension).
//		*/
//		inline void processDimension(Image<internals::RiSet>& ri, size_t dim, coord_t dir, Image<vector<internals::ActiveSpheresItem> >& activeSpheres, Image<int32_t>& result, bool showProgressInfo)
//		{
//			// Determine count of pixels to process
//			Vec3c reducedDimensions = getReducedDimensions(ri.dimensions(), dim);
//			coord_t rowCount = reducedDimensions.x * reducedDimensions.y * reducedDimensions.z;
//
//			activeSpheres.checkSize(reducedDimensions);
//
//			size_t counter = 0;
//#pragma omp parallel if(!omp_in_parallel() && result.pixelCount() > PARALLELIZATION_THRESHOLD)
//			{
//				// Temporary buffer
//				Image<internals::RiSet> centers(ri.dimension(dim));
//
//				// Determine start points of pixel rows and process each row
//#pragma omp for schedule(dynamic)
//				for (coord_t n = 0; n < rowCount; n++)
//				{
//					Vec3c start = indexToCoords(n, reducedDimensions);
//
//					// Make a copy of the current row as we update the row in the first pass but need
//					// the original data in the second pass.
//					Vec3c pos = start;
//					for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
//					{
//						centers(x) = ri(pos);
//						ri(pos).clear();
//					}
//
//					vector<internals::ActiveSpheresItem>& as = activeSpheres(start);
//
//					if (dim < ri.dimensionality() - 1)
//					{
//						singlePass(centers, ri, start, dim, dir, as);
//					}
//					else
//					{
//						singlePassFinal(centers, result, start, dim, dir, as);
//					}
//
//					showThreadProgress(counter, rowCount, showProgressInfo);
//				}
//			}
//		}
//
//		/*
//		Adds x to all xi values in all items of all lists in the given active spheres image.
//		*/
//		void addxi(Image<vector<internals::ActiveSpheresItem> >& activeSpheres, int32_t x)
//		{
//			for (coord_t n = 0; n < activeSpheres.pixelCount(); n++)
//			{
//				vector<internals::ActiveSpheresItem>& v = activeSpheres(n);
//				for (internals::ActiveSpheresItem& item : v)
//				{
//					item.xi += x;
//				}
//			}
//		}
//
//
//		/*
//		Copies squared sphere radius data from original centers image to ri image.
//		*/
//		inline void prepare(const Image<int32_t>& centers2, Image<internals::RiSet>& ri, const AABox<coord_t>& bounds)
//		{
//			ri.ensureSize(bounds.size());
//
//			for (coord_t z = 0; z < ri.depth(); z++)
//			{
//				for (coord_t y = 0; y < ri.height(); y++)
//				{
//					for (coord_t x = 0; x < ri.width(); x++)
//					{
//						Vec3c p(x, y, z);
//						ri(p).clear();
//						int32_t R2 = centers2(bounds.minc + p);
//						if (R2 > 0)
//							ri(p).push_back({ R2, R2 });
//					}
//				}
//			}
//			//for (coord_t n = 0; n < centers2.pixelCount(); n++)
//			//{
//			//	ri(n).clear();
//			//	int32_t R2 = centers2(n);
//			//	if (R2 > 0)
//			//	{
//			//		ri(n).push_back({ R2, R2 });
//			//	}
//			//}
//		}
//
//
//
//		struct Block
//		{
//			/*
//			Bounds of this block in image coordinates.
//			*/
//			AABox<coord_t> bounds;
//
//			/*
//			Initial active spheres list for each pixel row in each coordinate direction.
//			*/
//			array< Image<vector<internals::ActiveSpheresItem> >, 3 > activeSpheres;
//
//			/**
//			Checks if active spheres are ready up to given dimensionality.
//			*/
//			bool activeSpheresReady(size_t dimensionality) const
//			{
//				for (size_t n = 0; n < dimensionality; n++)
//					if (activeSpheres[n].pixelCount() <= 1)
//						return false;
//
//				return true;
//			}
//
//			/*
//			Resets active spheres images to not available-state.
//			*/
//			void resetActiveSpheres()
//			{
//				for (size_t n = 0; n < activeSpheres.size(); n++)
//					activeSpheres[n].ensureSize(1, 1, 1);
//			}
//
//			Block()
//			{
//
//			}
//
//			Block(const AABox<coord_t>& box) :
//				bounds(box)
//			{
//			}
//
//			Block(const Block& block) :
//				bounds(block.bounds)
//			{
//				for (size_t n = 0; n < activeSpheres.size(); n++)
//					setValue(activeSpheres[n], block.activeSpheres[n]);
//			}
//
//			Block &operator =(const Block & r)
//			{
//				bounds = r.bounds;
//				for (size_t n = 0; n < activeSpheres.size(); n++)
//					setValue(activeSpheres[n], r.activeSpheres[n]);
//				return *this;
//			}
//		};
//
//		void subdivide(const Vec3c& dmap2Dimensions, const Vec3c& blockSize, Image<Block>& blocks)
//		{
//			// Create blocks
//			Vec3c blockCount = ceil(Vec3d(dmap2Dimensions).componentwiseDivide(Vec3d(blockSize)));
//			blocks.ensureSize(blockCount);
//			for (coord_t z = 0; z < blockCount.z; z++)
//			{
//				for (coord_t y = 0; y < blockCount.y; y++)
//				{
//					for (coord_t x = 0; x < blockCount.x; x++)
//					{
//						Vec3c pos = blockSize.componentwiseMultiply(Vec3c(x, y, z));
//						Vec3c right = pos + blockSize;
//						clamp(right, Vec3c(0, 0, 0), dmap2Dimensions);
//						blocks(x, y, z).bounds = AABox<coord_t>(pos, right);
//					}
//				}
//			}
//		}
//
//		void thickmap2Blocks(Image<int32_t>& dmap2, const Vec3c& blockSize, double* extraBytes, bool showProgressInfo)
//		{
//			showProgressInfo = false;
//
//			internals::buildCircleLookup(max(dmap2));
//
//			Image<int32_t> result(dmap2.dimensions());
//
//			Image<Block> blocks;
//			subdivide(dmap2.dimensions(), blockSize, blocks);
//
//			// Find out which direction sets we must process
//			std::vector<Vec3c> directions;
//			if (dmap2.dimensionality() == 3)
//			{
//				directions = {
//					Vec3c(1, 1, 1),
//					Vec3c(1, 1, -1),
//					Vec3c(1, -1, 1),
//					Vec3c(1, -1, -1),
//					Vec3c(-1, 1, 1),
//					Vec3c(-1, 1, -1),
//					Vec3c(-1, -1, 1),
//					Vec3c(-1, -1, -1)
//				};
//			}
//			else if (dmap2.dimensionality() <= 2)
//			{
//				directions = {
//					Vec3c(1, 1, 0),
//					Vec3c(1, -1, 0),
//					Vec3c(-1, 1, 0),
//					Vec3c(-1, -1, 0)
//				};
//			}
//			else
//			{
//				// Enable dimensionality-independent determination of processing directions here
//				// if the system supports more than 3 dimensions later.
//				throw ITLException("Unsupported dimensionality.");
//			}
//
//			Image<internals::RiSet> ri;
//
//			// Process each direction separately, updating result image as we go.
//			for (coord_t dn = 0; dn < (coord_t)directions.size(); dn++)
//			{
//				Vec3c dirs = directions[dn];
//
//				// Initialize activeSpheres images to empty for those blocks that are in the edge of the image
//				// in processing direction, and undefined to all other blocks.
//				for (coord_t m = 0; m < (coord_t)blocks.pixelCount(); m++)
//				{
//					Block& b = blocks(m);
//
//					// First reset the images to "not available" state.
//					for (size_t dim = 0; dim < dmap2.dimensionality(); dim++)
//					{
//						if ((b.bounds.minc[dim] <= 0 && dirs[dim] == 1) ||
//							(b.bounds.maxc[dim] >= dmap2.dimension(dim) && dirs[dim] == -1))
//							b.activeSpheres[dim].ensureSize(getReducedDimensions(b.bounds.size(), dim));
//						else
//							b.activeSpheres[dim].ensureSize(1, 1, 1);
//					}
//				}
//
//				// Process blocks that have all activeSpheres images available.
//				bool changed;
//				do
//				{
//					changed = false;
//					for (coord_t z = 0; z < blocks.depth(); z++)
//					{
//						for (coord_t y = 0; y < blocks.height(); y++)
//						{
//							for (coord_t x = 0; x < blocks.width(); x++)
//							{
//								Block& b = blocks(x, y, z);
//
//								if (b.activeSpheresReady(dmap2.dimensionality()))
//								{
//									//cout << "Processing block " << Vec3c(x, y, z) << endl;
//
//									changed = true;
//
//									prepare(dmap2, ri, b.bounds);
//
//									Image<int32_t> resultBlock(b.bounds.size());
//									crop(result, resultBlock, b.bounds.position());
//
//									for (size_t n = 0; n < dirs.size(); n++)
//									{
//										if (dirs[n] == 0)
//											break;
//
//										processDimension(ri, n, dirs[n], b.activeSpheres[n], resultBlock, showProgressInfo);
//
//										// The final activeSpheres image of this block is initial activeSpheres image of the next block
//										// in the processing direction.
//										Vec3c nextBlockPos(x, y, z);
//										nextBlockPos[n] += dirs[n];
//
//										if (blocks.isInImage(nextBlockPos))
//										{
//											Block& bn = blocks(nextBlockPos);
//
//											if (bn.activeSpheres[n].pixelCount() > 1)
//												throw logic_error("Initial active spheres list of a block is already set.");
//
//											// Convert active sphere positions to the coordinates of the next block
//											if (dirs[n] > 0)
//												addxi(b.activeSpheres[n], -(int32_t)b.bounds.size()[n]);
//											else
//												addxi(b.activeSpheres[n], (int32_t)bn.bounds.size()[n]);
//
//											setValue(bn.activeSpheres[n], b.activeSpheres[n]);
//										}
//									}
//
//									setValue(result, resultBlock, b.bounds.position());
//
//									// Reset active spheres of this block to not ready state.
//									b.resetActiveSpheres();
//								}
//							}
//						}
//					}
//				} while (changed);
//
//			}
//
//			setValue(dmap2, result);
//
//		}
//
//	}
//
//
//
//	namespace propagation
//	{
//		struct SphereItem
//		{
//			Vec3sc center;
//			int32_t r2;
//		};
//
//		class Comparer
//		{
//		public:
//
//			bool operator()(const SphereItem& l, const SphereItem& r)
//			{
//				if (l.r2 != r.r2)
//				{
//					return l.r2 < r.r2;
//				}
//				else if (l.center.z != r.center.z)
//				{
//					return l.center.z < r.center.z;
//				}
//				else if (l.center.y != r.center.y)
//				{
//					return l.center.y < r.center.y;
//				}
//				else if (l.center.x != r.center.x)
//				{
//					return l.center.x < r.center.x;
//				}
//				else
//				{
//					return false;
//				}
//			}
//		};
//
//		//typedef set<SphereItem, Comparer> SphereSet;
//		typedef vector<SphereItem> SphereSet;
//
//		namespace internals
//		{
//			/**
//			Tests if sphere contains a point.
//			*/
//			bool contains(const SphereItem& sphere, const Vec3sc& p)
//			{
//				return (sphere.center - p).normSquared<int32_t>() < sphere.r2;
//			}
//
//			//void addSpheres(SphereSet& spheres, const Image<SphereSet>& activeSpheres, const Vec3sc& p, const Vec3sc& dir)
//			//{
//			//	if (activeSpheres.isInImage(p + dir))
//			//	{
//			//		const SphereSet& v = activeSpheres(p + dir);
//			//		for (const SphereItem& sphere : v)
//			//		{
//			//			if (contains(sphere, p))
//			//				spheres.insert(sphere);
//			//		}
//			//	}
//			//}
//
//			void addSpheres(SphereSet& dest, const Image<SphereSet>& activeSpheres, const Vec3sc& p, const Vec3sc& dir, SphereSet& temp)
//			{
//				Vec3sc pn = p + dir;
//				pn.z = 0;
//				if (activeSpheres.isInImage(pn))
//				{
//					const SphereSet& v = activeSpheres(pn);
//					set_union(v.begin(), v.end(), dest.begin(), dest.end(), std::back_inserter(temp), Comparer());
//					dest.swap(temp);
//					temp.clear();
//				}
//			}
//		}
//
//
//		class DirPPP
//		{
//		public:
//			Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims) const
//			{
//				return Vec3sc(x, y, z);
//			}
//
//			Vec3sc backNeighbourDir() const
//			{
//				return Vec3sc(-1, -1, -1);
//			}
//		};
//
//		class DirNPP
//		{
//		public:
//			Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims)
//			{
//				return Vec3sc((int32_t)dims.x - 1 - x, y, z);
//			}
//
//			Vec3sc backNeighbourDir() const
//			{
//				return Vec3sc(1, -1, -1);
//			}
//		};
//
//
//		/**
//		Iterates image diagonally.
//		@param dims Dimensions of the image.
//		@param dir Functor that transforms initial image coordinates such that the iteration direction is the desired one. Select one from DirXXX classes.
//		@param process Functor that is called for each iterated pixel. The signature must be process(const Vec3sc& position, const Vec3sc& backwardsDirection, int32_t iterationLevel). The function will be called from multiple threads simultaneously.
//		@param levelComplete Functor that is called after a diagonal is processed. The signature must be levelComplete(int32_t currLevel).
//		*/
//		template<typename dirfunc, typename processfunc, typename levelcompletefunc> void diagonalIteration(const Vec3c& dims, dirfunc dir, processfunc process, levelcompletefunc levelComplete, bool showProgressInfo = true)
//		{
//			int32_t w = (int32_t)dims.x;
//			int32_t h = (int32_t)dims.y;
//			int32_t d = (int32_t)dims.z;
//
//			int32_t maxLevel = w + h + d - 2;
//			size_t counter = 0;
//			for (int32_t start = 0; start < maxLevel; start++)
//			{
//#pragma omp parallel for if(!omp_in_parallel())
//				for (int32_t z = 0; z < d; z++)
//				{
//					int32_t x = start - z;
//					int32_t y = 0;
//
//					if (x >= w)
//					{
//						y += x - (w - 1);
//						x = (w - 1);
//					}
//
//					while (x >= 0 && y < h)
//					{
//						//if (!result.isInImage(x, y, z))
//						//	throw logic_error("Invalid program.");
//
//						process(dir(x, y, z, dims), dir.backNeighbourDir(), start);
//
//						x--;
//						y++;
//					}
//				}
//
//				levelComplete(start);
//
//				showThreadProgress(counter, maxLevel, showProgressInfo);
//			}
//		}
//
//		void thickmap2(Image<int32_t>& dmap2, double* extraBytes, bool showProgressInfo)
//		{
//			Image<int32_t> result(dmap2.dimensions());
//
//			Image<SphereSet> activeSpheres1(dmap2.width(), dmap2.height(), 1);
//			Image<SphereSet> activeSpheres2(dmap2.width(), dmap2.height(), 1);
//
//			Image<SphereSet>* currentAS = &activeSpheres1;
//			Image<SphereSet>* prevAS = &activeSpheres2;
//
//			// Allocate initial amount of memory
//			const size_t initialCapacity = 200;
//			vector<SphereSet> tempSets((size_t)omp_get_max_threads());
//			for (size_t n = 0; n < tempSets.size(); n++)
//			{
//				tempSets[n].reserve(initialCapacity);
//			}
//
//			for (coord_t n = 0; n < activeSpheres1.pixelCount(); n++)
//			{
//				activeSpheres1(n).reserve(initialCapacity);
//				activeSpheres2(n).reserve(initialCapacity);
//			}
//
//			size_t maxSize = 0;
//
//			diagonalIteration(result.dimensions(), DirPPP(), [&](const Vec3sc& p, const Vec3sc& backDir, int32_t level)
//			{
//				//result(p) = level;
//
//				SphereSet& tempSet = tempSets[omp_get_thread_num()];
//
//				SphereSet& as = (*currentAS)(p.x, p.y, 0);
//
//				// Get active spheres from previous neighbour points
//				//internals::addSpheres(as, *prevAS, p, Vec3sc(backDir.x, 0, 0));
//				//internals::addSpheres(as, *prevAS, p, Vec3sc(0, backDir.y, 0));
//				//internals::addSpheres(as, *prevAS, p, Vec3sc(0, 0, backDir.z));
//				internals::addSpheres(as, *prevAS, p, Vec3sc(backDir.x, 0, 0), tempSet);
//				internals::addSpheres(as, *prevAS, p, Vec3sc(0, backDir.y, 0), tempSet);
//				internals::addSpheres(as, *prevAS, p, Vec3sc(0, 0, backDir.z), tempSet);
//
//				// Remove unnecessary spheres
//				for (SphereItem& sphere : as)
//				{
//					if (internals::contains(sphere, p))
//						tempSet.push_back(sphere); // The order of the spheres does not change so no need to compare anything.
//				}
//				as.swap(tempSet);
//				tempSet.clear();
//
//
//				// Add possible current sphere
//				if (dmap2(p) != 0)
//				{
//					//as.insert({ p, dmap2(p) });
//					SphereItem item{ p, dmap2(p) };
//					as.insert(
//						std::upper_bound(as.begin(), as.end(), item, Comparer()),
//						item
//					);
//				}
//
//				// Choose max
//				if (as.size() > 0)
//					result(p) = std::max(result(p), as.rbegin()->r2);
//
//
//				maxSize = std::max(maxSize, as.size());
//			},
//				[&](int32_t level)
//			{
//				swap(currentAS, prevAS);
//				for (coord_t n = 0; n < currentAS->pixelCount(); n++)
//					if ((*currentAS)(n).size() > 0)
//						(*currentAS)(n).clear();
//			});
//
//
//			cout << "Max activeSpheres size = " << maxSize << endl;
//			//// Iterate diagonally:
//			//Image<SphereSet> activeSpheres(dmap2.dimensions());
//			//// x + y + z = n, where n is index of current iteration
//			//Vec3sc p0(0, 0, 0);
//			//int32_t maxn = (int32_t)dmap2.dimensions().sum() - 3;
//			//for (int32_t n = 0; n < maxn; n++)
//			//{
//			//	for (int32_t z = 0; z < dmap2.depth(); z++)
//			//	{
//			//		for (int32_t y = 0; y < dmap2.height(); y++)
//			//		{
//			//			for (int32_t x = 0; x < dmap2.width(); x++)
//			//			{
//			//				Vec3sc p(x, y, z);
//
//			//				result(p) = 0;
//
//			//				//SphereSet& as = activeSpheres(p);
//
//			//				//int32_t d = (p - p0).sum();
//
//			//				//if (d == n)
//			//				//{
//			//				//	// Get active spheres from previous neighbour points
//			//				//	internals::addSpheres(as, activeSpheres, p, -Vec3sc(1, 0, 0));
//			//				//	internals::addSpheres(as, activeSpheres, p, -Vec3sc(0, 1, 0));
//			//				//	internals::addSpheres(as, activeSpheres, p, -Vec3sc(0, 0, 1));
//			//				//	internals::addSpheres(as, activeSpheres, p, -Vec3sc(1, 1, 0));
//			//				//	internals::addSpheres(as, activeSpheres, p, -Vec3sc(0, 1, 1));
//			//				//	internals::addSpheres(as, activeSpheres, p, -Vec3sc(1, 0, 1));
//			//				//	internals::addSpheres(as, activeSpheres, p, -Vec3sc(1, 1, 1));
//
//			//				//	// Add possible current sphere
//			//				//	if (dmap2(p) != 0)
//			//				//	{
//			//				//		as.insert({ p, dmap2(p) });
//			//				//	}
//
//			//				//	// Choose max
//			//				//	if (as.size() > 0)
//			//				//		//result(p) = max_element(as.begin(), as.end(), [](const SphereItem& a, const SphereItem& b) { return a.r2 > b.r2; })->r2;
//			//				//		//result(p) = as.begin()->r2;
//			//				//		result(p) = as.rbegin()->r2;
//			//				//	else
//			//				//		result(p) = 0;
//			//				//}
//			//				//else if (d == n - 4)
//			//				//{
//			//				//	// Clear active spheres as this layer is not required anymore.
//			//				//	as.clear();
//			//				//	//as.shrink_to_fit();
//			//				//}
//
//			//				//if (result(p) == 0)
//			//				//{
//			//				//	if((p - p0).sum() == n)
//			//				//		result(p) = n;
//			//				//}
//			//			}
//			//		}
//			//	}
//
//			//	showProgress(n, maxn);
//			//}
//
//			setValue(dmap2, result);
//		}
//	}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//	namespace dimredblocks
//	{
//
//		struct Block
//		{
//			AABox<int32_t> bounds;
//
//			/**
//			Contains spheres that the block does not contain but that may touch the block.
//			*/
//			vector<Sphere2> spheres;
//
//			Block()
//			{
//			}
//
//			Block(const AABox<int32_t> box) : bounds(box)
//			{
//			}
//		};
//
//		/*
//		Checks which blocks do not contain sphere at (x, y, z) with radius sqrt(r2), but that the sphere touches.
//		Adds the sphere to spheres list of those blocks.
//		*/
//		void processSphere(const Vec3sc& pos, int32_t r2, vector<Block>& blocks)
//		{
//			for (Block& b : blocks)
//			{
//				if (!b.bounds.contains(pos))
//				{
//					int32_t dist2 = b.bounds.distanceSquared(pos);
//					if (dist2 <= r2 + 1) // NOTE: +1.0 so that floating point inaccuracy does not affect results.
//					{
//						// Convert position to block coordinates by subtracting position of the block.
//						b.spheres.push_back({ pos - b.bounds.minc, r2 });
//					}
//				}
//			}
//		}
//
//		/**
//		Sub-divides dmap2 into blocks that can be processed separately.
//		*/
//		vector<Block> subdivide(const Image<int32_t>& dmap2, const Vec3sc& blockSize)
//		{
//			// Create blocks
//			Vec3sc blockCount = Vec3sc(ceil(Vec3d(dmap2.dimensions()).componentwiseDivide(Vec3d(blockSize))));
//			vector<Block> blocks;
//			for (int32_t z = 0; z < blockCount.z; z++)
//			{
//				for (int32_t y = 0; y < blockCount.y; y++)
//				{
//					for (int32_t x = 0; x < blockCount.x; x++)
//					{
//						Vec3sc pos = blockSize.componentwiseMultiply(Vec3sc(x, y, z));
//						Vec3sc right = pos + blockSize;
//						clamp(right, Vec3sc(0, 0, 0), Vec3sc(dmap2.dimensions()));
//						blocks.push_back(Block(AABox<int32_t>(pos, right)));
//					}
//				}
//			}
//
//			// Find spheres that touch each block
//			size_t counter = 0;
//			// TODO: Test if enabling parallelization and adding locks to processSphere increases performance. Earlier it did not...
//			//#pragma omp parallel if(!omp_in_parallel() && result.pixelCount() > PARALLELIZATION_THRESHOLD)
//			for (int32_t z = 0; z < dmap2.depth(); z++)
//			{
//				for (int32_t y = 0; y < dmap2.height(); y++)
//				{
//					for (int32_t x = 0; x < dmap2.width(); x++)
//					{
//						int32_t r2 = dmap2(x, y, z);
//						if (r2 > 0)
//						{
//							processSphere(Vec3sc(x, y, z), r2, blocks);
//						}
//					}
//				}
//
//				showThreadProgress(counter, dmap2.depth());
//			}
//
//			//cout << blocks[0].bounds.minc << endl;
//			//cout << blocks[0].bounds.maxc << endl;
//			//cout << AABox<int32_t>(blocks[0].bounds.minc, blocks[0].bounds.maxc - Vec3sc(1, 1, 1)).distanceSquared(Vec3sc(199, 80, 0)) << endl;
//			//for (const auto& s : blocks[0].spheres)
//			//{
//			//	if (s.pos.x == 199 && s.r2 == 324)
//			//	{
//			//		cout << "***************ok****************" << endl;
//			//	}
//			//}
//
//			return blocks;
//		}
//
//
//
//
//		/*
//		@param centers Image containing only the row to be processed before first call to this method.
//		@param ri Full ri image that will be updated.
//		@param start Start point of the row to be processed.
//		@param dim Dimension that we are processing.
//		@param step +1 or -1 to indicate the direction of the pass.
//		*/
//		inline void singlePass(const Image<internals::RiSet>& centers, Image<internals::RiSet>& ri, const Vec3c& rowStart, coord_t dim, coord_t step, vector<internals::ActiveSpheresItem>& activeSpheres)
//		{
//			// Stores the spheres that have been encountered and that have not been passed yet.
//			// Stores the center point, original radius, and ri.
//			//vector<internals::ActiveSpheresItem> activeSpheres;
//			vector<internals::ActiveSpheresItem> activeSpheresTmp;
//			vector<internals::ActiveSpheresItem> Ctmp;
//
//			internals::RiSet resultTmp;
//			internals::RiSet resultTmp2;
//
//			//activeSpheres.reserve(40);
//			activeSpheresTmp.reserve(40);
//			Ctmp.reserve(40);
//			resultTmp.reserve(40);
//			resultTmp2.reserve(40);
//
//			// Set start point to the start or to the end of the current row.
//			Vec3c p = rowStart;
//			if (step < 0)
//				p[dim] += (ri.dimension(dim) - 1);
//
//			coord_t dimensionality = ri.dimensionality();
//
//			for (coord_t i = 0; i < ri.dimension(dim); i++, p[dim] += step)
//			{
//				int32_t x = (int32_t)p[dim];
//
//				// If there is one or more sphere centers at the current location, add them to the set of active spheres.
//				const internals::RiSet& C = centers(x);
//
//				// The C list is sorted by R and so is activeSpheres list.
//				// Instead of finding place for each item, convert each item in C list to activeSpheres format (add x coordinate)
//				// and then merge the two sorted lists to construct the new activeSpheres list.
//				Ctmp.clear();
//				for (auto& item : C)
//					Ctmp.push_back({ item.R2, item.ri2, x });
//
//				activeSpheresTmp.clear();
//				activeSpheresTmp.resize(C.size() + activeSpheres.size());
//				merge(Ctmp.begin(), Ctmp.end(), activeSpheres.begin(), activeSpheres.end(), activeSpheresTmp.begin(), internals::activeSpheresSorter);
//				swap(activeSpheres, activeSpheresTmp);
//
//
//				// Iterate through all active spheres and calculate radius for next dimension.
//				resultTmp.clear();
//				auto it = activeSpheres.begin();
//				while (it != activeSpheres.end())
//				{
//					int32_t Rorig2 = it->R2;
//					int32_t R2 = it->ri2;
//					int32_t cx = it->xi;
//
//					int32_t dx = abs(x - cx);
//
//					// Calculate new ri^2
//					int32_t rn2 = R2 - dx * dx;
//					if (rn2 > 0)
//					{
//						// Insert ry2 into the list (TODO: Don't insert if the item is going to be removed anyway?)
//						resultTmp.push_back({ Rorig2, rn2 });
//
//						// Remove possible duplicate Rorig and replace it by (Rorig, max(duplicates))
//						if (resultTmp.size() > 1 && Rorig2 == resultTmp[resultTmp.size() - 2].R2)
//						{
//							if (rn2 > resultTmp[resultTmp.size() - 2].ri2)
//								resultTmp[resultTmp.size() - 2].ri2 = rn2;
//
//							resultTmp.pop_back();
//						}
//
//						it++;
//					}
//					else
//					{
//						// ry is non-positive, i.e. dx >= R
//						// This sphere is not active anymore, so remove it from the list of active spheres.
//						it = activeSpheres.erase(it);
//					}
//				}
//
//				// Rebuild ri list from resultTmp. Don't include those items that are hidden by other items.
//				if (resultTmp.size() > 0)
//				{
//					// Add ri from the previous pass to the ri list and save the result to resultTmp2.
//					// Note that both resultTmp2 and rilist are sorted so we can just merge them.
//					internals::RiSet& rilist = ri(p);
//					resultTmp2.clear();
//					resultTmp2.resize(rilist.size() + resultTmp.size());
//					merge(rilist.begin(), rilist.end(), resultTmp.begin(), resultTmp.end(), resultTmp2.begin(), internals::origRSorter);
//
//
//					// NOTE: This is the basic version, the memory saving version is below.
//					//// Linear time algorithm for finding relevant ri (those not hidden by other items).
//					//// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//					//// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//					//rilist.clear();
//					//rilist.push_back(resultTmp2[0]);
//					//for (size_t n = 1; n < resultTmp2.size(); n++)
//					//{
//					//	int32_t currri2 = rilist[rilist.size() - 1].ri2;
//					//	int32_t newri2 = resultTmp2[n].ri2;
//					//	if (newri2 > currri2)
//					//	{
//					//		rilist.push_back(resultTmp2[n]);
//					//	}
//					//}
//
//					// Linear time algorithm for finding relevant ri (those not hidden by other items).
//					// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//					// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//					rilist.clear();
//					rilist.push_back(resultTmp2[0]);
//					//if (dim < ri.dimensionality() - 1) // In the last dimension only the item with maximum R is needed, but we have separate function for that
//					//{
//					for (size_t n = 1; n < resultTmp2.size(); n++)
//					{
//						int32_t currri2 = rilist[rilist.size() - 1].ri2;
//						int32_t newri2 = resultTmp2[n].ri2;
//
//						if (newri2 > currri2) // This is the basic condition that works always (but may include unnecessary items in the rilist)
//						{
//							if (dim == dimensionality - 2)
//							{
//								// In the second last dimension only really visible spans are needed as there's no next dimension whose ri we would calculate based on the spans.
//								if (largestIntWhoseSquareIsLessThan(newri2) > largestIntWhoseSquareIsLessThan(currri2))
//									rilist.push_back(resultTmp2[n]);
//							}
//							else if (dim == dimensionality - 3)
//							{
//								// In the third last dimension we know that only those spans are required that produce visible circles in the output.
//								if (!internals::doesDiscretizedCircle1FitInto2Cached(newri2, currri2))
//									rilist.push_back(resultTmp2[n]);
//							}
//							else
//							{
//								// Here we could insert test if discretized spheres fit into each other etc. etc.
//								rilist.push_back(resultTmp2[n]);
//							}
//						}
//					}
//					//}
//
//					// Make sure we don't use extraneous memory by storing empty items in each pixel of the ri image.
//					rilist.shrink_to_fit();
//				}
//			}
//		}
//
//		/*
//		Same than singlePass but generates result image instead of ri image, and is thus a little bit faster.
//		@param centers Image containing only the row to be processed before first call to this method.
//		@param Result image. Must be set to zero before first call to this method.
//		@param start Start point of the row to be processed.
//		@param dim Dimension that we are processing.
//		@param step +1 or -1 to indicate the direction of the pass.
//		*/
//		inline void singlePassFinal(const Image<internals::RiSet>& centers, Image<int32_t>& result, const Vec3c& rowStart, size_t dim, coord_t step, const vector<internals::ActiveSpheresItem>& initialActiveSpheres)
//		{
//			// Stores the spheres that have been encountered and that have not been passed yet.
//			// Stores the sphere with the largest R at the top of the priority queue.
//			priority_queue<internals::ActiveSpheresItem, vector<internals::ActiveSpheresItem>, std::function<decltype(internals::activeSpheresSorter2)> > activeSpheres(internals::activeSpheresSorter2);
//
//			for (const auto& item : initialActiveSpheres)
//				activeSpheres.push(item);
//
//
//			// Set start point to the start or end of the current row.
//			Vec3c p = rowStart;
//			if (step < 0)
//				p[dim] += (result.dimension(dim) - 1);
//
//			for (coord_t i = 0; i < result.dimension(dim); i++, p[dim] += step)
//			{
//				int32_t x = (int32_t)p[dim];
//
//				// If there is one or more sphere centers at the current location, add them to the set of active spheres.
//				const internals::RiSet& C = centers(x);
//				for (auto& item : C)
//				{
//					activeSpheres.push({ item.R2, item.ri2, x });
//				}
//
//				while (!activeSpheres.empty())
//				{
//					auto& item = activeSpheres.top();
//
//					int32_t Rorig2 = item.R2;
//					int32_t R2 = item.ri2;
//					int32_t cx = item.xi;
//					int32_t dx = abs(x - cx);
//
//					// Calculate new ri^2
//					int32_t rn2 = R2 - dx * dx;
//					if (rn2 > 0)
//					{
//						// Note that previous pass may have assigned larger value to the output.
//						if (Rorig2 > result(p))
//							result(p) = Rorig2;
//						break;
//					}
//					else
//					{
//						// ry is non-positive, i.e. dx >= R
//						// This sphere is not active anymore, so remove it from the list of active spheres.
//						activeSpheres.pop();
//					}
//
//				}
//
//			}
//		}
//
//		/**
//		Builds initial active spheres list for singlePass function, based on list of spheres that may contribute to the image but whose center point is not in the image.
//		*/
//		inline void createInitialActiveSpheres(Image<vector<internals::ActiveSpheresItem> >& initialActiveSpheres, vector<Sphere2>& extraSpheres, size_t dim, const Vec3c& blockDimensions, coord_t dir)
//		{
//			// Find all spheres that are out of bounds (on negative side)
//			vector<Sphere2>::iterator outOfBoundsBegin;
//
//			// Add each sphere to initial active spheres list
//			Vec3c reducedDimensions = blockDimensions;
//			reducedDimensions[dim] = 1;
//			initialActiveSpheres.ensureSize(reducedDimensions);
//
//
//			if (dim == 0)
//			{
//				// Process spheres that are out of bounds in x- but inside image in y- and z- directions.
//				if (dir >= 0)
//					outOfBoundsBegin = std::partition(extraSpheres.begin(), extraSpheres.end(), [&](const Sphere2& s) { return !(s.pos.x < 0 && 0 <= s.pos.y && s.pos.y < blockDimensions.y && 0 <= s.pos.z && s.pos.z < blockDimensions.z); });
//				else
//					outOfBoundsBegin = std::partition(extraSpheres.begin(), extraSpheres.end(), [&](const Sphere2& s) { return !(s.pos.x >= blockDimensions.x && 0 <= s.pos.y && s.pos.y < blockDimensions.y && 0 <= s.pos.z && s.pos.z < blockDimensions.z); });
//
//				for (auto it = outOfBoundsBegin; it != extraSpheres.end(); it++)
//				{
//					// The sphere needs to be added to only one active spheres list in order to propagate it correctly.
//
//					const Sphere2& sphere = *it;
//
//					internals::ActiveSpheresItem item;
//					item.R2 = sphere.r2;
//					item.ri2 = sphere.r2;
//					item.xi = sphere.pos.x;
//					initialActiveSpheres(0, sphere.pos.y, sphere.pos.z).push_back(item);
//				}
//			}
//			else if (dim == 1)
//			{
//				// Process spheres that are out of bounds in y- but inside image in z-direction.
//
//				if (dir >= 0)
//					outOfBoundsBegin = std::partition(extraSpheres.begin(), extraSpheres.end(), [&](const Sphere2& s) { return !(s.pos.y < 0 && 0 <= s.pos.z && s.pos.z < blockDimensions.z); });
//				else
//					outOfBoundsBegin = std::partition(extraSpheres.begin(), extraSpheres.end(), [&](const Sphere2& s) { return !(s.pos.y >= blockDimensions.y && 0 <= s.pos.z && s.pos.z < blockDimensions.z); });
//
//				for (auto it = outOfBoundsBegin; it != extraSpheres.end(); it++)
//				{
//					// The sphere needs to be added to x-directional row of active spheres lists, as
//					// the first iteration would have converted the sphere to a row of values.
//
//					const Sphere2& sphere = *it;
//
//					int32_t ri = largestIntWhoseSquareIsLessThan(sphere.r2);
//					int32_t xmin = sphere.pos.x - ri;
//					int32_t xmax = sphere.pos.x + ri;
//
//					clamp(xmin, 0, (int32_t)blockDimensions.x - 1);
//					clamp(xmax, 0, (int32_t)blockDimensions.x - 1);
//
//					for (int32_t x = xmin; x <= xmax; x++)
//					{
//						int32_t dx = x - sphere.pos.x;
//						int32_t ri2 = sphere.r2 - dx * dx;
//						if (ri2 > 0)
//						{
//							internals::ActiveSpheresItem item;
//							item.R2 = sphere.r2;
//							item.ri2 = ri2;
//							item.xi = sphere.pos.y;
//
//							vector<internals::ActiveSpheresItem>& list = initialActiveSpheres(Vec3sc(x, 0, sphere.pos.z));
//
//							//list.push_back(item);
//
//							if (list.size() > 0)
//							{
//								const internals::ActiveSpheresItem& top = list.front();
//								if (dir >= 0 &&
//									item.xi + ri < top.xi + largestIntWhoseSquareIsLessThan(top.ri2) &&
//									item.R2 < top.R2)
//								{
//									// Do not add
//								}
//								else if (dir < 0 &&
//									item.xi - ri > top.xi - largestIntWhoseSquareIsLessThan(top.ri2) &&
//									item.R2 < top.R2)
//								{
//									// Do not add
//								}
//								else
//								{
//									list.push_back(item);
//									push_heap(list.begin(), list.end(), internals::activeSpheresSorter2);
//								}
//							}
//							else
//							{
//								list.push_back(item);
//							}
//						}
//					}
//				}
//			}
//			else if (dim == 2)
//			{
//				// Process spheres that are out of bounds in z-direction.
//
//				if (dir >= 0)
//					outOfBoundsBegin = std::partition(extraSpheres.begin(), extraSpheres.end(), [&](const Sphere2& s) { return !(s.pos.z < 0); });
//				else
//					outOfBoundsBegin = std::partition(extraSpheres.begin(), extraSpheres.end(), [&](const Sphere2& s) { return !(s.pos.z >= blockDimensions.z); });
//
//				for (auto it = outOfBoundsBegin; it != extraSpheres.end(); it++)
//				{
//					// The sphere needs to be added to circular region of active spheres lists in xy-plane,
//					// as the first and second iterations would have converted it to a circle.
//
//					const Sphere2& sphere = *it;
//
//					int32_t ri = largestIntWhoseSquareIsLessThan(sphere.r2);
//					int32_t xmin = sphere.pos.x - ri;
//					int32_t xmax = sphere.pos.x + ri;
//					int32_t ymin = sphere.pos.y - ri;
//					int32_t ymax = sphere.pos.y + ri;
//
//					clamp(xmin, 0, (int32_t)blockDimensions.x - 1);
//					clamp(xmax, 0, (int32_t)blockDimensions.x - 1);
//					clamp(ymin, 0, (int32_t)blockDimensions.y - 1);
//					clamp(ymax, 0, (int32_t)blockDimensions.y - 1);
//
//#pragma omp parallel for if(ri > 3)
//					for (int32_t y = ymin; y <= ymax; y++)
//					{
//						int32_t dy = y - sphere.pos.y;
//						for (int32_t x = xmin; x <= xmax; x++)
//						{
//							int32_t dx = x - sphere.pos.x;
//							int32_t ri2 = sphere.r2 - dy * dy - dx * dx;
//							if (ri2 > 0)
//							{
//								internals::ActiveSpheresItem item;
//								item.R2 = sphere.r2;
//								item.ri2 = ri2;
//								item.xi = sphere.pos.z;
//								//initialActiveSpheres(Vec3sc(x, y, 0)).push_back(item);
//
//								vector<internals::ActiveSpheresItem>& list = initialActiveSpheres(Vec3sc(x, y, 0));
//
//								if (list.size() > 0)
//								{
//									const internals::ActiveSpheresItem& top = list.front();
//									if (dir >= 0 &&
//										item.xi + ri < top.xi + largestIntWhoseSquareIsLessThan(top.ri2) &&
//										item.R2 < top.R2)
//									{
//										// Do not add
//									}
//									else if (dir < 0 &&
//										item.xi - ri > top.xi - largestIntWhoseSquareIsLessThan(top.ri2) &&
//										item.R2 < top.R2)
//									{
//										// Do not add
//									}
//									else
//									{
//										list.push_back(item);
//										push_heap(list.begin(), list.end(), internals::activeSpheresSorter2);
//									}
//								}
//								else
//								{
//									list.push_back(item);
//								}
//
//							}
//						}
//					}
//				}
//			}
//			else
//			{
//				// This never happens unless images of larger dimension than 3 become supported.
//				throw logic_error("Block-by-block thickness map calculation is not implemented for dimensionalities greater than three.");
//			}
//
//			// Sort initial active spheres lists
//#pragma omp parallel for if(!omp_in_parallel())
//			for (coord_t n = 0; n < initialActiveSpheres.pixelCount(); n++)
//			{
//				auto& l = initialActiveSpheres(n);
//				std::sort(l.begin(), l.end(), internals::activeSpheresSorter);
//			}
//
//			// Erase spheres that have been processed already
//			extraSpheres.erase(outOfBoundsBegin, extraSpheres.end());
//		}
//
//		//inline double meanCount(const Image<vector<internals::ActiveSpheresItem> >& lists)
//		//{
//		//	double sum = 0.0;
//		//	for (coord_t n = 0; n < lists.pixelCount(); n++)
//		//	{
//		//		sum += lists(n).size();
//		//	}
//		//	return sum / lists.pixelCount();
//		//}
//
//		/*
//		Makes forward and backward pass in single dimension.
//		@param ri Image containing ri values from processing of last dimension.
//		@param dim Dimension to process.
//		*/
//		inline void processDimension(Image<internals::RiSet>& ri, size_t dim, Image<int32_t>& result, vector<Sphere2>& extraSpheres, bool showProgressInfo)
//		{
//			// Determine count of pixels to process
//			Vec3c reducedDimensions = ri.dimensions();
//			reducedDimensions[dim] = 1;
//			coord_t rowCount = reducedDimensions.x * reducedDimensions.y * reducedDimensions.z;
//
//			Image<vector<internals::ActiveSpheresItem> > initialActiveSpheresStart, initialActiveSpheresEnd;
//			createInitialActiveSpheres(initialActiveSpheresStart, extraSpheres, dim, ri.dimensions(), +1);
//			createInitialActiveSpheres(initialActiveSpheresEnd, extraSpheres, dim, ri.dimensions(), -1);
//
//			//cout << "Dim = " << dim << ": start count = " << meanCount(initialActiveSpheresStart) << ", end count = " << meanCount(initialActiveSpheresEnd) << endl;
//
//			size_t counter = 0;
//#pragma omp parallel if(!omp_in_parallel() && result.pixelCount() > PARALLELIZATION_THRESHOLD)
//			{
//				// Temporary buffer
//				Image<internals::RiSet> centers(ri.dimension(dim));
//
//				// Determine start points of pixel rows and process each row
//#pragma omp for schedule(dynamic)
//				for (coord_t n = 0; n < rowCount; n++)
//				{
//					Vec3c start = indexToCoords(n, reducedDimensions);
//
//					// Make a copy of the current row as we update the row in the first pass but need
//					// the original data in the second pass.
//					Vec3c pos = start;
//					for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
//					{
//						centers(x) = ri(pos);
//						ri(pos).clear();
//					}
//
//					if (dim < ri.dimensionality() - 1)
//					{
//						// Left to right pass
//						singlePass(centers, ri, start, dim, 1, initialActiveSpheresStart(start));
//
//						// Right to left pass
//						singlePass(centers, ri, start, dim, -1, initialActiveSpheresEnd(start));
//					}
//					else
//					{
//						// TODO: Test if it is better/faster to combine singlePass and singlePassFinal as in the ridge version
//
//						// Left to right pass
//						singlePassFinal(centers, result, start, dim, 1, initialActiveSpheresStart(start));
//
//						// Right to left pass
//						singlePassFinal(centers, result, start, dim, -1, initialActiveSpheresEnd(start));
//					}
//
//					showThreadProgress(counter, rowCount, showProgressInfo);
//				}
//			}
//		}
//
//
//		/*
//		Copies squared sphere radius data from original centers image to ri image.
//		Removes copied spheres from spheres list so that at output the spheres list contains only spheres that are out of bounds of the block.
//		*/
//		//void prepare(Image<internals::RiSet>& ri, vector<Sphere2>& spheres)
//		//{
//		//	auto inImageBegin = std::partition(spheres.begin(), spheres.end(), [&](const Sphere2& s) { return !ri.isInImage(s.pos); });
//
//		//	for (auto it = inImageBegin; it != spheres.end(); it++)
//		//	{
//		//		const Sphere2& sphere = *it;
//		//		ri(sphere.pos).push_back({ sphere.r2, sphere.r2 });
//		//	}
//
//		//	spheres.erase(inImageBegin, spheres.end());
//		//}
//
//		/*
//		Copies squared sphere radius data from original centers image to ri image.
//		*/
//		inline void prepare(const Image<int32_t>& centers2, Image<internals::RiSet>& ri)
//		{
//			ri.ensureSize(centers2);
//
//			for (coord_t n = 0; n < centers2.pixelCount(); n++)
//			{
//				int32_t R2 = centers2(n);
//				if (R2 > 0)
//				{
//					ri(n).push_back({ R2, R2 });
//				}
//			}
//		}
//
//
//		void thickmap2(Image<int32_t>& dmap2, vector<Sphere2>& extraSpheres, Vec3d* counts, bool showProgressInfo)
//		{
//			int32_t maxR2 = max(dmap2);
//
//			for (const Sphere2& s : extraSpheres)
//			{
//				if (s.r2 > maxR2)
//					maxR2 = s.r2;
//			}
//
//			internals::buildCircleLookup(maxR2);
//
//			Image<internals::RiSet> ri(dmap2.dimensions());
//			prepare(dmap2, ri);
//			setValue(dmap2, 0);
//			for (size_t n = 0; n < ri.dimensionality(); n++)
//			{
//				processDimension(ri, n, dmap2, extraSpheres, showProgressInfo);
//
//				if (counts)
//				{
//					double sum = 0;
//					for (coord_t n = 0; n < ri.pixelCount(); n++)
//					{
//						sum += ri(n).size();
//					}
//					(*counts)[n] = sum;
//				}
//			}
//		}
//
//		void thickmap2Blocks(Image<int32_t>& dmap2, const Vec3c& blockSize, double* extraBytes, bool showProgressInfo)
//		{
//			// Create blocks
//			vector<Block> blocks = subdivide(dmap2, Vec3sc(blockSize));
//
//			if (extraBytes)
//			{
//				*extraBytes = 0;
//				for (const Block& block : blocks)
//				{
//					(*extraBytes) += sizeof(Block) + block.spheres.size() * sizeof(Sphere2);
//				}
//			}
//
//
//			// Process each block
//			Vec3d counts;
//			for (Block& b : blocks)
//			{
//				Vec3c pos = Vec3c(b.bounds.position());
//
//				// Cut input block from the original
//				Image<int32_t> blockImage(Vec3c(b.bounds.size()));
//				crop(dmap2, blockImage, pos);
//
//				if (extraBytes)
//				{
//					Vec3d currentCounts;
//					thickmap2(blockImage, b.spheres, &currentCounts, showProgressInfo);
//					counts = componentwiseMax(counts, currentCounts);
//				}
//				else
//				{
//					thickmap2(blockImage, b.spheres, nullptr, showProgressInfo);
//				}
//
//				// Copy the result to the output
//#pragma omp parallel for if(blockImage.pixelCount() > PARALLELIZATION_THRESHOLD)
//				for (coord_t z = 0; z < blockImage.depth(); z++)
//				{
//					for (coord_t y = 0; y < blockImage.height(); y++)
//					{
//						for (coord_t x = 0; x < blockImage.width(); x++)
//						{
//							dmap2(pos + Vec3c(x, y, z)) = blockImage(x, y, z);
//						}
//					}
//				}
//			}
//
//			if (extraBytes)
//			{
//				(*extraBytes) += counts.max() * sizeof(internals::RiItem) + dmap2.pixelCount() * sizeof(internals::RiSet);
//			}
//		}
//
//	}
//
//
//
//
//	namespace internals
//	{
//		
//
//
//		///**
//		//Draws all spheres corresponding to squared distance map.
//		//Sphere radii are given by square roots of distance map values and centers by corresponding distance map points.
//		//*/
//		//template<typename out_t> void drawSpheres2(const Image<int32_t>& dmap2, Image<out_t>& result, out_t color, bool showProgressInfo = true)
//		//{
//		//	result.ensureSize(dmap2);
//
//		//	size_t counter = 0;
//		//	#pragma omp parallel if(!omp_in_parallel() && result.pixelCount() > PARALLELIZATION_THRESHOLD)
//		//	for (coord_t z = 0; z < result.depth(); z++)
//		//	{
//		//		for (coord_t y = 0; y < result.height(); y++)
//		//		{
//		//			for (coord_t x = 0; x < result.width(); x++)
//		//			{
//		//				int32_t r2 = dmap2(x, y, z);
//		//				internals::draw2(result, Vec3c(x, y, z), r2, color);
//		//			}
//		//		}
//		//		showThreadProgress(counter, result.depth(), showProgressInfo);
//		//	}
//		//}
//
//		/**
//		Draws sphere of radius sqrt(R2) to tmap to all those points that have nonzero value in centers image.
//		*/
//		inline void fillSpheres2(const Image<int32_t>& centers, Image<int32_t>& tmap, int32_t R2)
//		{
//			int32_t Rint = (int32_t)ceil(sqrt(R2));
//
//			// Create mask
//			coord_t size = 2 * Rint + 1;
//			Image<int32_t> xMask(size, size, 1, -1);
//#pragma omp parallel for if(xMask.pixelCount() > PARALLELIZATION_THRESHOLD)
//			for (int32_t z = 0; z < xMask.height(); z++)
//			{
//				int32_t dz = z - Rint;
//				for (int32_t y = 0; y < xMask.width(); y++)
//				{
//					int32_t dy = y - Rint;
//					int32_t dx2 = R2 - dy * dy - dz * dz;
//					if (dx2 >= 0)
//					{
//						xMask(y, z) = largestIntWhoseSquareIsLessThan(dx2);
//					}
//				}
//			}
//
//			// Fill all spheres
//			Vec3c maxCoords = tmap.dimensions() - Vec3c(1, 1, 1);
//#pragma omp parallel for if(centers.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
//			for (coord_t n = 0; n < centers.pixelCount(); n++)
//			{
//				Vec3c pos = tmap.getCoords(n);
//
//				Vec3c start0 = pos - Vec3c(Rint, Rint, Rint);
//				Vec3c end = pos + Vec3c(Rint, Rint, Rint);
//
//				Vec3c start = start0;
//				clamp(start, Vec3c(), maxCoords);
//				clamp(end, Vec3c(), maxCoords);
//
//				for (coord_t z = start.z; z <= end.z; z++)
//				{
//					for (coord_t y = start.y; y <= end.y; y++)
//					{
//						coord_t xr = xMask(y - start0.y, z - start0.z);
//
//						if (xr >= 0)
//						{
//							coord_t xStart = pos.x - xr;
//							if (xStart < 0)
//								xStart = 0;
//
//							coord_t xEnd = pos.x + xr;
//							if (xEnd > maxCoords.x)
//								xEnd = maxCoords.x;
//
//							size_t ind = tmap.getLinearIndex(xStart, y, z);
//							size_t endInd = tmap.getLinearIndex(xEnd, y, z);
//							while (ind <= endInd)
//							{
//								tmap(ind) = R2;
//								ind++;
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//
//
//
//
//	namespace thresholddecomposition
//	{
//
//		/**
//		Calculates squared local radius from squared distance map.
//		Uses threshold decomposition algorithm.
//		*/
//		void thickmap2(const Image<int32_t>& dmap2, Image<int32_t>& result, bool showProgressInfo)
//		{
//			result.mustNotBe(dmap2);
//			result.ensureSize(dmap2.dimensions());
//
//			// Find unique values
//			set<int32_t> values;
//			unique(dmap2, values);
//
//			Image<int32_t> layer(dmap2.dimensions());
//			Image<int32_t> tmp(dmap2.dimensions());
//			size_t counter = 0;
//			for (int32_t val : values)
//			{
//				if (val > 0)
//				{
//					// Extract image containing each value from the dmap
//					for (coord_t n = 0; n < dmap2.pixelCount(); n++)
//					{
//						layer(n) = dmap2(n) == val ? val : 0;
//						tmp(n) = 0;
//					}
//
//					// Does not work reliably! There is something wrong in reverseEuclideanDistanceTransform.
//					//reverseEuclideanDistanceTransform(layer, 0, false);
//					//multiply(layer, val);
//					//max(result, layer);
//
//					// Seems to work ok. Note that we can use drawSpheres2 instead of drawMaxSpheres2 as all the spheres in one layer have the same value.
//					//internals::drawSpheres2(layer, tmp, val, false);
//					//internals::fillSpheres2(layer, tmp, val);
//					//max(result, tmp);
//
//					// Seems to work ok. Note that as all the thickness values are the same, only maximum of one ri item is stored per pixel.
//					// TODO: Use Danielsson here?
//					//Vec3d counts;
//					//itl2::dimred::thickmap2(layer, &counts, false);
//					itl2::dimred::thickmap2(layer, (double*)0, false);
//					max(result, layer);
//				}
//
//
//				showProgress(counter, values.size());
//				counter++;
//			}
//		}
//	}
//
//
//
//
//
//
//
//
//
//
//	namespace spanmap
//	{
//		struct Span
//		{
//			/**
//			Center coordinate of this span.
//			*/
//			int32_t center;
//
//			/**
//			Original radius of the sphere that caused this span to be added.
//			*/
//			int32_t R2;
//
//			/**
//			Radius of this span
//			*/
//			int32_t r2;
//
//			int32_t start;
//			int32_t end;
//
//			//int32_t start() const
//			//{
//			//	return center - largestIntWhoseSquareIsLessThan(r2);
//			//}
//
//			//int32_t end() const
//			//{
//			//	return center + largestIntWhoseSquareIsLessThan(r2);
//			//}
//		};
//
//		//typedef priority_queue<Span, vector<Span>, std::function<decltype(spanSetComparer)> > SpanSet;
//		//bool spanSetComparer(const Span& a, const Span& b)
//		//{
//		//	if (a.start() != b.start())
//		//		return a.start() > b.start();
//		//	return a.R2 < b.R2;
//		//}
//
//
//		bool spanSetComparer(const Span& a, const Span& b)
//		{
//			if (a.start != b.start)
//				return a.start > b.start;
//			return a.R2 > b.R2;
//		}
//
//		bool activeSpanComparer(const Span& a, const Span& b)
//		{
//			if (a.R2 != b.R2)
//				return a.R2 < b.R2;			// Larger first
//			return a.start > b.start;	// Smaller first
//		}
//
//		typedef vector<Span> SpanSet;
//		//typedef set<Span, spanSetComparer> SpanSet;
//
//		size_t maxCount(const Image<SpanSet>& img)
//		{
//			size_t M = 0;
//			for (coord_t n = 0; n < img.pixelCount(); n++)
//				M = std::max(M, img(n).size());
//			return M;
//		}
//
//		size_t totalCount(const Image<SpanSet>& img)
//		{
//			size_t M = 0;
//			for (coord_t n = 0; n < img.pixelCount(); n++)
//				M += img(n).size();
//			return M;
//		}
//
//		void writeSize(const Image<SpanSet>& img, const string& filename)
//		{
//			Image<uint8_t> counts(img.dimensions());
//			for (coord_t n = 0; n < img.pixelCount(); n++)
//				counts(n) = (uint8_t)img(n).size();
//			raw::writed(counts, filename);
//		}
//
//		bool isOccluded(const SpanSet& s, const Span& span, coord_t invDim)
//		{
//			if (invDim == 1) // Second last dimension: next dimension needs only spans that are visible.
//			{
//				for (const Span& oldSpan : s)
//				{
//					if (oldSpan.R2 >= span.R2 &&
//						oldSpan.start <= span.start &&
//						oldSpan.end >= span.end)
//						return true;
//				}
//			}
//			else if (invDim == 2) // Third last dimension: next dimension needs only spans that form a visible circle.
//			{
//				for (const Span& oldSpan : s)
//				{
//					if (oldSpan.R2 >= span.R2 &&
//						oldSpan.start < span.start &&
//						oldSpan.end > span.end)
//						return true;
//				}
//			}
//			return false;
//		}
//
//		void processLine(SpanSet& s, int32_t coordMax, int32_t spanCenter, Image<SpanSet>& target, coord_t firstCoord, coord_t invDim)
//		{
//			SpanSet activeSpheres;
//			for (int32_t z = 0; z < coordMax; z++)
//			{
//				// First add all spans that start at this pixel to the activeSpheres queue.
//				while (!s.empty())
//				{
//					const Span& nextSpan = s.back();
//
//					int32_t d = nextSpan.center - z;
//					int32_t rn2 = nextSpan.r2 - d * d;
//
//					if (rn2 > 0)
//					{
//						activeSpheres.push_back(nextSpan);
//						s.pop_back();
//					}
//					else
//					{
//						break;
//					}
//				}
//
//				// Process all active spans
//				auto it = activeSpheres.begin();
//				while (it != activeSpheres.end())
//				{
//					const Span& activeSpan = *it;
//
//					int32_t d = activeSpan.center - z;
//					int32_t rn2 = activeSpan.r2 - d * d;
//
//					if (rn2 > 0)
//					{
//						// We are still inside the span.
//						// Calculate corresponding span in next dimension.
//
//						Span span;
//						span.R2 = activeSpan.R2;
//						span.r2 = rn2;
//						span.center = spanCenter;
//						span.start = span.center - largestIntWhoseSquareIsLessThan(span.r2);
//						span.end = span.center + largestIntWhoseSquareIsLessThan(span.r2);
//
//						// Ehka: Lisaa talle funktiolle parametri joka tallettaa ActiveSpheres-prioriteettijonon jokaiselle target-kuvan pikselille
//						// Pida ylla activeSpheres-jonoa kun spaneja lisataan listaan. Jata jaljelle vain ne spanit jotka ovat olleet aktiivisia
//						// TAI: korvaa target kokonaan activeSpheres-tyyppisilla prioriteettilistoilla, jossa paallimmaisena on kulloinkin aktiivinen span
//						// Joka tapauksessa avainasia on saada isOccluded nopeammaksi, nyt se kay lapi koko spanilistan.
//
//						auto& list = target(firstCoord, z);
//						if (!isOccluded(list, span, invDim))
//							list.push_back(span);
//						//target(firstCoord, z).push_back(span);
//						//auto& list = target(firstCoord, z);
//						//list.insert(std::lower_bound(list.begin(), list.end(), span, spanSetComparer), span);
//
//						it++;
//					}
//					else
//					{
//						// We have passed currently active span. Remove it from the active spans (and try next span).
//						it = activeSpheres.erase(it);
//					}
//				}
//
//				if (s.empty() && activeSpheres.empty())
//					break;
//			}
//		}
//
//
//		void processLineFinal(SpanSet& s, int32_t coordMax, int32_t spanCenter, Image<int32_t>& target, coord_t zz)
//		{
//			SpanSet activeSpheres;
//			for (int32_t y = 0; y < coordMax; y++)
//			{
//				// First add all spans that start at this pixel to the activeSpheres queue.
//				while (!s.empty())
//				{
//					const Span& nextSpan = s.back();
//
//					int32_t d = nextSpan.center - y;
//					int32_t rn2 = nextSpan.r2 - d * d;
//
//					if (rn2 > 0)
//					{
//						activeSpheres.push_back(nextSpan);
//						s.pop_back();
//					}
//					else
//					{
//						break;
//					}
//				}
//
//				// Process all active spans
//				auto it = activeSpheres.begin();
//				while (it != activeSpheres.end())
//				{
//					const Span& activeSpan = *it;
//
//					int32_t d = activeSpan.center - y;
//					int32_t rn2 = activeSpan.r2 - d * d;
//
//					if (rn2 > 0)
//					{
//						// We are still inside the span.
//						// Calculate corresponding span in next dimension.
//
//						Span span;
//						span.R2 = activeSpan.R2;
//						span.r2 = rn2;
//						span.center = spanCenter;
//						span.start = span.center - largestIntWhoseSquareIsLessThan(span.r2);
//						span.end = span.center + largestIntWhoseSquareIsLessThan(span.r2);
//
//						if (span.start < 0)
//							span.start = 0;
//						if (span.end > (int32_t)target.width() - 1)
//							span.end = (int32_t)target.width() - 1;
//
//
//						// Plot the span instead of saving
//						for (int32_t xx = span.start; xx <= span.end; xx++)
//						{
//							int32_t& p = target(xx, y, zz);
//							if (span.R2 > p)
//								p = span.R2;
//						}
//
//						it++;
//					}
//					else
//					{
//						// We have passed currently active span. Remove it from the active spans (and try next span).
//						it = activeSpheres.erase(it);
//					}
//				}
//
//				if (s.empty() && activeSpheres.empty())
//					break;
//			}
//		}
//
//
//		void thickmap2(Image<int32_t>& dmap2, double* extraBytes)
//		{
//			// Convert sphere centers in dmap to spans in z-direction
//			// ==> set of spans for each (x, y)
//			Image<SpanSet> spansz(dmap2.width(), dmap2.height());
//			for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//			{
//				for (int32_t y = 0; y < (int32_t)dmap2.height(); y++)
//				{
//					for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//					{
//						// Is there a sphere center at the current location?
//						int32_t R2 = dmap2(x, y, z);
//						if (R2 > 0)
//						{
//							// Insert corresponding span into the spans image.
//							Span span;
//							span.R2 = R2;
//							span.center = z;
//							span.r2 = R2;
//							span.start = span.center - largestIntWhoseSquareIsLessThan(span.r2);
//							span.end = span.center + largestIntWhoseSquareIsLessThan(span.r2);
//
//							spansz(x, y).push_back(span);
//						}
//					}
//				}
//			}
//
//			//cout << "spans max size = " << maxCount(spansz) << endl;
//			double cnt = (double)totalCount(spansz);
//			//cout << "spans total size = " << cnt / dmap2.pixelCount() << " x image pixel count" << endl;
//
//			if (extraBytes)
//				*extraBytes = cnt * sizeof(Span) + spansz.pixelCount() * sizeof(SpanSet);
//
//#pragma omp parallel for
//			for (coord_t n = 0; n < spansz.pixelCount(); n++)
//				sort(spansz(n).begin(), spansz(n).end(), spanSetComparer);
//
//			//writeSize(spansz, "./spanmap/spansz");
//
//			// Convert spans in z-direction to spans in y-direction
//			// ==> set of spans for each (x, z)
//			Image<SpanSet> spansy(dmap2.width(), dmap2.depth());
//
//#pragma omp parallel for schedule(dynamic)
//			for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//			{
//				for (int32_t y = 0; y < (int32_t)dmap2.height(); y++)
//				{
//					processLine(spansz(x, y), (int32_t)dmap2.depth(), y, spansy, x, 2);
//				}
//			}
//
//			//cout << "spans max size = " << maxCount(spansy) << endl;
//			cnt = (double)totalCount(spansy);
//			//cout << "spans total size = " << cnt / dmap2.pixelCount() << " x image pixel count" << endl;
//
//			if (extraBytes)
//				*extraBytes = std::max(*extraBytes, (cnt * sizeof(Span) + spansy.pixelCount() * sizeof(SpanSet)));
//
//#pragma omp parallel for
//			for (coord_t n = 0; n < spansy.pixelCount(); n++)
//				sort(spansy(n).begin(), spansy(n).end(), spanSetComparer);
//
//			//writeSize(spansy, "./spanmap/spansy");
//
//
//
//			if (true)
//			{
//				// Direct plotting
//#pragma omp parallel for schedule(dynamic)
//				for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//				{
//					for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//					{
//						processLineFinal(spansy(x, z), (int32_t)dmap2.height(), x, dmap2, z);
//					}
//				}
//			}
//			else
//			{
//				// Save spans, then convert to image
//
//				// Convert spans in y-direction to spans in x-direction
//				// ==> set of spans for each (z, y)
//				Image<SpanSet> spansx(dmap2.depth(), dmap2.height());
//
//#pragma omp parallel for schedule(dynamic)
//				for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//				{
//					for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//					{
//						processLine(spansy(x, z), (int32_t)dmap2.height(), x, spansx, z, 1);
//					}
//				}
//
//				//cout << "spans max size = " << maxCount(spansx) << endl;
//				cnt = (double)totalCount(spansx);
//				//cout << "spans total size = " << cnt / dmap2.pixelCount() << " x image pixel count" << endl;
//
//				if (extraBytes)
//					*extraBytes = std::max(*extraBytes, (cnt * sizeof(Span) + spansx.pixelCount() * sizeof(SpanSet)));
//
//#pragma omp parallel for
//				for (coord_t n = 0; n < spansx.pixelCount(); n++)
//					sort(spansx(n).begin(), spansx(n).end(), spanSetComparer);
//
//
//				//writeSize(spansx, "./spanmap/spansx");
//
//
//				// Convert spans to output image
//#pragma omp parallel for
//				for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//				{
//					for (int32_t y = 0; y < (int32_t)dmap2.height(); y++)
//					{
//						SpanSet& s = spansx(z, y);
//
//						priority_queue<Span, vector<Span>, std::function<decltype(activeSpanComparer)> > activeSpheres(activeSpanComparer);
//						for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//						{
//							// First add all spans that start at this pixel to the activeSpheres queue.
//							while (!s.empty())
//							{
//								const Span& nextSpan = s.back();
//
//								int32_t d = nextSpan.center - x;
//								int32_t rn2 = nextSpan.r2 - d * d;
//
//								if (rn2 > 0)
//								{
//									activeSpheres.push(nextSpan);
//									s.pop_back();
//								}
//								else
//								{
//									break;
//								}
//							}
//
//							// Find currently active span with largest R2 and assign output value.
//							while (!activeSpheres.empty())
//							{
//								const Span& activeSpan = activeSpheres.top();
//
//								int32_t d = activeSpan.center - x;
//								int32_t rn2 = activeSpan.r2 - d * d;
//
//								if (rn2 > 0)
//								{
//									dmap2(x, y, z) = activeSpan.R2;
//									break;
//								}
//								else
//								{
//									// We have passed currently active span. Remove it from the active spans (and try next span).
//									activeSpheres.pop();
//								}
//							}
//
//							if (s.empty() && activeSpheres.empty())
//								break;
//						}
//					}
//				}
//			}
//
//			// Initial 3D version that kind of works, but span ends are not symmetrical in all dimensions
//			//// Convert sphere centers in dmap to spans in z-direction
//			//// ==> set of spans for each (x, y)
//			//Image<SpanSet> spansz(dmap2.width(), dmap2.height());
//			//for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//			//{
//			//	for (int32_t y = 0; y < (int32_t)dmap2.height(); y++)
//			//	{
//			//		for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//			//		{
//			//			// Is there a sphere center at the current location?
//			//			int32_t R2 = dmap2(x, y, z);
//			//			if (R2 > 0)
//			//			{
//			//				// Insert corresponding span into the spans image.
//			//				Span s;
//			//				s.R2 = R2;
//			//				s.center = z;
//			//				s.r2 = R2;
//
//			//				spansz(x, y).push_back(s);
//			//			}
//			//		}
//			//	}
//			//}
//
//			//writeSize(spansz, "./spanmap/spansz");
//
//			//// Convert spans in z-direction to spans in y-direction
//			//// ==> set of spans for each (x, z)
//			//Image<SpanSet> spansy(dmap2.width(), dmap2.depth());
//
//			//for (int32_t y = 0; y < (int32_t)dmap2.height(); y++)
//			//{
//			//	for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//			//	{
//			//		//if (y == 54)
//			//		//	cout << "debug" << endl;
//
//			//		SpanSet& s = spansz(x, y);
//
//			//		SpanSet activeSpheres;
//			//		for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//			//		{
//			//			// First add all spans that start at this pixel to the activeSpheres queue.
//			//			while (!s.empty())
//			//			{
//			//				const Span& nextSpan = s.back();
//			//				if (nextSpan.start() <= z)
//			//				{
//			//					activeSpheres.push_back(nextSpan);
//			//					s.pop_back();
//			//				}
//			//				else
//			//				{
//			//					break;
//			//				}
//			//			}
//
//			//			// Process all active spans
//			//			auto it = activeSpheres.begin();
//			//			while(it != activeSpheres.end())
//			//			{
//			//				const Span& activeSpan = *it;
//
//			//				if (z <= activeSpan.end())
//			//				{
//			//					// We are still inside the span.
//			//					// Calculate corresponding span in next dimension.
//
//			//					int32_t dz = activeSpan.center - z;
//			//					int32_t rn2 = activeSpan.R2 - dz * dz;
//
//			//					Span span;
//			//					span.R2 = activeSpan.R2;
//			//					span.r2 = rn2;
//			//					span.center = y;
//
//			//					spansy(x, z).push_back(span);
//
//			//					it++;
//			//				}
//			//				else
//			//				{
//			//					// We have passed currently active span. Remove it from the active spans (and try next span).
//			//					it = activeSpheres.erase(it);
//			//				}
//			//			}
//
//			//			if (s.empty() && activeSpheres.empty())
//			//				break;
//			//		}
//			//	}
//			//}
//
//			//for (coord_t n = 0; n < spansy.pixelCount(); n++)
//			//	sort(spansy(n).begin(), spansy(n).end(), spanSetComparer);
//
//			//writeSize(spansy, "./spanmap/spansy");
//
//
//
//			//// Convert spans in y-direction to spans in x-direction
//			//// ==> set of spans for each (y, z)
//			//Image<SpanSet> spansx(dmap2.height(), dmap2.depth());
//
//			//for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//			//{
//			//	for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//			//	{
//			//		//if (y == 54)
//			//		//	cout << "debug" << endl;
//
//			//		SpanSet& s = spansy(x, z);
//
//			//		SpanSet activeSpheres;
//			//		for (int32_t y = 0; y < (int32_t)dmap2.height(); y++)
//			//		{
//			//			// First add all spans that start at this pixel to the activeSpheres queue.
//			//			while (!s.empty())
//			//			{
//			//				const Span& nextSpan = s.back();
//			//				if (nextSpan.start() <= y)
//			//				{
//			//					activeSpheres.push_back(nextSpan);
//			//					s.pop_back();
//			//				}
//			//				else
//			//				{
//			//					break;
//			//				}
//			//			}
//
//			//			// Process all active spans
//			//			auto it = activeSpheres.begin();
//			//			while (it != activeSpheres.end())
//			//			{
//			//				const Span& activeSpan = *it;
//
//			//				if (y <= activeSpan.end())
//			//				{
//			//					// We are still inside the span.
//			//					// Calculate corresponding span in next dimension.
//
//			//					int32_t dy = activeSpan.center - y;
//			//					int32_t rn2 = activeSpan.R2 - dy * dy;
//
//			//					Span span;
//			//					span.R2 = activeSpan.R2;
//			//					span.r2 = rn2;
//			//					span.center = x;
//
//			//					spansx(y, z).push_back(span);
//
//			//					it++;
//			//				}
//			//				else
//			//				{
//			//					// We have passed currently active span. Remove it from the active spans (and try next span).
//			//					it = activeSpheres.erase(it);
//			//				}
//			//			}
//
//			//			if (s.empty() && activeSpheres.empty())
//			//				break;
//			//		}
//			//	}
//			//}
//
//			//for (coord_t n = 0; n < spansx.pixelCount(); n++)
//			//	sort(spansx(n).begin(), spansx(n).end(), spanSetComparer);
//
//
//			//writeSize(spansx, "./spanmap/spansx");
//
//
//			//// Convert spans to output image
//			//for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//			//{
//			//	for (int32_t y = 0; y < (int32_t)dmap2.height(); y++)
//			//	{
//			//		//if (y == 54)
//			//		//	cout << "debug" << endl;
//
//			//		SpanSet& s = spansx(y, z);
//
//			//		priority_queue<Span, vector<Span>, std::function<decltype(activeSpanComparer)> > activeSpheres(activeSpanComparer);
//			//		for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//			//		{
//			//			// First add all spans that start at this pixel to the activeSpheres queue.
//			//			while (!s.empty())
//			//			{
//			//				const Span& nextSpan = s.back();
//			//				if (nextSpan.start() <= x)
//			//				{
//			//					activeSpheres.push(nextSpan);
//			//					s.pop_back();
//			//				}
//			//				else
//			//				{
//			//					break;
//			//				}
//			//			}
//
//			//			// Find currently active span with largest R2 and assign output value.
//			//			while (!activeSpheres.empty())
//			//			{
//			//				const Span& activeSpan = activeSpheres.top();
//
//			//				if (x <= activeSpan.end())
//			//				{
//			//					dmap2(x, y, z) = activeSpan.R2;
//			//					break;
//			//				}
//			//				else
//			//				{
//			//					// We have passed currently active span. Remove it from the active spans (and try next span).
//			//					activeSpheres.pop();
//			//				}
//			//			}
//
//			//			if (s.empty() && activeSpheres.empty())
//			//				break;
//			//		}
//			//	}
//			//}
//
//
//
//			// Initial 2D version that works
//			//Vec3c dims = dmap2.dimensions();
//			//dims[0] = 1;
//			//Image<SpanSet> spans(dims, SpanSet(spanSetComparer));
//
//			//// Convert sphere centers in dmap to spans in y-direction
//			//for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//			//{
//			//	for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//			//	{
//			//		for (int32_t y = 0; y < (int32_t)dmap2.height(); y++)
//			//		{
//			//			// Is there a sphere center at the current location?
//			//			int32_t R2 = dmap2(x, y, z);
//			//			if (R2 > 0)
//			//			{
//			//				// Insert corresponding spans into the spans image.
//			//				int32_t R = largestIntWhoseSquareIsLessThan(R2);
//			//				int32_t starty = std::max(0, y - R);
//			//				int32_t endy = std::min((int32_t)dmap2.height() - 1, y + R);
//			//				for (int32_t yy = starty; yy <= endy; yy++)
//			//				{
//			//					Span s;
//			//					s.R2 = R2;
//			//					s.center = x;
//			//					int32_t dy = yy - y;
//			//					s.r2 = R2 - dy * dy;
//			//					// TODO: Don't insert if the span is completely hidden by another span.
//			//					spans(0, yy, z).push(s);
//			//				}
//			//			}
//			//		}
//			//	}
//			//}
//
//			//// Convert spans to output image
//			//for (int32_t z = 0; z < (int32_t)dmap2.depth(); z++)
//			//{
//			//	for (int32_t y = 0; y < (int32_t)dmap2.height(); y++)
//			//	{
//			//		//if (y == 54)
//			//		//	cout << "debug" << endl;
//
//			//		SpanSet& s = spans(0, y, z);
//
//			//		priority_queue<Span, vector<Span>, std::function<decltype(activeSpanComparer)> > activeSpheres(activeSpanComparer);
//			//		for (int32_t x = 0; x < (int32_t)dmap2.width(); x++)
//			//		{
//			//			// First add all spans that start at this pixel to the activeSpheres queue.
//			//			while (!s.empty())
//			//			{
//			//				const Span& nextSpan = s.top();
//			//				if (nextSpan.start() <= x)
//			//				{
//			//					activeSpheres.push(nextSpan);
//			//					s.pop();
//			//				}
//			//				else
//			//				{
//			//					break;
//			//				}
//			//			}
//
//			//			// Find currently active span with largest R2 and assign output value.
//			//			while (!activeSpheres.empty())
//			//			{
//			//				const Span& activeSpan = activeSpheres.top();
//
//			//				if (x <= activeSpan.end())
//			//				{
//			//					// We are still inside the currently active span. Assign output value.
//			//					dmap2(x, y, z) = activeSpan.R2;
//			//					break;
//			//				}
//			//				else
//			//				{
//			//					// We have passed currently active span. Remove it from the active spans (and try next span).
//			//					activeSpheres.pop();
//			//				}
//			//			}
//
//			//			if (s.empty() && activeSpheres.empty())
//			//				break;
//			//		}
//			//	}
//			//}
//
//		}
//	}
//
//
//
//
//	namespace tests
//	{
//		void spanMap()
//		{
//			Image<uint8_t> geom(100, 100, 100);
//			//draw(geom, Sphere(Vec3f(50, 45, 50), 10.0f), (uint8_t)1);
//			//draw(geom, Sphere(Vec3f(58, 45, 50), 8.0f), (uint8_t)1);
//			//draw(geom, Sphere(Vec3f(45, 45, 50), 7.0f), (uint8_t)1);
//			generateSimpleGeometry(geom, 1237);
//
//			raw::writed(geom, "./spanmap/geom");
//
//			Image<int32_t> dmap2;
//			itl2::distanceTransform2(geom, dmap2);
//
//			//// Only one line at y=50
//			//for (coord_t z = 0; z < geom.depth(); z++)
//			//{
//			//	for (coord_t y = 0; y < geom.height(); y++)
//			//	{
//			//		//if (y != 45 || z != 50)
//			//		//if (z != 50)
//			//		if(!(y == 45 && (z == 50 || z == 51)))
//			//		{
//			//			for (coord_t x = 0; x < geom.width(); x++)
//			//				dmap2(x, y, z) = 0;
//			//		}
//			//	}
//			//}
//
//			//Image<int32_t> ridge2;
//			//itl2::centersOfLocallyMaximalSpheres(dmap2, ridge2);
//			//setValue(dmap2, ridge2);
//
//
//
//			//Image<float32_t> dmap(100, 100, 100);
//			////dmap(50, 50, 50) = 7;
//			////dmap(55, 50, 50) = 5;
//			////dmap(45, 50, 50) = 6;
//			////dmap(52, 50) = 2;
//			////dmap(55, 48) = 3;
//
//			////dmap(52, 50) = 6;
//
//			////dmap(50, 45, 50) = 10;
//			////dmap(58, 45, 50) = 8;
//			////dmap(45, 45, 50) = 7;
//
//			//dmap(50, 50, 50) = 7;
//			//dmap(50, 50, 52) = 8;
//			////dmap(50, 50, 75) = 6;
//
//			//// Convert dmap to dmap^2
//			//Image<int32_t> dmap2;
//			//{
//			//	Image<float32_t> tmp;
//			//	setValue(tmp, dmap);
//			//	multiply(tmp, tmp);
//			//	convert(tmp, dmap2);
//			//}
//
//			raw::writed(dmap2, "./spanmap/dmap2");
//
//			// Calculate local thickness using simple sphere plotting algorithm and distance map
//			//Image<int32_t> thicknessSimple;
//			//drawMaxSpheres2(dmap2, thicknessSimple);
//			//raw::writed(thicknessSimple, "./dimred2/thickness2_simple");
//
//
//			// Calculate local thickness using dimensionality reduction algorithm
//			Timer timer;
//			Image<int32_t> thickness2;
//			timer.start();
//			setValue(thickness2, dmap2);
//			itl2::spanmap::thickmap2(thickness2, 0);
//			timer.stop();
//			raw::writed(thickness2, "./spanmap/thickness2_spanmap");
//
//			cout << "Span map took " << timer.getTime() << " ms" << endl;
//
//
//			Image<int32_t> thickness2Simple;
//			timer.start();
//			setValue(thickness2Simple, dmap2);
//			itl2::optimized::thickmap2(thickness2Simple);
//			//itl2::standard::thickmap2(dmap2, thickness2Simple);
//			timer.stop();
//			raw::writed(thickness2Simple, "./spanmap/thickness2_simple");
//
//			cout << "Optimized took " << timer.getTime() << " ms" << endl;
//
//			testAssert(equals(thickness2, thickness2Simple), "spanmap gives wrong thickness");
//		}
//
//
//		void dimredBlocks2D()
//		{
//
//			Vec3sc c(15, 15, 0);
//			int32_t R = 20;
//			int32_t r = 10;
//
//			dimredblocks::Block block;
//			block.bounds = AABox<int32_t>::fromPosSize(Vec3sc(0, 0, 0), 2 * c + Vec3sc(1, 1, 1));
//			Image<int32_t> result(Vec3c(block.bounds.size()));
//
//			//block.spheres.push_back({ Vec3sc(8, 15, 0), 8 * 8 });
//			//block.spheres.push_back({ Vec3sc(-4, 10, 0), 10 * 10 });
//			//block.spheres.push_back({ Vec3sc(35, 13, 0), 9 * 9 });
//			//block.spheres.push_back({ Vec3sc(10, 35, 0), 9 * 9 });
//			//block.spheres.push_back({ Vec3sc(10, -5, 0), 11 * 11 });
//			//block.spheres.push_back({ Vec3sc(35, -5, 0), 12 * 12 });
//
//
//			block.spheres.push_back({ c + Vec3sc(-R, -R,  0), r * r });
//			block.spheres.push_back({ c + Vec3sc(0, -R,  0), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R, -R,  0), r * r });
//
//			block.spheres.push_back({ c + Vec3sc(-R,  0,  0), r * r });
//
//			//block.spheres.push_back({ c + Vec3sc(0,  0,  0), r * r });
//			result(c) = r * r;
//
//			block.spheres.push_back({ c + Vec3sc(+R,  0,  0), r * r });
//
//			block.spheres.push_back({ c + Vec3sc(-R, +R,  0), r * r });
//			block.spheres.push_back({ c + Vec3sc(0, +R,  0), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R, +R,  0), r * r });
//
//
//
//			dimredblocks::thickmap2(result, block.spheres, nullptr, true);
//
//			testAssert(block.spheres.size() <= 0, "not all extra spheres were consumed.");
//
//			raw::writed(result, "./dimredblocks/simple");
//		}
//
//		void dimredBlocks3D()
//		{
//			Vec3sc c(15, 15, 15);
//			int32_t R = 20;
//			int32_t r = 10;
//
//			dimredblocks::Block block;
//			block.bounds = AABox<int32_t>::fromPosSize(Vec3sc(0, 0, 0), 2 * c + Vec3sc(1, 1, 1));
//			Image<int32_t> result(Vec3c(block.bounds.size()));
//			block.spheres.push_back({ c + Vec3sc(-R, -R, -R), r * r });
//			block.spheres.push_back({ c + Vec3sc(0, -R, -R), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R, -R, -R), r * r });
//
//			block.spheres.push_back({ c + Vec3sc(-R,  0, -R), r * r });
//			block.spheres.push_back({ c + Vec3sc(0,  0, -R), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R,  0, -R), r * r });
//
//			block.spheres.push_back({ c + Vec3sc(-R, +R, -R), r * r });
//			block.spheres.push_back({ c + Vec3sc(0, +R, -R), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R, +R, -R), r * r });
//
//
//			block.spheres.push_back({ c + Vec3sc(-R, -R,  0), r * r });
//			block.spheres.push_back({ c + Vec3sc(0, -R,  0), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R, -R,  0), r * r });
//
//			block.spheres.push_back({ c + Vec3sc(-R,  0,  0), r * r });
//
//			//block.spheres.push_back({ c + Vec3sc( 0,  0,  0), r * r });
//			result(c) = r * r;
//
//			block.spheres.push_back({ c + Vec3sc(+R,  0,  0), r * r });
//
//			block.spheres.push_back({ c + Vec3sc(-R, +R,  0), r * r });
//			block.spheres.push_back({ c + Vec3sc(0, +R,  0), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R, +R,  0), r * r });
//
//			block.spheres.push_back({ c + Vec3sc(-R, -R, +R), r * r });
//			block.spheres.push_back({ c + Vec3sc(0, -R, +R), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R, -R, +R), r * r });
//
//			block.spheres.push_back({ c + Vec3sc(-R,  0, +R), r * r });
//			block.spheres.push_back({ c + Vec3sc(0,  0, +R), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R,  0, +R), r * r });
//
//			block.spheres.push_back({ c + Vec3sc(-R, +R, +R), r * r });
//			block.spheres.push_back({ c + Vec3sc(0, +R, +R), r * r });
//			block.spheres.push_back({ c + Vec3sc(+R, +R, +R), r * r });
//
//
//			dimredblocks::thickmap2(result, block.spheres, nullptr, true);
//
//			testAssert(block.spheres.size() <= 0, "not all extra spheres were consumed.");
//
//			raw::writed(result, "./dimredblocks/simple");
//		}
//	}
//
//
//
//}




// Old versions of some parts of dimred code:
///*
		//@param centers Image containing only the row to be processed before first call to this method.
		//@param ri Full ri image that will be updated.
		//@param start Start point of the row to be processed.
		//@param dim Dimension that we are processing.
		//@param step +1 or -1 to indicate the direction of the pass.
		//*/
		//inline void singlePass(const Image<internals::RiSet>& centers, Image<internals::RiSet>& ri, const Vec3c& rowStart, coord_t dim, coord_t step)
		//{
		//	// Stores the spheres that have been encountered and that have not been passed yet.
		//	// Stores the center point, original radius, and ri.
		//	vector<internals::ActiveSpheresItem> activeSpheres;
		//	vector<internals::ActiveSpheresItem> activeSpheresTmp;
		//	vector<internals::ActiveSpheresItem> Ctmp;

		//	internals::RiSet resultTmp;
		//	internals::RiSet resultTmp2;

		//	activeSpheres.reserve(40);
		//	activeSpheresTmp.reserve(40);
		//	Ctmp.reserve(40);
		//	resultTmp.reserve(40);
		//	resultTmp2.reserve(40);

		//	// Set start point to the start or to the end of the current row.
		//	Vec3c p = rowStart;
		//	if (step < 0)
		//		p[dim] += (ri.dimension(dim) - 1);

		//	coord_t dimensionality = ri.dimensionality();

		//	for (coord_t i = 0; i < ri.dimension(dim); i++, p[dim] += step)
		//	{
		//		int32_t x = (int32_t)p[dim];

		//		// If there is one or more sphere centers at the current location, add them to the set of active spheres.
		//		const internals::RiSet& C = centers(x);

		//		// The C list is sorted by R and so is activeSpheres list.
		//		// Instead of finding place for each item, convert each item in C list to activeSpheres format (add x coordinate)
		//		// and then merge the two sorted lists to construct the new activeSpheres list.
		//		Ctmp.clear();
		//		for (auto& item : C)
		//			Ctmp.push_back({ item.R2, item.ri2, x });

		//		activeSpheresTmp.clear();
		//		activeSpheresTmp.resize(C.size() + activeSpheres.size());
		//		merge(Ctmp.begin(), Ctmp.end(), activeSpheres.begin(), activeSpheres.end(), activeSpheresTmp.begin(), internals::activeSpheresSorter);
		//		swap(activeSpheres, activeSpheresTmp);


		//		// Iterate through all active spheres and calculate radius for next dimension.
		//		resultTmp.clear();
		//		auto it = activeSpheres.begin();
		//		while (it != activeSpheres.end())
		//		{
		//			int32_t Rorig2 = it->R2;
		//			int32_t R2 = it->ri2;
		//			int32_t cx = it->xi;

		//			int32_t dx = abs(x - cx);

		//			// Calculate new ri^2
		//			int32_t rn2 = R2 - dx * dx;
		//			if (rn2 > 0)
		//			{
		//				// Insert ry2 into the list (TODO: Don't insert if the item is going to be removed anyway?)
		//				resultTmp.push_back({ Rorig2, rn2 });

		//				// Remove possible duplicate Rorig and replace it by (Rorig, max(duplicates))
		//				if (resultTmp.size() > 1 && Rorig2 == resultTmp[resultTmp.size() - 2].R2)
		//				{
		//					if (rn2 > resultTmp[resultTmp.size() - 2].ri2)
		//						resultTmp[resultTmp.size() - 2].ri2 = rn2;

		//					resultTmp.pop_back();
		//				}

		//				it++;
		//			}
		//			else
		//			{
		//				// ry is non-positive, i.e. dx >= R
		//				// This sphere is not active anymore, so remove it from the list of active spheres.
		//				it = activeSpheres.erase(it);
		//			}
		//		}

		//		// Rebuild ri list from resultTmp. Don't include those items that are hidden by other items.
		//		if (resultTmp.size() > 0)
		//		{
		//			// Add ri from the previous pass to the ri list and save the result to resultTmp2.
		//			// Note that both resultTmp2 and rilist are sorted so we can just merge them.
		//			internals::RiSet& rilist = ri(p);
		//			resultTmp2.clear();
		//			resultTmp2.resize(rilist.size() + resultTmp.size());
		//			merge(rilist.begin(), rilist.end(), resultTmp.begin(), resultTmp.end(), resultTmp2.begin(), internals::origRSorter);


		//			// NOTE: This is the basic version, the memory saving version is below.
		//			//// Linear time algorithm for finding relevant ri (those not hidden by other items).
		//			//// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
		//			//// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
		//			//rilist.clear();
		//			//rilist.push_back(resultTmp2[0]);
		//			//for (size_t n = 1; n < resultTmp2.size(); n++)
		//			//{
		//			//	int32_t currri2 = rilist[rilist.size() - 1].ri2;
		//			//	int32_t newri2 = resultTmp2[n].ri2;
		//			//	if (newri2 > currri2)
		//			//	{
		//			//		rilist.push_back(resultTmp2[n]);
		//			//	}
		//			//}

		//			// Linear time algorithm for finding relevant ri (those not hidden by other items).
		//			// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
		//			// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
		//			rilist.clear();
		//			rilist.push_back(resultTmp2[0]);
		//			//if (dim < ri.dimensionality() - 1) // In the last dimension only the item with maximum R is needed, but we have separate function for that
		//			//{
		//			for (size_t n = 1; n < resultTmp2.size(); n++)
		//			{
		//				int32_t currri2 = rilist[rilist.size() - 1].ri2;
		//				int32_t newri2 = resultTmp2[n].ri2;

		//				if (newri2 > currri2) // This is the basic condition that works always (but may include unnecessary items in the rilist)
		//				{
		//					if (dim == dimensionality - 2)
		//					{
		//						// In the second last dimension only really visible spans are needed as there's no next dimension whose ri we would calculate based on the spans.
		//						if(largestIntWhoseSquareIsLessThan(newri2) > largestIntWhoseSquareIsLessThan(currri2))
		//							rilist.push_back(resultTmp2[n]);
		//					}
		//					else if (dim == dimensionality - 3)
		//					{
		//						// In the third last dimension we know that only those spans are required that produce visible circles in the output.
		//						if (!internals::doesDiscretizedCircle1FitInto2Cached(newri2, currri2))
		//							rilist.push_back(resultTmp2[n]);
		//					}
		//					else
		//					{
		//						// Here we could insert test if discretized spheres fit into each other etc. etc.
		//						rilist.push_back(resultTmp2[n]);
		//					}
		//				}
		//			}
		//			//}

		//			// Make sure we don't use extraneous memory by storing empty items in each pixel of the ri image.
		//			rilist.shrink_to_fit();
		//		}
		//	}
		//}

		///*
		//Same than singlePass but generates result image instead of ri image, and is thus a little bit faster.
		//@param centers Image containing only the row to be processed before first call to this method.
		//@param Result image. Must be set to zero before first call to this method.
		//@param start Start point of the row to be processed.
		//@param dim Dimension that we are processing.
		//@param step +1 or -1 to indicate the direction of the pass.
		//*/
		//inline void singlePassFinal(const Image<internals::RiSet>& centers, Image<int32_t>& result, const Vec3c& rowStart, size_t dim, coord_t step)
		//{
		//	// Stores the spheres that have been encountered and that have not been passed yet.
		//	// Stores the sphere with the largest R at the top of the priority queue.
		//	priority_queue<internals::ActiveSpheresItem, vector<internals::ActiveSpheresItem>, std::function<decltype(internals::activeSpheresSorter2)> > activeSpheres(internals::activeSpheresSorter2);


		//	// Set start point to the start or end of the current row.
		//	Vec3c p = rowStart;
		//	if (step < 0)
		//		p[dim] += (result.dimension(dim) - 1);

		//	for (coord_t i = 0; i < result.dimension(dim); i++, p[dim] += step)
		//	{
		//		int32_t x = (int32_t)p[dim];

		//		// If there is one or more sphere centers at the current location, add them to the set of active spheres.
		//		const internals::RiSet& C = centers(x);
		//		for (auto& item : C)
		//		{
		//			activeSpheres.push({ item.R2, item.ri2, x });
		//		}

		//		while (!activeSpheres.empty())
		//		{
		//			auto& item = activeSpheres.top();

		//			int32_t Rorig2 = item.R2;
		//			int32_t R2 = item.ri2;
		//			int32_t cx = item.xi;
		//			int32_t dx = abs(x - cx);

		//			// Calculate new ri^2
		//			int32_t rn2 = R2 - dx * dx;
		//			if (rn2 > 0)
		//			{
		//				// Note that previous pass may have assigned larger value to the output.
		//				if (Rorig2 > result(p))
		//					result(p) = Rorig2;
		//				break;
		//			}
		//			else
		//			{
		//				// ry is non-positive, i.e. dx >= R
		//				// This sphere is not active anymore, so remove it from the list of active spheres.
		//				activeSpheres.pop();
		//			}

		//		}

		//	}
		//}

		// This version is simplified from singlePass(...) function. This is slower than the above one.
		//inline void singlePassFinal(const Image<RiSet>& centers, /*Image<RiSet>& ri,*/Image<int32_t>& result, const Vec3c& rowStart, coord_t dim, coord_t step)
		//{
		//	// Stores the spheres that have been encountered and that have not been passed yet.
		//	// Stores the center point, original radius, and ri.
		//	vector<ActiveSpheresItem> activeSpheres;
		//	vector<ActiveSpheresItem> activeSpheresTmp;
		//	vector<ActiveSpheresItem> Ctmp;

		//	activeSpheres.reserve(40);
		//	activeSpheresTmp.reserve(40);
		//	Ctmp.reserve(40);

		//	// Set start point to the start or to the end of the current row.
		//	Vec3c p = rowStart;
		//	if (step < 0)
		//		p[dim] += (result.dimension(dim) - 1);

		//	coord_t dimensionality = result.dimensionality();

		//	for (coord_t i = 0; i < result.dimension(dim); i++, p[dim] += step)
		//	{
		//		int32_t x = (int32_t)p[dim];

		//		// If there is one or more sphere centers at the current location, add them to the set of active spheres.
		//		const RiSet& C = centers(x);

		//		// The C list is sorted by R and so is activeSpheres list.
		//		// Instead of finding place for each item, convert each item in C list to activeSpheres format (add x coordinate)
		//		// and then merge the two sorted lists to construct the new activeSpheres list.
		//		Ctmp.clear();
		//		for (auto& item : C)
		//			Ctmp.push_back({ item.R2, item.ri2, x });

		//		activeSpheresTmp.clear();
		//		activeSpheresTmp.resize(C.size() + activeSpheres.size());
		//		merge(Ctmp.begin(), Ctmp.end(), activeSpheres.begin(), activeSpheres.end(), activeSpheresTmp.begin(), activeSpheresSorter);
		//		swap(activeSpheres, activeSpheresTmp);


		//		// Iterate through all active spheres and calculate radius for next dimension.
		//		auto it = activeSpheres.begin();
		//		while (it != activeSpheres.end())
		//		{
		//			int32_t Rorig2 = it->R2;
		//			int32_t R2 = it->ri2;
		//			int32_t cx = it->xi;

		//			int32_t dx = abs(x - cx);

		//			// Calculate new ri^2
		//			int32_t rn2 = R2 - dx * dx;
		//			if (rn2 > 0)
		//			{
		//				if (Rorig2 > result(p))
		//					result(p) = Rorig2;

		//				//it++;
		//				break;
		//			}
		//			else
		//			{
		//				// ry is non-positive, i.e. dx >= R
		//				// This sphere is not active anymore, so remove it from the list of active spheres.
		//				it = activeSpheres.erase(it);
		//			}
		//		}
		//	}
		//}


		//inline bool ridgeTestComparerri(const internals::RiItem& a, const internals::RiItem& b)
		//{
		//	if (a.ri2 != b.ri2)
		//		return a.ri2 < b.ri2;
		//	return a.R2 < b.R2;
		//}

		///*
		//Makes forward and backward pass in single dimension.
		//@param ri Image containing ri values from processing of last dimension.
		//@param dim Dimension to process.
		//*/
		//inline void processDimension(Image<internals::RiSet>& ri, size_t dim, Image<int32_t>& result, bool showProgressInfo)
		//{
		//	// Determine count of pixels to process
		//	Vec3c reducedDimensions = internals::getReducedDimensions(ri.dimensions(), dim);
		//	coord_t rowCount = reducedDimensions.x * reducedDimensions.y * reducedDimensions.z;

		//	size_t counter = 0;
		//	#pragma omp parallel if(!omp_in_parallel() && result.pixelCount() > PARALLELIZATION_THRESHOLD)
		//	{
		//		// Temporary buffer
		//		Image<internals::RiSet> centers(ri.dimension(dim));
		//		vector<internals::ActiveSpheresItem> activeSpheres;
		//		activeSpheres.reserve(40);

		//		// Determine start points of pixel rows and process each row
		//		#pragma omp for schedule(dynamic)
		//		for (coord_t n = 0; n < rowCount; n++)
		//		{
		//			Vec3c start = indexToCoords(n, reducedDimensions);

		//			// Make a copy of the current row as we update the row in the forward pass but need
		//			// the original data in the backward pass.
		//			Vec3c pos = start;
		//			for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
		//			{
		//				centers(x) = ri(pos);
		//				ri(pos).clear();
		//			}

		//			if (dim < ri.dimensionality() - 1)
		//			{
		//				// Left to right pass
		//				activeSpheres.clear();
		//				singlePass(centers, ri, start, dim, 1, activeSpheres);

		//				// Right to left pass
		//				activeSpheres.clear();
		//				singlePass(centers, ri, start, dim, -1, activeSpheres);

		//				// TODO: This has not been tested except in 'ridge' version.
		//				//if (redtOnly)
		//				//{
		//				//	// Remove from ri everything except the biggest ri value.
		//				//	pos = start;
		//				//	for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
		//				//	{
		//				//		RiSet& list = ri(pos);
		//				//		if (list.size() > 1)
		//				//		{
		//				//			RiItem item = *max_element(list.begin(), list.end(), ridgeTestComparerri);
		//				//			list.clear();
		//				//			list.push_back(item);
		//				//			list.shrink_to_fit();
		//				//		}
		//				//	}
		//				//}
		//			}
		//			else
		//			{
		//				// TODO: Test if it is better/faster to combine singlePass and singlePassFinal as in the ridge version

		//				// Left to right pass
		//				singlePassFinal(centers, result, start, dim, 1);

		//				// Right to left pass
		//				singlePassFinal(centers, result, start, dim, -1);
		//			}

		//			showThreadProgress(counter, rowCount, showProgressInfo);
		//		}
		//	}
		//}

		///*
		//Copies squared sphere radius data from original centers image to ri image.
		//*/
		//inline void prepare(const Image<int32_t>& centers2, Image<internals::RiSet>& ri)
		//{
		//	ri.ensureSize(centers2);

		//	for (coord_t n = 0; n < centers2.pixelCount(); n++)
		//	{
		//		int32_t R2 = centers2(n);
		//		if (R2 > 0)
		//		{
		//			ri(n).push_back({ R2, R2 });
		//		}
		//	}
		//}







// dimred 'diagonal' tests
//
//namespace internals
//{
//
//	class DirPPP
//	{
//	public:
//		Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims) const
//		{
//			return Vec3sc(x, y, z);
//		}
//
//		Vec3sc backNeighbourDir() const
//		{
//			return Vec3sc(-1, -1, -1);
//		}
//	};
//
//	class DirNPP
//	{
//	public:
//		Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims)
//		{
//			return Vec3sc((int32_t)dims.x - 1 - x, y, z);
//		}
//
//		Vec3sc backNeighbourDir() const
//		{
//			return Vec3sc(1, -1, -1);
//		}
//	};
//
//	class DirPNP
//	{
//	public:
//		Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims)
//		{
//			return Vec3sc(x, (int32_t)dims.y - 1 - y, z);
//		}
//
//		Vec3sc backNeighbourDir() const
//		{
//			return Vec3sc(-1, 1, -1);
//		}
//	};
//
//	class DirPPN
//	{
//	public:
//		Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims)
//		{
//			return Vec3sc(x, y, (int32_t)dims.z - 1 - z);
//		}
//
//		Vec3sc backNeighbourDir() const
//		{
//			return Vec3sc(-1, -1, 1);
//		}
//	};
//
//
//
//	class DirNNP
//	{
//	public:
//		Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims)
//		{
//			return Vec3sc((int32_t)dims.x - 1 - x, (int32_t)dims.y - 1 - y, z);
//		}
//
//		Vec3sc backNeighbourDir() const
//		{
//			return Vec3sc(1, 1, -1);
//		}
//	};
//
//	class DirPNN
//	{
//	public:
//		Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims)
//		{
//			return Vec3sc(x, (int32_t)dims.y - 1 - y, (int32_t)dims.z - 1 - z);
//		}
//
//		Vec3sc backNeighbourDir() const
//		{
//			return Vec3sc(-1, 1, 1);
//		}
//	};
//
//	class DirNPN
//	{
//	public:
//		Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims)
//		{
//			return Vec3sc((int32_t)dims.x - 1 - x, y, (int32_t)dims.z - 1 - z);
//		}
//
//		Vec3sc backNeighbourDir() const
//		{
//			return Vec3sc(1, -1, 1);
//		}
//	};
//
//
//	class DirNNN
//	{
//	public:
//		Vec3sc operator()(int32_t x, int32_t y, int32_t z, const Vec3c& dims)
//		{
//			return Vec3sc((int32_t)dims.x - 1 - x, (int32_t)dims.y - 1 - y, (int32_t)dims.z - 1 - z);
//		}
//
//		Vec3sc backNeighbourDir() const
//		{
//			return Vec3sc(1, 1, 1);
//		}
//	};
//
//
//
//
//	/**
//	Iterates image diagonally.
//	@param dims Dimensions of the image.
//	@param dir Functor that transforms initial image coordinates such that the iteration direction is the desired one. Select one from DirXXX classes.
//	@param functions Object that defines three functions.
//
//	functions.initThread(int threadIndex)
//	is called from each thread before processing.
//
//	functions.processPoint(const Vec3sc& position, const Vec3sc& backwardsDirection, int32_t iterationLevel).
//	if called from multiple threads simultaneously to process each point.
//
//	functions.levelComplete(int32_t start)
//	is called from the main thread after each diagonal level has been processed.
//	*/
//	template<typename dirfunc, typename functionset> void diagonalIteration(const Vec3c& dims, dirfunc dir, functionset& functions, bool showProgressInfo = true)
//	{
//		int32_t w = (int32_t)dims.x;
//		int32_t h = (int32_t)dims.y;
//		int32_t d = (int32_t)dims.z;
//
//		int32_t maxLevel = w + h + d - 2;
//		size_t counter = 0;
//		for (int32_t start = 0; start < maxLevel; start++)
//		{
//#pragma omp parallel if(!omp_in_parallel())
//			{
//#pragma omp single
//				{
//					functions.initThreads(omp_get_num_threads());
//				}
//
//#pragma omp for schedule(dynamic)
//				for (int32_t z = 0; z < d; z++)
//				{
//					int32_t x = start - z;
//					int32_t y = 0;
//
//					if (x >= w)
//					{
//						y += x - (w - 1);
//						x = (w - 1);
//					}
//
//					while (x >= 0 && y < h)
//					{
//						//if (!result.isInImage(x, y, z))
//						//	throw logic_error("Invalid program.");
//
//						functions.processPoint(dir(x, y, z, dims), dir.backNeighbourDir(), start, omp_get_thread_num());
//
//						x--;
//						y++;
//					}
//				}
//
//			}
//
//			functions.levelComplete(start);
//
//			showThreadProgress(counter, maxLevel, showProgressInfo);
//		}
//	}
//
//	/**
//	@param p Point we are currently processing.
//	@param dim Dimension we are currently processing.
//	@param rilist At input list of center points from previous dimension, at output center points for next dimension.
//	@param activeSpheres Active spheres list from previous iteration.
//	@param dimensionality Dimensionality of the image.
//	*/
//	void step(int32_t x, size_t dim, internals::RiSet& rilist, vector<internals::ActiveSpheresItem>& activeSpheres, size_t dimensionality,
//		vector<internals::ActiveSpheresItem> activeSpheresTmp, vector<internals::ActiveSpheresItem> Ctmp, internals::RiSet resultTmp, internals::RiSet resultTmp2)
//
//	{
//
//
//		activeSpheresTmp.clear();
//		Ctmp.clear();
//		resultTmp.clear();
//
//		// If there is one or more sphere centers at the current location, add them to the set of active spheres.
//		// The rilist list is sorted by R and so is activeSpheres list.
//		// Instead of finding place for each item, convert each item in rilist list to activeSpheres format (add x coordinate)
//		// and then merge the two sorted lists to construct the new activeSpheres list.
//		for (auto& item : rilist)
//			Ctmp.push_back({ item.R2, item.ri2, x });
//
//		activeSpheresTmp.resize(rilist.size() + activeSpheres.size());
//		merge(Ctmp.begin(), Ctmp.end(), activeSpheres.begin(), activeSpheres.end(), activeSpheresTmp.begin(), internals::activeSpheresSorter);
//		swap(activeSpheres, activeSpheresTmp);
//
//
//		// Iterate through all active spheres and calculate radius for next dimension.
//		auto it = activeSpheres.begin();
//		while (it != activeSpheres.end())
//		{
//			int32_t Rorig2 = it->R2;
//			int32_t R2 = it->ri2;
//			int32_t cx = it->xi;
//
//			int32_t dx = abs(x - cx);
//
//			// Calculate new ri^2
//			int32_t rn2 = R2 - dx * dx;
//			if (rn2 > 0)
//			{
//				// Insert ry2 into the list (TODO: Don't insert if the item is going to be removed anyway?)
//				resultTmp.push_back({ Rorig2, rn2 });
//
//				// Remove possible duplicate Rorig and replace it by (Rorig, max(duplicates))
//				if (resultTmp.size() > 1 && Rorig2 == resultTmp[resultTmp.size() - 2].R2)
//				{
//					if (rn2 > resultTmp[resultTmp.size() - 2].ri2)
//						resultTmp[resultTmp.size() - 2].ri2 = rn2;
//
//					resultTmp.pop_back();
//				}
//
//				it++;
//			}
//			else
//			{
//				// ry is non-positive, i.e. dx >= R
//				// This sphere is not active anymore, so remove it from the list of active spheres.
//				it = activeSpheres.erase(it);
//			}
//		}
//
//		// Rebuild ri list from resultTmp. Don't include those items that are hidden by other items.
//		rilist.clear();
//		if (dim < dimensionality - 1)
//		{
//			if (resultTmp.size() > 0)
//			{
//				// Add ri from the previous pass to the ri list and save the result to resultTmp2.
//				// Note that both resultTmp2 and rilist are sorted so we can just merge them.
//				//internals::RiSet& rilist = ri(p);
//				resultTmp2.clear();
//				resultTmp2.resize(rilist.size() + resultTmp.size());
//				merge(rilist.begin(), rilist.end(), resultTmp.begin(), resultTmp.end(), resultTmp2.begin(), internals::origRSorter);
//
//
//				// NOTE: This is the basic version, the memory saving version is below.
//				//// Linear time algorithm for finding relevant ri (those not hidden by other items).
//				//// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//				//// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//				//rilist.clear();
//				//rilist.push_back(resultTmp2[0]);
//				//for (size_t n = 1; n < resultTmp2.size(); n++)
//				//{
//				//	int32_t currri2 = rilist[rilist.size() - 1].ri2;
//				//	int32_t newri2 = resultTmp2[n].ri2;
//				//	if (newri2 > currri2)
//				//	{
//				//		rilist.push_back(resultTmp2[n]);
//				//	}
//				//}
//
//				// Linear time algorithm for finding relevant ri (those not hidden by other items).
//				// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//				// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//				//rilist.clear();
//				rilist.push_back(resultTmp2[0]);
//				//if (dim < ri.dimensionality() - 1) // In the last dimension only the item with maximum R is needed, but we have separate function for that
//				//{
//				for (size_t n = 1; n < resultTmp2.size(); n++)
//				{
//					int32_t currri2 = rilist[rilist.size() - 1].ri2;
//					int32_t newri2 = resultTmp2[n].ri2;
//
//					if (newri2 > currri2) // This is the basic condition that works always (but may include unnecessary items in the rilist)
//					{
//						if (dim == dimensionality - 2)
//						{
//							// In the second last dimension only really visible spans are needed as there's no next dimension whose ri we would calculate based on the spans.
//							if (largestIntWhoseSquareIsLessThan(newri2) > largestIntWhoseSquareIsLessThan(currri2))
//								rilist.push_back(resultTmp2[n]);
//						}
//						else if (dim == dimensionality - 3)
//						{
//							// In the third last dimension we know that only those spans are required that produce visible circles in the output.
//							if (!internals::doesDiscretizedCircle1FitInto2Cached(newri2, currri2))
//								rilist.push_back(resultTmp2[n]);
//						}
//						else
//						{
//							// Here we could insert test if discretized spheres fit into each other etc. etc.
//							rilist.push_back(resultTmp2[n]);
//						}
//					}
//				}
//				//}
//
//				// Make sure we don't use extraneous memory by storing empty items in each pixel of the ri image.
//				rilist.shrink_to_fit();
//			}
//		}
//	}
//
//
//	vector<internals::ActiveSpheresItem>& getSphereList(Image<vector<internals::ActiveSpheresItem> >& activeSpheres, Vec3sc p, size_t dim)
//	{
//		p[dim] = 0;
//		return activeSpheres(p);
//	}
//
//
//	//inline bool origRSorter2(const RiItem& a, const RiItem& b)
//	//{
//	//	if (a.R2 != b.R2)
//	//		return a.R2 < b.R2;
//	//	return a.ri2 < b.ri2;
//	//}
//
//
//	struct DiagStepper
//	{
//		const Image<int32_t>& dmap2;
//		Image<int32_t>& result;
//		bool showProgressInfo;
//
//		Image<vector<internals::ActiveSpheresItem> >& asXImg;
//		Image<vector<internals::ActiveSpheresItem> >& asYImg;
//		Image<vector<internals::ActiveSpheresItem> >& asZImg;
//
//
//		DiagStepper(const Image<int32_t>& dmap2, Image<int32_t>& result,
//			Image<vector<internals::ActiveSpheresItem> >& asX, Image<vector<internals::ActiveSpheresItem> >& asY, Image<vector<internals::ActiveSpheresItem> >& asZ,
//			bool showProgressInfo) :
//			dmap2(dmap2),
//			result(result),
//			asXImg(asX),
//			asYImg(asY),
//			asZImg(asZ),
//			showProgressInfo(showProgressInfo)
//		{
//			asXImg.ensureSize(1, dmap2.height(), dmap2.depth());
//			asYImg.ensureSize(dmap2.width(), 1, dmap2.depth());
//			asZImg.ensureSize(dmap2.width(), dmap2.height(), 1);
//		}
//
//		void reset()
//		{
//			setValue(asXImg, vector<internals::ActiveSpheresItem>());
//			setValue(asYImg, vector<internals::ActiveSpheresItem>());
//			setValue(asZImg, vector<internals::ActiveSpheresItem>());
//		}
//
//
//		/**
//		Adds padding to list items. Does not have any effect really... but might be useful elsewhere.
//		*/
//		template<typename T> struct padder
//		{
//			T item;
//			uint8_t padding[256 - sizeof(T)];
//
//			padder(T& item) : item(item)
//			{
//
//			}
//
//			operator T()
//			{
//				return item;
//			}
//		};
//
//		/**
//		Temporary variables for each thread.
//		*/
//		vector<padder<vector<internals::ActiveSpheresItem> > > activeSpheresTmp;
//		vector<padder<vector<internals::ActiveSpheresItem> > > Ctmp;
//		vector<padder<internals::RiSet> > resultTmp;
//		vector<padder<internals::RiSet> > resultTmp2;
//		vector<padder<internals::RiSet> > Cs;
//
//		void initThreads(int threadCount)
//		{
//			if (activeSpheresTmp.size() < threadCount)
//			{
//				for (int n = 0; n < threadCount; n++)
//				{
//					activeSpheresTmp.push_back(vector<internals::ActiveSpheresItem>());
//					activeSpheresTmp.reserve(40);
//					Ctmp.push_back(vector<internals::ActiveSpheresItem>());
//					Ctmp.reserve(40);
//					resultTmp.push_back(internals::RiSet());
//					resultTmp.reserve(40);
//					resultTmp2.push_back(internals::RiSet());
//					resultTmp2.reserve(40);
//					Cs.push_back(internals::RiSet());
//					Cs.reserve(40);
//				}
//			}
//		}
//
//		void processPoint(const Vec3sc& p, const Vec3sc& backDir, int32_t level, int threadIndex)
//		{
//			auto& C = Cs[threadIndex].item;
//			C.clear();
//
//			int32_t R2 = dmap2(p);
//			if (R2 > 0)
//				C.push_back({ R2, R2 });
//
//			auto& tmp1 = activeSpheresTmp[threadIndex];
//			auto& tmp2 = Ctmp[threadIndex];
//			auto& tmp3 = resultTmp[threadIndex];
//			auto& tmp4 = resultTmp2[threadIndex];
//
//			vector<internals::ActiveSpheresItem>& asX = internals::getSphereList(asXImg, p, 0);
//			step(p[0], 0, C, asX, dmap2.dimensionality(), tmp1, tmp2, tmp3, tmp4);
//
//			vector<internals::ActiveSpheresItem>& asY = internals::getSphereList(asYImg, p, 1);
//			step(p[1], 1, C, asY, dmap2.dimensionality(), tmp1, tmp2, tmp3, tmp4);
//
//			vector<internals::ActiveSpheresItem>& asZ = internals::getSphereList(asZImg, p, 2);
//			step(p[2], 2, C, asZ, dmap2.dimensionality(), tmp1, tmp2, tmp3, tmp4);
//
//			auto it = max_element(asZ.begin(), asZ.end(), internals::activeSpheresSorter2);
//			if (it != asZ.end())
//			{
//				int32_t R2 = it->R2;
//				if (R2 > result(p))
//					result(p) = R2;
//			}
//
//			//auto it = max_element(C.begin(), C.end(), internals::origRSorter2);
//			//if (it != C.end())
//			//{
//			//	int32_t R2 = it->R2;
//			//	if (R2 > result(p))
//			//		result(p) = R2;
//			//}
//
//		}
//
//		void levelComplete(int32_t level)
//		{
//
//		}
//	};
//}
//
//
//namespace dimreddiagonal
//{
//
//	/**
//	Calculates local thickness without storing ri values for each dimension.
//	Uses diagonal iteration.
//	The problem in this version is that it does the same processing many times as
//	it can only do 1/8:th of sphere in one pass through the image.
//	I.e., first it propagates diagonally from (0, 0, 0) to (w, h, d) and stores the result.
//	This process consists (implicitly) of iterating over x, y, and z once.
//	Then, it propagates diagonally from (w, 0, 0) to (0, h, d), and stores the result.
//	This includes iterating over -x, y, and z once.
//	Now we have already iterated over y and z in the positive direction two times!
//	*/
//	void thickmap2(const Image<int32_t>& dmap2, Image<int32_t>& result, bool showProgressInfo)
//	{
//		internals::buildCircleLookup(max(dmap2));
//
//		result.ensureSize(dmap2.dimensions());
//		setValue(result, 0);
//
//		Image<vector<internals::ActiveSpheresItem> > asXImg;
//		Image<vector<internals::ActiveSpheresItem> > asYImg;
//		Image<vector<internals::ActiveSpheresItem> > asZImg;
//
//		internals::DiagStepper stepper(dmap2, result, asXImg, asYImg, asZImg, showProgressInfo);
//
//		internals::diagonalIteration(result.dimensions(), internals::DirPPP(),
//			stepper,
//			showProgressInfo);
//
//		stepper.reset();
//		internals::diagonalIteration(result.dimensions(), internals::DirNPP(),
//			stepper,
//			showProgressInfo);
//
//
//		stepper.reset();
//		internals::diagonalIteration(result.dimensions(), internals::DirPNP(),
//			stepper,
//			showProgressInfo);
//
//		stepper.reset();
//		internals::diagonalIteration(result.dimensions(), internals::DirNNP(),
//			stepper,
//			showProgressInfo);
//
//
//		internals::diagonalIteration(result.dimensions(), internals::DirPPN(),
//			stepper,
//			showProgressInfo);
//
//		stepper.reset();
//		internals::diagonalIteration(result.dimensions(), internals::DirNPN(),
//			stepper,
//			showProgressInfo);
//
//
//		stepper.reset();
//		internals::diagonalIteration(result.dimensions(), internals::DirPNN(),
//			stepper,
//			showProgressInfo);
//
//		stepper.reset();
//		internals::diagonalIteration(result.dimensions(), internals::DirNNN(),
//			stepper,
//			showProgressInfo);
//	}
//}




// Optimized, sorts in inverse order and does not use bit mask
//namespace internals
//{
//
//	inline void fillSpheres(const vector<tuple<int32_t, coord_t> >& points, coord_t starti, coord_t endi, Image<int32_t>& tmap)
//	{
//		int32_t R2 = get<0>(points[starti]);
//		int32_t Rint = (int32_t)ceil(sqrt(R2));
//
//		// Create mask
//		coord_t size = 2 * Rint + 1;
//		Image<int32_t> xMask(size, size, 1, -1);
//#pragma omp parallel for if(xMask.pixelCount() > PARALLELIZATION_THRESHOLD)
//		for (int32_t z = 0; z < xMask.height(); z++)
//		{
//			int32_t dz = z - Rint;
//			for (int32_t y = 0; y < xMask.width(); y++)
//			{
//				int32_t dy = y - Rint;
//				int32_t dx2 = R2 - dy * dy - dz * dz;
//				if (dx2 >= 0)
//				{
//					xMask(y, z) = largestIntWhoseSquareIsLessThan(dx2);
//				}
//			}
//		}
//
//		// Fill all spheres
//		Vec3c maxCoords = tmap.dimensions() - Vec3c(1, 1, 1);
//		bool centersParallel = endi - starti >= omp_get_max_threads();
//#pragma omp parallel for if(centersParallel)
//		for (coord_t n = starti; n < endi; n++)
//		{
//			Vec3c pos = tmap.getCoords(get<1>(points[n]));
//
//			Vec3c start0 = pos - Vec3c(Rint, Rint, Rint);
//			Vec3c end = pos + Vec3c(Rint, Rint, Rint);
//
//			Vec3c start = start0;
//			clamp(start, Vec3c(), maxCoords);
//			clamp(end, Vec3c(), maxCoords);
//
//#pragma omp parallel for if(!omp_in_parallel() && (end.z - start.z + 1) > omp_get_max_threads())
//			for (coord_t z = start.z; z <= end.z; z++)
//			{
//				for (coord_t y = start.y; y <= end.y; y++)
//				{
//					coord_t xr = xMask(y - start0.y, z - start0.z);
//
//					if (xr >= 0)
//					{
//						coord_t xStart = pos.x - xr;
//						if (xStart < 0)
//							xStart = 0;
//
//						coord_t xEnd = pos.x + xr;
//						if (xEnd > maxCoords.x)
//							xEnd = maxCoords.x;
//
//						coord_t startInd = (coord_t)tmap.getLinearIndex(xStart, y, z);
//						coord_t endInd = startInd + (xEnd - xStart);
//						for (coord_t ind = startInd; ind <= endInd; ind++)
//							tmap(ind) = R2;
//					}
//				}
//			}
//		}
//	}
//}
//
//namespace optimizedinverse
//{
//	/**
//	Calculates squared local radius from squared distance map
//	Uses standard Hildebrand & Ruegsegger algorithm (optimized version).
//	Plots maximal spheres corresponding to squared distance map, larger distance values replacing smaller ones.
//	*/
//	void thickmap2(Image<int32_t>& dmap2, double* extraBytes, bool showProgressInfo)
//	{
//		// In order to parallelize by sorting method, first find radii and locations (linear indices to save some memory) of non-zero pixels.
//		vector<tuple<int32_t, coord_t> > centers;
//		for (coord_t n = 0; n < dmap2.pixelCount(); n++)
//		{
//			if (dmap2(n) > 0)
//				centers.push_back(make_tuple(dmap2(n), n));
//		}
//
//		// Sort according to sphere diameter
//#if defined(_WIN32)
//		sort(execution::par_unseq, centers.begin(), centers.end(), internals::ridgeSorter2);
//#elif defined(__linux__)
//		__gnu_parallel::sort(centers.begin(), centers.end(), internals::ridgeSorter2);
//#endif
//
//		if (extraBytes)
//		{
//			(*extraBytes) = (double)(centers.size() * (sizeof(int32_t) + sizeof(coord_t)));
//		}
//
//		coord_t previ = 0;
//		for (coord_t n = 0; n < (coord_t)centers.size(); n++)
//		{
//			int32_t prevR = get<0>(centers[previ]);
//			int32_t r = get<0>(centers[n]);
//
//			if (r != prevR)
//			{
//				//cout << "Filling range " << previ << " - " << n << endl;
//				//cout << "Filling " << (n - previ) << " points" << endl;
//
//				// Draw spheres from previ to n-1
//				internals::fillSpheres(centers, previ, n, dmap2);
//				previ = n;
//			}
//
//			showProgress(n, centers.size(), showProgressInfo);
//		}
//
//		// Draw the final span of spheres from previ to end of ridgePoints list.
//		internals::fillSpheres(centers, previ, (coord_t)centers.size(), dmap2);
//	}
//}












// dimred super version that uses separate ri storage item struct optimized for space, and
// a array implementation for each pixel that uses minimal amount of space (pointer + size = 8 + 2 bytes = 10 bytes).

//namespace internals
//{
//#pragma pack(push, 1)
//
//	/**
//	Data entry that must be saved for each pixel for each sphere that must be considered in that pixel.
//	This data structure should be as small as possible.
//	*/
//	struct RiStorageItem
//	{
//		int16_t srcX, srcY;
//
//		// DEBUG
//		//int32_t R2, ri2;
//	};
//
//	/**
//	Simple array that consumes as little memory as possible.
//	*/
//	template<typename item_t, typename count_t = uint16_t> class SimpleVector
//	{
//	private:
//		unique_ptr<item_t[]> pData;
//		count_t count;
//
//	public:
//		SimpleVector(const vector<item_t>& items) : pData(nullptr), count(0)
//		{
//			setValues(items);
//		}
//
//		SimpleVector() : pData(nullptr), count(0)
//		{
//		}
//
//		SimpleVector(const SimpleVector& r) : pData(nullptr), count(r.count)
//		{
//			if (count > 0)
//			{
//				pData = make_unique<item_t[]>(count);
//				for (count_t n = 0; n < count; n++)
//					pData[n] = r.pData[n];
//			}
//		}
//
//		SimpleVector& operator=(SimpleVector other)
//		{
//			swap(count, other.count);
//			swap(pData, other.pData);
//
//			return *this;
//		}
//
//		/**
//		Sets values in this array to those in the items vector.
//		*/
//		void setValues(const vector<item_t>& items)
//		{
//			//delete[] pData;
//
//			if (items.size() <= 0)
//			{
//				pData = nullptr;
//				count = 0;
//			}
//			else
//			{
//				if (items.size() >= numeric_limits<count_t>::max())
//					throw invalid_argument("The source array has too many elements.");
//
//				count = (count_t)items.size();
//				//pData = new item_t[count];
//				pData = make_unique<item_t[]>(count);
//				for (count_t n = 0; n < count; n++)
//					pData[n] = items[n];
//			}
//		}
//
//		/**
//		Populates this array with single item.
//		*/
//		void setValues(item_t& item)
//		{
//			//delete[] pData;
//			count = 1;
//			//pData = new item_t[1];
//			pData = make_unique<item_t[]>(1);
//			pData[0] = item;
//		}
//
//		/**
//		Gets count of items in this array.
//		*/
//		count_t size() const
//		{
//			return count;
//		}
//
//		item_t& operator[](count_t index)
//		{
//			return pData[index];
//		}
//
//		const item_t& operator[](count_t index) const
//		{
//			return pData[index];
//		}
//	};
//
//
//	/*
//	ri image pixel type. Stores original radius and current ri for each sphere.
//	*/
//	//typedef vector<RiStorageItem> RiStorageSet;
//	typedef SimpleVector<RiStorageItem> RiStorageSet;
//
//
//	struct RiSuperItem
//	{
//		/*
//		Squared radius of the sphere that initiated creation of this extent.
//		*/
//		int32_t R2;
//		/*
//		Squared extent of the sphere in the next dimension.
//		*/
//		int32_t ri2;
//
//
//		int16_t srcX;
//		int16_t srcY;
//	};
//
//	struct ActiveSpheresSuperItem
//	{
//		/*
//		Squared radius of the sphere that initiated creation of this extent.
//		*/
//		int32_t R2;
//		/*
//		Squared extent of the sphere in the current dimension.
//		*/
//		int32_t ri2;
//		/*
//		Center point of the sphere in the current dimension.
//		*/
//		int32_t xi;
//
//		int16_t srcX;
//		int16_t srcY;
//	};
//
//#pragma pack(pop)
//
//	/*
//	ri image pixel type. Stores original radius and current ri for each sphere.
//	*/
//	typedef vector<RiSuperItem> RiSuperSet;
//
//	/*
//	Comparison functions for comparing RiItem and ActiveSpheresItem.
//	*/
//	inline bool origRSuperSorter(const RiSuperItem& a, const RiSuperItem& b)
//	{
//		if (a.R2 != b.R2)
//			return a.R2 > b.R2;
//		return a.ri2 > b.ri2;
//	}
//
//	inline bool activeSpheresSuperSorter(const ActiveSpheresSuperItem& a, const ActiveSpheresSuperItem& b)
//	{
//		if (a.R2 != b.R2)
//			return a.R2 > b.R2;
//		return a.ri2 > b.ri2;
//	}
//
//	/*
//	Sorts active spheres in inverse order
//	*/
//	inline bool activeSpheresSuperSorter2(const ActiveSpheresSuperItem& a, const ActiveSpheresSuperItem& b)
//	{
//		if (a.R2 != b.R2)
//			return a.R2 < b.R2;
//		return a.ri2 < b.ri2;
//	}
//
//	/*
//	@param centers Image containing only the row to be processed before first call to this method.
//	@param ri Full ri image that will be updated.
//	@param rowStart Start point of the pixel row to be processed.
//	@param dim Dimension that we are processing.
//	@param step +1 or -1 to indicate the direction of the pass.
//	@param activeSpheres List containing initial active spheres for the row. At exit, contains list of active spheres after processing the row.
//	*/
//	inline void singlePassSuper(const vector<internals::RiSuperSet>& centers, vector<internals::RiSuperSet>& ri, const Vec3c& dimensions, const Vec3c& rowStart, coord_t dim, coord_t step, vector<internals::ActiveSpheresSuperItem>& activeSpheres)
//	{
//		// Stores the spheres that have been encountered and that have not been passed yet.
//		// Stores the center point, original radius, and ri.
//		//vector<internals::ActiveSpheresItem> activeSpheres;
//		vector<internals::ActiveSpheresSuperItem> activeSpheresTmp;
//		vector<internals::ActiveSpheresSuperItem> Ctmp;
//
//		internals::RiSuperSet resultTmp;
//		internals::RiSuperSet resultTmp2;
//
//		//activeSpheres.reserve(40);
//		activeSpheresTmp.reserve(40);
//		Ctmp.reserve(40);
//		resultTmp.reserve(40);
//		resultTmp2.reserve(40);
//
//		// Set start point to the start or to the end of the current row.
//		Vec3c p = rowStart;
//		if (step < 0)
//			p[dim] += (dimensions[dim] - 1);
//
//		coord_t dimensionality = getDimensionality(dimensions);
//
//		for (coord_t i = 0; i < dimensions[dim]; i++, p[dim] += step)
//		{
//			int32_t x = (int32_t)p[dim];
//
//			// If there is one or more sphere centers at the current location, add them to the set of active spheres.
//			const internals::RiSuperSet& C = centers[x];
//
//			// The C list is sorted by R and so is activeSpheres list.
//			// Instead of finding place for each item, convert each item in C list to activeSpheres format (add x coordinate)
//			// and then merge the two sorted lists to construct the new activeSpheres list.
//			Ctmp.clear();
//			for (auto& item : C)
//				Ctmp.push_back({ item.R2, item.ri2, x, item.srcX, item.srcY });
//
//			activeSpheresTmp.clear();
//			activeSpheresTmp.resize(C.size() + activeSpheres.size());
//			merge(Ctmp.begin(), Ctmp.end(), activeSpheres.begin(), activeSpheres.end(), activeSpheresTmp.begin(), internals::activeSpheresSuperSorter);
//			swap(activeSpheres, activeSpheresTmp);
//
//
//			// Iterate through all active spheres and calculate radius for next dimension.
//			resultTmp.clear();
//			auto it = activeSpheres.begin();
//			while (it != activeSpheres.end())
//			{
//				int32_t Rorig2 = it->R2;
//				int32_t R2 = it->ri2;
//				int32_t cx = it->xi;
//
//				int16_t origcx = it->srcX;
//				int16_t origcy = it->srcY;
//
//				int32_t dx = abs(x - cx);
//
//				// Calculate new ri^2
//				int32_t rn2 = R2 - dx * dx;
//				if (rn2 > 0)
//				{
//					// Insert ry2 into the list, but don't insert duplicates.
//					if (resultTmp.size() > 0 && Rorig2 == resultTmp[resultTmp.size() - 1].R2)
//					{
//						// This is a duplicate R2 entry. Use the entry with the larger ri.
//						if (resultTmp[resultTmp.size() - 1].ri2 < rn2)
//							resultTmp[resultTmp.size() - 1] = { Rorig2, rn2, origcx, origcy };
//					}
//					else
//					{
//						resultTmp.push_back({ Rorig2, rn2, origcx, origcy });
//					}
//
//					it++;
//				}
//				else
//				{
//					// ry is non-positive, i.e. dx >= R
//					// This sphere is not active anymore, so remove it from the list of active spheres.
//					it = activeSpheres.erase(it);
//				}
//			}
//
//			// Rebuild ri list from resultTmp. Don't include those items that are hidden by other items.
//			if (resultTmp.size() > 0)
//			{
//				// Add ri from the previous pass to the ri list and save the result to resultTmp2.
//				// Note that both resultTmp2 and rilist are sorted so we can just merge them.
//				internals::RiSuperSet& rilist = ri[x];
//				resultTmp2.clear();
//				resultTmp2.resize(rilist.size() + resultTmp.size());
//				merge(rilist.begin(), rilist.end(), resultTmp.begin(), resultTmp.end(), resultTmp2.begin(), internals::origRSuperSorter);
//
//
//				// NOTE: This is the basic version, the memory saving version is below.
//				//// Linear time algorithm for finding relevant ri (those not hidden by other items).
//				//// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//				//// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//				//rilist.clear();
//				//rilist.push_back(resultTmp2[0]);
//				//for (size_t n = 1; n < resultTmp2.size(); n++)
//				//{
//				//	int32_t currri2 = rilist[rilist.size() - 1].ri2;
//				//	int32_t newri2 = resultTmp2[n].ri2;
//				//	if (newri2 > currri2)
//				//	{
//				//		rilist.push_back(resultTmp2[n]);
//				//	}
//				//}
//
//				// Linear time algorithm for finding relevant ri (those not hidden by other items).
//				// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//				// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//				rilist.clear();
//				rilist.push_back(resultTmp2[0]);
//				//if (dim < ri.dimensionality() - 1) // In the last dimension only the item with maximum R is needed, but we have separate function for that
//				//{
//				for (size_t n = 1; n < resultTmp2.size(); n++)
//				{
//					int32_t currri2 = rilist[rilist.size() - 1].ri2;
//					int32_t newri2 = resultTmp2[n].ri2;
//
//					if (newri2 > currri2) // This is the basic condition that works always (but may include unnecessary items in the rilist)
//					{
//						if (dim == dimensionality - 2) // 3 - 2 == 1 == 2nd dimension
//						{
//							// In the second last dimension only really visible spans are needed as there's no next dimension whose ri we would calculate based on the spans.
//							if (largestIntWhoseSquareIsLessThan(newri2) > largestIntWhoseSquareIsLessThan(currri2))
//								rilist.push_back(resultTmp2[n]);
//						}
//						else if (dim == dimensionality - 3) // 3 - 3 == 0 == 1st dimension
//						{
//							// In the third last dimension we know that only those spans are required that produce visible circles in the output.
//							if (!internals::doesDiscretizedCircle1FitInto2Cached(newri2, currri2))
//								rilist.push_back(resultTmp2[n]);
//						}
//						else
//						{
//							// Here we could insert test if discretized spheres fit into each other etc. etc.
//							rilist.push_back(resultTmp2[n]);
//						}
//					}
//				}
//				//}
//
//				// Make sure we don't use extraneous memory by storing empty items in each pixel of the ri image.
//				rilist.shrink_to_fit();
//
//				//debugCheckRiSetRedundancy(rilist, dim);
//			}
//		}
//	}
//
//	/*
//	Same than singlePass but generates result image instead of ri image, and is thus a little bit faster.
//	@param centers Image containing only the row to be processed before first call to this method.
//	@param Result image. Pixels must be set to zero before first call to this method.
//	@param rowStart Start point of the pixel row to be processed.
//	@param dim Dimension that we are processing.
//	@param step +1 or -1 to indicate the direction of the pass.
//	@param initialActiveSpheres Initial active spheres list. Contains final active spheres list at output. Set to nullptr to assume empty list.
//	*/
//	inline void singlePassFinalSuper(const vector<internals::RiSuperSet>& centers, Image<int32_t>& result, const Vec3c& rowStart, size_t dim, coord_t step, vector<internals::ActiveSpheresSuperItem>* initialActiveSpheres)
//	{
//		// Stores the spheres that have been encountered and that have not been passed yet.
//		// Stores the sphere with the largest R at the top of the priority queue.
//		priority_queue<internals::ActiveSpheresSuperItem, vector<internals::ActiveSpheresSuperItem>, std::function<decltype(internals::activeSpheresSuperSorter2)> > activeSpheres(internals::activeSpheresSuperSorter2);
//
//		if (initialActiveSpheres)
//		{
//			for (const auto& item : *initialActiveSpheres)
//				activeSpheres.push(item);
//		}
//
//
//		// Set start point to the start or end of the current row.
//		Vec3c p = rowStart;
//		if (step < 0)
//			p[dim] += (result.dimension(dim) - 1);
//
//		for (coord_t i = 0; i < result.dimension(dim); i++, p[dim] += step)
//		{
//			int32_t x = (int32_t)p[dim];
//
//			// If there is one or more sphere centers at the current location, add them to the set of active spheres.
//			const internals::RiSuperSet& C = centers[x];
//			for (auto& item : C)
//			{
//				activeSpheres.push({ item.R2, item.ri2, x });
//			}
//
//			while (!activeSpheres.empty())
//			{
//				auto& item = activeSpheres.top();
//
//				int32_t Rorig2 = item.R2;
//				int32_t R2 = item.ri2;
//				int32_t cx = item.xi;
//				int32_t dx = abs(x - cx);
//
//				// Calculate new ri^2
//				int32_t rn2 = R2 - dx * dx;
//				if (rn2 > 0)
//				{
//					// Note that previous pass may have assigned larger value to the output.
//					if (Rorig2 > result(p))
//						result(p) = Rorig2;
//					break;
//				}
//				else
//				{
//					// ry is non-positive, i.e. dx >= R
//					// This sphere is not active anymore, so remove it from the list of active spheres.
//					activeSpheres.pop();
//				}
//
//			}
//
//		}
//
//		// Put stuff back to initialActiveSpheres list as it will be the initial active spheres list for possible adjacent block.
//		if (initialActiveSpheres)
//		{
//			initialActiveSpheres->clear();
//			while (!activeSpheres.empty())
//			{
//				initialActiveSpheres->push_back(activeSpheres.top());
//				activeSpheres.pop();
//			}
//			reverse(initialActiveSpheres->begin(), initialActiveSpheres->end());
//			// TODO: How to avoid sorting here?
//			sort(initialActiveSpheres->begin(), initialActiveSpheres->end(), internals::activeSpheresSuperSorter);
//			//if (!is_sorted(initialActiveSpheres.begin(), initialActiveSpheres.end(), internals::activeSpheresSorter))
//			//{
//			//	cout << "Not sorted" << endl;
//			//	sort(initialActiveSpheres.begin(), initialActiveSpheres.end(), internals::activeSpheresSorter);
//			//}
//		}
//	}
//
//
//	RiSuperItem toRiItem(const RiStorageItem& si, const Vec3c& p, const Image<int32_t>& dmap2)
//	{
//		RiSuperItem out;
//		int32_t R2 = dmap2(si.srcX, si.srcY, p.z);
//		int32_t dx = (int32_t)p.x - si.srcX;
//		int32_t dy = (int32_t)p.y - si.srcY;
//		int32_t ri2 = R2 - dx * dx - dy * dy;
//
//		//if (R2 != si.R2 || ri2 != si.ri2)
//		//{
//		//	cout << "Error" << endl;
//		//}
//
//		out.R2 = R2;
//		out.ri2 = ri2;
//
//		// DEBUG
//		//out.R2 = si.R2;
//		//out.ri2 = si.ri2;
//
//
//		out.srcX = si.srcX;
//		out.srcY = si.srcY;
//		return out;
//	}
//
//	RiStorageItem toRiStorageItem(const RiSuperItem& item)
//	{
//		RiStorageItem out;
//		out.srcX = item.srcX;
//		out.srcY = item.srcY;
//
//		// DEBUG
//		//out.R2 = item.R2;
//		//out.ri2 = item.ri2;
//
//		return out;
//	}
//
//	void toRiSet(const RiStorageSet& in, RiSuperSet& out, const Vec3c& p, const Image<int32_t>& dmap2)
//	{
//		out.clear();
//		//for (auto& item : in)
//		//	out.push_back(toRiItem(item, p, dmap2));
//		for (uint16_t n = 0; n < in.size(); n++)
//			out.push_back(toRiItem(in[n], p, dmap2));
//	}
//
//	void toStorageSet(const RiSuperSet& in, RiStorageSet& out)
//	{
//		//out.clear();
//		//for (auto& item : in)
//		//	out.push_back(toRiStorageItem(item));
//		vector<RiStorageItem> temp;
//		for (auto& item : in)
//			temp.push_back(toRiStorageItem(item));
//		out.setValues(temp);
//	}
//
//	/*
//	Makes one pass over image in specific dimension and direction.
//	@param ri Image containing ri values from processing of previous dimension.
//	@param dim Dimension to process.
//	@param doForwardPass Set to true if forward pass should be made.
//	@param activeSpheresStart Initial active spheres list for each pixel row for forward pass. Set to nullptr to assume empty lists.
//	@param doBackwardPass Set to true if backward pass should be made.
//	@param activeSpheresEnd Initial active spheres list for each pixel row for backward pass. Set to nullptr to assume empty lists.
//	@param result Result image.
//	@param showProgressInfo Set to true to show progress indicator.
//	*/
//	inline void processDimensionSuper(Image<internals::RiStorageSet>& ri, size_t dim,
//		bool doForwardPass, Image<vector<internals::ActiveSpheresSuperItem> >* activeSpheresStart,
//		bool doBackwardPass, Image<vector<internals::ActiveSpheresSuperItem> >* activeSpheresEnd,
//		const Image<int32_t>& dmap2,
//		Image<int32_t>& result, bool showProgressInfo)
//	{
//		// Determine count of pixels to process
//		Vec3c reducedDimensions = internals::getReducedDimensions(ri.dimensions(), dim);
//		coord_t rowCount = reducedDimensions.x * reducedDimensions.y * reducedDimensions.z;
//
//		if (doForwardPass && activeSpheresStart)
//			activeSpheresStart->checkSize(reducedDimensions);
//
//		if (doBackwardPass && activeSpheresEnd)
//			activeSpheresEnd->checkSize(reducedDimensions);
//
//
//		bool isFinalPass = !(dim < ri.dimensionality() - 1);
//
//		size_t counter = 0;
//#pragma omp parallel if(!omp_in_parallel() && ri.pixelCount() > PARALLELIZATION_THRESHOLD)
//		{
//			// Temporary buffer
//			vector<internals::RiSuperSet> inRow(ri.dimension(dim));
//			vector<internals::RiSuperSet> outRow(ri.dimension(dim));
//			vector<internals::ActiveSpheresSuperItem> activeSpheres;
//			activeSpheres.reserve(40);
//
//			// Determine start points of pixel rows and process each row
//#pragma omp for schedule(dynamic)
//			for (coord_t n = 0; n < rowCount; n++)
//			{
//				Vec3c start = indexToCoords(n, reducedDimensions);
//
//				// Make a copy of the current row as we update the row in the forward pass but need
//				// the original data in the backward pass.
//				Vec3c pos = start;
//				for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
//				{
//					toRiSet(ri(pos), inRow[x], pos, dmap2);
//					outRow[x].clear();
//					//ri(pos).clear();
//				}
//
//
//				if (!isFinalPass)
//				{
//					if (doForwardPass)
//					{
//						// Use either temporary activeSpheres or element from activeSpheresStart
//						if (activeSpheresStart)
//						{
//							singlePassSuper(inRow, outRow, ri.dimensions(), start, dim, 1, (*activeSpheresStart)(start));
//						}
//						else
//						{
//							activeSpheres.clear();
//							singlePassSuper(inRow, outRow, ri.dimensions(), start, dim, 1, activeSpheres);
//						}
//					}
//					if (doBackwardPass)
//					{
//						// Use either temporary activeSpheres or element from activeSpheresEnd
//						if (activeSpheresEnd)
//						{
//							singlePassSuper(inRow, outRow, ri.dimensions(), start, dim, -1, (*activeSpheresEnd)(start));
//						}
//						else
//						{
//							activeSpheres.clear();
//							singlePassSuper(inRow, outRow, ri.dimensions(), start, dim, -1, activeSpheres);
//						}
//					}
//
//					// Copy data back to storage
//					pos = start;
//					for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
//					{
//						toStorageSet(outRow[x], ri(pos));
//					}
//				}
//				else
//				{
//					if (doForwardPass)
//					{
//						if (activeSpheresStart)
//							singlePassFinalSuper(inRow, result, start, dim, 1, &(*activeSpheresStart)(start));
//						else
//							singlePassFinalSuper(inRow, result, start, dim, 1, nullptr);
//					}
//					if (doBackwardPass)
//					{
//						if (activeSpheresEnd)
//							singlePassFinalSuper(inRow, result, start, dim, -1, &(*activeSpheresEnd)(start));
//						else
//							singlePassFinalSuper(inRow, result, start, dim, -1, nullptr);
//					}
//				}
//
//				showThreadProgress(counter, rowCount, showProgressInfo);
//			}
//		}
//	}
//
//
//	/*
//	Copies squared sphere radius data from original centers image to ri image.
//	@param centers2 Original distance ridge.
//	@param ri Temporary image that is to be initialized.
//	@param bounds If the ri covers only a block of centers2, a box defining the block. Pass AABox from [0, 0, 0] to centers2.dimensions() to prepare whole image.
//	*/
//	inline void prepareSuper(const Image<int32_t>& centers2, Image<internals::RiStorageSet>& ri, const AABox<coord_t>& bounds)
//	{
//		if (centers2.dimensions().max() >= numeric_limits<int16_t>::max())
//			throw ITLException(string("Linear image size exceeds ") + toString(numeric_limits<int16_t>::max()) + " pixels. This algorithm is not configured for that big images.");
//
//		ri.ensureSize(bounds.size());
//
//		for (coord_t z = 0; z < ri.depth(); z++)
//		{
//			for (coord_t y = 0; y < ri.height(); y++)
//			{
//				for (coord_t x = 0; x < ri.width(); x++)
//				{
//					//Vec3c p(x, y, z);
//					//ri(p).clear();
//					//int32_t R2 = centers2(bounds.minc + p);
//					//if (R2 > 0)
//					//	// DEBUG
//					//	//ri(p).push_back({ (int16_t)x, (int16_t)y, R2, R2 });
//					//	ri(p).push_back({ (int16_t)x, (int16_t)y });
//
//					Vec3c p(x, y, z);
//					int32_t R2 = centers2(bounds.minc + p);
//					if (R2 > 0)
//					{
//						ri(p).setValues(RiStorageItem{ (int16_t)x, (int16_t)y });
//					}
//				}
//			}
//		}
//	}
//}
//
//
//
//namespace dimredsuper
//{
//
//	void thickmap2(const Image<int32_t>& dmap2, Image<int32_t>& tmap2, Vec3d* counts, bool showProgressInfo)
//	{
//		tmap2.mustNotBe(dmap2);
//
//		internals::buildCircleLookup(max(dmap2));
//
//		tmap2.ensureSize(dmap2);
//
//		Image<internals::RiStorageSet> ri(dmap2.dimensions());
//		internals::prepareSuper(dmap2, ri, AABox<coord_t>(Vec3c(0, 0, 0), dmap2.dimensions()));
//		setValue(tmap2, 0);
//		for (size_t n = 0; n < ri.dimensionality(); n++)
//		{
//			internals::processDimensionSuper(ri, n, true, nullptr, true, nullptr, dmap2, tmap2, showProgressInfo);
//
//			//size_t M = 0;
//			//for (coord_t n = 0; n < ri.pixelCount(); n++)
//			//{
//			//	M = std::max(ri(n).size(), M);
//			//}
//			//cout << "max ricount = " << M << endl;
//
//			if (counts)
//			{
//				double sum = 0;
//				for (coord_t n = 0; n < ri.pixelCount(); n++)
//				{
//					sum += ri(n).size();
//				}
//				(*counts)[n] = sum;
//			}
//		}
//	}
//
//	void thickmap2(const Image<int32_t>& dmap2, Image<int32_t>& tmap2, double* extraBytes, Vec3d* riCounts, bool showProgressInfo)
//	{
//		Vec3d tmp;
//		if (extraBytes && !riCounts)
//		{
//			riCounts = &tmp;
//		}
//		thickmap2(dmap2, tmap2, riCounts, showProgressInfo);
//		if (extraBytes)
//		{
//			(*extraBytes) = riCounts->max() * sizeof(internals::RiStorageItem) + dmap2.pixelCount() * sizeof(internals::RiStorageSet) + tmap2.pixelCount() * sizeof(int32_t);
//		}
//	}
//
//}






// Original dimred version. There's nothing wrong in this version but dimred super is faster and consumes less memory!
//
//namespace internals
//{
//	struct RiItem
//	{
//		/*
//		Squared radius of the sphere that initiated creation of this extent.
//		*/
//		int32_t R2;
//		/*
//		Squared extent of the sphere in the next dimension.
//		*/
//		int32_t ri2;
//	};
//
//	struct ActiveSpheresItem
//	{
//		/*
//		Squared radius of the sphere that initiated creation of this extent.
//		*/
//		int32_t R2;
//		/*
//		Squared extent of the sphere in the current dimension.
//		*/
//		int32_t ri2;
//		/*
//		Center point of the sphere in the current dimension.
//		*/
//		int32_t xi;
//	};
//
//	/*
//	ri image pixel type. Stores original radius and current ri for each sphere.
//	*/
//	typedef vector<RiItem> RiSet;
//
//	/*
//	Comparison functions for comparing RiItem and ActiveSpheresItem.
//	*/
//	inline bool origRSorter(const RiItem& a, const RiItem& b)
//	{
//		if (a.R2 != b.R2)
//			return a.R2 > b.R2;
//		return a.ri2 > b.ri2;
//	}
//
//	inline bool activeSpheresSorter(const ActiveSpheresItem& a, const ActiveSpheresItem& b)
//	{
//		if (a.R2 != b.R2)
//			return a.R2 > b.R2;
//		return a.ri2 > b.ri2;
//	}
//
//	/*
//	Sorts active spheres in inverse order
//	*/
//	inline bool activeSpheresSorter2(const ActiveSpheresItem& a, const ActiveSpheresItem& b)
//	{
//		if (a.R2 != b.R2)
//			return a.R2 < b.R2;
//		return a.ri2 < b.ri2;
//	}
//
//
//	///*
//	//This can be used the find redundant extents in a RiSet.
//	//*/
//	//void debugCheckRiSetRedundancy(const RiSet& set, coord_t dim)
//	//{
//	//	for (size_t i = 0; i < set.size(); i++)
//	//	{
//	//		int32_t Ri2 = set[i].R2;
//	//		int32_t ri2 = set[i].ri2;
//	//		for (size_t j = 0; j < set.size(); j++)
//	//		{
//	//			if (i != j)
//	//			{
//	//				int32_t Rj2 = set[j].R2;
//	//				int32_t rj2 = set[j].ri2;
//
//	//				if (dim == 0)
//	//				{
//	//					if (Ri2 <= Rj2 && doesDiscretizedCircle1FitInto2Cached(ri2, rj2))
//	//						cout << "Redundant circle found" << endl;
//	//				}
//	//				else if (dim == 1)
//	//				{
//	//					if(Ri2 <= Rj2 && largestIntWhoseSquareIsLessThan(ri2) <= largestIntWhoseSquareIsLessThan(rj2))
//	//						cout << "Redundant span found" << endl;
//	//				}
//	//			}
//	//		}
//	//	}
//	//}
//
//
//	/*
//	@param centers Image containing only the row to be processed before first call to this method.
//	@param ri Full ri image that will be updated.
//	@param rowStart Start point of the pixel row to be processed.
//	@param dim Dimension that we are processing.
//	@param step +1 or -1 to indicate the direction of the pass.
//	@param activeSpheres List containing initial active spheres for the row. At exit, contains list of active spheres after processing the row.
//	*/
//	inline void singlePass(const Image<internals::RiSet>& centers, Image<internals::RiSet>& ri, const Vec3c& rowStart, coord_t dim, coord_t step, vector<internals::ActiveSpheresItem>& activeSpheres)
//	{
//		// Stores the spheres that have been encountered and that have not been passed yet.
//		// Stores the center point, original radius, and ri.
//		//vector<internals::ActiveSpheresItem> activeSpheres;
//		vector<internals::ActiveSpheresItem> activeSpheresTmp;
//		vector<internals::ActiveSpheresItem> Ctmp;
//
//		internals::RiSet resultTmp;
//		internals::RiSet resultTmp2;
//
//		//activeSpheres.reserve(40);
//		activeSpheresTmp.reserve(40);
//		Ctmp.reserve(40);
//		resultTmp.reserve(40);
//		resultTmp2.reserve(40);
//
//		// Set start point to the start or to the end of the current row.
//		Vec3c p = rowStart;
//		if (step < 0)
//			p[dim] += (ri.dimension(dim) - 1);
//
//		coord_t dimensionality = ri.dimensionality();
//
//		for (coord_t i = 0; i < ri.dimension(dim); i++, p[dim] += step)
//		{
//			int32_t x = (int32_t)p[dim];
//
//			// If there is one or more sphere centers at the current location, add them to the set of active spheres.
//			const internals::RiSet& C = centers(x);
//
//			// The C list is sorted by R and so is activeSpheres list.
//			// Instead of finding place for each item, convert each item in C list to activeSpheres format (add x coordinate)
//			// and then merge the two sorted lists to construct the new activeSpheres list.
//			Ctmp.clear();
//			for (auto& item : C)
//				Ctmp.push_back({ item.R2, item.ri2, x });
//
//			activeSpheresTmp.clear();
//			activeSpheresTmp.resize(C.size() + activeSpheres.size());
//			merge(Ctmp.begin(), Ctmp.end(), activeSpheres.begin(), activeSpheres.end(), activeSpheresTmp.begin(), internals::activeSpheresSorter);
//			swap(activeSpheres, activeSpheresTmp);
//
//
//			// Iterate through all active spheres and calculate radius for next dimension.
//			resultTmp.clear();
//			auto it = activeSpheres.begin();
//			while (it != activeSpheres.end())
//			{
//				int32_t Rorig2 = it->R2;
//				int32_t R2 = it->ri2;
//				int32_t cx = it->xi;
//
//				int32_t dx = abs(x - cx);
//
//				// Calculate new ri^2
//				int32_t rn2 = R2 - dx * dx;
//				if (rn2 > 0)
//				{
//					// Insert ry2 into the list (TODO: Don't insert if the item is going to be removed anyway?)
//					resultTmp.push_back({ Rorig2, rn2 });
//
//					// Remove possible duplicate Rorig and replace it by (Rorig, max(duplicates))
//					if (resultTmp.size() > 1 && Rorig2 == resultTmp[resultTmp.size() - 2].R2)
//					{
//						if (rn2 > resultTmp[resultTmp.size() - 2].ri2)
//							resultTmp[resultTmp.size() - 2].ri2 = rn2;
//
//						resultTmp.pop_back();
//					}
//
//					it++;
//				}
//				else
//				{
//					// ry is non-positive, i.e. dx >= R
//					// This sphere is not active anymore, so remove it from the list of active spheres.
//					it = activeSpheres.erase(it);
//				}
//			}
//
//			// Rebuild ri list from resultTmp. Don't include those items that are hidden by other items.
//			if (resultTmp.size() > 0)
//			{
//				// Add ri from the previous pass to the ri list and save the result to resultTmp2.
//				// Note that both resultTmp2 and rilist are sorted so we can just merge them.
//				internals::RiSet& rilist = ri(p);
//				resultTmp2.clear();
//				resultTmp2.resize(rilist.size() + resultTmp.size());
//				merge(rilist.begin(), rilist.end(), resultTmp.begin(), resultTmp.end(), resultTmp2.begin(), internals::origRSorter);
//
//
//				// NOTE: This is the basic version, the memory saving version is below.
//				//// Linear time algorithm for finding relevant ri (those not hidden by other items).
//				//// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//				//// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//				//rilist.clear();
//				//rilist.push_back(resultTmp2[0]);
//				//for (size_t n = 1; n < resultTmp2.size(); n++)
//				//{
//				//	int32_t currri2 = rilist[rilist.size() - 1].ri2;
//				//	int32_t newri2 = resultTmp2[n].ri2;
//				//	if (newri2 > currri2)
//				//	{
//				//		rilist.push_back(resultTmp2[n]);
//				//	}
//				//}
//
//				// Linear time algorithm for finding relevant ri (those not hidden by other items).
//				// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
//				// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
//				rilist.clear();
//				rilist.push_back(resultTmp2[0]);
//				//if (dim < ri.dimensionality() - 1) // In the last dimension only the item with maximum R is needed, but we have separate function for that
//				//{
//				for (size_t n = 1; n < resultTmp2.size(); n++)
//				{
//					int32_t currri2 = rilist[rilist.size() - 1].ri2;
//					int32_t newri2 = resultTmp2[n].ri2;
//
//					if (newri2 > currri2) // This is the basic condition that works always (but may include unnecessary items in the rilist)
//					{
//						if (dim == dimensionality - 2) // 3 - 2 == 1 == 2nd dimension
//						{
//							// In the second last dimension only really visible spans are needed as there's no next dimension whose ri we would calculate based on the spans.
//							if (largestIntWhoseSquareIsLessThan(newri2) > largestIntWhoseSquareIsLessThan(currri2))
//								rilist.push_back(resultTmp2[n]);
//						}
//						else if (dim == dimensionality - 3) // 3 - 3 == 0 == 1st dimension
//						{
//							// In the third last dimension we know that only those spans are required that produce visible circles in the output.
//							if (!internals::doesDiscretizedCircle1FitInto2Cached(newri2, currri2))
//								rilist.push_back(resultTmp2[n]);
//						}
//						else
//						{
//							// Here we could insert test if discretized spheres fit into each other etc. etc.
//							rilist.push_back(resultTmp2[n]);
//						}
//					}
//				}
//				//}
//
//				// Make sure we don't use extraneous memory by storing empty items in each pixel of the ri image.
//				rilist.shrink_to_fit();
//
//				//debugCheckRiSetRedundancy(rilist, dim);
//			}
//		}
//	}
//
//	/*
//	Same than singlePass but generates result image instead of ri image, and is thus a little bit faster.
//	@param centers Image containing only the row to be processed before first call to this method.
//	@param Result image. Pixels must be set to zero before first call to this method.
//	@param rowStart Start point of the pixel row to be processed.
//	@param dim Dimension that we are processing.
//	@param step +1 or -1 to indicate the direction of the pass.
//	@param initialActiveSpheres Initial active spheres list. Contains final active spheres list at output. Set to nullptr to assume empty list.
//	*/
//	inline void singlePassFinal(const Image<internals::RiSet>& centers, Image<int32_t>& result, const Vec3c& rowStart, size_t dim, coord_t step, vector<internals::ActiveSpheresItem>* initialActiveSpheres)
//	{
//		// Stores the spheres that have been encountered and that have not been passed yet.
//		// Stores the sphere with the largest R at the top of the priority queue.
//		priority_queue<internals::ActiveSpheresItem, vector<internals::ActiveSpheresItem>, std::function<decltype(internals::activeSpheresSorter2)> > activeSpheres(internals::activeSpheresSorter2);
//
//		if (initialActiveSpheres)
//		{
//			for (const auto& item : *initialActiveSpheres)
//				activeSpheres.push(item);
//		}
//
//
//		// Set start point to the start or end of the current row.
//		Vec3c p = rowStart;
//		if (step < 0)
//			p[dim] += (result.dimension(dim) - 1);
//
//		for (coord_t i = 0; i < result.dimension(dim); i++, p[dim] += step)
//		{
//			int32_t x = (int32_t)p[dim];
//
//			// If there is one or more sphere centers at the current location, add them to the set of active spheres.
//			const internals::RiSet& C = centers(x);
//			for (auto& item : C)
//			{
//				activeSpheres.push({ item.R2, item.ri2, x });
//			}
//
//			while (!activeSpheres.empty())
//			{
//				auto& item = activeSpheres.top();
//
//				int32_t Rorig2 = item.R2;
//				int32_t R2 = item.ri2;
//				int32_t cx = item.xi;
//				int32_t dx = abs(x - cx);
//
//				// Calculate new ri^2
//				int32_t rn2 = R2 - dx * dx;
//				if (rn2 > 0)
//				{
//					// Note that previous pass may have assigned larger value to the output.
//					if (Rorig2 > result(p))
//						result(p) = Rorig2;
//					break;
//				}
//				else
//				{
//					// ry is non-positive, i.e. dx >= R
//					// This sphere is not active anymore, so remove it from the list of active spheres.
//					activeSpheres.pop();
//				}
//
//			}
//
//		}
//
//		// Put stuff back to initialActiveSpheres list as it will be the initial active spheres list for possible adjacent block.
//		if (initialActiveSpheres)
//		{
//			initialActiveSpheres->clear();
//			while (!activeSpheres.empty())
//			{
//				initialActiveSpheres->push_back(activeSpheres.top());
//				activeSpheres.pop();
//			}
//			reverse(initialActiveSpheres->begin(), initialActiveSpheres->end());
//			// TODO: How to avoid sorting here?
//			sort(initialActiveSpheres->begin(), initialActiveSpheres->end(), internals::activeSpheresSorter);
//			//if (!is_sorted(initialActiveSpheres.begin(), initialActiveSpheres.end(), internals::activeSpheresSorter))
//			//{
//			//	cout << "Not sorted" << endl;
//			//	sort(initialActiveSpheres.begin(), initialActiveSpheres.end(), internals::activeSpheresSorter);
//			//}
//		}
//	}
//
//	/*
//	Makes one pass over image in specific dimension and direction.
//	@param ri Image containing ri values from processing of previous dimension.
//	@param dim Dimension to process.
//	@param doForwardPass Set to true if forward pass should be made.
//	@param activeSpheresStart Initial active spheres list for each pixel row for forward pass. Set to nullptr to assume empty lists.
//	@param doBackwardPass Set to true if backward pass should be made.
//	@param activeSpheresEnd Initial active spheres list for each pixel row for backward pass. Set to nullptr to assume empty lists.
//	@param result Result image.
//	@param showProgressInfo Set to true to show progress indicator.
//	*/
//	inline void processDimension(Image<internals::RiSet>& ri, size_t dim,
//		bool doForwardPass, Image<vector<internals::ActiveSpheresItem> >* activeSpheresStart,
//		bool doBackwardPass, Image<vector<internals::ActiveSpheresItem> >* activeSpheresEnd,
//		Image<int32_t>& result, bool showProgressInfo)
//	{
//		// Determine count of pixels to process
//		Vec3c reducedDimensions = internals::getReducedDimensions(ri.dimensions(), dim);
//		coord_t rowCount = reducedDimensions.x * reducedDimensions.y * reducedDimensions.z;
//
//		if (doForwardPass && activeSpheresStart)
//			activeSpheresStart->checkSize(reducedDimensions);
//
//		if (doBackwardPass && activeSpheresEnd)
//			activeSpheresEnd->checkSize(reducedDimensions);
//
//
//		bool isFinalPass = !(dim < ri.dimensionality() - 1);
//
//		size_t counter = 0;
//#pragma omp parallel if(!omp_in_parallel() && ri.pixelCount() > PARALLELIZATION_THRESHOLD)
//		{
//			// Temporary buffer
//			Image<internals::RiSet> centers(ri.dimension(dim));
//			vector<internals::ActiveSpheresItem> activeSpheres;
//			activeSpheres.reserve(40);
//
//			// Determine start points of pixel rows and process each row
//#pragma omp for schedule(dynamic)
//			for (coord_t n = 0; n < rowCount; n++)
//			{
//				Vec3c start = indexToCoords(n, reducedDimensions);
//
//				// Make a copy of the current row as we update the row in the forward pass but need
//				// the original data in the backward pass.
//				Vec3c pos = start;
//				for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
//				{
//					centers(x) = ri(pos);
//					//if(!isFinalPass) //TODO: Did commenting this if break something? I suppose not.
//					ri(pos).clear();
//				}
//
//
//				if (!isFinalPass)
//				{
//					if (doForwardPass)
//					{
//						// Use either temporary activeSpheres or element from activeSpheresStart
//						if (activeSpheresStart)
//						{
//							singlePass(centers, ri, start, dim, 1, (*activeSpheresStart)(start));
//						}
//						else
//						{
//							activeSpheres.clear();
//							singlePass(centers, ri, start, dim, 1, activeSpheres);
//						}
//					}
//					if (doBackwardPass)
//					{
//						// Use either temporary activeSpheres or element from activeSpheresEnd
//						if (activeSpheresEnd)
//						{
//							singlePass(centers, ri, start, dim, -1, (*activeSpheresEnd)(start));
//						}
//						else
//						{
//							activeSpheres.clear();
//							singlePass(centers, ri, start, dim, -1, activeSpheres);
//						}
//					}
//				}
//				else
//				{
//					if (doForwardPass)
//					{
//						if (activeSpheresStart)
//							singlePassFinal(centers, result, start, dim, 1, &(*activeSpheresStart)(start));
//						else
//							singlePassFinal(centers, result, start, dim, 1, nullptr);
//					}
//					if (doBackwardPass)
//					{
//						if (activeSpheresEnd)
//							singlePassFinal(centers, result, start, dim, -1, &(*activeSpheresEnd)(start));
//						else
//							singlePassFinal(centers, result, start, dim, -1, nullptr);
//					}
//				}
//
//				showThreadProgress(counter, rowCount, showProgressInfo);
//			}
//		}
//	}
//
//	/*
//	Copies squared sphere radius data from original centers image to ri image.
//	@param centers2 Original distance ridge.
//	@param ri Temporary image that is to be initialized.
//	@param bounds If the ri covers only a block of centers2, a box defining the block. Pass AABox from [0, 0, 0] to centers2.dimensions() to prepare whole image.
//	*/
//	inline void prepare(const Image<int32_t>& centers2, Image<internals::RiSet>& ri, const AABox<coord_t>& bounds)
//	{
//		ri.ensureSize(bounds.size());
//
//		for (coord_t z = 0; z < ri.depth(); z++)
//		{
//			for (coord_t y = 0; y < ri.height(); y++)
//			{
//				for (coord_t x = 0; x < ri.width(); x++)
//				{
//					Vec3c p(x, y, z);
//					ri(p).clear();
//					int32_t R2 = centers2(bounds.minc + p);
//					if (R2 > 0)
//						ri(p).push_back({ R2, R2 });
//				}
//			}
//		}
//	}
//}
//
//
//namespace dimred
//{
//
//	void thickmap2(Image<int32_t>& dmap2, Vec3d* counts, bool showProgressInfo)
//	{
//		internals::buildCircleLookup(max(dmap2));
//
//		Image<internals::RiSet> ri(dmap2.dimensions());
//		internals::prepare(dmap2, ri, AABox<coord_t>(Vec3c(0, 0, 0), dmap2.dimensions()));
//		setValue(dmap2, 0);
//		for (size_t n = 0; n < ri.dimensionality(); n++)
//		{
//			internals::processDimension(ri, n, true, nullptr, true, nullptr, dmap2, showProgressInfo);
//
//			//size_t M = 0;
//			//for (coord_t n = 0; n < ri.pixelCount(); n++)
//			//{
//			//	M = std::max(ri(n).size(), M);
//			//}
//			//cout << "max ricount = " << M << endl;
//
//			if (counts)
//			{
//				double sum = 0;
//				for (coord_t n = 0; n < ri.pixelCount(); n++)
//				{
//					sum += ri(n).size();
//				}
//				(*counts)[n] = sum;
//			}
//		}
//	}
//
//	void thickmap2(Image<int32_t>& dmap2, double* extraBytes, Vec3d* riCounts, bool showProgressInfo)
//	{
//		Vec3d tmp;
//		if (extraBytes && !riCounts)
//		{
//			riCounts = &tmp;
//		}
//		thickmap2(dmap2, riCounts, showProgressInfo);
//		if (extraBytes)
//		{
//			(*extraBytes) = riCounts->max() * sizeof(internals::RiItem) + dmap2.pixelCount() * sizeof(internals::RiSet);
//		}
//	}
//
//}







// This is a version of dimredsmartblocks2 that is compatible with standard dimred

//namespace dimredsmartblocks2
//{
//
//	/**
//	Describes one block in block-wise processing.
//	*/
//	struct Block
//	{
//		/**
//		Bounds of this block in image coordinates.
//		*/
//		AABox<coord_t> bounds;
//
//	private:
//		/**
//		Initial active spheres list for each pixel row in each coordinate direction.
//		Access as activeSpheres[dimension][0/1]
//		where 0 = start and 1 = end
//		*/
//		array<array<tuple<bool, Image<vector<internals::ActiveSpheresItem> > >, 2>, 3 > activeSpheres1, activeSpheres2;
//
//		/**
//		Pointers to current active spheres for reading and writing.
//		Used because we have to be able to read the old active spheres while updating.
//		*/
//		array<array<tuple<bool, Image<vector<internals::ActiveSpheresItem> > >, 2>, 3 > *activeSpheresRead, *activeSpheresWrite;
//
//
//		void clearAll(Image<vector<internals::ActiveSpheresItem> >& img)
//		{
//			for (coord_t n = 0; n < img.pixelCount(); n++)
//				img(n).clear();
//		}
//
//
//		size_t totalCount(const Image<vector<internals::ActiveSpheresItem> >& img)
//		{
//			size_t count = 0;
//			for (coord_t n = 0; n < img.pixelCount(); n++)
//				count += img(n).size();
//			return count;
//		}
//
//	public:
//
//		/**
//		Calculates approximation of amount of memory allocated to this object in bytes.
//		*/
//		size_t memSize()
//		{
//			size_t count = 0;
//
//			count += totalCount(get<1>((activeSpheres1[0][0])));
//			count += totalCount(get<1>((activeSpheres1[0][1])));
//			count += totalCount(get<1>((activeSpheres1[1][0])));
//			count += totalCount(get<1>((activeSpheres1[1][1])));
//			count += totalCount(get<1>((activeSpheres1[2][0])));
//			count += totalCount(get<1>((activeSpheres1[2][1])));
//
//			count += totalCount(get<1>((activeSpheres2[0][0])));
//			count += totalCount(get<1>((activeSpheres2[0][1])));
//			count += totalCount(get<1>((activeSpheres2[1][0])));
//			count += totalCount(get<1>((activeSpheres2[1][1])));
//			count += totalCount(get<1>((activeSpheres2[2][0])));
//			count += totalCount(get<1>((activeSpheres2[2][1])));
//
//			return count * sizeof(internals::ActiveSpheresItem);
//		}
//
//		/**
//		Used to initialize activeSpheres.
//		*/
//		void setActiveSpheres(size_t dim, size_t dir, const Vec3c& size)
//		{
//			{
//				auto& tuple = (*activeSpheresRead)[dim][dir];
//				get<1>(tuple).ensureSize(size);
//				get<0>(tuple) = true;
//			}
//
//			{
//				auto& tuple = (*activeSpheresWrite)[dim][dir];
//				get<1>(tuple).ensureSize(size);
//				get<0>(tuple) = true;
//			}
//		}
//
//		/*
//		Gets bool that indicates if active spheres for given dimension and direction are set,
//		and the corresponding active spheres image.
//		*/
//		tuple<bool, Image<vector<internals::ActiveSpheresItem> > >& getActiveSpheresSource(size_t dim, size_t dir)
//		{
//			return (*activeSpheresRead)[dim][dir];
//		}
//
//		/**
//		Copies values in the given activeSpheres image to active spheres of given dimension and direction.
//		*/
//		void setActiveSpheres(size_t dim, size_t dir, const Image<vector<internals::ActiveSpheresItem> >& activeSpheres)
//		{
//			//if (!activeSpheres.sizeEquals(getReducedDimensions(bounds.size(), dim)))
//			//	throw logic_error("Bad activeSpheres image size.");
//
//			auto& tuple = (*activeSpheresWrite)[dim][dir];
//
//			setValue(get<1>(tuple), activeSpheres);
//			get<0>(tuple) = true;
//		}
//
//		/**
//		Called after all blocks in single iteration have been processed.
//		Discards old active spheres list and uses the newly written lists as active spheres for next iteration.
//		*/
//		void commitChanges()
//		{
//			swap(activeSpheresRead, activeSpheresWrite);
//
//			for (size_t n = 0; n < (*activeSpheresWrite).size(); n++)
//			{
//				for (size_t m = 0; m < 2; m++)
//				{
//					//get<0>((*activeSpheresRead)[n][m]) = get<0>((*activeSpheresWrite)[n][m]);
//					//setValue(get<1>((*activeSpheresRead)[n][m]), get<1>((*activeSpheresWrite)[n][m]));
//
//					clearAll(get<1>((*activeSpheresWrite)[n][m]));
//
//				}
//			}
//		}
//
//		/**
//		Constructor
//		*/
//		Block()
//		{
//			activeSpheresRead = &activeSpheres1;
//			activeSpheresWrite = &activeSpheres2;
//		}
//
//		/**
//		Constructor.
//		*/
//		Block(const AABox<coord_t>& box) :
//			bounds(box)
//		{
//			activeSpheresRead = &activeSpheres1;
//			activeSpheresWrite = &activeSpheres2;
//		}
//
//
//		/**
//		Copy constructor.
//		TODO: This is a bit hacky. Better way would be to have operator= in Image, but that kind of breaks the mindset behind it...
//		The pointers would need to be handled separately anyway.
//		*/
//		Block(const Block& block)
//		{
//			*this = block;
//		}
//
//		/**
//		Assignment operator.
//		*/
//		Block &operator =(const Block& block)
//		{
//			// NOTE: Swap template would not work here (directly) because of the array members that cannot be swapped because
//			// there is no operator= for the Image class.
//			// NOTE: Self-assignment should be ok.
//
//			bounds = block.bounds;
//
//			for (size_t n = 0; n < block.activeSpheres1.size(); n++)
//			{
//				for (size_t m = 0; m < 2; m++)
//				{
//					get<0>(activeSpheres1[n][m]) = get<0>(block.activeSpheres1[n][m]);
//					setValue(get<1>(activeSpheres1[n][m]), get<1>(block.activeSpheres1[n][m]));
//
//					get<0>(activeSpheres2[n][m]) = get<0>(block.activeSpheres2[n][m]);
//					setValue(get<1>(activeSpheres2[n][m]), get<1>(block.activeSpheres2[n][m]));
//				}
//			}
//
//			if (block.activeSpheresRead == &block.activeSpheres1)
//				activeSpheresRead = &activeSpheres1;
//			else
//				activeSpheresRead = &activeSpheres2;
//
//			if (block.activeSpheresWrite == &block.activeSpheres1)
//				activeSpheresWrite = &activeSpheres1;
//			else
//				activeSpheresWrite = &activeSpheres2;
//
//			return *this;
//		}
//	};
//
//	/**
//	Builds subdivision of distance ridge into blocks of given size.
//	Does not actually copy pixel data but creates array/image of Block structures.
//	*/
//	void subdivide(const Vec3c& dmap2Dimensions, const Vec3c& blockSize, Image<Block>& blocks)
//	{
//		// Create blocks
//		Vec3c blockCount = ceil(Vec3d(dmap2Dimensions).componentwiseDivide(Vec3d(blockSize)));
//		blocks.ensureSize(blockCount);
//		for (coord_t z = 0; z < blockCount.z; z++)
//		{
//			for (coord_t y = 0; y < blockCount.y; y++)
//			{
//				for (coord_t x = 0; x < blockCount.x; x++)
//				{
//					Vec3c pos = blockSize.componentwiseMultiply(Vec3c(x, y, z));
//					Vec3c right = pos + blockSize;
//					clamp(right, Vec3c(0, 0, 0), dmap2Dimensions);
//					blocks(x, y, z).bounds = AABox<coord_t>(pos, right);
//				}
//			}
//		}
//	}
//
//
//	/*
//	Adds x to all xi values in all items of all lists in the given active spheres image.
//	*/
//	void addxi(Image<vector<internals::ActiveSpheresItem> >& activeSpheres, int32_t x)
//	{
//		for (coord_t n = 0; n < activeSpheres.pixelCount(); n++)
//		{
//			vector<internals::ActiveSpheresItem>& v = activeSpheres(n);
//			for (internals::ActiveSpheresItem& item : v)
//			{
//				item.xi += x;
//			}
//		}
//	}
//
//
//	/**
//	Sets initial active spheres of neighbouring image based on final active spheres of current image.
//	@param pos Position of current image.
//	@param blocks Blocks image.
//	@param dim Dimension where the neighbours are considered. Pass 0 to indicate next or previous block in x-direction, 1 in y-direction etc.
//	@param dirn Index of neighbour: 0 for next, 1 for previous.
//	*/
//	void propagate(const Vec3c& pos, Image<Block>& blocks, coord_t dim, coord_t dirn, Image<vector<internals::ActiveSpheresItem> >& activeSpheres)
//	{
//		coord_t dirs[2] = { 1, -1 };
//		Vec3c nextBlockPos = pos;
//		nextBlockPos[dim] += dirs[dirn];
//
//		if (blocks.isInImage(nextBlockPos))
//		{
//			Block& b = blocks(pos);
//			Block& bn = blocks(nextBlockPos);
//
//			// Convert active sphere positions to the coordinates of the next block
//			if (dirs[dirn] > 0)
//				addxi(activeSpheres, -(int32_t)b.bounds.size()[dim]);
//			else
//				addxi(activeSpheres, (int32_t)bn.bounds.size()[dim]);
//
//			bn.setActiveSpheres(dim, dirn, activeSpheres);
//		}
//	}
//
//	/**
//	Calculate local thickness in one block and update initial active spheres of neighouring blocks.
//	@param pos Position of the block.
//	@param blocks Blocks image.
//	@param dmap2 Original distance ridge.
//	@param result Full result image.
//	@param showProgressInfo Indicates if progress indicator should be shown.
//	@param counts Set to non-null value to calculate total count of elements in ri image in each iteration.
//	*/
//	void processBlock(const Vec3c& pos, Image<Block>& blocks, const Image<int32_t>& dmap2, Image<int32_t>& result, bool showProgressInfo, Vec3d* counts = nullptr)
//	{
//		Block& b = blocks(pos);
//
//		// Create ri of correct size and load data from dmap2 to it.
//		Image<internals::RiSet> ri;
//		prepare(dmap2, ri, b.bounds);
//
//		// Create output image.
//		Image<int32_t> resultBlock(b.bounds.size());
//		crop(result, resultBlock, b.bounds.position());
//
//		for (size_t n = 0; n < dmap2.dimensionality(); n++)
//		{
//			// Check if initial active spheres image is available.
//			auto& tupleStart = b.getActiveSpheresSource(n, 0);
//			auto& tupleEnd = b.getActiveSpheresSource(n, 1);
//			bool startAvailable = get<0>(tupleStart);
//			bool endAvailable = get<0>(tupleEnd);
//			if (startAvailable || endAvailable)
//			{
//				// Initialize active spheres lists from stored data.
//				Image<vector<internals::ActiveSpheresItem> > activeSpheresStart, activeSpheresEnd;
//				setValue(activeSpheresStart, get<1>(tupleStart));
//				setValue(activeSpheresEnd, get<1>(tupleEnd));
//
//				// Process
//				processDimension(ri, n, startAvailable, &activeSpheresStart, endAvailable, &activeSpheresEnd, resultBlock, showProgressInfo);
//
//				// The final activeSpheres image of this block is initial activeSpheres image of the next block
//				// in the processing direction.
//				if (startAvailable)
//					propagate(pos, blocks, n, 0, activeSpheresStart);
//				if (endAvailable)
//					propagate(pos, blocks, n, 1, activeSpheresEnd);
//			}
//
//			if (counts)
//			{
//				double sum = 0;
//				for (coord_t n = 0; n < ri.pixelCount(); n++)
//				{
//					sum += ri(n).size();
//				}
//				(*counts)[n] = sum;
//			}
//		}
//
//		setValue(result, resultBlock, b.bounds.position());
//	}
//
//	/**
//	Calculates city-block distance between p1 and p2.
//	*/
//	coord_t blockDistance(const Vec3c& p1, const Vec3c& p2)
//	{
//		Vec3c diff = p1 - p2;
//		return diff.abs().sum();
//	}
//
//	/**
//	Finds positions of blocks that should be processed in iteration 'it'.
//	Returns list of positions of all block whose city-block distance to some corner of the image is 'it'.
//	*/
//	vector<Vec3c> findReadyBlocks(Image<Block>& blocks, size_t it)
//	{
//		coord_t w = blocks.width() - 1;
//		coord_t h = blocks.height() - 1;
//		coord_t d = blocks.depth() - 1;
//
//		vector<Vec3c> result;
//		for (coord_t z = 0; z < blocks.depth(); z++)
//		{
//			for (coord_t y = 0; y < blocks.height(); y++)
//			{
//				for (coord_t x = 0; x < blocks.width(); x++)
//				{
//					Vec3c p(x, y, z);
//					array<coord_t, 8> ds;
//					ds[0] = blockDistance(p, Vec3c(0, 0, 0));
//					ds[1] = blockDistance(p, Vec3c(w, 0, 0));
//					ds[2] = blockDistance(p, Vec3c(w, h, 0));
//					ds[3] = blockDistance(p, Vec3c(0, h, 0));
//					ds[4] = blockDistance(p, Vec3c(0, 0, d));
//					ds[5] = blockDistance(p, Vec3c(w, 0, d));
//					ds[6] = blockDistance(p, Vec3c(w, h, d));
//					ds[7] = blockDistance(p, Vec3c(0, h, d));
//					if (std::find(ds.begin(), ds.end(), it) != ds.end())
//						result.push_back(p);
//				}
//			}
//		}
//
//		return result;
//	}
//
//	//void saveBlockStateCheckImage(const Image<Block>& blocks, coord_t it)
//	//{
//	//	Image<uint8_t> readyState(blocks.dimensions());
//	//	for (coord_t n = 0; n < blocks.pixelCount(); n++)
//	//	{
//	//		const Block& b = blocks(n);
//	//		//readyState(n) = b.readyCount();
//	//		//readyState(n) = b.processable() ? 255 : 0;
//	//		readyState(n) = pixelRound<uint8_t>(b.round);
//	//	}
//
//	//	raw::writed(readyState, "./dimredsmartblocks/debug_readystate_" + toString(it));
//	//}
//
//	/*
//	Perform one iteration on the blocks,
//	return count of blocks that had to be processed.
//	*/
//	size_t iteration(Image<Block>& blocks, const Image<int32_t>& dmap2, Image<int32_t>& result, size_t it, int reportingIt, bool showProgressInfo, double* maxRiCount = nullptr)
//	{
//		vector<Vec3c> toProcess = findReadyBlocks(blocks, it);
//
//		if (toProcess.size() <= 0)
//			return 0;
//
//		//cout << "Iteration " << reportingIt << ", processing " << toProcess.size() << " blocks..." << endl;
//
//		//for (Vec3c pos : toProcess)
//		//	blocks(pos).round = reportingIt;
//		//saveBlockStateCheckImage(blocks, reportingIt);
//
//		if (!maxRiCount)
//		{
//			for (Vec3c pos : toProcess)
//				processBlock(pos, blocks, dmap2, result, false);
//		}
//		else
//		{
//			// This version calculates the maximum number of ri elements in memory at the same time.
//			for (Vec3c pos : toProcess)
//			{
//				Vec3d counts;
//				processBlock(pos, blocks, dmap2, result, false, &counts);
//				*maxRiCount = std::max(*maxRiCount, (double)counts.max());
//			}
//		}
//
//		for (coord_t n = 0; n < blocks.pixelCount(); n++)
//			blocks(n).commitChanges();
//
//		return toProcess.size();
//	}
//
//	void thickmap2Blocks(const Image<int32_t>& dmap2, Image<int32_t>& result, const Vec3c& blockSize, double* extraBytes, double* diskBytes, bool showProgressInfo)
//	{
//		dmap2.mustNotBe(result);
//
//		internals::buildCircleLookup(max(dmap2));
//
//		result.ensureSize(dmap2);
//
//		Image<Block> blocks;
//		subdivide(dmap2.dimensions(), blockSize, blocks);
//
//		// Initialize activeSpheres images to empty for those blocks that are in the edge of the image.
//		for (size_t skip = 0; skip < dmap2.dimensionality(); skip++)
//		{
//			Vec3c reducedsize = blocks.dimensions();
//			reducedsize[skip] = 1;
//
//			for (coord_t z = 0; z < reducedsize.z; z++)
//			{
//				for (coord_t y = 0; y < reducedsize.y; y++)
//				{
//					for (coord_t x = 0; x < reducedsize.x; x++)
//					{
//						Vec3c coords(x, y, z);
//
//						coords[skip] = 0;
//						// Start block is blocks(coords)
//						{
//							Block& b = blocks(coords);
//							b.setActiveSpheres(skip, 0, internals::getReducedDimensions(b.bounds.size(), skip));
//						}
//
//						coords[skip] = blocks.dimension(skip) - 1;
//						// End block is blocks(coords)
//						{
//							Block& b = blocks(coords);
//							b.setActiveSpheres(skip, 1, internals::getReducedDimensions(b.bounds.size(), skip));
//						}
//					}
//				}
//			}
//		}
//
//		// Iterate until no changes are made
//		size_t totalProcessed = 0;
//		coord_t it = 0;
//		double maxri = 0;
//		while (true)
//		{
//			size_t count;
//			if (!extraBytes)
//				count = iteration(blocks, dmap2, result, it, (int)it, showProgressInfo);
//			else
//				count = iteration(blocks, dmap2, result, it, (int)it, showProgressInfo, &maxri);
//
//			if (count <= 0)
//				break;
//			totalProcessed += count;
//			it++;
//		}
//
//
//		//cout << "Processed each block on average " << ((double)totalProcessed / (double)blocks.pixelCount()) << " times." << endl;
//		//cout << "Could process on average " << ((double)totalProcessed / (double)it) << " blocks in parallel." << endl;
//
//		if (extraBytes)
//		{
//			// Calculate amount of memory required to be in memory at once.
//
//			*extraBytes = 0;
//
//			// Maximum size of all block data structures
//			size_t maxBlockSize = 0;
//			for (coord_t n = 0; n < blocks.pixelCount(); n++)
//			{
//				maxBlockSize = std::max(maxBlockSize, blocks(n).memSize());
//			}
//			*extraBytes += maxBlockSize;
//
//			// One block result image allocated at once
//			*extraBytes += blockSize.x * blockSize.y * blockSize.z * sizeof(int32_t);
//
//			// Arrays for holding ri values in one block
//			*extraBytes += blockSize.x * blockSize.y * blockSize.z * sizeof(internals::RiSet);
//
//			// Maximum size of ri's in memory
//			*extraBytes += maxri * sizeof(internals::RiItem);
//
//
//			/* This calculates size of memory used by block-wise processing when all the blocks are in memory at the same time. Not very relevant.
//			// Separate result image
//			*extraBytes = (double)result.pixelCount() * sizeof(int32_t);
//
//			// Size of all block data structures
//			for (coord_t n = 0; n < blocks.pixelCount(); n++)
//			{
//				*extraBytes += blocks(n).memSize();
//			}
//
//			// One block result image allocated at once
//			*extraBytes += blockSize.x * blockSize.y * blockSize.z * sizeof(int32_t);
//
//			// Array holding ri values
//			*extraBytes += blockSize.x * blockSize.y * blockSize.z * sizeof(vector<internals::RiItem>);
//
//			// Maximum size of ri's in memory
//			*extraBytes += maxri * sizeof(internals::RiItem);
//			*/
//		}
//
//		if (diskBytes)
//		{
//			// Blocks data structures must be saved to disk between processing of each block.
//			size_t totalBlockSize = 0;
//			for (coord_t n = 0; n < blocks.pixelCount(); n++)
//			{
//				totalBlockSize += blocks(n).memSize();
//			}
//
//			*diskBytes = (double)totalBlockSize;
//		}
//	}
//}