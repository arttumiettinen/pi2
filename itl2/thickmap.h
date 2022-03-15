#pragma once

#include <vector>
#include <queue>
#include <functional>

#if defined(__linux__)
	#include <parallel/algorithm>
#elif defined(_WIN32)
	#include <execution>
#elif defined(__APPLE__)
	#include <algorithm>
#endif


#include "image.h"
#include "math/vec3.h"
#include "math/aabox.h"
#include "danielsson.h"

namespace itl2
{
	namespace internals
	{
		/**
		Tests if discretized circle whose squared radius is r1Square fit completely into discretized circle whose squared radius is r2Square.
		*/
		bool doesDiscretizedCircle1FitInto2(int32_t r1Square, int32_t r2Square);

		/**
		Builds lookup table that doesDiscretizedCircle1FitInto2Cached uses.
		@param maxrSquare Maximum squared radius found in the image.
		*/
		void buildCircleLookup(int32_t maxrSquare);

		/**
		Lookup table for doesDiscretizedCircle1FitInto2Cached.
		Element circleLookup[r2] stores the maximal squared radius of a circle that fits into a circle of squared radius r2.
		*/
		extern std::vector<int32_t> circleLookup;

		/**
		Tests if discretized circle whose squared radius is r1Square fit completely into discretized circle whose squared radius is r2Square.
		Assumes buildCircleLookup has been called with argument whose value is greater than or equal to r2Square.
		*/
		inline bool doesDiscretizedCircle1FitInto2Cached(int32_t r1Square, int32_t r2Square)
		{
			return r1Square <= circleLookup[r2Square];
		}
	}


	namespace internals
	{
#pragma pack(push, 1)

		/**
		Data entry that must be saved for each pixel for each sphere that must be considered in that pixel.
		This data structure should be as small as possible.
		*/
		struct RiStorageItem
		{
			/**
			x and y coordinate of the center of the sphere that this item corresponds to.
			*/
			int16_t srcX, srcY;

			bool operator==(const RiStorageItem& r) const
			{
				return srcX == r.srcX && srcY == r.srcY;
			}

			bool operator!=(const RiStorageItem& rhs) const
			{
				return !(*this == rhs);
			}
		};

		/**
		Simple array that consumes as little memory as possible.
		Uses small array optimization, i.e. if all the items fill into the space taken by a pointer,
		does not allocate any memory but stores the data in the place of the pointer.
		Consequently, for small arrays, requires sizeof(count_t) + sizeof(item_t*) bytes of memory.
		*/
		template<typename item_t, typename count_t = uint16_t> class SmallVector
		{
		private:
			/**
			Stores count of items.
			*/
			count_t count;

			/**
			Element count where the vector changes from local storage to heap-allocated storage.
			*/
			static constexpr size_t ARR_SIZE = sizeof(item_t*) / sizeof(item_t);

			static_assert(ARR_SIZE > 0, "TODO: SmallVector class is not designed to work if item size is greater than pointer size.");

			/**
			Data store, either pointer to heap-allocated memory or array of items.
			*/
			union
			{
				item_t* pData;
				item_t arr[ARR_SIZE];
			};

			/**
			Frees memory associated with this vector.
			*/
			void deleteData()
			{
				if (!isUsingOptimization())
				{
					delete[] pData;
					pData = nullptr;
				}
				else
				{
					for (size_t n = 0; n < count; n++)
						arr[n].~item_t();
				}
				count = 0;
			}

		public:
			/**
			Constructor
			*/
			SmallVector(const std::vector<item_t>& items) : pData(nullptr), count(0)
			{
				setValues(items);
			}

			/**
			Constructor
			*/
			SmallVector() : pData(nullptr), count(0)
			{
			}

			/**
			Constructor
			*/
			SmallVector(const SmallVector& other) : pData(nullptr), count(other.count)
			{
				if (isUsingOptimization())
				{
					for (size_t n = 0; n < count; n++)
						arr[n] = other.arr[n];
				}
				else
				{
					pData = new item_t[count];
					for (count_t n = 0; n < count; n++)
						pData[n] = other.pData[n];
				}
			}

			/**
			Assignment
			*/
			SmallVector& operator=(SmallVector other)
			{
				deleteData();
				std::swap(count, other.count);
				if (isUsingOptimization())
				{
					for (size_t n = 0; n < count; n++)
						arr[n] = other.arr[n];
				}
				else
				{
					std::swap(pData, other.pData);
				}

				return *this;
			}

			/**
			Destructor
			*/
			~SmallVector()
			{
				deleteData();
			}

			/**
			Sets values in this array to those in the given vector.
			*/
			void setValues(const std::vector<item_t>& items)
			{
				deleteData();

				count = (count_t)items.size();

				if (isUsingOptimization())
				{
					for (size_t n = 0; n < count; n++)
						arr[n] = items[n];
				}
				else
				{
					if (items.size() >= std::numeric_limits<count_t>::max())
						throw std::invalid_argument("The source array has too many elements.");

					pData = new item_t[count];
					for (count_t n = 0; n < count; n++)
						pData[n] = items[n];
				}
			}

			/**
			Populates this array with single item.
			*/
			void setValues(const item_t& item)
			{
				deleteData();
				count = 1;
				arr[0] = item;
			}

			/**
			Gets count of items in this array.
			*/
			count_t size() const
			{
				return count;
			}

			/**
			Gets maximum element count that fits into this vector without additional memory allocation.
			*/
			static constexpr size_t SmallVectorSize()
			{
				return ARR_SIZE;
			}

			/**
			Gets a value indicating if the small vector optimization is in use.
			*/
			bool isUsingOptimization() const
			{
				return count <= ARR_SIZE;
			}

			/**
			Returns amount of memory allocated by this object.
			*/
			size_t memSize() const
			{
				if (isUsingOptimization())
					return sizeof(SmallVector);
				else
					return sizeof(SmallVector) + count * sizeof(item_t);
			}

			/**
			Returns reference to element at given index.
			*/
			item_t& operator[](count_t index)
			{
				return const_cast<item_t&>(const_cast<const SmallVector<item_t, count_t>*>(this)->operator[](index));
			}

			/**
			Returns reference to element at given index.
			*/
			const item_t& operator[](count_t index) const
			{
				if (isUsingOptimization())
					return arr[index];
				else
					return pData[index];
			}

			/**
			Equality operator
			*/
			bool operator==(const SmallVector& rhs) const
			{
				if (size() != rhs.size())
					return false;

				for (count_t n = 0; n < count; n++)
				{
					if ((*this)[n] != rhs[n])
						return false;
				}

				return true;
			}

			bool operator!=(const SmallVector& rhs) const
			{
				return !(*this == rhs);
			}
		};


		/*
		ri image pixel type. Stores original radius and current ri for each sphere.
		*/
		typedef SmallVector<RiStorageItem> RiStorageSet;


		template<typename pixel_t> struct RiSuperItem
		{
			/*
			Squared radius of the sphere that initiated creation of this extent.
			*/
			pixel_t R2;

			/*
			Squared extent of the sphere in the next dimension.
			*/
			pixel_t ri2;

			/**
			x and y coordinate of the center of the sphere that this item corresponds to.
			*/
			int16_t srcX;
			int16_t srcY;
		};

		template<typename pixel_t> struct ActiveSpheresSuperItem
		{
			using signed_t = typename NumberUtils<pixel_t>::SignedType;

			/*
			Squared radius of the sphere that initiated creation of this extent.
			*/
			pixel_t R2;

			/*
			Squared extent of the sphere in the current dimension.
			*/
			pixel_t ri2;

			/*
			Center point of the sphere in the current dimension.
			This must be signed in order to support block-wise processing.
			*/
			signed_t xi;

			/**
			x and y coordinate of the center of the sphere that this item corresponds to.
			*/
			int16_t srcX;
			int16_t srcY;
		};

#pragma pack(pop)

		/*
		ri image pixel type. Stores original radius and current ri for each sphere.
		*/
		template<typename pixel_t>
		using RiSuperSet = std::vector<RiSuperItem<pixel_t> >;

		/*
		Comparison functions for comparing RiItem and ActiveSpheresItem.
		*/
		template<typename pixel_t> bool origRSuperSorter(const RiSuperItem<pixel_t>& a, const RiSuperItem<pixel_t>& b)
		{
			if (a.R2 != b.R2)
				return a.R2 > b.R2;
			return a.ri2 > b.ri2;
		}

		template<typename pixel_t> bool activeSpheresSuperSorter(const ActiveSpheresSuperItem<pixel_t>& a, const ActiveSpheresSuperItem<pixel_t>& b)
		{
			if (a.R2 != b.R2)
				return a.R2 > b.R2;
			return a.ri2 > b.ri2;
		}

		/*
		Sorts active spheres in inverse order
		*/
		template<typename pixel_t> bool activeSpheresSuperSorter2(const ActiveSpheresSuperItem<pixel_t>& a, const ActiveSpheresSuperItem<pixel_t>& b)
		{
			if (a.R2 != b.R2)
				return a.R2 < b.R2;
			return a.ri2 < b.ri2;
		}



		/*
		Sets size[dim] = 0 and returns the result.
		*/
		inline Vec3c getReducedDimensions(Vec3c size, size_t dim)
		{
			size[dim] = 1;
			return size;
		}

		/*
		@param centers Image containing only the row to be processed before first call to this method.
		@param ri Full ri image that will be updated.
		@param rowStart Start point of the pixel row to be processed.
		@param dim Dimension that we are processing.
		@param step +1 or -1 to indicate the direction of the pass.
		@param activeSpheres List containing initial active spheres for the row. At exit, contains list of active spheres after processing the row.
		*/
		template<typename pixel_t> void singlePassSuper(const std::vector<internals::RiSuperSet<pixel_t> >& centers, std::vector<internals::RiSuperSet<pixel_t> >& ri, size_t dimensionality, const Vec3c& dimensions, const Vec3c& rowStart, coord_t dim, coord_t step)
		{
			using signed_t = typename NumberUtils<pixel_t>::SignedType;

			// Stores the spheres that have been encountered and that have not been passed yet.
			// Stores the center point, original radius, and ri.
			std::vector<internals::ActiveSpheresSuperItem<pixel_t> > activeSpheres;
			std::vector<internals::ActiveSpheresSuperItem<pixel_t> > activeSpheresTmp;
			std::vector<internals::ActiveSpheresSuperItem<pixel_t> > Ctmp;

			internals::RiSuperSet<pixel_t> resultTmp;
			internals::RiSuperSet<pixel_t> resultTmp2;

			activeSpheres.reserve(40);
			activeSpheresTmp.reserve(40);
			Ctmp.reserve(40);
			resultTmp.reserve(40);
			resultTmp2.reserve(40);

			// Set start point to the start or to the end of the current row.
			Vec3c p = rowStart;
			if (step < 0)
				p[dim] += (dimensions[dim] - 1);

			//coord_t dimensionality = getDimensionality(dimensions); // This does not work if processing block that has size 1 as last dimension.

			for (coord_t i = 0; i < dimensions[dim]; i++, p[dim] += step)
			{
				signed_t x = (signed_t)p[dim];

				// If there is one or more sphere centers at the current location, add them to the set of active spheres.
				const internals::RiSuperSet<pixel_t>& C = centers[p[dim]];

				// The C list is sorted by R and so is activeSpheres list.
				// Instead of finding place for each item, convert each item in C list to activeSpheres format (add x coordinate)
				// and then merge the two sorted lists to construct the new activeSpheres list.
				Ctmp.clear();
				for (auto& item : C)
					Ctmp.push_back({ item.R2, item.ri2, x, item.srcX, item.srcY });

				activeSpheresTmp.clear();
				activeSpheresTmp.resize(C.size() + activeSpheres.size());
				merge(Ctmp.begin(), Ctmp.end(), activeSpheres.begin(), activeSpheres.end(), activeSpheresTmp.begin(), internals::activeSpheresSuperSorter<pixel_t>);
				swap(activeSpheres, activeSpheresTmp);


				// Iterate through all active spheres and calculate radius for next dimension.
				resultTmp.clear();
				auto it = activeSpheres.begin();
				while (it != activeSpheres.end())
				{
					pixel_t Rorig2 = it->R2;
					pixel_t R2 = it->ri2;
					signed_t cx = it->xi;

					int16_t origcx = it->srcX;
					int16_t origcy = it->srcY;

					signed_t dx = abs(x - cx);

					// Calculate new ri^2
					signed_t rn2 = (signed_t)R2 - dx * dx;
					if (rn2 > 0)
					{
						// Insert ry2 into the list, but don't insert duplicates.
						if (resultTmp.size() > 0 && Rorig2 == resultTmp[resultTmp.size() - 1].R2)
						{
							// This is a duplicate R2 entry. Use the entry with the larger ri.
							if (resultTmp[resultTmp.size() - 1].ri2 < (pixel_t)rn2)
								resultTmp[resultTmp.size() - 1] = { Rorig2, (pixel_t)rn2, origcx, origcy };
						}
						else
						{
							resultTmp.push_back({ Rorig2, (pixel_t)rn2, origcx, origcy });
						}

						it++;
					}
					else
					{
						// ry is non-positive, i.e. dx >= R
						// This sphere is not active anymore, so remove it from the list of active spheres.
						it = activeSpheres.erase(it);
					}
				}

				// Rebuild ri list from resultTmp. Don't include those items that are hidden by other items.
				if (resultTmp.size() > 0)
				{
					// Add ri from the previous pass to the ri list and save the result to resultTmp2.
					// Note that both resultTmp2 and rilist are sorted so we can just merge them.
					internals::RiSuperSet<pixel_t>& rilist = ri[p[dim]];
					resultTmp2.clear();
					resultTmp2.resize(rilist.size() + resultTmp.size());
					merge(rilist.begin(), rilist.end(), resultTmp.begin(), resultTmp.end(), resultTmp2.begin(), internals::origRSuperSorter<pixel_t>);


					// NOTE: This is the basic version, the memory saving version is below.
					//// Linear time algorithm for finding relevant ri (those not hidden by other items).
					//// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
					//// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
					//rilist.clear();
					//rilist.push_back(resultTmp2[0]);
					//for (size_t n = 1; n < resultTmp2.size(); n++)
					//{
					//	int32_t currri2 = rilist[rilist.size() - 1].ri2;
					//	int32_t newri2 = resultTmp2[n].ri2;
					//	if (newri2 > currri2)
					//	{
					//		rilist.push_back(resultTmp2[n]);
					//	}
					//}

					// Linear time algorithm for finding relevant ri (those not hidden by other items).
					// Requires that ri list is sorted according to R in descending order and then by ri in descending order.
					// Non-hidden items satisfy R > R_prev || ri > ri_prev. The first condition is always true as we sort the list by R.
					rilist.clear();
					rilist.push_back(resultTmp2[0]);
					//if (dim < ri.dimensionality() - 1) // In the last dimension only the item with maximum R is needed, but we have separate function for that
					//{
					for (size_t n = 1; n < resultTmp2.size(); n++)
					{
						pixel_t currri2 = rilist[rilist.size() - 1].ri2;
						pixel_t newri2 = resultTmp2[n].ri2;

						if (newri2 > currri2) // This is the basic condition that works always (but may include unnecessary items in the rilist)
						{
							if (dim == dimensionality - 2) // 3 - 2 == 1 == 2nd dimension
							{
								// In the second last dimension only really visible spans are needed as there's no next dimension whose ri we would calculate based on the spans.
								if (largestIntWhoseSquareIsLessThan(newri2) > largestIntWhoseSquareIsLessThan(currri2))
									rilist.push_back(resultTmp2[n]);
							}
							else if (dim == dimensionality - 3) // 3 - 3 == 0 == 1st dimension
							{
								// In the third last dimension we know that only those spans are required that produce visible circles in the output.
								if (!internals::doesDiscretizedCircle1FitInto2Cached((int32_t)newri2, (int32_t)currri2))
									rilist.push_back(resultTmp2[n]);
							}
							else
							{
								// Here we could insert test if discretized spheres fit into each other etc. etc.
								rilist.push_back(resultTmp2[n]);
							}
						}
					}
					//}

					// Make sure we don't use extraneous memory by storing empty items in each pixel of the ri image.
					rilist.shrink_to_fit();

					//debugCheckRiSetRedundancy(rilist, dim);
				}
			}
		}

		/*
		Same than singlePass but generates result image instead of ri image, and is thus a little bit faster.
		@param centers Image containing only the row to be processed before first call to this method.
		@param Result image. Pixels must be set to zero before first call to this method.
		@param rowStart Start point of the pixel row to be processed.
		@param dim Dimension that we are processing.
		@param step +1 or -1 to indicate the direction of the pass.
		@param initialActiveSpheres Initial active spheres list. Contains final active spheres list at output. Set to nullptr to assume empty list.
		*/
		template<typename pixel_t> void singlePassFinalSuper(const std::vector<internals::RiSuperSet<pixel_t> >& centers, Image<pixel_t>& result, const Vec3c& rowStart, size_t dim, coord_t step)
		{
			using signed_t = typename NumberUtils<pixel_t>::SignedType;

			// Stores the spheres that have been encountered and that have not been passed yet.
			// Stores the sphere with the largest R at the top of the priority queue.
			std::priority_queue<internals::ActiveSpheresSuperItem<pixel_t>, std::vector<internals::ActiveSpheresSuperItem<pixel_t> >, std::function<decltype(internals::activeSpheresSuperSorter2<pixel_t>)> > activeSpheres(internals::activeSpheresSuperSorter2<pixel_t>);

			// Set start point to the start or end of the current row.
			Vec3c p = rowStart;
			if (step < 0)
				p[dim] += (result.dimension(dim) - 1);

			for (coord_t i = 0; i < result.dimension(dim); i++, p[dim] += step)
			{
				signed_t x = (signed_t)p[dim];

				// If there is one or more sphere centers at the current location, add them to the set of active spheres.
				const internals::RiSuperSet<pixel_t>& C = centers[p[dim]];
				for (auto& item : C)
				{
					activeSpheres.push({ item.R2, item.ri2, x });
				}

				while (!activeSpheres.empty())
				{
					auto& item = activeSpheres.top();

					pixel_t Rorig2 = item.R2;
					pixel_t R2 = item.ri2;
					signed_t cx = item.xi;

					signed_t dx = abs(x - cx);

					// Calculate new ri^2
					signed_t rn2 = R2 - dx * dx;
					if (rn2 > 0)
					{
						// Note that previous pass may have assigned larger value to the output.
						if (Rorig2 > result(p))
							result(p) = Rorig2;
						break;
					}
					else
					{
						// ry is non-positive, i.e. dx >= R
						// This sphere is not active anymore, so remove it from the list of active spheres.
						activeSpheres.pop();
					}

				}

			}
		}

		//size_t readsFromBlock = 0;
		//size_t readsFromFull = 0;

		/**
		Converts RiStorageItem to RiSuperItem.
		@param p Position (in the block) where the item is taken from.
		*/
		template<typename pixel_t> RiSuperItem<pixel_t> toRiItem(const RiStorageItem& si, const Vec3c& p, const Image<pixel_t>& dmap2, const Vec3sc& blockPos, const Image<pixel_t>& dmap2Full)
		{
			RiSuperItem<pixel_t> out;

			pixel_t R2;

			// This works but reads always from the full image:
			//R2 = dmap2Full(si.srcX, si.srcY, p.z + currBlock.position().z);
			//int32_t dx = (int32_t)p.x - (si.srcX - currBlock.position().x);
			//int32_t dy = (int32_t)p.y - (si.srcY - currBlock.position().y);

			// This version reads from block when possible, and from the full image when not.
			Vec3sc pBlock(si.srcX - blockPos.x, si.srcY - blockPos.y, (int32_t)p.z);
			if (dmap2.isInImage(pBlock))
			{
				R2 = dmap2(pBlock);
				//readsFromBlock++;
			}
			else
			{
				R2 = dmap2Full(si.srcX, si.srcY, (int32_t)p.z + blockPos.z);
				//readsFromFull++;
			}
			int32_t dx = (int32_t)p.x - pBlock.x;
			int32_t dy = (int32_t)p.y - pBlock.y;

			// This is the original non-block supporting version
			//int32_t dx = (int32_t)p.x - si.srcX;
			//int32_t dy = (int32_t)p.y - si.srcY;

			pixel_t ri2 = R2 - dx * dx - dy * dy;

			//if (R2 != si.R2 || ri2 != si.ri2)
			//{
			//	cout << "Error" << endl;
			//}

			out.R2 = R2;
			out.ri2 = ri2;

			out.srcX = si.srcX;
			out.srcY = si.srcY;
			return out;
		}

		/**
		Converts RiSuperItem to RiStorageItem
		*/
		template<typename pixel_t> RiStorageItem toRiStorageItem(const RiSuperItem<pixel_t>& item)
		{
			RiStorageItem out;
			out.srcX = item.srcX;
			out.srcY = item.srcY;
			return out;
		}

		/**
		Converts RiStorageSet to RiSuperSet.
		@param in Source set
		@param p Position (in the block) where the in set is taken from.
		*/
		template<typename pixel_t> void toRiSet(const RiStorageSet& in, RiSuperSet<pixel_t>& out, const Vec3c& p, const Image<pixel_t>& dmap2, const Vec3sc& blockPos, const Image<pixel_t>& dmap2Full)
		{
			out.clear();
			for (uint16_t n = 0; n < in.size(); n++)
				out.push_back(toRiItem(in[n], p, dmap2, blockPos, dmap2Full));
		}

		/**
		Converts RiSuperSet to RiStorageSet.
		*/
		template<typename pixel_t> void toStorageSet(const RiSuperSet<pixel_t>& in, RiStorageSet& out)
		{
			std::vector<RiStorageItem> temp;
			for (auto& item : in)
				temp.push_back(toRiStorageItem(item));
			out.setValues(temp);
		}

		/*
		Makes one pass over image in specific dimension and direction.
		@param ri Image containing ri values from processing of previous dimension.
		@param dim Dimension to process.
		@param result Result image.
		@param showProgressInfo Set to true to show progress indicator.
		*/
		template<typename pixel_t> void processDimensionSuper(Image<internals::RiStorageSet>& ri, size_t dim,
			const Image<pixel_t>& dmap2,
			Image<pixel_t>& result, bool showProgressInfo,
			const AABox<coord_t>& currBlock, const Image<pixel_t>& dmap2Full)
		{
			// Determine count of pixels to process
			Vec3c reducedDimensions = internals::getReducedDimensions(ri.dimensions(), dim);
			coord_t rowCount = reducedDimensions.x * reducedDimensions.y * reducedDimensions.z;


			Vec3sc blockPos(currBlock.position());

			bool isFinalPass = !(dim < dmap2Full.dimensionality() - 1);

			size_t counter = 0;
#pragma omp parallel if(!omp_in_parallel() && ri.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				// Temporary buffer
				std::vector<internals::RiSuperSet<pixel_t> > inRow(ri.dimension(dim));
				std::vector<internals::RiSuperSet<pixel_t> > outRow(ri.dimension(dim));

				// Determine start points of pixel rows and process each row
#pragma omp for schedule(dynamic)
				for (coord_t n = 0; n < rowCount; n++)
				{
					Vec3c start = indexToCoords(n, reducedDimensions);

					// Make a copy of the current row as we update the row in the forward pass but need
					// the original data in the backward pass.
					Vec3c pos = start;
					for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
					{
						toRiSet(ri(pos), inRow[x], pos, dmap2, blockPos, dmap2Full);
						outRow[x].clear();
						//ri(pos).clear();
					}


					if (!isFinalPass)
					{
						singlePassSuper(inRow, outRow, dmap2Full.dimensionality(), ri.dimensions(), start, dim, 1);
						singlePassSuper(inRow, outRow, dmap2Full.dimensionality(), ri.dimensions(), start, dim, -1);

						// Copy data back to storage
						pos = start;
						for (coord_t x = 0; x < ri.dimension(dim); x++, pos[dim]++)
						{
							toStorageSet(outRow[x], ri(pos));
						}
					}
					else
					{
						singlePassFinalSuper<pixel_t>(inRow, result, start, dim, 1);
						singlePassFinalSuper<pixel_t>(inRow, result, start, dim, -1);
					}

					showThreadProgress(counter, rowCount, showProgressInfo);
				}
			}
		}


		/*
		Copies squared sphere radius data from original centers image to ri image.
		@param centers2 Original distance ridge.
		@param ri Temporary image that is to be initialized.
		@param bounds If the ri covers only a block of centers2, a box defining the block. Pass AABox from [0, 0, 0] to centers2.dimensions() to prepare whole image.
		*/
		template<typename pixel_t> void prepareSuper(const Image<pixel_t>& centers2, Image<internals::RiStorageSet>& ri, const AABox<coord_t>& bounds)
		{
			if (centers2.dimensions().max() >= std::numeric_limits<int16_t>::max())
				throw ITLException(string("Linear image size exceeds ") + toString(std::numeric_limits<int16_t>::max()) + " pixels. This implementation is not configured for that big images.");

			ri.ensureSize(bounds.size());

			for (coord_t z = 0; z < ri.depth(); z++)
			{
				for (coord_t y = 0; y < ri.height(); y++)
				{
					for (coord_t x = 0; x < ri.width(); x++)
					{
						Vec3c p(x, y, z);
						pixel_t R2 = centers2(bounds.minc + p);
						if (R2 > 0)
						{
							ri(p).setValues(RiStorageItem{ (int16_t)(x + bounds.position().x), (int16_t)(y + bounds.position().y) });
						}
					}
				}
			}
		}
	}


	namespace dimredsuper
	{
		/*
		Calculate squared local radius from squared distance map.
		@param dmap2 Squared distance map.
		@param tmap2 At output, squared radius map.
		@param extraBytes Approximation of memory used in addition to the input image.
		*/
		template<typename pixel_t> void thickmap2(const Image<pixel_t>& dmap2, Image<pixel_t>& tmap2, double* extraBytes = nullptr, Vec3d* riCounts = nullptr, bool showProgressInfo = true)
		{
			tmap2.mustNotBe(dmap2);

			pixel_t M = max(dmap2);
			if (intuitive::ge(M, std::numeric_limits<int32_t>::max()))
				throw ITLException("The squared distance map contains too large values. (This error is easily avoidable by changing buildCircleLookup functionality.)");
			internals::buildCircleLookup((int32_t)M);

			tmap2.ensureSize(dmap2);

			Image<internals::RiStorageSet> ri(dmap2.dimensions());
			internals::prepareSuper(dmap2, ri, AABox<coord_t>(Vec3c(0, 0, 0), dmap2.dimensions()));
			setValue(tmap2, 0);
			double totalRiMem = 0;
			for (size_t n = 0; n < ri.dimensionality(); n++)
			{
				internals::processDimensionSuper<pixel_t>(ri, n, dmap2, tmap2, showProgressInfo, AABox<coord_t>(Vec3c(0, 0, 0), dmap2.dimensions()), dmap2);

				if (extraBytes || riCounts)
				{
					double memSum = 0;
					double sum = 0;

					for (coord_t nn = 0; nn < ri.pixelCount(); nn++)
					{
						memSum += ri(nn).memSize();
						sum += ri(nn).size();
					}

					if (memSum > totalRiMem)
						totalRiMem = memSum;

					if (riCounts)
						(*riCounts)[n] = sum;
				}
			}

			if (extraBytes)
			{
				*extraBytes = totalRiMem + tmap2.pixelCount() * tmap2.pixelSize();
			}
		}
	}








	namespace internals
	{


		template<typename pixel_t> bool ridgeSorter(const std::tuple<pixel_t, coord_t>& a, const std::tuple<pixel_t, coord_t>& b)
		{
			return std::get<0>(a) > std::get<0>(b);
		}

		inline bool isBitSet(uint8_t* bitarray, size_t index)
		{
			return ((bitarray[index / 8] >> (index & 0x7)) & 1U) != 0;
		}

		inline void clearBitSafe(uint8_t* bitarray, size_t index)
		{
			uint8_t* ptr = &bitarray[index / 8];
			uint8_t i = ~(uint8_t)(1UL << (index & 0x7));
#pragma omp atomic
			*ptr &= i;
		}


		template<typename pixel_t> void fillSpheres(const std::vector<std::tuple<pixel_t, coord_t> >& points, coord_t starti, coord_t endi, Image<pixel_t>& tmap, Image<uint8_t>& resultBitMask)
		{
			// NOTE: This function supports only R2 values that can be converted to int32

			pixel_t R2orig = std::get<0>(points[starti]);
			int32_t R2 = (int32_t)R2orig;
			int32_t Rint = (int32_t)ceil(sqrt(R2));

			// Create mask
			coord_t size = 2 * Rint + 1;
			Image<int32_t> xMask(size, size, 1, -1);
#pragma omp parallel for if(xMask.pixelCount() > PARALLELIZATION_THRESHOLD)
			for (int32_t z = 0; z < xMask.height(); z++)
			{
				int32_t dz = z - Rint;
				for (int32_t y = 0; y < xMask.width(); y++)
				{
					int32_t dy = y - Rint;
					int32_t dx2 = R2 - dy * dy - dz * dz;
					if (dx2 >= 0)
					{
						xMask(y, z) = largestIntWhoseSquareIsLessThan(dx2);
					}
				}
			}

			// Fill all spheres
			Vec3c maxCoords = tmap.dimensions() - Vec3c(1, 1, 1);
			bool centersParallel = endi - starti >= omp_get_max_threads();
#pragma omp parallel for if(centersParallel)
			for (coord_t n = starti; n < endi; n++)
			{
				Vec3c pos = tmap.getCoords(std::get<1>(points[n]));

				Vec3c start0 = pos - Vec3c(Rint, Rint, Rint);
				Vec3c end = pos + Vec3c(Rint, Rint, Rint);

				Vec3c start = start0;
				clamp(start, Vec3c(), maxCoords);
				clamp(end, Vec3c(), maxCoords);

#pragma omp parallel for if(!omp_in_parallel() && (end.z - start.z + 1) > omp_get_max_threads())
				for (coord_t z = start.z; z <= end.z; z++)
				{
					for (coord_t y = start.y; y <= end.y; y++)
					{
						coord_t xr = xMask(y - start0.y, z - start0.z);

						if (xr >= 0)
						{
							coord_t xStart = pos.x - xr;
							if (xStart < 0)
								xStart = 0;

							coord_t xEnd = pos.x + xr;
							if (xEnd > maxCoords.x)
								xEnd = maxCoords.x;

							// Loop over bytes in the result bit mask
							coord_t biInd = (coord_t)resultBitMask.getLinearIndex(0, y, z);
							uint8_t* bitRowStart = &resultBitMask(biInd);
							for (coord_t bi = xStart / 8; bi <= xEnd / 8; bi++)
							{
								if (resultBitMask(biInd + bi) != 0) // No pragma omp atomic - possible erroneous intermediate value read can only be != 0 and in that case we check the individual bits again below.
								{
									coord_t x0 = std::max(8 * bi, xStart);
									coord_t x1 = std::min(x0 + 8, xEnd + 1);
									coord_t ind = (coord_t)tmap.getLinearIndex(x0, y, z);
									for (coord_t x = x0; x < x1; x++)
									{
										if (isBitSet(bitRowStart, x))
										{
											tmap(ind) = R2orig;

											clearBitSafe(bitRowStart, x);
										}
										ind++;
									}
								}
							}
						}
					}
				}
			}
		}
	}



	namespace optimized
	{
		/**
		Calculates squared local radius from squared distance map
		Uses standard Hildebrand & Ruegsegger algorithm.
		Plots maximal spheres corresponding to squared distance map, larger distance values replacing smaller ones.
		@param extraBytes Approximation of memory used in addition to the input image.
		*/
		template<typename pixel_t> void thickmap2(Image<pixel_t>& dmap2, double* extraBytes = nullptr, bool showProgressInfo = true)
		{
			// In order to parallelize by sorting method, first find radii and locations (linear indices to save some memory) of non-zero pixels.
			std::vector<std::tuple<pixel_t, coord_t> > centers;
			for (coord_t n = 0; n < dmap2.pixelCount(); n++)
			{
				if (dmap2(n) > 0)
					centers.push_back(std::make_tuple(dmap2(n), n));
			}

			// Sort according to sphere diameter
#if defined(_WIN32)
			std::sort(std::execution::par_unseq, centers.begin(), centers.end(), internals::ridgeSorter<pixel_t>);
#elif defined(__linux__)
			__gnu_parallel::sort(centers.begin(), centers.end(), internals::ridgeSorter<pixel_t>);
#elif defined(__APPLE__)
			std::sort(centers.begin(), centers.end(), internals::ridgeSorter<pixel_t>);
#endif

			// Bitmask where filled pixels are marked. This enables testing 8 filled pixels at once
			// (and never read from output)
			Image<uint8_t> bitmask(dmap2.width() / 8 + 1, dmap2.height(), dmap2.depth(), 255);

			if (extraBytes)
			{
				(*extraBytes) = (double)(centers.size() * sizeof(std::tuple<int32_t, coord_t>) + bitmask.pixelCount() * sizeof(uint8_t));
			}

			coord_t previ = 0;
			for (coord_t n = 0; n < (coord_t)centers.size(); n++)
			{
				pixel_t prevR = std::get<0>(centers[previ]);
				pixel_t r = std::get<0>(centers[n]);

				if (r != prevR)
				{
					// Draw spheres from previ to n-1
					internals::fillSpheres<pixel_t>(centers, previ, n, dmap2, bitmask);
					previ = n;
				}

				showProgress(n, centers.size(), showProgressInfo);
			}

			// Draw the final span of spheres from previ to end of ridgePoints list.
			internals::fillSpheres<pixel_t>(centers, previ, (coord_t)centers.size(), dmap2, bitmask);
		}
	}



	namespace standard
	{
		/**
		Calculates squared local radius from squared distance map
		Uses standard Hildebrand & Ruegsegger algorithm without optimizations.
		Plots maximal spheres corresponding to squared distance map, larger distance values replacing smaller ones.
		*/
		void thickmap2(const Image<int32_t>& dmap2, Image<int32_t>& result, bool showProgressInfo = true);
	}


	/**
	Rounds squared distance map such that (non-squared) distance values are integers.
	*/
	template<typename in_t> void roundDistanceRidge2(Image<in_t>& ridge2)
	{
#pragma omp parallel for if(ridge2.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t n = 0; n < ridge2.pixelCount(); n++)
		{
			double v = (double)ridge2(n);
			v = std::round(sqrt(v));
			ridge2(n) = pixelRound<in_t>(v * v);
		}
	}


	/**
	Convert squared local radius to thickness map (takes square root and multiplies by 2).
	*/
	template<typename in_t, typename out_t> void finalizeThickmap(const Image<in_t>& thickmap2, Image<out_t>& thickmap)
	{
		thickmap.ensureSize(thickmap2);

#pragma omp parallel for if(thickmap2.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t n = 0; n < thickmap.pixelCount(); n++)
			thickmap(n) = pixelRound<out_t>(2.0 * sqrt(thickmap2(n)));
	}




	/**
	These utilities are used in distributed processing.
	*/
	namespace internals
	{

		/**
		Writes ri image block to file.
		*/
		void writeRiBlock(Image<internals::RiStorageSet>& ri, const string& indexFilePrefix, uint16_t blockIndex, const Vec3c& filePosition, const Vec3c& fileDimensions);

		/**
		Reads ri image block from file.
		*/
		void readRiBlock(Image<internals::RiStorageSet>& ri, std::string indexFilePrefix, const Vec3c& start);
	}


	/**
	Calculates thickness map from geometry image.
	@param input Binary image containing the geometry. The image must be of data type that can hold large enough squared distance values (e.g. uint32 is often suitable). This image is used as temporary storage space and its contents are replaced by (rounded) distance ridge.
	@param output Output image that will store the thickness map. For exact results, output data type should be floating-point.
	@param roundedTmap Set to true to round thickness values to nearest integer. Rounding increases performance.
	@param bgVal Value corresponding to background.
	*/
	template<typename pixel_t, typename out_t> void thicknessMap(Image<pixel_t>& input, Image<out_t>& output, bool roundedTmap = false, pixel_t bgVal = (pixel_t)0)
	{
		{
			Image<pixel_t> temp(input.dimensions());
			distanceTransform2(input, input, nullptr, bgVal);
			centersOfLocallyMaximalSpheres(input, temp);
			if (roundedTmap)
				roundDistanceRidge2(temp);
			dimredsuper::thickmap2(temp, input);
		} // temp image is deleted from memory here, and output image is created below if it does not exist.
		finalizeThickmap(input, output);
	}


	
	namespace tests
	{
		void thickmapsEquality();
		void thickmapScaling();
		void dimred2D();
		void dimred3D();
		void discretizedCircles();
		void indexItem();
		void readWriteRi();
		void testDataForFiji();
		void thickmapRounding();
		void thickmapScaling();
	}


}
