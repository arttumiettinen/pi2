#pragma once

#include <vector>

#include "image.h"
#include "math/vec3.h"
#include "neighbourhood.h"
#include "fft.h"
#include "transform.h"
#include "inpaint.h"
#include "io/raw.h"
#include "projections.h"

using math::Vec3;
using math::Vec3d;
using std::vector;

namespace itl2
{
	
	//namespace internals
	//{
	//	inline void clampedBlockRadius(const Vec3c& p, const Vec3c& blockRadius, const Vec3c& dimensions, Vec3c& pos, Vec3c& r)
	//	{
	//		Vec3c m = p - blockRadius;
	//		Vec3c M = p + blockRadius;
	//		clamp(m, Vec3c(0, 0, 0), dimensions);
	//		clamp(M, Vec3c(0, 0, 0), dimensions);
	//		r = (M - m - Vec3c(1, 1, 1)) / 2;
	//		pos = m + r;
	//	}
	//}

	namespace internals
	{
		/*
		Block matcher.
		@param reference Reference image.
		@param deformed Deformed image.
		@param blockRadius Radius of matching block. The value is also maximum change in defPoint that can be found.
		@param refPoint Point in the reference image.
		@param defPoint Point in the deformed image. Input value is used as initial guess of the shift. On output, will contain the shift that was estimated.
		@param accuracy Stores a measure of the accuracy of the matching.
		*/
		template<typename ref_t, typename def_t> void blockMatchOnePoint(const Image<ref_t>& reference, const Image<def_t>& deformed, const Vec3c& blockRadius, const Vec3c& refPoint, Vec3d& defPoint, double& accuracy, size_t binningSize = 1)
		{
			Vec3c r = blockRadius;
			for (size_t n = reference.dimensionality(); n < 3; n++)
				r[n] = 0;

			Vec3c blockSize = 2 * r + Vec3c(1, 1, 1);

			Vec3c defPointRounded = round(defPoint);

			Image<float32_t> refBlock;
			Image<float32_t> defBlock;
			if (binningSize > 1)
			{
				Image<float32_t> refBlockOrig(blockSize);
				Image<float32_t> defBlockOrig(blockSize);

				getNeighbourhood(reference, refPoint, r, refBlockOrig, Zero);
				getNeighbourhood(deformed, defPointRounded, r, defBlockOrig, Zero);

				//raw::writed(refBlockOrig, "./block_ref_nobinning");
				//raw::writed(defBlockOrig, "./block_def_nobinning");

				maskedBinning(refBlockOrig, refBlock, binningSize, (float32_t)0, (float32_t)0, false);
				maskedBinning(defBlockOrig, defBlock, binningSize, (float32_t)0, (float32_t)0, false);
			}
			else
			{
				refBlock.ensureSize(blockSize);
				defBlock.ensureSize(blockSize);
				getNeighbourhood(reference, refPoint, r, refBlock, Zero);
				getNeighbourhood(deformed, defPointRounded, r, defBlock, Zero);
			}

			// For testing
			//raw::writed(refBlock, "./block_ref_noinpaint");
			//raw::writed(defBlock, "./block_def_noinpaint");


			// Set zeros to nearest non-zero value. This has effect particularly in the edges and corners of non-rectangular images.
			inpaintNearest(refBlock);
			inpaintNearest(defBlock);
			//inpaintGarcia(refBlock);
			//inpaintGarcia(defBlock);

			// For testing
			//raw::writed(refBlock, "./block_ref");
			//raw::writed(defBlock, "./block_def");

			Vec3d shift = phaseCorrelation(refBlock, defBlock, r / binningSize, accuracy);
			shift *= (double)binningSize;

			defPoint = Vec3d(defPointRounded) - shift;
		}

		/*
		Block matcher that first block matches with low resolution and then improves the result by block matching with full resolution.
		@param reference Reference image.
		@param deformed Deformed image.
		@param coarseBlockRadius Radius of matching block in the coarse matching phase. The value is also maximum change in defPoint that can be found in this phase.
		@param coarseBinning Binning applied to the coarse block to lower resolution. If set to 1, no fine resolution matching is done.
		@param fineBlockRadius Radius of matching block in the full-resolution matching phase. The value is also maximum change that can be registered on top of the result of the coarse matching phase.
		@param refPoint Point in the reference image.
		@param defPoint Point in the deformed image. Input value is used as initial guess of the shift. On output, will contain the shift that was estimated.
		@param accuracy Stores a measure of the accuracy of the matching.
		*/
		template<typename ref_t, typename def_t> void blockMatchOnePointMultires(const Image<ref_t>& reference, const Image<def_t>& deformed, const Vec3c& coarseBlockRadius, size_t coarseBinning, const Vec3c& fineBlockRadius, size_t fineBinning, const Vec3c& refPoint, Vec3d& defPoint, double& accuracy)
		{
			blockMatchOnePoint(reference, deformed, coarseBlockRadius, refPoint, defPoint, accuracy, coarseBinning);

			if(accuracy > 0 && coarseBinning > fineBinning)
				blockMatchOnePoint(reference, deformed, fineBlockRadius, refPoint, defPoint, accuracy, fineBinning);
		}
	}

	/*
	Phase correlation-based block matching.
	@param reference Reference image. Zeros are considered empty values.
	@param deformed Deformed image. Zeros are considered empty values.
	@param refPoints Points in reference image whose locations in deformed image should be determined.
	@param defPoints Locations of points in deformed image corresponding to reference points. As input this vector contains initial guess of the point locations.
	@param searchRadius Radius of search region around initial guess.
	@param compRadius Radius of a block of reference image extracted at each calculation point for comparing between reference and deformed.
	*/
	template<typename ref_t, typename def_t> void blockMatch(const Image<ref_t>& reference, const Image<def_t>& deformed, const vector<Vec3c>& refPoints, vector<Vec3d>& defPoints, vector<double>& accuracy, const Vec3c& blockRadius)
	{
		while (defPoints.size() < refPoints.size())
			defPoints.push_back(Vec3d());

		while (accuracy.size() < refPoints.size())
			accuracy.push_back(0);

		//Vec3c blockSize = 2 * blockRadius + Vec3c(1, 1, 1);

		size_t counter = 0;
		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t n = 0; n < (coord_t)refPoints.size(); n++)
		{
			internals::blockMatchOnePoint(reference, deformed, blockRadius, refPoints[n], defPoints[n], accuracy[n]);
			//Vec3c refPoint = refPoints[n];
			//Vec3c defPointRounded = round(defPoints[n]);

			//Image<float32_t> refBlock(blockSize);
			//Image<float32_t> defBlock(blockSize);

			//getNeighbourhoodClamp(reference, refPoint, blockRadius, refBlock);
			//getNeighbourhoodClamp(deformed, defPointRounded, blockRadius, defBlock);

			//// Prefetch next cubes. NOTE: This is not correct if the next block is processed by some other thread, but
			//// that should not be a big problem.
			//if (n < (coord_t)refPoints.size() - 1)
			//{
			//	reference.prefetch(refPoints[n + 1] - blockRadius, refPoints[n + 1] + blockRadius);
			//	deformed.prefetch(round(defPoints[n + 1]) - blockRadius, round(defPoints[n + 1]) + blockRadius);
			//}
			//
			//// Set zeros to nearest non-zero value. This has effect particularly in the edges and corners of non-rectangular images.
			//inpaintZerosNearest(refBlock);
			//inpaintZerosNearest(defBlock);

			//// For testing
			////raw::writed(refBlock, "./block_ref");
			////raw::writed(defBlock, "./block_def");

			//double goodness;
			//Vec3d shift = phaseCorrelation(refBlock, defBlock, blockSize / 2, goodness);

			//accuracy[n] = goodness;
			//defPoints[n] -= shift;

			showThreadProgress(counter, refPoints.size());
		}
	}

	/*
	1D point grid. Stores min, max and step between points. First point is at 'first'.
	*/
	template<typename T> class PointGrid1D
	{
	public:
		T first;
		T maximum;
		T step;

		PointGrid1D(T first, T maximum, T step) :
			first(first),
			maximum(maximum),
			step(step)
		{
		}

		PointGrid1D(): PointGrid1D(0, 0, 0)
		{

		}

		T operator()(coord_t index) const
		{
			//return (T)::round((double)first + (double)index * step);
			return first + index * step;
		}

		/*
		Calculates step between points in the grid.
		*/
		//double step() const
		//{
		//	return (double)(last - first) / (double)(count - 1);
		//}

		/*
		Gets point position at given index.
		No bounds checking is made.
		*/
		double getPosition(double index) const
		{
			//return first + index / (count - 1) * (last - first);
			return first + index * step;
		}

		/*
		Get index of point at given location.
		*/
		double getIndex(double point) const
		{
			//return 0 + (val - first) / (last - first) * (count - 1);
			return (point - first) / step;
		}

		/*
		Tests if the given position is inside the grid (between first and last grid points)
		*/
		template<typename real_t> bool contains(real_t val) const
		{
			//return first <= (double)val && (double)val < last;
			return first <= (double)val && (double)val <= first + (pointCount() - 1) * step;
		}

		/*
		Gets count of points in this grid.
		*/
		coord_t pointCount() const
		{
			//return count;
			return (coord_t)::floor(getIndex((double)maximum)) + 1;
		}
	};

	/*
	Represents a grid of points.
	*/
	template<typename T> class PointGrid3D
	{
	public:
		PointGrid1D<T> xg;
		PointGrid1D<T> yg;
		PointGrid1D<T> zg;
	
		PointGrid3D() : PointGrid3D(PointGrid1D<T>(), PointGrid1D<T>(), PointGrid1D<T>())
		{

		}
		
		PointGrid3D(const PointGrid1D<T>& xGrid, const PointGrid1D<T>& yGrid, const PointGrid1D<T>& zGrid) :
			xg(xGrid), yg(yGrid), zg(zGrid)
		{

		}

		/*
		Returns coordinates of point whose index is given.
		*/
		Vec3<T> operator()(coord_t x, coord_t y, coord_t z) const
		{
			return Vec3<T>(xg(x), yg(y), zg(z));
		}

		Vec3c pointCounts() const
		{
			return Vec3c((coord_t)xg.pointCount(), (coord_t)yg.pointCount(), (coord_t)zg.pointCount());
		}

		/*
		Gets count of points in the grid.
		*/
		size_t pointCount() const
		{
			return xg.pointCount() * yg.pointCount() * zg.pointCount();
		}

		/*
		Tests if the given point is inside the grid.
		*/
		template<typename real_t> bool contains(const Vec3<real_t>& x) const
		{
			return xg.contains(x.x) && yg.contains(x.y) && zg.contains(x.z);
		}

		/*
		Calculates locations of all points in this grid and adds them to the given list.
		*/
		void getAllPoints(vector<Vec3c>& out) const
		{
			for (coord_t zi = 0; zi < zg.pointCount(); zi++)
			{
				for (coord_t yi = 0; yi < yg.pointCount(); yi++)
				{
					for (coord_t xi = 0; xi < xg.pointCount(); xi++)
					{
						out.push_back((*this)(xi, yi, zi));
					}
				}
			}
		}
	};

	///*
	//Block matching for point grid and image output.
	//*/
	//template<typename ref_t, typename def_t> void blockMatch(const Image<ref_t>& reference, const Image<def_t>& deformed, const PointGrid3D<coord_t>& refGrid, Image<Vec3d>& defPoints, Image<float32_t>& accuracy, const Vec3c& blockRadius)
	//{
	//	vector<Vec3c> refs;
	//	vector<Vec3d> defs;
	//	vector<double> accs;

	//	refs.reserve(refGrid.pointCount());
	//	refGrid.getAllPoints(refs);

	//	accuracy.ensureSize(refGrid.pointCounts());
	//	accs.reserve(refGrid.pointCount());

	//	defPoints.ensureSize(refGrid.pointCounts());
	//	defs.reserve(defPoints.pixelCount());
	//	for (coord_t n = 0; n < defPoints.pixelCount(); n++)
	//		defs.push_back(defPoints(n));

	//	// TODO: If conversion between point list and image causes performance issues, migrate grid and image output in actual blockMatch code.
	//	blockMatch(reference, deformed, refs, defs, accs, blockRadius);

	//	size_t n = 0;
	//	for (coord_t z = 0; z < defPoints.depth(); z++)
	//	{
	//		for (coord_t y = 0; y < defPoints.height(); y++)
	//		{
	//			for (coord_t x = 0; x < defPoints.width(); x++)
	//			{
	//				defPoints(x, y, z) = defs[n];
	//				accuracy(x, y, z) = (float32_t)accs[n];
	//				n++;
	//			}
	//		}
	//	}
	//}

	/*
	Block matching for point grid and image output.
	*/
	template<typename ref_t, typename def_t> void blockMatch(const Image<ref_t>& reference, const Image<def_t>& deformed, const PointGrid3D<coord_t>& refGrid, Image<Vec3d>& defPoints, Image<float32_t>& accuracy, const Vec3c& blockRadius)
	{
		accuracy.ensureSize(refGrid.pointCounts());
		defPoints.ensureSize(refGrid.pointCounts());

		size_t counter = 0;
		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t z = 0; z < defPoints.depth(); z++)
		{
			for (coord_t y = 0; y < defPoints.height(); y++)
			{
				//reference.prefetch(Vec3c(0, y, z) - blockRadius, Vec3c(reference.width(), y, z) + blockRadius);
				//deformed.prefetch(round(defPoints(0, y, z)) - 2 * blockRadius, round(defPoints(defPoints.width(), y, z)) + 2 * blockRadius);

				for (coord_t x = 0; x < defPoints.width(); x++)
				{
					Vec3c refPoint = refGrid(x, y, z);
					Vec3d defPoint = defPoints(x, y, z);
					double gof;

					internals::blockMatchOnePoint(reference, deformed, blockRadius, refPoint, defPoint, gof);

					defPoints(x, y, z) = defPoint;
					accuracy(x, y, z) = (float32_t)gof;
				}

				showThreadProgress(counter, defPoints.depth() * defPoints.height());
			}
		}
	}

	/*
	Block matching for point grid and image output, loads images only partially.
	*/
	template<typename ref_t, typename def_t> void blockMatchPartialLoad(const string& referenceFile, const Vec3c& refImageDimensions, const string& deformedFile, const Vec3c& defImageDimensions, const PointGrid3D<coord_t>& refGrid, Image<Vec3d>& defPoints, Image<float32_t>& accuracy,
		const Vec3c& coarseBlockRadius, size_t coarseBinning,
		const Vec3c& fineBlockRadius, size_t fineBinning,
		bool normalize, double& normFact)
	{
		accuracy.ensureSize(refGrid.pointCounts());
		defPoints.ensureSize(refGrid.pointCounts());

		Vec3c defStart = round(defPoints(0, 0, 0)) - coarseBlockRadius;
		Vec3c defEnd = round(defPoints(defPoints.dimensions() - Vec3c(1, 1, 1))) + coarseBlockRadius + Vec3c(1, 1, 1);
		Vec3c refStart = refGrid(0, 0, 0) - coarseBlockRadius;
		Vec3c refEnd = refGrid(refGrid.xg.pointCount() - 1, refGrid.yg.pointCount() - 1, refGrid.zg.pointCount() - 1) + coarseBlockRadius + Vec3c(1, 1, 1);

		clamp(defStart, Vec3c(0, 0, 0), defImageDimensions);
		clamp(defEnd, Vec3c(0, 0, 0), defImageDimensions);
		clamp(refStart, Vec3c(0, 0, 0), refImageDimensions);
		clamp(refEnd, Vec3c(0, 0, 0), refImageDimensions);

		// TODO: debug
		//while (defEnd.z - defStart.z < refEnd.z - refStart.z)
		//	defEnd.z++;

		Image<ref_t> referenceBlock(refEnd - refStart);
		Image<def_t> deformedBlock(defEnd - defStart);

		cout << "Loading block of reference image, size = " << (referenceBlock.pixelCount() * sizeof(ref_t) / (1024 * 1024)) << " MiB." << endl;
		cout << "Loading block of deformed image, size  = " << (deformedBlock.pixelCount() * sizeof(ref_t) / (1024 * 1024)) << " MiB." << endl;

		raw::readBlock(referenceBlock, referenceFile, refImageDimensions, refStart, true);
		raw::readBlock(deformedBlock, deformedFile, defImageDimensions, defStart, true);

		// Calculate normalization factor for gray values
		normFact = 1;
		double meanRef = maskedmean(referenceBlock, (ref_t)0);
		double meanDef = maskedmean(deformedBlock, (def_t)0);
		//normFact = 1 / meanDef * meanRef;
		normFact = meanRef - meanDef;
		if (normalize)
		{
			//if (!NumberUtils<double>::equals(normFact, 0))
			//	multiply(deformedBlock, normFact);
			maskedAdd(deformedBlock, normFact, (def_t)0);
		}

		
		// Initial translation test using MIP matching technique
		//Vec3d mipTranslation = mipMatch(referenceBlock, deformedBlock);
		//for (int i = 0; i < 3; i++)
		//	if (-initialShift[i] > reference.dimension(i) / 2)
		//		mipTranslation[i] = reference.dimension(i) + mipTranslation[i];
		//cout << "Initial translation = " << mipTranslation << endl;


		size_t counter = 0;
		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t z = 0; z < defPoints.depth(); z++)
		{
			for (coord_t y = 0; y < defPoints.height(); y++)
			{
				for (coord_t x = 0; x < defPoints.width(); x++)
				{
					Vec3c refPoint = refGrid(x, y, z) - refStart;
					Vec3d defPoint = defPoints(x, y, z) - Vec3d(defStart);
					double gof;

					//internals::blockMatchOnePoint(referenceBlock, deformedBlock, blockRadius, refPoint, defPoint, gof);
					internals::blockMatchOnePointMultires(referenceBlock, deformedBlock, coarseBlockRadius, coarseBinning, fineBlockRadius, fineBinning, refPoint, defPoint, gof);

					defPoints(x, y, z) = defPoint + Vec3d(defStart);
					accuracy(x, y, z) = (float32_t)gof;
				}

				showThreadProgress(counter, defPoints.depth() * defPoints.height());
			}
		}
	}

	/*
	Finds bad displacement values and replaces them with nans.
	*/
	void filterDisplacements(const PointGrid3D<coord_t>& refGrid, Image<Vec3d>& defPoints, Image<float32_t>& accuracy, size_t filterRadius = 5, float32_t threshold = 3);

	/*
	Calculates translation between two 3D images by Maximum Intensity projecting them in two planes
	and by using phase correlation on the projections.
	*/
	template<typename pixel_t> Vec3d mipMatch(const Image<pixel_t>& ref, const Image<pixel_t>& def)
	{
		// XY
		Image<float32_t> refP, defP;
		max(ref, 0, refP);
		max(def, 0, defP);

		raw::writed(refP, "./registration/mipmatch_xy_ref");
		raw::writed(defP, "./registration/mipmatch_xy_def");

		double goodness;
		Vec3c maxXY;
		maxXY.x = refP.width() / 2;
		maxXY.y = refP.height() / 2;
		maxXY.z = 0;
		Vec3d shiftXY = phaseCorrelation(refP, defP, maxXY, goodness);
		
		raw::writed(refP, "./registration/mipmatch_xy_correlation");

		// YZ
		max(ref, 1, refP);
		max(def, 1, defP);

		raw::writed(refP, "./registration/mipmatch_yz_ref");
		raw::writed(defP, "./registration/mipmatch_yz_def");

		Vec3c maxYZ;
		maxYZ.x = refP.width() / 2;
		maxYZ.y = refP.height() / 2;
		maxYZ.z = 0;
		Vec3d shiftYZ = phaseCorrelation(refP, defP, maxYZ, goodness);

		raw::writed(refP, "./registration/mipmatch_yz_correlation");

		// Total shift
		Vec3d shift(shiftYZ.x, shiftXY.y, -(shiftXY.x - shiftYZ.y) / 2.0);

		return shift;
	}

	namespace internals
	{
		///*
		//Calculates coordinates of a point in deformed image corresponding to a point in reference coordinates.
		//*/
		//inline Vec3d projectPointToDeformed(const Vec3d& xRef, const vector<Vec3d>& refPoints, const vector<Vec3d>& defPoints)
		//{
		//	// Inverse distance interpolation for p == 4
		//	//constexpr double p = 4;
		//	Vec3d sum;
		//	double wsum = 0;
		//	for (size_t n = 0; n < refPoints.size(); n++)
		//	{
		//		double dist2 = (xRef - refPoints[n]).normSquared();

		//		if (dist2 < 1e-14)
		//		{
		//			return defPoints[n];
		//		}

		//		double w = (1 / (dist2 * dist2));

		//		sum += w * defPoints[n];
		//		wsum += w;
		//	}

		//	return sum / wsum;

		//	//// Inverse distance interpolation
		//	//constexpr double p = 3.25;
		//	//Vec3d sum;
		//	//double wsum = 0;
		//	//for (size_t n = 0; n < refPoints.size(); n++)
		//	//{
		//	//	double dist = (xRef - refPoints[n]).norm();

		//	//	if (dist < 1e-7)
		//	//	{
		//	//		return defPoints[n];
		//	//	}

		//	//	double w = (1 / pow(dist, p));

		//	//	sum += w * defPoints[n];
		//	//	wsum += w;
		//	//}

		//	//return sum / wsum;
		//}


		/*
		Calculates coordinates of a point in deformed image corresponding to a point in reference coordinates.
		*/
		template<typename real_t> Vec3<real_t> projectPointToDeformed(const Vec3<real_t>& xRef, const PointGrid3D<coord_t>& refPoints, const Image<Vec3<real_t> >& defPoints, const Interpolator<Vec3<real_t>, Vec3<real_t>, real_t>& interpolator)
		{
			real_t fx = (real_t)refPoints.xg.getIndex(xRef.x);
			real_t fy = (real_t)refPoints.yg.getIndex(xRef.y);
			real_t fz = (real_t)refPoints.zg.getIndex(xRef.z);

			//return linearInterpolationClamp<Vec3d, Vec3d>(defPoints, fx, fy, fz);
			//return LinearInterpolator<Vec3d, Vec3d, double, Vec3d>(Nearest)(defPoints, fx, fy, fz);
			return interpolator(defPoints, fx, fy, fz);
		}
	}


	/*
	Reverses deformation so that deformed image becomes similar to the original image.
	@param deformed Deformed image.
	@param pullback Result image. This image will contain deformed image reversed to the coordinates of the original, non-deformed, image. Size of this image must be set by the caller.
	@param refPoints Points in the reference image whose locations in the deformed image have been be determined.
	@param defPoints Locations of points in deformed image corresponding to reference points.
	*/
	template<typename def_t, typename result_t> void reverseDeformation(const Image<def_t>& deformed, Image<result_t>& pullback, const PointGrid3D<coord_t>& refGrid, const Image<Vec3d>& defPoints, const Interpolator<result_t, def_t, double>& interpolator = LinearInterpolator<result_t, def_t, double, double>(Nearest))
	{
		if (refGrid.pointCount() != defPoints.pixelCount())
			throw ITLException("refGrid and defPoints must have equal number of elements.");

		// TODO: determine region of pullback that must be processed

		Image<Vec3d> shifts(defPoints.dimensions());
		for (coord_t z = 0; z < defPoints.depth(); z += 1)
		{
			for (coord_t y = 0; y < defPoints.height(); y += 1)
			{
				for (coord_t x = 0; x < defPoints.width(); x += 1)
				{
					shifts(x, y, z) = defPoints(x, y, z) - Vec3d(refGrid(x, y, z));
				}
			}
		}

		LinearInterpolator<Vec3d, Vec3d, double, Vec3d> shiftInterpolator(Nearest);

		size_t counter = 0;
#pragma omp parallel for if(pullback.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < pullback.depth(); z += 1)
		{
			for (coord_t y = 0; y < pullback.height(); y += 1)
			{
				for (coord_t x = 0; x < pullback.width(); x += 1)
				{
					Vec3d xRef((double)x, (double)y, (double)z);
					Vec3d shift = internals::projectPointToDeformed(xRef, refGrid, shifts, shiftInterpolator);
					
					Vec3d xDef = xRef + shift;
					result_t val = interpolator(deformed, xDef);
					pullback(x, y, z) = val;
				}
			}

			showThreadProgress(counter, pullback.depth());
		}
	}

	///*
	//Reverses deformation so that deformed image becomes similar to the original image.
	//@param deformed Deformed image.
	//@param pullback Result image. This image will contain deformed image reversed to the coordinates of the original, non-deformed, image. Size of this image must be set by the caller.
	//@param refPoints Points in the reference image whose locations in the deformed image have been be determined.
	//@param defPoints Locations of points in deformed image corresponding to reference points.
	//*/
	//template<typename def_t, typename result_t> void reverseDeformation(const Image<def_t>& deformed, Image<result_t>& pullback, const vector<Vec3d>& refPoints, const vector<Vec3d>& defPoints)
	//{
	//	if (refPoints.size() != defPoints.size())
	//		throw ITLException("refPoints and defPoints must have same size.");

	//	// TODO: determine region of pullback that must be processed

	//	vector<Vec3d> shifts;
	//	shifts.reserve(refPoints.size());
	//	for (coord_t n = 0; n < refPoints.size(); n++)
	//		//shifts.push_back(refPoints[n] - defPoints[n]);
	//		shifts.push_back(defPoints[n] - refPoints[n]);

	//	//Image<float32_t> shiftX(pullback.dimensions());

	//	size_t counter = 0;
	//	#pragma omp parallel for if(pullback.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
	//	for (coord_t z = 0; z < pullback.depth(); z += 1)
	//	//for (coord_t z = 256; z < 257; z += 1)
	//	{
	//		for (coord_t y = 0; y < pullback.height(); y += 1)
	//		{
	//			for (coord_t x = 0; x < pullback.width(); x += 1)
	//			{
	//				Vec3d xRef(x, y, z);
	//				Vec3d shift = internals::projectPointToDeformed(xRef, refPoints, shifts);
	//				Vec3d xDef = xRef + shift;
	//				result_t val = pixelRound<result_t>(linearInterpolation(deformed, xDef));
	//				pullback(x, y, z) = val;

	//				//shiftX(x, y, z) = shift.x;
	//			}
	//		}

	//		showThreadProgress(counter, pullback.depth());
	//	}

	//	//raw::writed(shiftX, "./registration/shiftX");
	//}

	/*
	Writes a result of blockmatch operation to disk.
	*/
	void writeBlockMatchResult(const string& filenamePrefix, const PointGrid3D<coord_t>& refPoints, const Image<Vec3d>& defPoints, const Image<float32_t>& gof, double normFact = 1);

	/*
	Reads a result of blockmatch operation from disk.
	*/
	void readBlockMatchResult(const string& filenamePrefix, PointGrid3D<coord_t>& refPoints, Image<Vec3d>& defPoints, Image<float32_t>& gof, double& normFact);

	/*
	Used to read PointGrid1D from file.
	*/
	PointGrid1D<coord_t> readPointGrid1D(ifstream& in);

	namespace tests
	{
		void blockMatch1();
		void blockMatch2Match();
		void blockMatch2Pullback();
		void blockMatch3MatchNormal();
		void blockMatch3PullbackNormal();
		void blockMatch3MatchPartialLoad();
		void blockMatch3PullbackPartialLoad();
		void blockMatch4Match();
		void blockMatch4Pullback();
		void blockMatch5Match();
		void blockMatch5Pullback();
		void blockMatch6MatchPartialLoad();
		void blockMatch7MatchPartialLoad();
		void mipMatch();
	}
}