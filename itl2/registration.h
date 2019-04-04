#pragma once

#include <vector>

#include "image.h"
#include "math/vec3.h"
#include "neighbourhood.h"
#include "fft.h"
#include "transform.h"
#include "inpaint.h"
#include "projections.h"
#include "io/io.h"

using math::Vec3;
using math::Vec3d;
using std::vector;

namespace itl2
{

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

				getNeighbourhood(reference, refPoint, r, refBlockOrig, BoundaryCondition::Zero);
				getNeighbourhood(deformed, defPointRounded, r, defBlockOrig, BoundaryCondition::Zero);

				maskedBinning(refBlockOrig, refBlock, binningSize, (float32_t)0, (float32_t)0, false);
				maskedBinning(defBlockOrig, defBlock, binningSize, (float32_t)0, (float32_t)0, false);
			}
			else
			{
				refBlock.ensureSize(blockSize);
				defBlock.ensureSize(blockSize);
				getNeighbourhood(reference, refPoint, r, refBlock, BoundaryCondition::Zero);
				getNeighbourhood(deformed, defPointRounded, r, defBlock, BoundaryCondition::Zero);
			}

			// Set zeros to nearest non-zero value. This has effect particularly in the edges and corners of non-rectangular images.
			inpaintNearest(refBlock);
			inpaintNearest(defBlock);

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

		size_t counter = 0;
		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t n = 0; n < (coord_t)refPoints.size(); n++)
		{
			internals::blockMatchOnePoint(reference, deformed, blockRadius, refPoints[n], defPoints[n], accuracy[n]);

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
			return first + index * step;
		}

		/*
		Gets point position at given index.
		No bounds checking is made.
		*/
		double getPosition(double index) const
		{
			return first + index * step;
		}

		/*
		Get index of point at given location.
		*/
		double getIndex(double point) const
		{
			return (point - first) / step;
		}

		/*
		Tests if the given position is inside the grid (between first and last grid points)
		*/
		template<typename real_t> bool contains(real_t val) const
		{
			return first <= (double)val && (double)val <= first + (pointCount() - 1) * step;
		}

		/*
		Gets count of points in this grid.
		*/
		coord_t pointCount() const
		{
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
	template<typename ref_t, typename def_t> void blockMatchPartialLoad(const string& referenceFile, const string& deformedFile, const PointGrid3D<coord_t>& refGrid, Image<Vec3d>& defPoints, Image<float32_t>& accuracy,
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

		ImageDataType refDt, defDt;
		Vec3c refImageDimensions, defImageDimensions;
		io::getInfo(referenceFile, refImageDimensions, refDt);
		io::getInfo(deformedFile, defImageDimensions, defDt);

		clamp(defStart, Vec3c(0, 0, 0), defImageDimensions);
		clamp(defEnd, Vec3c(0, 0, 0), defImageDimensions);
		clamp(refStart, Vec3c(0, 0, 0), refImageDimensions);
		clamp(refEnd, Vec3c(0, 0, 0), refImageDimensions);

		Image<ref_t> referenceBlock(refEnd - refStart);
		Image<def_t> deformedBlock(defEnd - defStart);

		cout << "Loading block of reference image, size = " << (referenceBlock.pixelCount() * sizeof(ref_t) / (1024 * 1024)) << " MiB." << endl;
		io::readBlock(referenceBlock, referenceFile, refStart, true);

		cout << "Loading block of deformed image, size  = " << (deformedBlock.pixelCount() * sizeof(ref_t) / (1024 * 1024)) << " MiB." << endl;
		io::readBlock(deformedBlock, deformedFile, defStart, true);

		// Calculate normalization factor for gray values
		normFact = 1;
		double meanRef = maskedmean(referenceBlock, (ref_t)0);
		double meanDef = maskedmean(deformedBlock, (def_t)0);
		normFact = meanRef - meanDef;
		
		if (math::isnan(normFact))
			normFact = 1;

		if (normalize)
		{
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


	namespace internals
	{
		/*
		Calculates coordinates of a point in deformed image corresponding to a point in reference coordinates.
		*/
		template<typename real_t> Vec3<real_t> projectPointToDeformed(const Vec3<real_t>& xRef, const PointGrid3D<coord_t>& refPoints, const Image<Vec3<real_t> >& defPoints, const Interpolator<Vec3<real_t>, Vec3<real_t>, real_t>& interpolator)
		{
			real_t fx = (real_t)refPoints.xg.getIndex(xRef.x);
			real_t fy = (real_t)refPoints.yg.getIndex(xRef.y);
			real_t fz = (real_t)refPoints.zg.getIndex(xRef.z);

			return interpolator(defPoints, fx, fy, fz);
		}
	}


	/*
	Reverses deformation so that a deformed image becomes similar to the original image.
	The reference image points are expressed as a rectangular grid, allowing fast interpolation operations
	in the reverse deformation direction.

	[
	NOTE:
	Generally, for fast interpolation operations the rectangular grid must be in the coordinates of the 'target', i.e.
	- if the grid is in the coordinates of the reference image, then it is fast to transform deformed image to reference image.
	- if the grid is in the coordinates of the deformed image, then it is fast to transform reference image to deformed image.
	]

	@param deformed Deformed image.
	@param pullback Result image. This image will contain deformed image reversed to the coordinates of the original, non-deformed, image. Size of this image must be set by the caller.
	@param refPoints Points in the reference image whose locations in the deformed image have been be determined.
	@param defPoints Locations of points in the deformed image corresponding to the reference points.
	*/
	template<typename def_t, typename result_t> void reverseDeformation(const Image<def_t>& deformed, Image<result_t>& pullback, const PointGrid3D<coord_t>& refGrid, const Image<Vec3d>& defPoints, const Interpolator<result_t, def_t, double>& interpolator = LinearInterpolator<result_t, def_t, double, double>(BoundaryCondition::Nearest))
	{
		deformed.mustNotBe(pullback);

		if (refGrid.pointCounts() != defPoints.dimensions())
			throw ITLException("Arguments refGrid and defPoints must have the same dimensions.");

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

		LinearInterpolator<Vec3d, Vec3d, double, Vec3d> shiftInterpolator(BoundaryCondition::Nearest);

		size_t counter = 0;
		#pragma omp parallel for if(pullback.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < pullback.depth(); z++)
		{
			for (coord_t y = 0; y < pullback.height(); y++)
			{
				for (coord_t x = 0; x < pullback.width(); x++)
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
		void mipMatch();
	}
}