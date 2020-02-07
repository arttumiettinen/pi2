//#pragma once
//
//#include <vector>
//
//#include "image.h"
//#include "math/vec3.h"
//#include "sphere.h"
//
//namespace itl2
//{
//
//	namespace thresholddecomposition
//	{
//		/**
//		Calculates squared local radius from squared distance map.
//		Uses threshold decomposition algorithm.
//		*/
//		void thickmap2(const Image<int32_t>& dmap2, Image<int32_t>& result, bool showProgressInfo = true);
//	}
//
//	namespace dimredblocks
//	{
//		/**
//		Calculate squared local radius from squared distance map.
//		@param dmap2 Squared distance map at input, squared local thickness map at output.
//		@param extraSpheres List of spheres that may contribute to the image but whose location is outside of the image. The vector is emptied during processing.
//		@param counts If nonzero, outputs total count of ri values after processing of each dimension.
//		*/
//		void thickmap2(Image<int32_t>& dmap2, std::vector<Sphere2>& extraSpheres, Vec3d* counts, bool showProgressInfo = true);
//
//		/**
//		Calculate squared local radius from squared distance map, but divide the input image into blocks of size blocSize and process each block separately.
//		@param dmap2 Squared distance map at input, squared local thickness map at output.
//		@param blockSize Size of one calculation block.
//		@param extraBytes Approximation of memory used in addition to the input image.
//		*/
//		void thickmap2Blocks(Image<int32_t>& dmap2, const Vec3c& blockSize, double* extraBytes, bool showProgressInfo = true);
//	}
//
//	namespace dimredsmartblocks
//	{
//		/*
//		Calculate squared local radius from squared distance map.
//		@param dmap2 Squared distance map at input, squared local thickness map at output.
//		@param counts If nonzero, outputs total count of ri values after processing of each dimension.
//		*/
//		void thickmap2(Image<int32_t>& dmap2, Vec3d* counts, bool showProgressInfo = true);
//
//		/*
//		Calculate squared local radius from squared distance map.
//		@param dmap2 Squared distance map at input, squared local thickness map at output.
//		@param extraBytes Approximation of memory used in addition to the input image.
//		*/
//		//void thickmap2(Image<int32_t>& dmap2, double* extraBytes = nullptr, bool showProgressInfo = true);
//
//		/**
//		Calculate squared local radius from squared distance map, but divide the input image into blocks of size blocSize and process each block separately.
//		@param dmap2 Squared distance map at input, squared local thickness map at output.
//		@param blockSize Size of one calculation block.
//		@param extraBytes Approximation of memory used in addition to the input image.
//		*/
//		void thickmap2Blocks(Image<int32_t>& dmap2, const Vec3c& blockSize, double* extraBytes, bool showProgressInfo = true);
//	}
//



	//namespace dimred
	//{
	//	/*
	//	Calculate squared local radius from squared distance map.
	//	@param dmap2 Squared distance map at input, squared local thickness map at output.
	//	@param counts If nonzero, outputs total count of ri values after processing of each dimension.
	//	*/
	//	void thickmap2(Image<int32_t>& dmap2, Vec3d* counts, bool showProgressInfo = true);

	//	/*
	//	Calculate squared local radius from squared distance map.
	//	@param dmap2 Squared distance map at input, squared local thickness map at output.
	//	@param extraBytes Approximation of memory used in addition to the input image.
	//	*/
	//	void thickmap2(Image<int32_t>& dmap2, double* extraBytes = nullptr, Vec3d* counts = nullptr, bool showProgressInfo = true);
	//}

//}