#pragma once

// TODO: This is direct port from itl2 and needs further cleaning to comply with style of itl2.

#include "image.h"
#include "pointprocess.h"
#include "timer.h"
#include "io/raw.h"
#include "misc.h"
#include "math/vectoroperations.h"

#include <vector>
#include <string>

namespace itl2
{
	/**
	Length measures for path opening.
	*/
	enum class LengthType
	{
		/**
		Distance to all neighbours is 1.
		*/
		Ones,
		/**
		Use chamfer mask. The final result must be divided by 3.
		*/
		Chamfer,
		/**
		Use real distance. Should be used only if length image type is a floating point type.
		*/
		RealDistance
	};

	namespace internals
	{
		/**
		Builds list of strides to neighbours in the main direction.
		3D version.
		*/
		inline void buildNeighbours(std::vector<Vec3c>& nbs, const Vec3c& mainDir, size_t dimensionality)
		{
			nbs.clear();

			coord_t mx = mainDir.x;
			coord_t my = mainDir.y;
			coord_t mz = mainDir.z;
			for (coord_t x = -1; x <= 1; x++)
			{
				for (coord_t y = -1; y <= 1; y++)
				{
					for (coord_t z = -1; z <= 1; z++)
					{
						// Skip zero vector
						if (x == 0 && y == 0 && z == 0)
							continue;

						// Skip neighbours with z != 0 if 2D version.
						if (dimensionality == 2 && z != 0)
							continue;

						if (abs(x - mainDir.x) <= 1 && abs(y - mainDir.y) <= 1 && abs(z - mainDir.z) <= 1)
						{
							if ((x == mx && x != 0) ||
								(y == my && y != 0) ||
								(z == mz && z != 0))
							{
								nbs.push_back(Vec3c(x, y, z));
							}
						}
					}
				}
			}

		}

		template<typename LT> LT calcDistance(const Vec3c& nb, LengthType type)
		{
			if (type == LengthType::Ones)
			{
				return 1;
			}
			else if (type == LengthType::Chamfer)
			{
				double nsq = nb.normSquared<double>();
				if (NumberUtils<double>::equals(nsq, 1.0))
					return 3;
				if (NumberUtils<double>::equals(nsq, 2.0))
					return 4;
				if (NumberUtils<double>::equals(nsq, 3.0))
					return 5;
			}
			else if (type == LengthType::RealDistance)
			{
				return (LT)nb.norm();
			}

			throw ITLException("Invalid length type.");
		}

		/**
		Calculate stride to neighbour in 3D case.
		*/
		inline coord_t calcStride(const Vec3c& nb, const Vec3c& s)
		{
			return nb.x * s.x + nb.y * s.y + nb.z * s.z;
		}

		/**
		Calculate strides to neighbours and neighbour distances in 3D case.
		*/
		template<typename LENGTH_TYPE> void calcStrides(std::vector<coord_t>& nplus, std::vector<coord_t>& nminus, std::vector<LENGTH_TYPE>& nds, const std::vector<Vec3c>& nbs, const Vec3c& s, LengthType lengthType)
		{
			nplus.clear();
			nminus.clear();
			nds.clear();

			for (size_t n = 0; n < nbs.size(); n++)
			{
				coord_t x = nbs[n].x;
				coord_t y = nbs[n].y;
				coord_t z = nbs[n].z;

				nds.push_back(calcDistance<LENGTH_TYPE>(Vec3c(x, y, z), lengthType));
				nplus.push_back(calcStride(Vec3c(x, y, z), s));
			}

			nminus = -nplus;
		}

		


		///**
		//Performs one path length scan over the image.
		//LT: Data type of length images.
		//PT: Data type of geometry image.
		//@param img Geometry image. Pixels == 1 are the ROI and pixels == 0 are background. No other pixel values should be present in the image.
		//@param nf, nb Lists of strides to forward and backward neighbours.
		//@param nfc, nbc Stride to the neighbour that is in the main orientation, in forward and backward direction for nplusc and nminusc, respectively.
		//@param nds List of distances to neighbours whose positions are given by the strides. The list must be the same for both the forward and backward neighbours.
		//@param ndc Distance to neighbour corresponding to the main orientation.
		//@param lambda Resulting lambda image.
		//@param lambdac Resulting lambdac image.
		//@param mainDir The main traversal direction.
		//*/
		//template<typename LT, typename PT> void longestPathOneScanC2(BasicImage<PT>& img, const vector<coord_t>& nf, coord_t nfc, const vector<coord_t>& nb, coord_t nbc, const vector<LT>& nds, LT ndc, BasicImage<LT>& lambda, BasicImage<LT>& lambdac, const Vec3c& mainDir)
		//{
		//	Timer timer;
		//	timer.start();

		//	// Init lambdas to maximum possible value
		//	pointprocess<LT, LT, basicfilters::setValue>(lambda, numeric_limits<LT>::max());
		//	pointprocess<LT, LT, basicfilters::setValue>(lambdac, numeric_limits<LT>::max());

		//	// length is zero where img is zero (=outside structures)
		//	pointprocess<LT, PT, basicfilters::multiply<LT, PT> >(lambda, img);
		//	pointprocess<LT, PT, basicfilters::multiply<LT, PT> >(lambdac, img);



		//	coord_t w = (coord_t)img.getDimension(0);
		//	coord_t h = (coord_t)img.getDimension(1);
		//	coord_t d = (coord_t)img.getDimension(2);

		//	// Calculate minimum and maximum level
		//	vector<coord_t> l;
		//	l.push_back(mainDir.dott(Vec3c(0, 0, 0)));
		//	l.push_back(mainDir.dott(Vec3c(w - 1, 0, 0)));
		//	l.push_back(mainDir.dott(Vec3c(0, h - 1, 0)));
		//	l.push_back(mainDir.dott(Vec3c(w - 1, h - 1, 0)));
		//	l.push_back(mainDir.dott(Vec3c(0, 0, d - 1)));
		//	l.push_back(mainDir.dott(Vec3c(w - 1, 0, d - 1)));
		//	l.push_back(mainDir.dott(Vec3c(0, h - 1, d - 1)));
		//	l.push_back(mainDir.dott(Vec3c(w - 1, h - 1, d - 1)));
		//	coord_t lmin = mathutils::min(l);
		//	coord_t lmax = mathutils::max(l);


		//	size_t pixelsAtCurrentLevel = 0;

		//	Timer roundTimer;
		//	for (coord_t currLevel = lmin; currLevel <= lmax; currLevel++)
		//	{
		//		roundTimer.start();

		//		pixelsAtCurrentLevel = 0;

		//		for (coord_t z = 1; z < d - 1; z++)
		//		{
		//			coord_t orderz = mainDir.z * z;
		//			for (coord_t y = 1; y < h - 1; y++)
		//			{
		//				coord_t orderyz = mainDir.y * y + orderz;

		//				coord_t order = orderyz + mainDir.x; // As x starts from 1, add one mainDir.x to the order here.
		//				for (coord_t x = 1; x < w - 1; x++, order += mainDir.x) // Optimized version, this does not waste time for multiplication
		//				//for(coord_t x = 1; x < w-1; x++)
		//				{
		//					// Process only pixels at the current level.
		//					//coord_t order = mainDir.x * x + mainDir.y * y + mainDir.z * z;
		//					//coord_t order = mainDir.x * x + orderyz;
		//					if (order == currLevel)
		//					{
		//						size_t q = (size_t)z * (size_t)w * (size_t)h + (size_t)y * (size_t)w + (size_t)x;

		//						// Process only pixels that don't have length value yet.
		//						if (lambda.getPixel(q) >= numeric_limits<LT>::max())
		//						{

		//							// update constrained length from normal length
		//							LT l = lambda.getPixel(q + nbc);
		//							//if(l >= numeric_limits<LT>::max())
		//							//	throw "1. Impossible, l == maxval";

		//							l += ndc;
		//							lambdac.setPixel(q, l);


		//							// update normal length from constrained length

		//							// l = max of lambdas at nb neighbours
		//							// i.e. maximum path length in neighbours in 'backward' direction + 1
		//							l = 0;
		//							for (size_t i = 0; i < nb.size(); i++)
		//							{
		//								LT ll = lambdac.getPixel(q + nb[i]);
		//								//if(ll >= numeric_limits<LT>::max())
		//								//	throw "2. Impossible, l == maxval";
		//								ll += nds[i];
		//								//l = max(l, ll);   // Take max as we are searching for the longest path, not the shortest.
		//								if (ll > l)
		//									l = ll;
		//							}

		//							//if(l >= numeric_limits<LT>::max())
		//							//	throw "3. Impossible, l == maxval";

		//							lambda.setPixel(q, l);
		//						}


		//						pixelsAtCurrentLevel++;
		//					}
		//				}
		//			}
		//		}


		//		roundTimer.stop();

		//		cout << "Level " << currLevel << ", processing time " << roundTimer.getSeconds() << " s, processed " << pixelsAtCurrentLevel << " pixels       \r" << flush;
		//	}


		//	timer.stop();
		//	cout << endl << "Operation took " << timer.getSeconds() << " s" << endl;

		//}

		/**
		Performs one path length scan over the image.
		LT: Data type of length images.
		PT: Data type of geometry image.
		@param img Geometry image. Pixels == 1 are the ROI and pixels == 0 are background. No other pixel values should be present in the image.
		@param nf, nb Lists of strides to forward and backward neighbours.
		@param nfc, nbc Stride to the neighbour that is in the main orientation, in forward and backward direction for nplusc and nminusc, respectively.
		@param nds List of distances to neighbours whose positions are given by the strides. The list must be the same for both the forward and backward neighbours.
		@param ndc Distance to neighbour corresponding to the main orientation.
		@param lambda Resulting lambda image.
		@param lambdac Resulting lambdac image.
		@param mainDir The main traversal direction.
		*/
		template<typename LT, typename PT> void longestPathOneScanC2_parallel(Image<PT>& img, const std::vector<coord_t>& nf, coord_t nfc, const std::vector<coord_t>& nb, coord_t nbc, const std::vector<LT>& nds, LT ndc, Image<LT>& lambda, Image<LT>& lambdac, const Vec3c& mainDir)
		{
			Timer timer;
			timer.start();

			// Init lambdas to maximum possible value
			setValue(lambda, std::numeric_limits<LT>::max());
			setValue(lambdac, std::numeric_limits<LT>::max());

			// length is zero where img is zero (=outside structures)
			multiply(lambda, img);
			multiply(lambdac, img);
			

			coord_t w = (coord_t)img.dimension(0);
			coord_t h = (coord_t)img.dimension(1);
			coord_t d = (coord_t)img.dimension(2);

			// Calculate minimum and maximum level
			std::vector<coord_t> l;
			l.push_back(mainDir.dot<coord_t>(Vec3c(0, 0, 0)));
			l.push_back(mainDir.dot<coord_t>(Vec3c(w - 1, 0, 0)));
			l.push_back(mainDir.dot<coord_t>(Vec3c(0, h - 1, 0)));
			l.push_back(mainDir.dot<coord_t>(Vec3c(w - 1, h - 1, 0)));
			if (img.dimensionality() >= 3)
			{
				l.push_back(mainDir.dot<coord_t>(Vec3c(0, 0, d - 1)));
				l.push_back(mainDir.dot<coord_t>(Vec3c(w - 1, 0, d - 1)));
				l.push_back(mainDir.dot<coord_t>(Vec3c(0, h - 1, d - 1)));
				l.push_back(mainDir.dot<coord_t>(Vec3c(w - 1, h - 1, d - 1)));
			}
			coord_t lmin = min(l);
			coord_t lmax = max(l);


			coord_t minz = 0;
			coord_t maxz = 1;
			if (img.dimensionality() >= 3)
			{
				minz = 1;
				maxz = d - 1;
			}

			//size_t pixelsAtCurrentLevel = 0;

			// TODO: This can be done without traversing all image pixels for each level. See fastminmaxfilters.h, function lineOp(...).

			Timer roundTimer;
			for (coord_t currLevel = lmin; currLevel <= lmax; currLevel++)
			{
				roundTimer.start();

				//pixelsAtCurrentLevel = 0;

#pragma omp parallel for
				for (coord_t z = minz; z < maxz; z++)
				{
					coord_t orderz = mainDir.z * z;
					for (coord_t y = 1; y < h - 1; y++)
					{
						coord_t orderyz = mainDir.y * y + orderz;

						coord_t order = orderyz + mainDir.x; // As x starts from 1, add one mainDir.x to the order here.
						for (coord_t x = 1; x < w - 1; x++, order += mainDir.x) // Optimized version, this does not waste time for multiplication
						//for(coord_t x = 1; x < w-1; x++)
						{
							// Process only pixels at the current level.
							//coord_t order = mainDir.x * x + mainDir.y * y + mainDir.z * z;
							//coord_t order = mainDir.x * x + orderyz;
							if (order == currLevel)
							{
								size_t q = (size_t)z * (size_t)w * (size_t)h + (size_t)y * (size_t)w + (size_t)x;

								// Process only pixels that don't have length value yet.
								if (lambda(q) >= std::numeric_limits<LT>::max())
								{

									// update constrained length from normal length
									LT l = lambda(q + nbc);
									//if(l >= numeric_limits<LT>::max())
									//	throw "1. Impossible, l == maxval";

									l += ndc;
									lambdac(q) = l;


									// update normal length from constrained length

									// l = max of lambdas at nb neighbours
									// i.e. maximum path length in neighbours in 'backward' direction + 1
									l = 0;
									for (size_t i = 0; i < nb.size(); i++)
									{
										LT ll = lambdac(q + nb[i]);
										//if(ll >= numeric_limits<LT>::max())
										//	throw "2. Impossible, l == maxval";
										ll += nds[i];
										//l = max(l, ll);   // Take max as we are searching for the longest path, not the shortest.
										if (ll > l)
											l = ll;
									}

									//if(l >= numeric_limits<LT>::max())
									//	throw "3. Impossible, l == maxval";

									lambda(q) = l;
								}


								//pixelsAtCurrentLevel++;
							}
						}
					}
				}


				roundTimer.stop();

				//cout << "Level " << currLevel << ", processing time " << roundTimer.getSeconds() << " s, processed " << pixelsAtCurrentLevel << " pixels       \r" << flush;
				std::cout << "Level " << currLevel << ", processing time " << roundTimer.getSeconds() << " s       \r" << std::flush;
			}


			timer.stop();
			std::cout << std::endl << "Operation took " << timer.getSeconds() << " s" << std::endl;

		}
	}

	/**
	Longest path calculation for 3D binary image. Memory saver version. Uses "2nd generation" algorithm that does not need the queue in Cris' algorithm.
	@param img Original binary geometry image. The edges of this image will be set to zeroes.
	@param lengths Output image that will contain path length for each pixel.
	@param tempName Name stem for temporary files.
	*/
	template<typename PIXEL_TYPE> void pathLength2Binary3dNormalOrChamferMemorySave(Image<PIXEL_TYPE>& img, Image<float32_t>& lengths, LengthType lengthType = LengthType::Ones, const std::string& tempName = "TEMP")
	{
		if (img.dimensionality() != 3 && img.dimensionality() != 2)
			throw ITLException("This method supports only 2- and 3-dimensional images.");

		// Constrain paths to the image
		setEdges<PIXEL_TYPE>(img, 0);

		// Make sure that the original image contains only zero and one. (so that it can be used as temporary space, too)
		threshold(img, 0);

		// List of main directions
		std::vector<Vec3c> mainDirs;
		mainDirs.push_back(Vec3c(1, -1, 0));
		mainDirs.push_back(Vec3c(1, 0, 0));
		mainDirs.push_back(Vec3c(1, 1, 0));
		mainDirs.push_back(Vec3c(0, 1, 0));

		if (img.dimensionality() >= 3)
		{
			mainDirs.push_back(Vec3c(1, -1, -1));
			mainDirs.push_back(Vec3c(1, -1, 1));
			mainDirs.push_back(Vec3c(1, 0, -1));
			mainDirs.push_back(Vec3c(1, 0, 1));
			mainDirs.push_back(Vec3c(1, 1, -1));
			mainDirs.push_back(Vec3c(1, 1, 1));
			mainDirs.push_back(Vec3c(0, 1, -1));
			mainDirs.push_back(Vec3c(0, 1, 1));
			mainDirs.push_back(Vec3c(0, 0, 1));
		}

		// Stride vector
		coord_t s1 = 1;
		coord_t s2 = (coord_t)img.dimension(0);
		coord_t s3 = (coord_t)img.dimension(0) * (coord_t)img.dimension(1);
		Vec3c s(s1, s2, s3);

		std::cout << "Calculating lambdas..." << std::endl;
		std::vector<std::string> lambdaplusFilenames;
		std::vector<std::string> lambdapluscFilenames;
		std::vector<std::string> lambdaminusFilenames;
		std::vector<std::string> lambdaminuscFilenames;
		for (size_t n = 0; n < mainDirs.size(); n++)
		{
			std::cout << (n + 1) << " / " << mainDirs.size() << std::endl;
			Vec3c mainDir = mainDirs[n];

			// Build neighbour list
			std::vector<coord_t> nplus, nminus;
			std::vector<uint16_t> nds;
			nplus.reserve(8);
			nminus.reserve(8);
			nds.reserve(8);

			std::vector<Vec3c> nbs;
			internals::buildNeighbours(nbs, mainDir, img.dimensionality());
			internals::calcStrides(nplus, nminus, nds, nbs, s, lengthType);

			// Run the algorithm
			coord_t mainStride = internals::calcStride(mainDir, s);
			uint16_t mainDistance = internals::calcDistance<uint16_t>(mainDir, lengthType);


			Image<uint16_t> lambda(img.dimensions());	// lambdaplus(x) = Length of maximal path starting at x.
			Image<uint16_t> lambdac(img.dimensions());	// lambdaplusc(x) = Length of maximal path starting at x for which the previous step is in the main direction.

			// Skip existing files

			std::string nameminus = tempName + toString(n) + toString("-");
			//nameminus = raw::concatDimensions(nameminus, lambda.dimensions());
			std::string nameminusc = tempName + toString(n) + toString("-c");
			//nameminusc = raw::concatDimensions(nameminusc, lambdac.dimensions());

			Vec3c dimensions;
			ImageDataType dt;
			std::string reason;
			if (!raw::getInfo(nameminus, dimensions, dt, reason) || !raw::getInfo(nameminusc, dimensions, dt, reason))
			{
				std::cout << "Calculating " << nameminus << " and " << nameminusc << std::endl;

				internals::longestPathOneScanC2_parallel<uint16_t, PIXEL_TYPE>(img, nplus, mainStride, nminus, -mainStride, nds, mainDistance, lambda, lambdac, mainDir);
				raw::writed(lambda, nameminus);
				raw::writed(lambdac, nameminusc);
			}
			else
			{
				std::cout << "Skipping " << nameminus << " and " << nameminusc << std::endl;
			}

			lambdaminusFilenames.push_back(nameminus);
			lambdaminuscFilenames.push_back(nameminusc);



			string nameplus = tempName + toString(n) + toString("+");
			//nameplus = raw::concatDimensions(nameplus, lambda.dimensions());
			string nameplusc = tempName + toString(n) + toString("+c");
			//nameplusc = raw::concatDimensions(nameplusc, lambdac.dimensions());

			if (!raw::getInfo(nameplus, dimensions, dt, reason) || !raw::getInfo(nameplusc, dimensions, dt, reason))
			{
				std::cout << "Calculating " << nameplus << " and " << nameplusc << std::endl;

				internals::longestPathOneScanC2_parallel<uint16_t, PIXEL_TYPE>(img, nminus, -mainStride, nplus, mainStride, nds, mainDistance, lambda, lambdac, -mainDir);
				raw::writed(lambda, nameplus);
				raw::writed(lambdac, nameplusc);
			}
			else
			{
				std::cout << "Skipping " << nameplus << " and " << nameplusc << std::endl;
			}
			lambdaplusFilenames.push_back(nameplus);
			lambdapluscFilenames.push_back(nameplusc);
		}


		// Build final lengths
		std::cout << "Building lengths..." << std::endl;
		lengths.ensureSize(img);
		for (size_t k = 0; k < mainDirs.size(); k++)
		{
			std::cout << (k + 1) << " / " << mainDirs.size() << std::endl;

			Image<uint16_t> lambdaplus(lambdaplusFilenames[k], true, img.dimensions());
			Image<uint16_t> lambdaplusc(lambdapluscFilenames[k], true, img.dimensions());
			Image<uint16_t> lambdaminus(lambdaminusFilenames[k], true, img.dimensions());
			Image<uint16_t> lambdaminusc(lambdaminuscFilenames[k], true, img.dimensions());

			// Calculate lengths
			{
				ProgressIndicator prog(img.pixelCount());
				#pragma omp parallel for
				for (coord_t p = 0; p < img.pixelCount(); p++)
				{
					// If the data type works like integer (max value + 1 = min value), this won't work.
					//LT l1 = lambdaplus.getPixel(p) + lambdaminusc.getPixel(p) - 1; // Length of path through the pixel = length of path starting at the pixel + length of path ending at the pixel - 1
					//LT l2 = lambdaplusc.getPixel(p) + lambdaminus.getPixel(p) - 1;
					//LT l = max(l1, l2);
					//if(lengths.getPixel(p) < l)
					//	lengths.setPixel(p, l);  // If we have length, then it's the final value for this pixel, but only set it if it is larger than the old length at that pixel.

					uint16_t lp = lambdaplus(p);
					uint16_t lpc = lambdaplusc(p);
					uint16_t lm = lambdaminus(p);
					uint16_t lmc = lambdaminusc(p);

					if (lp < std::numeric_limits<uint16_t>::max() &&
						lpc < std::numeric_limits<uint16_t>::max() &&
						lm < std::numeric_limits<uint16_t>::max() &&
						lmc < std::numeric_limits<uint16_t>::max())
					{
						uint16_t l = std::max(lp + lmc, lm + lpc);
						if (l > 0)
							l = l - 1;
						if (lengths(p) < l)
							lengths(p) = l;
					}

					prog.step();
				}
			}
		}

		// Delete temp files
		// TODO: This does not work in Windows as the memory mappings above prevent the deletion for some reason.
		for (const auto& f : lambdaminusFilenames)
			deleteFile(f);
		for (const auto& f : lambdaminuscFilenames)
			deleteFile(f);
		for (const auto& f : lambdaplusFilenames)
			deleteFile(f);
		for (const auto& f : lambdapluscFilenames)
			deleteFile(f);

		if (lengthType == LengthType::Chamfer)
			divide(lengths, (float32_t)3.0);

	}


	namespace tests
	{
		void pathopening2d();
		void pathopening();
	}
}