#pragma once

#include "image.h"
#include "particleanalysis.h"
#include "getslice.h"
#include "generation.h"
#include "transform.h"

namespace itl2
{
	namespace internals
	{
		inline bool isOnSliceEdge(std::vector<Vec3sc> points, coord_t w, coord_t h)
		{
			for (size_t n = 0; n < points.size(); n++)
			{
				Vec3sc& p = points[n];
				if (p.x <= 0 || p.y <= 0 || p.x >= w - 1 || p.y >= h - 1)
					return true;
			}

			return false;
		}
	}

	/**
	Analyze cross-sections of fibres.
	@param original Original binary image.
	@param oriEnergy Image corresponding to the 'energy' output of cylinderorientation command.
	@param oriPhi The azimuthal angle of local fibre orientation. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.
	@param oriTheta The polar angle of local fibre orientation. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.
	@param pLength Pointer to image containing local length for each fibre pixel. If null, length won't be analyzed. If not null, length values are read from this image and their per-slice statistics are added to the end of results array in columns 'length min', 'length max', 'lengh mean', and 'length mode'.
	@param analyzers Particle analyzers to apply to each fibre cross-section. In addition to the analyzer results, a few other quantities will be added to the results array: (X3, Y3, Z3) corresponds to a location in the original image that is inside the cross-section; (CX3, CY3, CZ3) corresponds to the centroid of the cross-section transformed into the coordinates of the original image.
	@param results Results table for the particle analysis.
	@param sliceRadius Radius of each slice. Slice width and height will be 2 x sliceRadius + 1
	@param sliceCount Count of slices to extract and analyze.
	@param randseed Seed for random number generator.
	@param pSlices If not 0, the generated cross-sectional slices (extracted from the original image) will be stored into this image.
	@param pLengthSlices If not 0, the generated cross-sectional length slices (extracted from the length image) will be stored into this image.
	@param pVisualization If not 0, the locations of the extracted slices will be drawn into this image.
	@param fillColor Temporary color. This must be a non-zero color other that what is used to mark the fibres in the original image.
	*/
	template<typename pixel_t> void csa(const Image<pixel_t>& original,
		const Image<float32_t>& oriEnergy, const Image<float32_t>& oriPhi, const Image<float32_t>& oriTheta,
		const Image<float32_t>* pLength,
		const AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results,
		coord_t sliceRadius,
		coord_t sliceCount,
		size_t randseed = 123,
		Image<pixel_t>* pSlices = 0,
		Image<float32_t>* pLengthSlices = 0,
		Image<pixel_t>* pVisualization = 0,
		pixel_t fillColor = internals::SpecialColors<pixel_t>::fillColor())
	{
		if (sliceCount < 1)
			throw ITLException("Slice count is too small.");

		oriEnergy.checkSize(oriPhi);
		oriEnergy.checkSize(oriTheta);

		coord_t sliceSize = 2 * sliceRadius + 1;

		double oriScale = (double)oriEnergy.width() / (double)original.width();
		coord_t oriSliceSize = pixelRound<coord_t>(oriScale * sliceSize);

		double lengthScale = 1.0;
		coord_t lengthSliceSize = 1;
		if (pLength)
		{
			lengthScale = (double)pLength->width() / (double)original.width();
			lengthSliceSize = pixelRound<coord_t>(lengthScale * sliceSize);
		}

		Image<pixel_t> slice(sliceSize, sliceSize);
		Image<float32_t> lengthSlice(lengthSliceSize, lengthSliceSize);

		if (pSlices)
			pSlices->ensureSize(sliceSize, sliceSize, sliceCount);
		if(pLengthSlices)
			pLengthSlices->ensureSize(sliceSize, sliceSize, sliceCount);
		if (pVisualization)
			pVisualization->ensureSize(original);

		srand((unsigned int)randseed);

		std::vector<Vec3sc> filledPoints;
		filledPoints.reserve(100);
		
		results.headers() = analyzers.headers();

		// Coordinates of a point that is in the slice in 3D
		results.headers().push_back("X3 [pixel]");
		results.headers().push_back("Y3 [pixel]");
		results.headers().push_back("Z3 [pixel]");

		size_t cxi, cyi, l1i, l2i, alphai, bsi;
		try
		{
			cxi = results.getColumnIndex("CX [pixel]");
			cyi = results.getColumnIndex("CY [pixel]");
			l1i = results.getColumnIndex("l1 [pixel]");
			l2i = results.getColumnIndex("l2 [pixel]");
			alphai = results.getColumnIndex("alpha [rad]");
			bsi = results.getColumnIndex("bounding scale [1]");

			// Position of the slice centroid in 3D
			results.headers().push_back("CX3 [pixel]");
			results.headers().push_back("CY3 [pixel]");
			results.headers().push_back("CZ3 [pixel]");
		}
		catch (ITLException)
		{
			// Not all columns required for semi-axis visualization are found.
			// Disable visualization
			cxi = std::numeric_limits<size_t>::max();
		}

		if (pLength)
		{
			results.headers().push_back("length min");
			results.headers().push_back("length max");
			results.headers().push_back("length mean");
			results.headers().push_back("length mode");
		}

		std::vector<float32_t> lSamples;
		lSamples.reserve(100);

		coord_t totalTrials = 0;
		size_t inBackgroundCount = 0;
		size_t noOrientationCount = 0;
		size_t tooSmallCount = 0;
		size_t touchesOriginalEdgeCount = 0;
		size_t touchesSliceEdgeCount = 0;

		NearestNeighbourInterpolator<pixel_t, pixel_t, double> interp(BoundaryCondition::Zero);

		{
			
			ProgressIndicator prog(sliceCount);
			for (coord_t n = 0; n < sliceCount; )
			{
				totalTrials++;
				Vec3c pos(randc(original.width()), randc(original.height()), randc(original.depth()));
				
				if (original.isInImage(pos) && original(pos) != 0)
				{
					Vec3c oriPos = itl2::round(Vec3d(pos) * oriScale);
					if (oriEnergy.isInImage(oriPos) && oriEnergy(oriPos) != 0)
					{
						double polar = oriTheta(oriPos);
						double azimuthal = oriPhi(oriPos);
						Vec3d dir = toCartesian(1.0, azimuthal, polar);

						bool touchesEdge = false;
						getSlice<pixel_t>(original, Vec3d(pos), dir, slice, &touchesEdge, interp);

						if (!touchesEdge)
						{

							filledPoints.clear();
							itl2::floodfillSingleThreaded(slice, Vec3c(sliceRadius, sliceRadius, 0), fillColor, fillColor, Connectivity::AllNeighbours, nullptr, &filledPoints);

							if (filledPoints.size() > 3)
							{
								if (!internals::isOnSliceEdge(filledPoints, sliceSize, sliceSize))
								{
									std::vector<double> resultLine;
									analyzers.analyze(filledPoints, resultLine);

									resultLine.push_back((double)pos.x);
									resultLine.push_back((double)pos.y);
									resultLine.push_back((double)pos.z);

									if (cxi != std::numeric_limits<size_t>::max())
									{
										// Visualize semi-axes of the particle
										double cx = resultLine[cxi];
										double cy = resultLine[cyi];
										double l1 = resultLine[l1i];
										double l2 = resultLine[l2i];
										double alpha = resultLine[alphai];
										double b = resultLine[bsi];
										l1 = b * l1;
										l2 = b * l2;
										Vec3d sdir(cos(alpha), sin(alpha), 0);
										Vec3d c(cx, cy, 0);
										Vec3d perp = sdir.rotate(Vec3d(0, 0, 1), PI / 2);
										Vec3d end1 = c + sdir * l1;
										Vec3d end2 = c + perp * l2;
										draw<pixel_t, double>(slice, Line<double>(c, end1), std::numeric_limits<pixel_t>::max() - 1);
										draw<pixel_t, double>(slice, Line<double>(c, end2), std::numeric_limits<pixel_t>::max() - 1);


										Matrix3x3d rot = internals::sliceRotationMatrix(dir);
										Vec3d p(cx - (double)sliceRadius, cy - (double)sliceRadius, 0);
										Vec3d cx3 = rot * p + Vec3d(pos);

										resultLine.push_back(cx3.x);
										resultLine.push_back(cx3.y);
										resultLine.push_back(cx3.z);
									}

									if (pLength)
									{
										Vec3d lengthPos = Vec3d(pos) * lengthScale;
										clamp(lengthPos, Vec3d(0, 0, 0), Vec3d(pLength->dimensions()) - Vec3d(1, 1, 1));
										getSlice(*pLength, lengthPos, dir, lengthSlice);

										lSamples.clear();
										for (size_t nn = 0; nn < filledPoints.size(); nn++)
										{
											Vec3c p = round(Vec3d(filledPoints[nn]) * lengthScale);
											if (lengthSlice.isInImage(p))
												lSamples.push_back((float32_t)(lengthSlice(p) / lengthScale));
										}

										float32_t lmin = min(lSamples);
										float32_t lmax = max(lSamples);
										float32_t lmean = mean(lSamples);
										float32_t lmode = mode(lSamples);

										resultLine.push_back(lmin);
										resultLine.push_back(lmax);
										resultLine.push_back(lmean);
										resultLine.push_back(lmode);

										if (pLengthSlices)
											copyValues(*pLengthSlices, lengthSlice, Vec3c(0, 0, n));
									}

									results.push_back(resultLine);

									if (pSlices)
										copyValues(*pSlices, slice, Vec3c(0, 0, n));
									
									if(pVisualization)
										drawSlice<pixel_t>(*pVisualization, Vec3d(pos), dir, Vec2d((double)slice.width(), (double)slice.height()), std::numeric_limits<pixel_t>::max() / 2);

									n++;
									prog.step();
								}
								else
								{
									touchesSliceEdgeCount++;
								}
							}
							else
							{
								tooSmallCount++;
							}
						}
						else
						{
							touchesOriginalEdgeCount++;
						}
					}
					else
					{
						noOrientationCount++;
					}
				}
				else
				{
					inBackgroundCount++;
				}
			}
		}

		std::cout << "Total points tested: " << totalTrials << std::endl
				<< "of which " << std::endl
				<< "in background: " << inBackgroundCount << " (" << std::setprecision(2) << ((double)inBackgroundCount / totalTrials * 100) << " %)" << std::endl
				<< "too small: " << tooSmallCount << " (" << std::setprecision(2) << ((double)tooSmallCount  / totalTrials * 100) << " %)" << std::endl
				<< "touches edge of original: " << touchesOriginalEdgeCount << " (" << std::setprecision(2) << ((double)touchesOriginalEdgeCount / totalTrials * 100) << " %)" << std::endl
				<< "touches edge of slice: " << touchesSliceEdgeCount << " (" << std::setprecision(2) << ((double)touchesSliceEdgeCount / totalTrials * 100) << " %)" << std::endl
				<< "no orientation available: " << noOrientationCount << " (" << std::setprecision(2) << ((double)noOrientationCount / totalTrials * 100) << " %)" << std::endl;
	}


	namespace tests
	{
		void csa();
	}
}