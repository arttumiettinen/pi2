#pragma once

#include "image.h"
#include "particleanalysis.h"
#include "getslice.h"
#include "transform.h"
#include "generation.h"
#include "structure.h"

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
	@param oriOrig Binary image that has been used to calculate orientation (phi and theta).
	@param phi The azimuthal angle of local fibre orientation. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.
	@param theta The polar angle of local fibre orientation. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.
	@param pLength Pointer to image containing local length for each fibre pixel. If null, length won't be analyzed. If not null, length values are read from this image and their per-slice statistics are added to the end of results array in columns 'length min', 'length max', 'lengh mean', and 'length mode'.
	@param analyzers Particle analyzers to apply to eahch fibre cross-section.
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
		const Image<pixel_t>& oriOrig, const Image<float32_t>& phi, const Image<float32_t>& theta,
		const Image<float32_t>* pLength,
		const AnalyzerSet<Vec3sc, pixel_t>& analyzers, Results& results,
		coord_t sliceRadius,
		coord_t sliceCount,
		int randseed = 123,
		Image<pixel_t>* pSlices = 0,
		Image<pixel_t>* pLengthSlices = 0,
		Image<pixel_t>* pVisualization = 0,
		pixel_t fillColor = internals::SpecialColors<pixel_t>::fillColor())
	{
		if (sliceCount < 1)
			throw ITLException("Slice count is too small.");

		oriOrig.checkSize(phi);
		oriOrig.checkSize(theta);

		coord_t sliceSize = 2 * sliceRadius + 1;

		double oriScale = (double)oriOrig.width() / (double)original.width();
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

		if (pLength)
		{
			results.headers().push_back("length min");
			results.headers().push_back("length max");
			results.headers().push_back("length mean");
			results.headers().push_back("length mode");
		}

		size_t cxi, cyi, l1i, l2i, alphai, bsi;
		try
		{
			cxi = results.getColumnIndex("CX [pixel]");
			cyi = results.getColumnIndex("CY [pixel]");
			l1i = results.getColumnIndex("l1 [pixel]");
			l2i = results.getColumnIndex("l2 [pixel]");
			alphai = results.getColumnIndex("alpha [rad]");
			bsi = results.getColumnIndex("bounding scale");
		}
		catch (ITLException)
		{
			// Not all columns required for semi-axis visualization are found.
			// Disable visualization
			cxi = std::numeric_limits<size_t>::max();
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
					if (oriOrig.isInImage(oriPos) && oriOrig(oriPos) != 0)
					{
						double polar = theta(oriPos);
						double azimuthal = phi(oriPos);
						Vec3d dir = toCartesian(1.0, azimuthal, polar);

						bool touchesEdge = false;
						getSlice<pixel_t>(original, Vec3d(pos), dir, slice, &touchesEdge, interp);

						if (!touchesEdge)
						{

							filledPoints.clear();
							itl2::floodfill(slice, Vec3c(sliceRadius, sliceRadius, 0), fillColor, fillColor, Connectivity::AllNeighbours, nullptr, &filledPoints);

							if (filledPoints.size() > 3)
							{
								if (!internals::isOnSliceEdge(filledPoints, sliceSize, sliceSize))
								{
									std::vector<double> resultLine;
									analyzers.analyze(filledPoints, resultLine);

									if (pLength)
									{
										Vec3d lengthPos = Vec3d(pos) * lengthScale;
										clamp(lengthPos, Vec3d(0, 0, 0), Vec3d(pLength->dimensions()) - Vec3d(1, 1, 1));
										getSlice(*pLength, lengthPos, dir, lengthSlice);

										lSamples.clear();
										for (size_t n = 0; n < filledPoints.size(); n++)
										{
											Vec3c p = round(Vec3d(filledPoints[n]) * lengthScale);
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
										Vec3d dir(cos(alpha), sin(alpha), 0);
										Vec3d c(cx, cy, 0);
										Vec3d perp = dir.rotate(Vec3d(0, 0, 1), PI / 2);
										Vec3d end1 = c + dir * l1;
										Vec3d end2 = c + perp * l2;
										draw<pixel_t, double>(slice, Line<double>(c, end1), std::numeric_limits<pixel_t>::max() - 1);
										draw<pixel_t, double>(slice, Line<double>(c, end2), std::numeric_limits<pixel_t>::max() - 1);
									}

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

	/**
	Creates list of particle analyzers suitable for csa function.
	*/
	template<typename pixel_t> AnalyzerSet<Vec3sc, pixel_t> allCrossSectionAnalyzers()
	{
		AnalyzerSet<Vec3sc, pixel_t> analyzers;
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::Coordinates2D<Vec3sc, pixel_t>()));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::Volume<Vec3sc, pixel_t>()));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::PCA2D<Vec3sc, pixel_t>()));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::ConvexHull2D<Vec3sc, pixel_t>()));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::BoundingBox2D<Vec3sc, pixel_t>()));

		return analyzers;
	}

	namespace tests
	{
		void csa()
		{
			// Create original
			Image<uint8_t> orig(100, 100, 100);
			draw(orig, Capsule<double>(Vec3d(10, 10, 10), Vec3d(70, 70, 70), 5), (uint8_t)255);
			draw(orig, Capsule<double>(Vec3d(20, 80, 15), Vec3d(50, 50, 70), 8), (uint8_t)255);

			// Calculate orientation
			Image<float32_t> energy, phi, theta;
			//convert(orig, energy);
			//cylinderOrientation(energy, phi, theta, 3, 3);
			raw::read(energy, "./csa/energy");
			raw::read(phi, "./csa/phi");
			raw::read(theta, "./csa/theta");

			// Calculate length
			// TODO

			// Calculate CSA
			Results results;
			Image<uint8_t> slices;
			Image<uint8_t> vis;
			convert(orig, vis);
			itl2::csa<uint8_t>(orig, orig, phi, theta, nullptr, allCrossSectionAnalyzers<uint8_t>(), results, 20, 300, 123, &slices, nullptr, &vis);

			// Save results
			raw::writed(orig, "./csa/orig");
			raw::writed(energy, "./csa/energy");
			raw::writed(phi, "./csa/phi");
			raw::writed(theta, "./csa/theta");
			raw::writed(slices, "./csa/slices");
			raw::writed(vis, "./csa/vis");
		}
	}
}