#pragma once

#include <tuple>
#include <fstream>
#include <memory>
#include <map>

#include "math/mathutils.h"
#include "math/matrix3x3.h"
#include "pointprocess.h"
#include "transform.h"
#include "registration.h"
#include "conversions.h"
#include "io/io.h"

namespace itl2
{
	/*
	Encapsulates transformation between world coordinates and image.
	*/
	template<typename real_t> struct Transformation
	{
		/*
		Rotation matrix
		*/
		Matrix3x3<real_t> R;

		/*
		Inverse of rotation matrix
		*/
		Matrix3x3<real_t> Rinv;

		/*
		Scaling factor
		*/
		real_t a;

		/*
		Translation
		*/
		Vec3<real_t> c;

		/*
		Normalization factor for gray values.
		*/
		real_t normFactor;

		/*
		The same quantities for all parent images + reference grid and shifts.
		*/
		std::vector<real_t> parenta;
		std::vector<Vec3<real_t> > parentc;
		std::vector<Matrix3x3<real_t> > parentR;
		std::vector<Matrix3x3<real_t> > parentRinv;
		std::vector<PointGrid3D<coord_t> > parentRefGrid;
		std::vector<std::shared_ptr<Image<Vec3d> > > parentShifts;

		/*
		Interpolator object.
		*/
		LinearInterpolator<Vec3d, Vec3d, double, Vec3d> shiftInterpolator;


		Transformation(const Vec3<real_t>& c, real_t a, const Matrix3x3<real_t>& R, real_t normFact) :
			a(a),
			c(c),
			R(R),
			normFactor(normFact),
			shiftInterpolator(BoundaryCondition::Nearest)
		{
			if (!R.inverse(Rinv))
				throw ITLException("Rotation matrix is not invertible.");
		}


		/*
		Converts from world coordinates to image coordinates.
		Accounts for elastic deformations.
		*/
		Vec3<real_t> toImage(const Vec3<real_t>& x) const
		{
			/*
			take a point in world coordinates,
			convert to parent scan coordinates,
			if point is inside displacement field region,
				interpolate(u, v, w),
				calculate point in child scan coordinates
			else
				convert world point to child scan coordinates using similarity transformation
			*/
			
			
			Vec3<real_t> p(0, 0, 0);
			real_t count = 0;
			
			for (size_t n = 0; n < parenta.size(); n++)
			{
				Vec3<real_t> xP = 1 / parenta[n] * (parentRinv[n] * x) + parentc[n];

				if (parentRefGrid[n].contains(xP))
				{
					Vec3d shift = internals::projectPointToDeformed(Vec3d(xP), parentRefGrid[n], *parentShifts[n], shiftInterpolator);
					p += xP + Vec3<real_t>(shift);
					count++;
				}
			}

			if (count > 0)
				return p / count;
			

			// No parent scan displacement field contains the point, so use the similarity transformation.
			return 1 / a * (Rinv * x) + c;
		}

		/*
		Converts from image coordinates to world coordinates, but does not take elastic deformations into account.
		*/
		Vec3<real_t> toWorldApprox(const Vec3<real_t>& x) const 
		{
			return a * (R * (x - c));
		}

		/*
		Add a displacement field to this transformation.
		*/
		void addField(const Vec3<real_t>& c, real_t a, const Matrix3x3<real_t>& R, PointGrid3D<coord_t>& refGrid, std::shared_ptr<Image<Vec3d> > shifts)
		{
			Matrix3x3<real_t> Rinv;
			if (!R.inverse(Rinv))
				throw ITLException("Parent scan rotation matrix is not invertible.");

			parentc.push_back(c);
			parenta.push_back(a);
			parentR.push_back(R);
			parentRinv.push_back(Rinv);
			parentRefGrid.push_back(refGrid);
			parentShifts.push_back(shifts);
		}

		void zeroThirdDimension()
		{
			for (size_t n = 0; n < parentShifts.size(); n++)
			{
				Image<Vec3d>& shifts = *(parentShifts[n].get());

				for (coord_t n = 0; n < shifts.pixelCount(); n++)
				{
					shifts(n).z = 0;
				}
			}
		}

		/*
		Inpaints NaN values in displacement fields.
		*/
		void fixMissingValues()
		{
			for (size_t n = 0; n < parentShifts.size(); n++)
			{
				Image<Vec3d>& shifts = *(parentShifts[n].get());

				// Convert shifts to three separate images, convert NaN values to FLAG values.
				Image<double> u(shifts.dimensions());
				Image<double> v(shifts.dimensions());
				Image<double> w(shifts.dimensions());

				constexpr double FLAG = std::numeric_limits<double>::max();
				for (coord_t n = 0; n < shifts.pixelCount(); n++)
				{
					Vec3d U = shifts(n);
					u(n) = std::isnan(U.x) ? FLAG : U.x;
					v(n) = std::isnan(U.y) ? FLAG : U.y;
					w(n) = std::isnan(U.z) ? FLAG : U.z;
				}

				// Inpaint
				inpaintGarcia(u, FLAG, true);
				inpaintGarcia(v, FLAG, true);
				inpaintGarcia(w, FLAG, true);

				// Convert shifts back to vector image
				for (coord_t n = 0; n < shifts.pixelCount(); n++)
				{
					Vec3d U(u(n), v(n), w(n));
					shifts(n) = U;
				}
			}
		}
	};


	namespace internals
	{
		/*
		Tests if the given source bounds overlap with output image bounds.
		*/
		inline bool overlapsWithOutput(const Vec3c& cc, const Vec3c& cd, const Vec3c& outPos, const Vec3c& outSize)
		{
			coord_t xmin = cc.x;
			coord_t ymin = cc.y;
			coord_t zmin = cc.z;
			coord_t xmax = cd.x;
			coord_t ymax = cd.y;
			coord_t zmax = cd.z;

			// TODO: Add some padding so that deformations caused by the local deformation field do not result in cut image.

			// Return if this source image is not in the region of target image we want to process
			if (xmin > outPos.x + outSize.x ||
				ymin > outPos.y + outSize.y ||
				zmin > outPos.z + outSize.z ||
				xmax < outPos.x ||
				ymax < outPos.y ||
				zmax < outPos.z)
				return false;

			return true;
		}


		/*
		Reads a line from stream and returns it converted to the given type.
		*/
		template<typename out_t> out_t fromString(std::ifstream& in)
		{
			string line;
			getline(in, line);
			return itl2::fromString<out_t>(line);
		}

		template<typename real_t> std::tuple<Vec3<real_t>, real_t, Matrix3x3<real_t>, real_t, real_t, real_t> readcaR(std::ifstream& in)
		{
			// First line is a
			real_t a = fromString<real_t>(in);

			// Next three lines are components of c
			Vec3<real_t> c;
			c.x = fromString<real_t>(in);
			c.y = fromString<real_t>(in);
			c.z = fromString<real_t>(in);

			// Next three lines contain components of R

			std::stringstream ss1, ss2, ss3;
			std::string line;

			std::getline(in, line);
			ss1 << line;
			std::getline(in, line);
			ss2 << line;
			std::getline(in, line);
			ss3 << line;

			std::getline(ss1, line, ' ');
			real_t r00 = itl2::fromString<real_t>(line);
			std::getline(ss1, line, ' ');
			real_t r01 = itl2::fromString<real_t>(line);
			std::getline(ss1, line, ' ');
			real_t r02 = itl2::fromString<real_t>(line);

			std::getline(ss2, line, ' ');
			real_t r10 = itl2::fromString<real_t>(line);
			std::getline(ss2, line, ' ');
			real_t r11 = itl2::fromString<real_t>(line);
			std::getline(ss2, line, ' ');
			real_t r12 = itl2::fromString<real_t>(line);

			std::getline(ss3, line, ' ');
			real_t r20 = itl2::fromString<real_t>(line);
			std::getline(ss3, line, ' ');
			real_t r21 = itl2::fromString<real_t>(line);
			std::getline(ss3, line, ' ');
			real_t r22 = itl2::fromString<real_t>(line);

			Matrix3x3<real_t> R(r00, r01, r02,
				r10, r11, r12,
				r20, r21, r22);

			// Then normalization factor
			real_t normFact = fromString<real_t>(in);
			real_t normFactStd = fromString<real_t>(in);
			real_t meanDef = fromString<real_t>(in);

			return std::make_tuple(c, a, R, normFact, normFactStd, meanDef);
		}

		/*
		Reads refpoints saved by TransformationVer2::save.
		*/
		template<typename real_t> void readRefPoints(const std::string& prefix, PointGrid3D<coord_t>& refPoints, real_t& normFact, real_t& normFactStd, real_t& meanDef)
		{
			std::ifstream in(prefix + "_refpoints.txt");

			if (!in.good())
				throw ITLException(string("File not found: ") + prefix + "_refpoints.txt");

			PointGrid1D<coord_t> g1 = readPointGrid1D(in);
			PointGrid1D<coord_t> g2 = readPointGrid1D(in);
			PointGrid1D<coord_t> g3 = readPointGrid1D(in);
			refPoints = PointGrid3D<coord_t>(g1, g2, g3);

			normFact = fromString<real_t>(in);
			normFactStd = fromString<real_t>(in);
			meanDef = fromString<real_t>(in);
		}

		/*
		Reads shifts saved by TransformationVer2::save.
		*/
		template<typename real_t> void readShifts(const std::string& prefix, const PointGrid3D<coord_t>& refPoints, Image<Vec3<real_t> >& shifts)
		{
			shifts.init(refPoints.pointCounts());
			raw::read(shifts, prefix + "_shifts_" + toString(refPoints.xg.pointCount()) + "x" + toString(refPoints.yg.pointCount()) + "x" + toString(refPoints.zg.pointCount()) + ".raw");
		}

		template<typename real_t> struct TransformationVer2
		{
			/*
			Local to world rotation matrix
			*/
			Matrix3x3<real_t> R;

			/*
			Inverse of rotation matrix (world to local rotation)
			*/
			Matrix3x3<real_t> Rinv;

			/*
			Scaling factor
			*/
			real_t a = 1;

			/*
			Translation
			*/
			Vec3<real_t> c;

			/*
			Normalization factors for gray values.
			*/
			real_t normFactor = 0, normFactorStd = 1, meanDef = 0;

			/*
			World to parent refpoints, world to parent shifts, parent to me refpoints, parent to me shifts.
			*/
			std::vector<std::tuple<PointGrid3D<coord_t>, std::shared_ptr<Image<Vec3<real_t> > >, PointGrid3D<coord_t>, std::shared_ptr<Image<Vec3<real_t> > > > > parentToMeGrids;

			/*
			World to local displacement field.
			*/
			PointGrid3D<coord_t> worldGrid;
			std::shared_ptr<Image<Vec3<real_t> > > worldShifts;

			/*
			Converts (repoints, defpoints) pair to (refpoints, shifts) pair.
			*/
			void convertToShifts(const PointGrid3D<coord_t>& refPoints, Image<Vec3<real_t> >& defPoints)
			{
				for (coord_t z = 0; z < defPoints.depth(); z++)
				{
					for (coord_t y = 0; y < defPoints.height(); y++)
					{
						for (coord_t x = 0; x < defPoints.width(); x++)
						{
							defPoints(x, y, z) = defPoints(x, y, z) - Vec3<real_t>(refPoints(x, y, z));
						}
					}
				}
			}

			/*
			Reads transformation from file saved by elastic_stitcher_3D_2.py script.
			*/
			void readFromFile(const std::string& filename)
			{
				std::ifstream in(filename);

				// Read c, a, R, and normalization factors
				std::tuple<Vec3<real_t>, real_t, Matrix3x3<real_t>, real_t, real_t, real_t> caR = readcaR<real_t>(in);

				c = std::get<0>(caR);
				a = std::get<1>(caR);
				R = std::get<2>(caR);
				normFactor = std::get<3>(caR);
				normFactorStd = std::get<4>(caR);
				meanDef = std::get<5>(caR);

				R.inverse(Rinv);

				// Next line contains count of neighbours
				coord_t nbCount = round(fromString<double>(in));

				for (coord_t n = 0; n < nbCount; n++)
				{
					// World to parent local transformation name
					// NOTE: This transformation is saved in (refpoints, shifts) format!
					std::string wpPrefix;
					std::getline(in, wpPrefix);

					PointGrid3D<coord_t> wpRefPoints;
					std::shared_ptr<Image<Vec3<real_t> > > wpShifts(new Image<Vec3<real_t> >());
					real_t dummy, dummy2, dummy3;
					readRefPoints(wpPrefix, wpRefPoints, dummy, dummy2, dummy3);
					readShifts<real_t>(wpPrefix, wpRefPoints, *wpShifts);

					// Parent to me transformation name
					// NOTE: This transformation is saved in (refpoints, defpoints) format, so we convert defpoints to shifts!
					std::string pmPrefix;
					getline(in, pmPrefix);

					PointGrid3D<coord_t> pmRefPoints;
					Image<Vec3d> tmp;
					Image<float32_t> gof;
					std::shared_ptr<Image<Vec3<real_t> > > pmDefPoints(new Image<Vec3<real_t> >());
					double ddummy2, ddummy3, ddummy4;
					readBlockMatchResult(pmPrefix, pmRefPoints, tmp, gof, ddummy2, ddummy3, ddummy4);
					convert(tmp, *pmDefPoints);
					tmp.deleteData();
					convertToShifts(pmRefPoints, *pmDefPoints);

					parentToMeGrids.push_back(std::make_tuple(wpRefPoints, wpShifts, pmRefPoints, pmDefPoints));
				}
			}


			void inpaintNanShifts(Image<Vec3<real_t> >& shifts, float32_t tolerance = 0, bool indicateProgress = false)
			{
				// Convert shifts to three separate images, convert NaN values to FLAG values.
				Image<double> u(shifts.dimensions());
				Image<double> v(shifts.dimensions());
				Image<double> w(shifts.dimensions());

				constexpr double FLAG = std::numeric_limits<double>::max();
				for (coord_t n = 0; n < shifts.pixelCount(); n++)
				{
					Vec3<real_t> U = shifts(n);
					u(n) = std::isnan(U.x) ? FLAG : U.x;
					v(n) = std::isnan(U.y) ? FLAG : U.y;
					w(n) = std::isnan(U.z) ? FLAG : U.z;
				}

				// Inpaint
				inpaintGarcia(u, FLAG, indicateProgress, tolerance);
				inpaintGarcia(v, FLAG, indicateProgress, tolerance);
				inpaintGarcia(w, FLAG, indicateProgress, tolerance);

				// Convert shifts back to vector image
				for (coord_t n = 0; n < shifts.pixelCount(); n++)
				{
					Vec3d U(u(n), v(n), w(n));
					shifts(n).x = (real_t)U.x;
					shifts(n).y = (real_t)U.y;
					shifts(n).z = (real_t)U.z;
				}
			}

			/*
			Inpaints NaN values in parent to me displacement fields.
			*/
			void fixMissingValues()
			{
				for (size_t n = 0; n < parentToMeGrids.size(); n++)
				{
					Image<Vec3<real_t> >& shifts = *std::get<3>(parentToMeGrids[n]);

					inpaintNanShifts(shifts);
				}
			}

			/*
			Converts from image coordinates to world coordinates, but does not take local deformations into account.
			*/
			Vec3<real_t> toWorldApprox(const Vec3<real_t>& x) const
			{
				return a * (R * (x - c));
			}

			/*
			Calculates bounds of image of size srcDimensions in world coordinates, discards local deformations.
			*/
			void boundsInWorldCoordinates(const Vec3c& srcDimensions, Vec3c& cc, Vec3c& cd)
			{
				real_t w = (real_t)srcDimensions.x;
				real_t h = (real_t)srcDimensions.y;
				real_t d = (real_t)srcDimensions.z;

				// Calculate bounds of the input image in world coordinates.
				Vec3c c1 = floor(toWorldApprox(Vec3<real_t>(0, 0, 0)));
				Vec3c c2 = floor(toWorldApprox(Vec3<real_t>(w, 0, 0)));
				Vec3c c3 = floor(toWorldApprox(Vec3<real_t>(w, h, 0)));
				Vec3c c4 = floor(toWorldApprox(Vec3<real_t>(0, h, 0)));
				Vec3c c5 = floor(toWorldApprox(Vec3<real_t>(0, 0, d)));
				Vec3c c6 = floor(toWorldApprox(Vec3<real_t>(w, 0, d)));
				Vec3c c7 = floor(toWorldApprox(Vec3<real_t>(w, h, d)));
				Vec3c c8 = floor(toWorldApprox(Vec3<real_t>(0, h, d)));

				cc = c1;
				cc = min(cc, c2);
				cc = min(cc, c3);
				cc = min(cc, c4);
				cc = min(cc, c5);
				cc = min(cc, c6);
				cc = min(cc, c7);
				cc = min(cc, c8);

				cd = c1;
				cd = max(cd, c2);
				cd = max(cd, c3);
				cd = max(cd, c4);
				cd = max(cd, c5);
				cd = max(cd, c6);
				cd = max(cd, c7);
				cd = max(cd, c8);
			}

			/*
			Zeros z-component of world to local shifts.
			*/
			void zeroThirdDimension()
			{
				for (coord_t n = 0; n < worldShifts->pixelCount(); n++)
				{
					(*worldShifts)(n).z = 0;
				}
			}

			enum class DefAvgMode
			{
				/*
				Take the first encountered shift.
				*/
				TakeOldest,
				/*
				Take average of all points. This seems to cause artefacts if some of the shifts are bad.
				*/
				Average,
				/*
				Take shift with image that we have most overlap with.
				*/
				MostOverlap
			};

			/*
			Populate world to local transformation from (c, a, R) and parent-to-me transformations.
			*/
			void convertToWorld(Vec3c imageDimensions, bool allowLocalShifts)
			{
				//LinearInterpolator<Vec3<real_t>, Vec3<real_t>, real_t, Vec3<real_t> > shiftInterpolator = LinearInterpolator<Vec3<real_t>, Vec3<real_t>, real_t, Vec3<real_t> >(BoundaryCondition::Nearest);
				CubicInterpolator<Vec3<real_t>, Vec3<real_t>, real_t, Vec3<real_t> > shiftInterpolator = CubicInterpolator<Vec3<real_t>, Vec3<real_t>, real_t, Vec3<real_t> >(BoundaryCondition::Nearest);

				// Set to true to average shifts where there are multiple overlapping parent images.
				// Set to false to take shifts from the first (oldest) encountered parent image.
				// NOT averaging seems to make sharper images. Averaging seems to work better for some not so well matching images.
				// TODO: Make this a parameter?
				//constexpr bool average = false;
                //constexpr bool average = true;
				constexpr DefAvgMode avgMode = DefAvgMode::MostOverlap;

				imageDimensions = max(imageDimensions, Vec3c(1, 1, 1));

				Vec3c cc;
				Vec3c cd;
				boundsInWorldCoordinates(imageDimensions, cc, cd);

				// TODO: Make this step parameter or calculate it from grid step used in local phase correlation
				constexpr size_t step = 10;

				// Add step to the end so that we really fill the whole image region.
				cd.x += step;
				cd.y += step;
				if (imageDimensions.z > 1)
					cd.z += step;

				worldGrid = PointGrid3D<coord_t>(PointGrid1D<coord_t>(cc.x, cd.x, step), PointGrid1D<coord_t>(cc.y, cd.y, step), PointGrid1D<coord_t>(cc.z, cd.z, step));
				worldShifts.reset(new Image<Vec3<real_t> >(worldGrid.pointCounts()));

				// Fill shifts with FLAG value.
				const Vec3<real_t> FLAG(std::numeric_limits<real_t>::max(), std::numeric_limits<real_t>::max(), std::numeric_limits<real_t>::max());
				for (coord_t n = 0; n < worldShifts->pixelCount(); n++)
					(*worldShifts)(n) = FLAG;

				if (allowLocalShifts)
				{
					// Fill shifts with values derived from parent to me displacement fields, where available.
					size_t counter = 0;
#pragma omp parallel for if(worldShifts->pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
					for (coord_t z = 0; z < worldShifts->depth(); z++)
					{
						for (coord_t y = 0; y < worldShifts->height(); y++)
						{
							for (coord_t x = 0; x < worldShifts->width(); x++)
							{
								Vec3<real_t> worldP = Vec3<real_t>(worldGrid(x, y, z));

								// Convert worldP to parent coordinates using parent world to me grid,
								// resulting in parentP.
								// Convert parentP to own coordinates using parent to me fields. 
								// Average over all parents

								Vec3<real_t> avgP(0, 0, 0);
								real_t count = 0;

								for (size_t k = 0; k < parentToMeGrids.size(); k++)
								{

									PointGrid3D<coord_t>& worldToParentRefPoints = std::get<0>(parentToMeGrids[k]);
									Image<Vec3<real_t> >& worldToParentShifts = *std::get<1>(parentToMeGrids[k]);

									PointGrid3D<coord_t>& parentToMeRefPoints = std::get<2>(parentToMeGrids[k]);
									Image<Vec3<real_t> >& parentToMeShifts = *std::get<3>(parentToMeGrids[k]);


									if (worldToParentRefPoints.contains(worldP))
									{
										// Convert from world to parent coordinates
										Vec3<real_t> parentP = worldP + internals::projectPointToDeformed<real_t>(worldP, worldToParentRefPoints, worldToParentShifts, shiftInterpolator);

										if (parentToMeRefPoints.contains(parentP))
										{
											Vec3<real_t> myP = parentP + internals::projectPointToDeformed<real_t>(parentP, parentToMeRefPoints, parentToMeShifts, shiftInterpolator);
											switch (avgMode)
											{
											case DefAvgMode::Average:
												// Average all points
												avgP += myP;
												count++;
												break;
											case DefAvgMode::TakeOldest:
												// Take the first point encountered.
												if (count <= 0)
												{
													avgP = myP;
													count = 1;
												}
												break;
											case DefAvgMode::MostOverlap:
												// Take the point with neighbour with largest overlapping region
												if (parentToMeRefPoints.pointCount() > count)
												{
													avgP = myP;
													count = (real_t)parentToMeRefPoints.pointCount();
												}
												break;
											default:
												throw ITLException("Unsupported avg mode.");
											}
											//if (average || (!average && count <= 0))
											//{
											//	Vec3<real_t> myP = parentP + internals::projectPointToDeformed<real_t>(parentP, parentToMeRefPoints, parentToMeShifts, shiftInterpolator);
											//	avgP += myP;
											//	count++;
											//}
										}
										else
										{
											// Parent to me displacement field does not contain the point.
										}
									}
									else
									{
										// World to parent transformation is not valid at this point.
										// The world point does not overlap with parent image.
									}
								}

								if (count > 0)
								{
									switch (avgMode)
									{
									case DefAvgMode::Average:
										avgP /= count;
										break;
									case DefAvgMode::TakeOldest:
										break;
									case DefAvgMode::MostOverlap:
										break;
									default:
										throw ITLException("Unsupported avg mode.");
									}

									(*worldShifts)(x, y, z) = avgP - worldP;
								}
								//if (count > 0)
								//{
								//	avgP /= count;

								//	(*worldShifts)(x, y, z) = avgP - worldP;
								//}


							}
						}

						showThreadProgress(counter, worldShifts->depth());
					}

					// Set shifts near to valid values to nan to create a border of nans around values set from parent to me transformations.
					coord_t r = std::max<coord_t>(1, pixelRound<coord_t>(0.1 * worldShifts->dimensions().min()));
#pragma omp parallel for if(worldShifts->pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
					for (coord_t z = 0; z < worldShifts->depth(); z++)
					{
						for (coord_t y = 0; y < worldShifts->height(); y++)
						{
							for (coord_t x = 0; x < worldShifts->width(); x++)
							{
								Vec3<real_t> v = (*worldShifts)(x, y, z);
								if (v != FLAG && !std::isnan(v.x))
								{

									coord_t zmin = std::max<coord_t>(0, z - r);
									coord_t ymin = std::max<coord_t>(0, y - r);
									coord_t xmin = std::max<coord_t>(0, x - r);
									coord_t zmax = std::min<coord_t>(worldShifts->depth() - 1, z + r);
									coord_t ymax = std::min<coord_t>(worldShifts->height() - 1, y + r);
									coord_t xmax = std::min<coord_t>(worldShifts->width() - 1, x + r);
									for (coord_t zz = zmin; zz <= zmax; zz++)
									{
										for (coord_t yy = ymin; yy <= ymax; yy++)
										{
											for (coord_t xx = xmin; xx <= xmax; xx++)
											{
												// This should not cause erroneous output although it is a race condition
												Vec3<real_t>& p = (*worldShifts)(xx, yy, zz);
												if (p == FLAG)
													p = Vec3<real_t>(std::numeric_limits<real_t>::signaling_NaN(), std::numeric_limits<real_t>::signaling_NaN(), std::numeric_limits<real_t>::signaling_NaN());
											}
										}
									}
								}
							}
						}
					}
				}

				// Fill unset shifts with similarity transformation
#pragma omp parallel for if(worldShifts->pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t z = 0; z < worldShifts->depth(); z++)
				{
					for (coord_t y = 0; y < worldShifts->height(); y++)
					{
						for (coord_t x = 0; x < worldShifts->width(); x++)
						{
							if ((*worldShifts)(x, y, z) == FLAG)
							{
								Vec3<real_t> worldP = Vec3<real_t>(worldGrid(x, y, z));
								Vec3<real_t> simP = 1 / a * (Rinv * worldP) + c;
								(*worldShifts)(x, y, z) = simP - worldP;
							}
						}
					}
				}

				// Inpaint nan values.
				inpaintNanShifts(*worldShifts, 0.5, true);

				if (imageDimensions.z <= 1)
					zeroThirdDimension();
			}

			void save(const string& prefix)
			{
				raw::writed(*worldShifts, prefix + "_shifts");

				std::ofstream out;
				out.open(prefix + "_refpoints.txt");
				out << worldGrid.xg.first << ", " << worldGrid.xg.maximum << ", " << worldGrid.xg.step << std::endl;
				out << worldGrid.yg.first << ", " << worldGrid.yg.maximum << ", " << worldGrid.yg.step << std::endl;
				out << worldGrid.zg.first << ", " << worldGrid.zg.maximum << ", " << worldGrid.zg.step << std::endl;
				out << normFactor << std::endl;
				out << normFactorStd << std::endl;
				out << meanDef << std::endl;
			}

		};

		///*
		//Stitch src image and the corresponding transformation to output image, and update weight image.
		//This function only processes region starting at outPos and having the size of output image.
		//*/
		//template<typename pixel_t, typename output_t, typename real_t> void stitchOneVer2(
		//	const Image<pixel_t>& src,
		//	const PointGrid3D<coord_t>& refPoints, const Image<Vec3<real_t> >& shifts, pixel_t normFactor,
		//	const Vec3c& outPos, Image<output_t>& output, Image<real_t>& weight,
		//	bool normalize)
		//{
		//	//const Interpolator<real_t, pixel_t, real_t>& interpolator = NearestNeighbourInterpolator<real_t, pixel_t, real_t>(BoundaryCondition::Zero);
		//	//const Interpolator<real_t, pixel_t, real_t>& interpolator = LinearInvalidValueInterpolator<real_t, pixel_t, real_t>(BoundaryCondition::Zero, 0, 0);
		//	//const Interpolator<real_t, pixel_t, real_t>& interpolator = CubicInvalidValueInterpolator<real_t, pixel_t, real_t>(BoundaryCondition::Zero, 0, 0);
		//	const Interpolator<real_t, pixel_t, real_t>& interpolator = CubicInterpolator<real_t, pixel_t, real_t>(BoundaryCondition::Zero);
		//	//const Interpolator<Vec3<real_t>, Vec3<real_t>, real_t>& shiftInterpolator = LinearInterpolator<Vec3<real_t>, Vec3<real_t>, real_t, Vec3<real_t> >(BoundaryCondition::Nearest);
		//	const Interpolator<Vec3<real_t>, Vec3<real_t>, real_t>& shiftInterpolator = CubicInterpolator<Vec3<real_t>, Vec3<real_t>, real_t, Vec3<real_t> >(BoundaryCondition::Nearest);

		//	Vec3c cc(refPoints.xg.first, refPoints.yg.first, refPoints.zg.first);
		//	Vec3c cd(refPoints.xg.maximum, refPoints.yg.maximum, refPoints.zg.maximum);

		//	coord_t xmin = cc.x;
		//	coord_t ymin = cc.y;
		//	coord_t zmin = cc.z;
		//	coord_t xmax = cd.x;
		//	coord_t ymax = cd.y;
		//	coord_t zmax = cd.z;

		//	// TODO: Add some padding so that deformations caused by the local deformation field do not result in cut image.

		//	xmin = std::max(xmin, outPos.x);
		//	ymin = std::max(ymin, outPos.y);
		//	zmin = std::max(zmin, outPos.z);

		//	xmax = std::min(xmax, outPos.x + output.width());
		//	ymax = std::min(ymax, outPos.y + output.height());
		//	zmax = std::min(zmax, outPos.z + output.depth());

		//	Vec3c srcDimensions = src.dimensions();
		//	coord_t s = srcDimensions.min();

		//	if (src.dimensionality() < 3)
		//	{
		//		zmin = 0;
		//		zmax = 1;
		//	}

		//	std::cout << "Transforming..." << endl;
		//	// Process all pixels in the relevant region of the target image and find source image value at each location.
		//	size_t counter = 0;
		//	#pragma omp parallel for if((zmax-zmin)*(ymax-ymin)*(xmax-xmin) > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		//	for (coord_t z = zmin; z < zmax; z++)
		//	{
		//		for (coord_t y = ymin; y < ymax; y++)
		//		{
		//			for (coord_t x = xmin; x < xmax; x++)
		//			{
		//				// X is position in the output image
		//				Vec3<real_t> X((real_t)x, (real_t)y, (real_t)z);

		//				// Convert X to p, position in the input image.
		//				Vec3<real_t> p = X + internals::projectPointToDeformed(X, refPoints, shifts, shiftInterpolator);

		//				// Convert p to pdot, position in the input block.
		//				//Vec3<real_t> pdot = p - Vec3<real_t>(srcBlockPos);
		//				Vec3<real_t> pdot = p; // No src block support

		//				if (src.isInImage(pdot))
		//				{
		//					real_t pix = interpolator(src, pdot);
		//					if (pix != 0) // Don't process pixels that could not be interpolated (are given background value)
		//					{
		//						if (normalize)
		//							pix += normFactor;

		//						real_t w1 = 2 * std::min(p.x, srcDimensions.x - 1 - p.x) / s;
		//						real_t w2 = 2 * std::min(p.y, srcDimensions.y - 1 - p.y) / s;
		//						real_t w3 = 2 * std::min(p.z, srcDimensions.z - 1 - p.z) / s;

		//						if (src.dimensionality() < 3)
		//							w3 = 1;

		//						real_t ww = w1 * w2 * w3;

		//						if (ww > 0)
		//						{
		//							// Finally transform to the coordinates of the region of the target image that we are processing
		//							coord_t xo = x - outPos.x;
		//							coord_t yo = y - outPos.y;
		//							coord_t zo = z - outPos.z;

		//							if (xo >= 0 && yo >= 0 && zo >= 0 && xo < output.width() && yo < output.height() && zo < output.depth())
		//							{
		//								output(xo, yo, zo) += pixelRound<output_t>(pix * ww);
		//								weight(xo, yo, zo) += ww;
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}

		//		showThreadProgress(counter, zmax - zmin);
		//	}
		//}

		/*
		Stitch src image and the corresponding transformation to output image, and update weight image.
		This function only processes region starting at outPos and having the size of output image.
		This version can also calculate standard deviation of overlapping images in the overlapping regions.
		After calling the method for all input images:
		- image mean does not need further processing.
		- image S must be divided by image weight and to get standard deviation, sqrt must be taken.
		*/
		template<typename pixel_t, typename real_t> void stitchOneVer3(
			const Image<pixel_t>& src,
			const PointGrid3D<coord_t>& refPoints, const Image<Vec3<real_t> >& shifts,
			real_t normFactor, real_t normFactorStd, real_t meanDef,
			const Vec3c& outPos, Image<real_t>& mean, Image<real_t>& weight, Image<real_t>* S,
			bool normalize)
		{
			//const Interpolator<real_t, pixel_t, real_t>& interpolator = NearestNeighbourInterpolator<real_t, pixel_t, real_t>(BoundaryCondition::Zero);
			//const Interpolator<real_t, pixel_t, real_t>& interpolator = LinearInvalidValueInterpolator<real_t, pixel_t, real_t>(BoundaryCondition::Zero, 0, 0);
			const Interpolator<real_t, pixel_t, real_t>& interpolator = CubicInvalidValueInterpolator<real_t, pixel_t, real_t>(BoundaryCondition::Zero, 0, 0);
			//const Interpolator<real_t, pixel_t, real_t>& interpolator = CubicInterpolator<real_t, pixel_t, real_t>(BoundaryCondition::Zero);

			// NOTE: Cubic interpolator will overshoot, linear is rough. Perhaps monotone cubic would be the best for shifts?
			//const Interpolator<Vec3<real_t>, Vec3<real_t>, real_t>& shiftInterpolator = LinearInterpolator<Vec3<real_t>, Vec3<real_t>, real_t, Vec3<real_t> >(BoundaryCondition::Nearest);
			const Interpolator<Vec3<real_t>, Vec3<real_t>, real_t>& shiftInterpolator = CubicInterpolator<Vec3<real_t>, Vec3<real_t>, real_t, Vec3<real_t> >(BoundaryCondition::Nearest);

			Vec3c cc(refPoints.xg.first, refPoints.yg.first, refPoints.zg.first);
			Vec3c cd(refPoints.xg.maximum, refPoints.yg.maximum, refPoints.zg.maximum);

			coord_t xmin = cc.x;
			coord_t ymin = cc.y;
			coord_t zmin = cc.z;
			coord_t xmax = cd.x;
			coord_t ymax = cd.y;
			coord_t zmax = cd.z;

			// TODO: Add some padding so that deformations caused by the local deformation field do not result in cut image.

			xmin = std::max(xmin, outPos.x);
			ymin = std::max(ymin, outPos.y);
			zmin = std::max(zmin, outPos.z);

			xmax = std::min(xmax, outPos.x + mean.width());
			ymax = std::min(ymax, outPos.y + mean.height());
			zmax = std::min(zmax, outPos.z + mean.depth());

			Vec3c srcDimensions = src.dimensions();
			coord_t s = srcDimensions.min();

			if (src.dimensionality() < 3)
			{
				zmin = 0;
				zmax = 1;
			}

			std::cout << "Transforming..." << std::endl;
			// Process all pixels in the relevant region of the target image and find source image value at each location.
			size_t counter = 0;
#pragma omp parallel for if((zmax-zmin)*(ymax-ymin)*(xmax-xmin) > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t z = zmin; z < zmax; z++)
			{
				for (coord_t y = ymin; y < ymax; y++)
				{
					for (coord_t x = xmin; x < xmax; x++)
					{
						// X is position in the output image
						Vec3<real_t> X((real_t)x, (real_t)y, (real_t)z);

						// Convert X to p, position in the input image.
						Vec3<real_t> p = X + internals::projectPointToDeformed(X, refPoints, shifts, shiftInterpolator);

						// Convert p to pdot, position in the input block.
						//Vec3<real_t> pdot = p - Vec3<real_t>(srcBlockPos);
						Vec3<real_t> pdot = p; // No src block support

						if (src.isInImage(pdot))
						{
							real_t pix = interpolator(src, pdot);
							if (pix != 0) // Don't process pixels that could not be interpolated (are given background value)
							{
								if (normalize)
									pix = (pix - meanDef) * normFactorStd + meanDef + normFactor;
									//pix += normFactor;

								real_t w1 = 2 * std::min(p.x, srcDimensions.x - 1 - p.x) / s;
								real_t w2 = 2 * std::min(p.y, srcDimensions.y - 1 - p.y) / s;
								real_t w3 = 2 * std::min(p.z, srcDimensions.z - 1 - p.z) / s;

								if (src.dimensionality() < 3)
									w3 = 1;

								real_t ww = w1 * w2 * w3;

								if (ww > 0)
								{
									// Finally transform to the coordinates of the region of the target image that we are processing
									coord_t xo = x - outPos.x;
									coord_t yo = y - outPos.y;
									coord_t zo = z - outPos.z;

									if (xo >= 0 && yo >= 0 && zo >= 0 && xo < mean.width() && yo < mean.height() && zo < mean.depth())
									{
										// Calculate mean and, if requested, standard deviation.
										// This uses algorithm from West, D. H. D. (1979). "Updating Mean and Variance Estimates: An Improved Method". 
										// See also https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance -> Weighted incremental algorithm

										real_t& wSum = weight(xo, yo, zo);
										wSum += ww;
										
										real_t& meanNew = mean(xo, yo, zo);
										real_t meanOld = meanNew;

										meanNew = meanOld + (ww / wSum) * (pix - meanOld);

										if (S)
											(*S)(xo, yo, zo) += ww * (pix - meanOld) * (pix - meanNew);

									}
								}
							}
						}
					}
				}

				showThreadProgress(counter, zmax - zmin);
			}
		}
	}


	/*
	Use to determine world to local transformation to one image.
	*/
	inline void determineWorldToLocal(const string& transfFile, const Vec3c& imageSize, const string& worldToLocalPrefix, bool allowLocalShifts)
	{
		internals::TransformationVer2<float32_t> transformation;
		transformation.readFromFile(transfFile);
		transformation.fixMissingValues();
		transformation.convertToWorld(imageSize, allowLocalShifts);
		transformation.save(worldToLocalPrefix);
	}

	//template<typename pixel_t> void stitchVer2(const string& indexFile, const Vec3c& outputPos, const Vec3c& outputSize, Image<pixel_t>& output, bool normalize)
	//{

	//	// Read index file
	//	std::cout << "Reading index file..." << std::endl;
	//	std::ifstream in(indexFile);
	//	std::vector<std::tuple<string, string> > images;
	//	while (in.good())
	//	{
	//		string imgFile;
	//		string wlFile;
	//		getline(in, imgFile);
	//		getline(in, wlFile);
	//		if(imgFile.length() > 0)
	//			images.push_back(std::make_tuple(imgFile, wlFile));
	//	}

	//	std::cout << "Stitching..." << std::endl;

	//	Image<float32_t> weight(outputSize);
	//	Image<float32_t> out(outputSize);

	//	// Process each image that is to be stitched
	//	for (auto it = images.begin(); it != images.end(); it++)
	//	{
	//		string imgFile = std::get<0>(*it);
	//		string wlPrefix = std::get<1>(*it);

	//		//raw::expandRawFilename(imgFile);
	//		std::cout << "Stitching " << imgFile << std::endl;

	//		Vec3c srcDimensions;
	//		ImageDataType dt;
	//		if (!io::getInfo(imgFile, srcDimensions, dt))
	//			throw ITLException("Unable to find dimensions of source image.");

	//		if(dt != imageDataType<pixel_t>())
	//			throw ITLException(string("Data type of input image ") + imgFile + " is " + toString(dt) + " but " + toString(imageDataType<pixel_t>()) + " was expected.");

	//		PointGrid3D<coord_t> refPoints;
	//		float32_t normFact;
	//		internals::readRefPoints(wlPrefix, refPoints, normFact);

	//		Vec3c cc(refPoints.xg.first, refPoints.yg.first, refPoints.zg.first);
	//		Vec3c cd(refPoints.xg.maximum, refPoints.yg.maximum, refPoints.zg.maximum);
	//		
	//		if (internals::overlapsWithOutput(cc, cd, outputPos, outputSize))
	//		{

	//			//// If source image is big, load it in blocks.
	//			//// Python script assumes 116 GB blocks for output (2 x float32_t), so we use 80 GB blocks (1 x uint16_t) for input.
	//			//Vec3c srcBlockSize(3500, 3500, 3500);
	//			//for (coord_t z = 0; z < srcDimensions.z; z += srcBlockSize.z)
	//			//{
	//			//	for (coord_t y = 0; y < srcDimensions.y; y += srcBlockSize.y)
	//			//	{
	//			//		for (coord_t x = 0; x < srcDimensions.x; x += srcBlockSize.x)
	//			//		{
	//			//			std::cout << "Processing source image block (" << x << ", " << y << ", " << z << ")..." << std::endl;

	//			//			// Take extra pixel layer so that interpolation succeeds at (x, y, z)
	//			//			Vec3c srcBlockPos = Vec3c(x, y, z) - Vec3c(1, 1, 1);
	//			//			Vec3c srcBlockEnd = srcBlockPos + srcBlockSize + 2 * Vec3c(1, 1, 1);

	//			//			clamp(srcBlockPos, Vec3c(0, 0, 0), srcDimensions);
	//			//			clamp(srcBlockEnd, Vec3c(0, 0, 0), srcDimensions);

	//			//			Image<pixel_t> srcBlock(srcBlockEnd - srcBlockPos);
	//			//			raw::readBlock(srcBlock, imgFile, srcDimensions, srcBlockPos);

	//			//			internals::stitchOne<pixel_t, float32_t, float32_t>(srcBlockPos, srcDimensions, srcBlock, transformation, outputPos, out, weight, normalize, cc, cd, true);
	//			//		}
	//			//	}
	//			//}

	//			Image<Vec3<float32_t> > shifts;
	//			internals::readShifts(wlPrefix, refPoints, shifts);

	//			Image<pixel_t> src(srcDimensions);
	//			//raw::read(src, imgFile);
	//			io::read(src, imgFile);
	//			internals::stitchOneVer2<pixel_t, float32_t, float32_t>(src, refPoints, shifts, pixelRound<pixel_t>(normFact), outputPos, out, weight, normalize);
	//		}
	//	}

	//	// Divide output by weight
	//	std::cout << "Final division..." << std::endl;
	//	divide(out, weight);

	//	// Weight is not needed anymore
	//	weight.deleteData();

	//	// TODO: This can be skipped (values assigned directly to output) if output data type is float32_t
	//	std::cout << "Final data type conversion..." << std::endl;
	//	convert(out, output);
	//}

	template<typename pixel_t> void stitchVer3(const string& indexFile, const Vec3c& outputPos, const Vec3c& outputSize, Image<pixel_t>& output, Image<pixel_t>* std, bool normalize)
	{

		// Read index file
		std::cout << "Reading index file..." << std::endl;
		std::ifstream in(indexFile);
		std::vector<std::tuple<string, string> > images;
		while (in.good())
		{
			string imgFile;
			string wlFile;
			getline(in, imgFile);
			getline(in, wlFile);
			if (imgFile.length() > 0)
				images.push_back(std::make_tuple(imgFile, wlFile));
		}

		std::cout << "Stitching..." << std::endl;

		Image<float32_t> weight(outputSize);
		Image<float32_t> out(outputSize);

		Image<float32_t> stdtmp(1);
		if (std)
			stdtmp.ensureSize(outputSize);

		// Process each image that is to be stitched
		for (auto it = images.begin(); it != images.end(); it++)
		{
			string imgFile = std::get<0>(*it);
			string wlPrefix = std::get<1>(*it);

			std::cout << "Stitching " << imgFile << std::endl;

			Vec3c srcDimensions;
			ImageDataType dt;
			string reason;
			if (!io::getInfo(imgFile, srcDimensions, dt, reason))
				throw ITLException(string("Unable to find dimensions of source image. ") + reason);

			if (dt != imageDataType<pixel_t>())
				throw ITLException(string("Data type of input image ") + imgFile + " is " + toString(dt) + " but " + toString(imageDataType<pixel_t>()) + " was expected.");

			PointGrid3D<coord_t> refPoints;
			float32_t normFact, normFactStd, meanDef;
			internals::readRefPoints(wlPrefix, refPoints, normFact, normFactStd, meanDef);

			Vec3c cc(refPoints.xg.first, refPoints.yg.first, refPoints.zg.first);
			Vec3c cd(refPoints.xg.maximum, refPoints.yg.maximum, refPoints.zg.maximum);

			if (internals::overlapsWithOutput(cc, cd, outputPos, outputSize))
			{

				//// If source image is big, load it in blocks.
				//// Python script assumes 116 GB blocks for output (2 x float32_t), so we use 80 GB blocks (1 x uint16_t) for input.
				//Vec3c srcBlockSize(3500, 3500, 3500);
				//for (coord_t z = 0; z < srcDimensions.z; z += srcBlockSize.z)
				//{
				//	for (coord_t y = 0; y < srcDimensions.y; y += srcBlockSize.y)
				//	{
				//		for (coord_t x = 0; x < srcDimensions.x; x += srcBlockSize.x)
				//		{
				//			std::cout << "Processing source image block (" << x << ", " << y << ", " << z << ")..." << std::endl;

				//			// Take extra pixel layer so that interpolation succeeds at (x, y, z)
				//			Vec3c srcBlockPos = Vec3c(x, y, z) - Vec3c(1, 1, 1);
				//			Vec3c srcBlockEnd = srcBlockPos + srcBlockSize + 2 * Vec3c(1, 1, 1);

				//			clamp(srcBlockPos, Vec3c(0, 0, 0), srcDimensions);
				//			clamp(srcBlockEnd, Vec3c(0, 0, 0), srcDimensions);

				//			Image<pixel_t> srcBlock(srcBlockEnd - srcBlockPos);
				//			raw::readBlock(srcBlock, imgFile, srcDimensions, srcBlockPos);

				//			internals::stitchOne<pixel_t, float32_t, float32_t>(srcBlockPos, srcDimensions, srcBlock, transformation, outputPos, out, weight, normalize, cc, cd, true);
				//		}
				//	}
				//}

				Image<Vec3<float32_t> > shifts;
				internals::readShifts(wlPrefix, refPoints, shifts);

				Image<pixel_t> src(srcDimensions);
				io::read(src, imgFile);
				internals::stitchOneVer3<pixel_t, float32_t>(src, refPoints, shifts, normFact, normFactStd, meanDef, outputPos, out, weight, std ? &stdtmp : nullptr, normalize);
			}
		}

		if (std)
		{
			// Divide S by weight
			std::cout << "Final division..." << std::endl;
			divide(stdtmp, weight);
			squareRoot(stdtmp);
		}

		// Weight is not needed anymore
		weight.deleteData();

		// TODO: This can be skipped (values assigned directly to output) if output data type is float32_t
		std::cout << "Final data type conversion..." << std::endl;
		convert(out, output);

		if (std)
		{
			convert(stdtmp, *std);
		}
	}
}
