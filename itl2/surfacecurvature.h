#pragma once

#include <vector>
#include <map>
#include <queue>
#include <set>

#include "image.h"
#include "math/vec3.h"
#include "io/raw.h"
#include "transform.h"
#include "interpolation.h"
#include "sphere.h"
#include "math/matrix.h"


namespace itl2
{
	

	namespace internals
	{
		/**
		Tests if pixel at given location is nonzero and if it has one or more zero 6-connected neighbours.
		*/
		template<typename pixel_t> bool isSurfacePixel6(const Image<pixel_t>& img, const Vec3c& p, BoundaryCondition bc)
		{
			if (img(p) == 0)
				return false;

			return getPixelSafe(img, p.x - 1, p.y, p.z, bc) == 0 ||
				getPixelSafe(img, p.x + 1, p.y, p.z, bc) == 0 ||
				getPixelSafe(img, p.x, p.y - 1, p.z, bc) == 0 ||
				getPixelSafe(img, p.x, p.y + 1, p.z, bc) == 0 ||
				getPixelSafe(img, p.x, p.y, p.z - 1, bc) == 0 ||
				getPixelSafe(img, p.x, p.y, p.z + 1, bc) == 0;
		}

		/**
		Finds neighbours of p that are surface pixels.
		*/
		template<typename pixel_t> std::vector<Vec3c> surfaceNeighbours(const Image<pixel_t>& img, const Vec3c& p, BoundaryCondition bc)
		{
			std::vector<Vec3c> nbs;

			for (coord_t z = -1; z <= 1; z++)
			{
				for (coord_t y = -1; y <= 1; y++)
				{
					for (coord_t x = -1; x <= 1; x++)
					{
						if (x != 0 || y != 0 || z != 0)
						{
							Vec3c pn = p + Vec3c(x, y, z);
							if (img.isInImage(pn) && isSurfacePixel6(img, pn, bc))
								nbs.push_back(pn);
						}
					}
				}
			}

			return nbs;
		}

		/**
		Finds surface normal (up to sign) for the given set of points.
		Returns tangents and normal, tuple(t1, t2, n).
		*/
		inline std::tuple<Vec3d, Vec3d, Vec3d> findFrame(const std::vector<Vec3d>& points)
		{
			// Mean position
			Vec3d pbar = mean<Vec3d, Vec3d, double>(points);

			Matrix3x3d CI;
			for (const auto& p : points)
			{
				Vec3d d = p - pbar;
				CI += Matrix3x3d::outer(d, d);
			}
			CI /= (double)points.size();

			Vec3d t1, t2, n;
			double l1, l2, ln;
			CI.eigsym(t1, t2, n, l1, l2, ln);


			return std::make_tuple(t1, t2, n);
		}
	}

	

	/**
	Calculates curvature of surfaces.
	Uses quadratic surface fitting algorithms in Petitjean - A Survey of Methods for Recovering Quadrics in Triangle Meshes.
	Surface normals are determined using principal component analysis of the covariance matrix of surface points near the center point.
	The surface normal orientation is chosen so that it points toward background voxels.
	The curvature is determined by transforming surface points near center point to a coordinate system where z-direction is
	parallel to surface normal, and then fitting a surface
	f(x, y) = a x^2 + b x y + c y^2 + d
	to the tranformed points. The curvature values and directions are calculated from the coefficients a, b and c.
	Directions are then transformed back to original coordinates.
	@param img Image containing the geometry. Nonzero pixels are assumed to be foreground.
	@param radius Radius of neighbourhood around surface point considered when determining curvature. Typically e.g. 5 gives good results.
	@param kappa1, kappa2 If not nullptr, largest and smallest principal curvature will be placed to these images.
	@param dir1, dir2 If not nullptr, directions of largest and smallest principal curvature will be filled to these images.
	@param bc Boundary condition.
	@param nonSurfaceValue Value that is used to fill the non-surface points in the kappa1 and kappa2 images.
	@param showProgressInfo Flag indicating whether or not progress information should be printed.
	*/
	template<typename pixel_t, typename out_t> void surfaceCurvature(const Image<pixel_t>& img, float32_t radius, Image<out_t>* kappa1 = nullptr, Image<out_t>* kappa2 = nullptr, Image<Vec3f>* dir1 = nullptr, Image<Vec3f>* dir2 = nullptr, BoundaryCondition bc = BoundaryCondition::Nearest, out_t nonSurfaceValue = 0, bool showProgressInfo = true)
	{
// Define this to save debug information for plotting in Matlab
//#define SAVE_DEBUG
#if defined(SAVE_DEBUG)
std::vector<Vec3f> vsurfacePoints;
std::vector<Vec3f> vsurfaceNormals;
std::vector<float32_t> vcurvature1, vcurvature2;
std::vector<Vec3f> vdir1, vdir2;
#endif
		// Initialize output images
		if (kappa1)
		{
			kappa1->mustNotBe(img);
			kappa1->ensureSize(img);
			setValue(*kappa1, nonSurfaceValue);
		}

		if (kappa2)
		{
			kappa2->mustNotBe(img);
			kappa2->ensureSize(img);
			setValue(*kappa2, nonSurfaceValue);
		}

		if (dir1)
		{
			dir1->mustNotBe(img);
			dir1->ensureSize(img);
		}

		if (dir2)
		{
			dir2->mustNotBe(img);
			dir2->ensureSize(img);
		}

		size_t counter = 0;
		#pragma omp parallel for if(img.pixelCount() >= PARALLELIZATION_THRESHOLD)
		for (coord_t z = 0; z < img.depth(); z++)
		{
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					Vec3c p0(x, y, z);
					if (internals::isSurfacePixel6(img, p0, bc))
					{
						// 1. Find surface points around p
						// -------------------------------

						// This version estimates surface normal on the fly and does not allow neighbours that indicate bending of the surface under itself.
						// tuple<point, geodesic distance>
						auto distanceComparer = [](const std::tuple<Vec3c, float32_t>& e1, const std::tuple<Vec3c, float32_t>& e2) { return std::get<1>(e1) > std::get<1>(e2); };
						std::priority_queue<std::tuple<Vec3c, float32_t>, std::vector<std::tuple<Vec3c, float32_t> >, std::reference_wrapper<decltype(distanceComparer)> > q(distanceComparer);
						q.push(std::make_tuple(p0, 0.0f));

						// Get neighbours that are near to the center point and the corresponding geodesic distances to the center point.
						std::map<Vec3c, float32_t, decltype(vecComparer<coord_t>)*> nbs(&vecComparer<coord_t>);
						std::set<Vec3c, decltype(vecComparer<coord_t>)*> backgroundNeighbours(&vecComparer<coord_t>);
						Vec3d normalEstimate;
						while (!q.empty())
						{
							std::priority_queue<std::tuple<Vec3c, float32_t>, std::vector<std::tuple<Vec3c, float32_t> >, std::reference_wrapper<decltype(distanceComparer)> > q2(distanceComparer);

							while (!q.empty())
							{
								auto item = q.top();
								Vec3c pi = std::get<0>(item);
								float32_t dist = std::get<1>(item);
								q.pop();

								if (nbs.count(pi) <= 0 || dist < nbs[pi])
								{
									nbs[pi] = dist;

									const std::vector<Vec3c>& snbs = internals::surfaceNeighbours(img, pi, bc);
									for (const Vec3c& t : snbs)
									{
										float32_t newDist = dist + (t - pi).norm();
										if (NumberUtils<float32_t>::lessThanOrEqual(newDist, radius))
										{
											if(normalEstimate == Vec3d(0, 0, 0) ||
												NumberUtils<double>::greaterThanOrEqual(projectToPlane(Vec3d(t - p0), normalEstimate).normSquared(), projectToPlane(Vec3d(pi - p0), normalEstimate).normSquared(), 1e-5)) // Only add if the new point is further away from the p0 point than its parent.
											q2.push(std::make_tuple(t, newDist));
										}
									}

								}
							}

							// Update surface normal estimate
							std::vector<Vec3d> points;
							points.reserve(nbs.size());
							for (auto const& p : nbs)
								points.push_back(Vec3d(p.first));

							if(points.size() >= 3)
								std::tie(std::ignore, std::ignore, normalEstimate) = internals::findFrame(points);

							// Note: q is empty in this phase!
							q.swap(q2);
						}

						//// This version does not adjust surface normal. Therefore it does not work correctly for details smaller than 2*radius.
						//// tuple<point, geodesic distance>
						//auto distanceComparer = [](const tuple<Vec3c, float32_t>& e1, const tuple<Vec3c, float32_t>& e2) { return get<1>(e1) > get<1>(e2); };
						//std::priority_queue<tuple<Vec3c, float32_t>, std::vector<tuple<Vec3c, float32_t> >, decltype(distanceComparer) > q(distanceComparer);
						//q.push(make_tuple(p0, 0.0f));

						//// Get neighbours that are near to the center point and the corresponding geodesic distances to the center point.
						//std::map<Vec3c, float32_t, decltype(vecComparer<coord_t>)*> nbs(&vecComparer<coord_t>);
						//std::set<Vec3c, decltype(vecComparer<coord_t>)*> backgroundNeighbours(&vecComparer<coord_t>);
						//while (!q.empty())
						//{
						//	auto item = q.top();
						//	Vec3c pi = get<0>(item);
						//	float32_t dist = get<1>(item);
						//	q.pop();

						//	if (nbs.count(pi) <= 0 || dist < nbs[pi])
						//	{
						//		nbs[pi] = dist;

						//		const std::vector<Vec3c>& snbs = internals::surfaceNeighbours(img, pi, bc);
						//		for (const Vec3c& t : snbs)
						//		{
						//			float32_t newDist = dist + (t - pi).norm();
						//			if (newDist <= radius)
						//				q.push(make_tuple(t, newDist));
						//		}

						//	}
						//}

						// 2. Find surface normal by PCA
						// -----------------------------

						std::vector<Vec3d> points;
						points.reserve(nbs.size());
						for (auto const& p : nbs)
							points.push_back(Vec3d(p.first));

						auto [t1, t2, n] = internals::findFrame(points);

						// Orient the normal such that majority of background points neighbouring p0
						// are below the surface defined by p0 and n.
						// I.e. the normal is oriented so that it points to the direction of background pixels when
						// looking from p0.
						// TODO: If this procedure fails, try determining the neighbour of p0 where the normal points to,
						// and if that is foreground, flip the normal.
						double backgroundCount = 0;
						double okCount = 0;
						for (coord_t z = -1; z <= 1; z++)
						{
							for (coord_t y = -1; y <= 1; y++)
							{
								for (coord_t x = -1; x <= 1; x++)
								{
									if (x != 0 || y != 0 || z != 0)
									{
										Vec3c pn = p0 + Vec3c(x, y, z);
										if (getPixelSafe(img, pn.x, pn.y, pn.z, bc) == 0)
										{
											backgroundCount++;
											Vec3d d = Vec3d(pn) - Vec3d(p0);
											if (d.dot(n) > 0)
												okCount++;
										}
									}
								}
							}
						}
						if (okCount < backgroundCount / 2)
							n *= -1;

						// This is the initial test.
						//// Store neighbouring background points
						//for (coord_t z = -1; z <= 1; z++)
						//{
						//	for (coord_t y = -1; y <= 1; y++)
						//	{
						//		for (coord_t x = -1; x <= 1; x++)
						//		{
						//			if (x != 0 || y != 0 || z != 0)
						//			{
						//				Vec3c pn = pi + Vec3c(x, y, z);
						//				if (getPixelSafe(img, pn.x, pn.y, pn.z, bc) == 0)
						//				{
						//					backgroundNeighbours.emplace(pn);
						//				}
						//			}
						//		}
						//	}
						//}

						//int above = 0;
						//for(const Vec3c pn : backgroundNeighbours)
						//{
						//	Vec3d d = Vec3d(pn) - Vec3d(p0);
						//	if (d.dot(n) > 0)
						//		above++;
						//}
						//if (above <= backgroundNeighbours.size() / 2)
						//	n *= -1;


						// 3. Transform point coordinates to [t1, t2, n] coordinate system
						//    and fit surface to the transformed points.
						// ----------------------------------------------------------------

						// Create matrix that rotates t1 to [1, 0, 0], t1 to [0, 1, 0] and n to [0, 0, 1].
						Matrix3x3d R(t1, t2, n);
						R.transpose();
						//Vec3d x = R * t1;
						//Vec3d y = R * t2;
						//Vec3d z = R * n;

						// Transform and calculate M and Z matrices in the same loop.
						// M = [xi^2, xi * yi, yi^2, 1]
						// Z = [zi]
						Matrix Z(points.size(), 1);
						Matrix M(points.size(), 4);

						for(size_t n = 0; n < points.size(); n++)
						{
							Vec3d pn = points[n];
							Vec3d pp = R * (pn - Vec3d(p0)); // p0 is at the origin.
							M(n, 0) = pp.x * pp.x;
							M(n, 1) = pp.x * pp.y;
							M(n, 2) = pp.y * pp.y;
							M(n, 3) = 1;
							Z(n, 0) = pp.z;
						}

						// Solve M * C = Z to find a, b, and c such that surface
						// a * x^2 + b * x * y + c * y^2 + d
						// fits best the points.
						double a, b, c, d;
						try
						{
							// Least-squares solution
							Matrix C = M.solve(Z);
							a = C(0);
							b = C(1);
							c = C(2);
							d = C(3);
						}
						catch (ITLException)
						{
							// The matrix is rank deficient. There are probably not enough points to fit the surface.
							a = std::numeric_limits<double>::signaling_NaN();
							b = a;
							c = a;
							d = a;
						}

						// Calculate minimum and maximum principal curvature and
						// their directions in the rotated coordinate system.
						double k1 = a + c + sqrt((a - c) * (a - c) + b * b);
						double k2 = a + c - sqrt((a - c) * (a - c) + b * b);
						double alpha = 0.5 * atan2(b, a - c);
						Vec3d k1h(cos(alpha), sin(alpha), 0);
						Vec3d k2h(-sin(alpha), cos(alpha), 0);

						// Transform directions back to world coordinates
						R.transpose();
						k1h = R * k1h;
						k2h = R * k2h;


						// Fill the curvature values to output images
						if (kappa1)
							(*kappa1)(p0) = pixelRound<out_t>(k1);
						if (kappa2)
							(*kappa2)(p0) = pixelRound<out_t>(k2);
						if (dir1)
							(*dir1)(p0) = pixelRound<Vec3f>(k1h);
						if (dir2)
							(*dir2)(p0) = pixelRound<Vec3f>(k2h);

#if defined(SAVE_DEBUG)
						// For debugging
						vsurfacePoints.push_back(Vec3f(p0));
						vsurfaceNormals.push_back(Vec3f(n));
						vcurvature1.push_back((float32_t)k1);
						vcurvature2.push_back((float32_t)k2);
						vdir1.push_back(Vec3f(k1h));
						vdir2.push_back(Vec3f(k2h));
#endif
					}
				}
			}

			showThreadProgress(counter, img.depth(), showProgressInfo);
		}

#if defined(SAVE_DEBUG)
		// Save for visualization
		{
			Image<float32_t> vi(3, vsurfacePoints.size());
			for (size_t n = 0; n < vsurfacePoints.size(); n++)
			{
				vi(0, n) = vsurfacePoints[n].x;
				vi(1, n) = vsurfacePoints[n].y;
				vi(2, n) = vsurfacePoints[n].z;
			}
			raw::writed(vi, "./surface_curvature/surface_points");
		}

		{
			Image<float32_t> dir(3, vsurfaceNormals.size());
			for (size_t n = 0; n < vsurfaceNormals.size(); n++)
			{
				dir(0, n) = vsurfaceNormals[n].x;
				dir(1, n) = vsurfaceNormals[n].y;
				dir(2, n) = vsurfaceNormals[n].z;
			}
			raw::writed(dir, "./surface_curvature/surface_normals");
		}

		{
			Image<float32_t> dir(3, vdir1.size());
			for (size_t n = 0; n < vdir1.size(); n++)
			{
				dir(0, n) = vdir1[n].x;
				dir(1, n) = vdir1[n].y;
				dir(2, n) = vdir1[n].z;
			}
			raw::writed(dir, "./surface_curvature/dir1");
		}

		{
			Image<float32_t> dir(3, vdir2.size());
			for (size_t n = 0; n < vdir2.size(); n++)
			{
				dir(0, n) = vdir2[n].x;
				dir(1, n) = vdir2[n].y;
				dir(2, n) = vdir2[n].z;
			}
			raw::writed(dir, "./surface_curvature/dir2");
		}

		{
			Image<float32_t> dir(1, vcurvature1.size());
			for (size_t n = 0; n < vcurvature1.size(); n++)
			{
				dir(n) = vcurvature1[n];
			}
			raw::writed(dir, "./surface_curvature/curvature1");
		}

		{
			Image<float32_t> dir(1, vcurvature2.size());
			for (size_t n = 0; n < vcurvature2.size(); n++)
			{
				dir(n) = vcurvature2[n];
			}
			raw::writed(dir, "./surface_curvature/curvature2");
		}
#endif
	}


	namespace tests
	{
		void surfaceCurvature();
	}
}
