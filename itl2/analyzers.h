#pragma once

// TODO: This file is mostly work in progress. Most of the required code exists in previous pi versions, but is not ported here yet.

#include <string>
#include <vector>
#include <memory>
#include <algorithm>

#include "math/matrix2x2.h"
#include "math/matrix3x3.h"
#include "image.h"
#include "resultstable.h"
#include "sphere.h"

#include "marchingcubes.h"
#include "ellipsoid.h"

#include "convexhull.h"
//#include "sphere.h"


namespace itl2
{
	/**
	Class that specifies one analyzer for particle analysis.
	@param POINT Storage type for point coordinate class. Should be indexable with [0] and have size() member.
	*/
	template<class POINT, class pixel_t> class Analyzer
	{
	public:
		/**
		* Virtual destructor
		*/
		virtual ~Analyzer()
		{
		}

		/**
		Get title of this analysis.
		*/
		virtual std::vector<string> getTitles() const = 0;

		/**
		Analyze the particles and return the result.
		*/
		virtual std::vector<double> analyze(const std::vector<POINT>& points) const = 0;

		/**
		Gets name of this analyzers.
		*/
		virtual std::string name() const = 0;

		/**
		Gets description of the analysis this analyzer does.
		*/
		virtual std::string description() const = 0;
	};



	/**
	Set of Analyzers.
	*/
	template<class POINT, typename pixel_t> class AnalyzerSet : public std::vector<std::shared_ptr<Analyzer<POINT, pixel_t> > >
	{
	public:

		/**
		Get list of column headers for all the analyzers in this vector.
		*/
		virtual Headers headers() const
		{
			Headers headers;
			headers.reserve(this->size());
			for (size_t n = 0; n < this->size(); n++)
			{
				std::vector<string> titles = (*this)[n]->getTitles();
				for (size_t m = 0; m < titles.size(); m++)
					headers.push_back(titles[m]);
			}

			return headers;
		}

		/**
		Find in index of element name in headers array.
		*/
		size_t getColumnIndex(const string& name)
		{
			return headers().getColumnIndex(name);
		}

		using std::vector<std::shared_ptr<Analyzer<POINT, pixel_t> > >::push_back;

		/**
		Shorthand for adding analyzers to the set.
		*/
		virtual void push_back(Analyzer<POINT, pixel_t>* pAnalyzer)
		{
			std::vector<std::shared_ptr<Analyzer<POINT, pixel_t> > >::push_back(std::shared_ptr<Analyzer<POINT, pixel_t> >(pAnalyzer));
		}

		/**
		Perform analysis using all the analyzers in this vector.
		@param particlePoints Points on the particle.
		@param resultLine Results will be added to this array.
		*/
		virtual void analyze(const std::vector<POINT>& particlePoints, std::vector<double>& resultLine) const
		{
			size_t s = std::vector<std::shared_ptr<Analyzer<POINT, pixel_t> > >::size();

			if (particlePoints.size() <= 0)
				throw ITLException("Particle analyzer received no input points.");

			resultLine.reserve(s);
			for (size_t n = 0; n < s; n++)
			{
				std::vector<double> currentResults = (*this)[n]->analyze(particlePoints);
				for (size_t m = 0; m < currentResults.size(); m++)
					resultLine.push_back(currentResults[m]);
			}
		}
	};


	/**
	Contains particle analyzers.
	*/
	namespace analyzers
	{
		///**
		//Returns the coordinates of the first point on the particle.
		//*/
		//template<class POINT, typename pixel_t> class PointCoords2D : public Analyzer<POINT, pixel_t>
		//{
		//public:
		//	virtual std::vector<string> getTitles() const
		//	{
		//		std::vector<string> labels;
		//		labels.push_back("X [pixel]");
		//		labels.push_back("Y [pixel]");
		//		return labels;
		//	}

		//	virtual std::vector<double> analyze(const std::vector<POINT >& points) const
		//	{
		//		std::vector<double> results;
		//		results.push_back(points[0][0]);
		//		results.push_back(points[0][1]);
		//		return results;
		//	}
		//};

		///**
		//Returns the color of the point where the particle was found.
		//*/
		//template<class POINT, typename pixel_t> class Color : public Analyzer<POINT, pixel_t>
		//{
		//public:
		//    virtual std::vector<string> getTitles() const
		//    {
		//		std::vector<string> labels;
		//		labels.push_back("Color");
		//			return labels;
		//    }

		//    virtual std::vector<double> analyze(pixel_t color, const std::vector<POINT >& points) const
		//    {
		//		std::vector<double> results;
		//		results.push_back((double)color);
		//		return results;
		//    }
		//};

		/**
		Returns user-specified value.
		*/
		template<class POINT, typename pixel_t> class ArbitraryValue : public Analyzer<POINT, pixel_t>
		{
		private:
			const string title;
			double value;
		public:
			ArbitraryValue(const string& title, const double value) : title(title), value(value)
			{
			}

			virtual std::vector<string> getTitles() const override
			{
				std::vector<string> labels;
				labels.push_back(title);
				return labels;
			}

			virtual std::vector<double> analyze(const std::vector<POINT >& points) const override
			{
				std::vector<double> results;
				results.push_back(value);
				return results;
			}

			void setValue(double newValue)
			{
				value = newValue;
			}

			virtual std::string name() const override
			{
				return "Value";
			}

			virtual std::string description() const override
			{
				return "Shows user-specified value that is the same for all the particles. Output column name is 'Value'.";
			}
		};

		/**
		Returns the coordinates of the first point on the particle.
		*/
		template<class POINT, typename pixel_t> class Coordinates3D : public Analyzer<POINT, pixel_t>
		{
		public:
			virtual std::vector<std::string> getTitles() const override
			{
				std::vector<std::string> labels;
				labels.push_back("X [pixel]");
				labels.push_back("Y [pixel]");
				labels.push_back("Z [pixel]");
				return labels;
			}

			virtual std::vector<double> analyze(const std::vector<POINT>& points) const override
			{
				// This is needed so that threaded and non-threaded particle analysis versions
				// give exactly the same results.
				// TODO: The comparer template is not compatible with all POINT arguments.
				auto pi = std::min_element(points.begin(), points.end(), vecComparer<int32_t>);
				POINT p = *pi;

				std::vector<double> results;
				results.push_back((double)p[0]);
				results.push_back((double)p[1]);
				results.push_back((double)p[2]);
				return results;
			}

			virtual std::string name() const override
			{
				return "coordinates";
			}

			virtual std::string description() const override
			{
				return "Shows (zero-based) coordinates of a pixel that is guaranteed to be inside the particle. Output column names are 'X', 'Y', and 'Z'.";
			}
		};


		/**
		Returns the coordinates of the first point on the particle.
		*/
		template<class POINT, typename pixel_t> class Coordinates2D : public Analyzer<POINT, pixel_t>
		{
		public:
			virtual std::vector<std::string> getTitles() const override
			{
				std::vector<std::string> labels;
				labels.push_back("X [pixel]");
				labels.push_back("Y [pixel]");
				return labels;
			}

			virtual std::vector<double> analyze(const std::vector<POINT>& points) const override
			{
				// This is needed so that threaded and non-threaded particle analysis versions
				// give exactly the same results.
				// TODO: The comparer template is not compatible with all POINT arguments.
				auto pi = std::min_element(points.begin(), points.end(), vecComparer<int32_t>);
				POINT p = *pi;

				std::vector<double> results;
				results.push_back((double)p[0]);
				results.push_back((double)p[1]);
				return results;
			}

			virtual std::string name() const override
			{
				return "coordinates2d";
			}

			virtual std::string description() const override
			{
				return "Shows (zero-based) coordinates of a pixel that is guaranteed to be inside the particle. Output column names are 'X' and 'Y'.";
			}
		};

		/**
		Calculates the volume of the particle.
		*/
		template<class POINT, typename pixel_t> class Volume : public Analyzer<POINT, pixel_t>
		{
		public:
			virtual std::vector<std::string> getTitles() const override
			{
				std::vector<std::string> labels;
				labels.push_back("Volume [pixel]");
				return labels;
			}

			virtual std::vector<double> analyze(const std::vector<POINT >& points) const override
			{
				std::vector<double> results;
				results.push_back((double)points.size());
				return results;
			}

			virtual std::string name() const override
			{
				return "volume";
			}

			virtual std::string description() const override
			{
				return "Shows total volume of the particle in pixels, or total area of the particle in the 2D case. Outputs one column with name 'Volume'.";
			}
		};

		/**
		Calculates the surface area of a 3D object using Marching Cubes.
		*/
		template<class POINT, typename pixel_t> class SurfaceArea3D : public Analyzer<POINT, pixel_t>
		{
		public:

			virtual std::vector<std::string> getTitles() const override
			{
				std::vector<std::string> labels;
				labels.push_back("Surface area [pixel^2]");
				return labels;
			}

			virtual std::vector<double> analyze(const std::vector<POINT>& points) const override
			{
				// Calculate bounding box for the points.
				POINT min;
				POINT max;
				for (size_t n = 0; n < points.size(); n++)
				{
					POINT p = points[n];
					if (p.x < min.x)
						min.x = p.x;
					if (p.y < min.y)
						min.y = p.y;
					if (p.z < min.z)
						min.z = p.z;

					if (p.x > max.x)
						max.x = p.x;
					if (p.y > max.y)
						max.y = p.y;
					if (p.z > max.z)
						max.z = p.z;
				}

				coord_t w = (coord_t)std::ceil(max.x) - (coord_t)std::floor(min.x) + 2;
				coord_t h = (coord_t)std::ceil(max.y) - (coord_t)std::floor(min.y) + 2;
				coord_t d = (coord_t)std::ceil(max.z) - (coord_t)std::floor(min.z) + 2;

				// Plot the points into a temporary image
				Image<uint8_t> block(w, h, d);

				for (size_t n = 0; n < points.size(); n++)
				{
					POINT p = points[n] - min + POINT(1, 1, 1);
					block(p) = 255;
				}

				// Calculate the area of the particle using Marching Cubes
				double area = getMarchingCubesArea<uint8_t>(block, 128);

				std::vector<double> results;
				results.push_back(area);
				return results;
			}

			virtual std::string name() const override
			{
				return "surfacearea";
			}

			virtual std::string description() const override
			{
				return "Shows surface area of particles determined with the Marching Cubes method. Outputs one column 'Surface area'.";
			}
		};

		/**
		Calculates measures of the convex hull of the particle.
		*/
        template<class POINT, typename pixel_t> class ConvexHull2D : public Analyzer<POINT, pixel_t>
        {
		private:
			/**
			Calculates area of convex 2D polygon.
			*/
			static double convexPolyArea(const std::vector<POINT>& v)
			{
				double area = 0;
				for(size_t n = 1; n < v.size()-1; n++)
				{
					POINT A = v[0];
					POINT B = v[n];
					POINT C = v[n + 1];
					area += abs((A.x-C.x)*(B.y-A.y)-(A.x-B.x)*(C.y-A.y));
				}
				area *= 0.5;

				return area;
			}

			/**
			Calculates edge length of polygon.
			*/
			static double polyPerimeter(const std::vector<POINT>& v)
			{
				double p = 0;
				for(size_t n = 0; n < v.size() - 1; n++)
					p += (v[n + 1] - v[n]).norm();
				p += (v[0] - v[v.size() - 1]).norm();
				return p;
			}

        public:
			virtual std::vector<std::string> getTitles() const override
			{
				std::vector<std::string> labels;
				labels.push_back("Convex area [pixel^2]");
				labels.push_back("Convex perimeter [pixel]");
				return labels;
			}

			virtual std::vector<double> analyze(const std::vector<POINT>& points) const override
			{
				std::vector<double> results;

				// Calculate convex hull from edge points so that area of convex hull is always >= nonconvex area.
				std::vector<POINT > edgePoints;
				for(size_t n = 0; n < points.size(); n++)
				{
					POINT p = points[n];
					int ix = (int)p.x;
					int iy = (int)p.y;

					edgePoints.push_back(POINT(ix, iy, 0));
					edgePoints.push_back(POINT(ix+1, iy, 0));
					edgePoints.push_back(POINT(ix, iy+1, 0));
					edgePoints.push_back(POINT(ix+1, iy+1, 0));
				}

				std::vector<POINT> hull;
				convexHull2D(edgePoints, hull);
				double area = convexPolyArea(hull);
				double perimeter = polyPerimeter(hull);

				results.push_back(area);
				results.push_back(perimeter);
				return results;
            }

			virtual std::string name() const override
			{
				return "convexhull2d";
			}

			virtual std::string description() const override
			{
				return "Shows area and perimeter of the convex hull of the particle. Outputs columns 'Convex area' and 'Convex perimeter'.";
			}
        };


		/**
		Tests whether the particle touches image edge.
		*/
		template<class POINT, typename pixel_t> class IsOnEdge : public Analyzer<POINT, pixel_t>
		{
		private:
			Vec3c dimensions;
		public:
			/**
			Constructor
			@param dimensions Dimensions of the image.
			*/
			IsOnEdge(const Vec3c& dimensions)
			{
				this->dimensions = dimensions;
			}

			virtual std::vector<std::string> getTitles() const override
			{
				std::vector<std::string> titles;
				titles.push_back("Is on edge [0/1]");
				return titles;
			}

			/**
			Tests if any of the given points is on the edge of image whose dimensions are also given.
			*/
			static bool isOnEdge(const std::vector<POINT>& points, const Vec3c& dimensions)
			{
				// Traverse all the points in the particle
				for (size_t n = 0; n < points.size(); n++)
				{
					// Test each dimension of the point for image edge status
					const POINT& p = points[n];
					for (size_t dimension = 0; dimension < p.size(); dimension++)
					{
						if (p[dimension] <= 0 || p[dimension] >= dimensions[dimension] - 1)
						{
							return true;
						}
					}
				}

				return false;
			}

			virtual std::vector<double> analyze(const std::vector<POINT>& points) const override
			{
				std::vector<double> results;
				results.push_back(0.0);

				if (isOnEdge(points, dimensions))
					results[0] = 1.0;

				return results;
			}

			virtual std::string name() const override
			{
				return "isonedge";
			}

			virtual std::string description() const override
			{
				return "Tests whether the particle touches image edge. If it does, returns 1 in a column 'Is on edge'; if it doesn't, returns zero.";
			}
		};

		/**
		Calculates principal components.
		*/
		template<class POINT, typename pixel_t> class PCA3D : public Analyzer<POINT, pixel_t>
		{
			virtual std::string name() const override
			{
				return "pca";
			}

			virtual std::string description() const override
			{
				return "Calculates orientation of the particle using principal component analysis. Outputs:\n"
					"Centroid of the particle in columns 'CX', 'CY', and 'CZ'.\n"
					"Meridional eccentricity in column 'e (meridional)'.\n"
					"Standard deviations of the projections of the particle points to the principal axes in columns 'l1', 'l2', and 'l3'.\n"
					"Orientations of the three principal axes of the particle, in spherical coordinates, in columns\n"
					"'phiN' and 'thetaN', where N is in 1...3, phiN is the azimuthal angle and thetaN is the polar angle.\n"
					"Maximum distance from the centroid to the edge of the particle in colum 'rmax'.\n"
					"Maximum diameter of the projection of the particle on each principal component in columns 'd1', 'd2', and 'd3'.\n"
					"Scaling factor for a bounding ellipsoid in column 'bounding scale'. An ellipsoid that bounds the\n"
					"particle and whose semi-axes correspond to the principal components has semi-axis lengths\n"
					"$b lN$, where $b$ is the bounding scale and $lN$ are the lengths of the principal axes.";
			}

			virtual std::vector<string> getTitles() const override
			{
				std::vector<string> titles;
				titles.push_back("CX [pixel]");
				titles.push_back("CY [pixel]");
				titles.push_back("CZ [pixel]");
				titles.push_back("e (meridional)");
				titles.push_back("l1 [pixel]");
				titles.push_back("l2 [pixel]");
				titles.push_back("l3 [pixel]");
				titles.push_back("phi1 [rad]");
				titles.push_back("theta1 [rad]");
				titles.push_back("phi2 [rad]");
				titles.push_back("theta2 [rad]");
				titles.push_back("phi3 [rad]");
				titles.push_back("theta3 [rad]");
				titles.push_back("rmax [pixel]");	// Overall maximum radius from center point
				titles.push_back("d1 [pixel]");     // Maximum diameter in direction of first principal component.
				titles.push_back("d2 [pixel]");     // Maximum diameter in direction of second principal component.
				titles.push_back("d3 [pixel]");     // Maximum diameter in direction of third principal component.
				titles.push_back("bounding scale [1]");
				return titles;
			}

			virtual std::vector<double> analyze(const std::vector<POINT>& points) const override
			{
				std::vector<double> results;


				Vec3d centroid = mean<POINT, Vec3d, double>(points);

				Matrix3x3d CI;
				for (const auto& p : points)
				{
					Vec3d d = Vec3d(p) - centroid;
					CI += Matrix3x3d::outer(d, d);
				}
				CI /= (double)points.size();

				Vec3d t1, t2, t3;
				double lambda1, lambda2, lambda3;
				CI.eigsym(t1, t2, t3, lambda1, lambda2, lambda3);


				if(lambda1 < lambda2 || lambda2 < lambda3 || lambda1 < lambda3)
					throw ITLException("Eigenvalues are unsorted.");

				double r1, phi1, theta1;
				double r2, phi2, theta2;
				double r3, phi3, theta3;

				toSpherical(t1.x, t1.y, t1.z, r1, phi1, theta1);
				toSpherical(t2.x, t2.y, t2.z, r2, phi2, theta2);
				toSpherical(t3.x, t3.y, t3.z, r3, phi3, theta3);
				
				double e = 0.0;
				if(abs(lambda1) > 1e-6)
					e = sqrt(1.0 - (lambda3 * lambda3) / (lambda1 * lambda1));

				double l1 = sqrt(lambda1);
				double l2 = sqrt(lambda2);
				double l3 = sqrt(lambda3);

				results.push_back(centroid.x);
				results.push_back(centroid.y);
				results.push_back(centroid.z);

				results.push_back(e);

				results.push_back(l1);
				results.push_back(l2);
				results.push_back(l3);

				results.push_back(phi1);
				results.push_back(theta1);
				
				results.push_back(phi2);
				results.push_back(theta2);
				
				results.push_back(phi3);
				results.push_back(theta3);


				// Calculate maximum radius, projections to principal axes, etc.
				double maxr = 0;
				t1.normalize();
				t2.normalize();
				t3.normalize();
				double maxd1 = 0;
				double mind1 = 0;
				double maxd2 = 0;
				double mind2 = 0;
				double maxd3 = 0;
				double mind3 = 0;
				double boundingScale = 0;
				for(size_t n = 0; n < points.size(); n++)
				{
					Vec3d p = Vec3d(points[n]) - centroid;
					double r = p.norm();
					maxr = std::max(maxr, r);

					r = t1.dot(p);
					if(r > maxd1)
						maxd1 = r;
                    if(r < mind1)
                        mind1 = r;

					r = t2.dot(p);
					if(r > maxd2)
						maxd2 = r;
                    if(r < mind2)
                        mind2 = r;

					r = t3.dot(p);
					if(r > maxd3)
						maxd3 = r;
                    if(r < mind3)
                        mind3 = r;


					double f = getEllipsoidFunctionValue(Vec3d(points[n]),
								centroid,
								l1, l2, l3,
								phi1, theta1,
								phi2, theta2,
								phi3, theta3);
					double scale = sqrt(f);
					if(!std::isnan(scale) && !std::isinf(scale))
						boundingScale = std::max(boundingScale, scale);
				}
				results.push_back(maxr);
				results.push_back(maxd1 - mind1);
				results.push_back(maxd2 - mind2);
				results.push_back(maxd3 - mind3);
				results.push_back(boundingScale);

				return results;
			}
		};

		/**
		Calculates principal components.
		*/
		template<class POINT, typename pixel_t> class PCA2D : public Analyzer<POINT, pixel_t>
		{
			virtual std::string name() const override
			{
				return "pca2d";
			}

			virtual std::string description() const override
			{
				return "Calculates orientation of the particle using principal component analysis. Outputs:\n"
					"Centroid of the particle in columns 'CX', and 'CY'.\n"
					"Eccentricity in column 'e'.\n"
					"Standard deviations of the projections of the particle points to the principal axes in columns 'l1', and 'l2'.\n"
					"Orientation of the first principal axis of the particle in column 'alpha'.\n"
					"The alpha value is angle in radians between positive $x$-axis and particle orientation.\n"
					"Maximum distance from the centroid to the edge of the particle in colum 'rmax'.\n"
					"Maximum diameter of the projection of the particle on each principal component in columns 'd1', and 'd2'.\n"
					"Scaling factor for a bounding ellipse in column 'bounding scale'. An ellipse that bounds the\n"
					"particle and whose semi-axes correspond to the principal components, has semi-axis lengths\n"
					"$b lN$, where $b$ is the bounding scale and $lN$ are the lengths of the principal axes.";
			}

			virtual std::vector<string> getTitles() const override
			{
				std::vector<string> titles;
				titles.push_back("CX [pixel]");
				titles.push_back("CY [pixel]");
				titles.push_back("e");
				titles.push_back("l1 [pixel]");
				titles.push_back("l2 [pixel]");
				titles.push_back("alpha [rad]");
				titles.push_back("rmax [pixel]");	// Overall maximum radius from center point
				titles.push_back("d1 [pixel]");     // Maximum diameter in direction of first principal component.
				titles.push_back("d2 [pixel]");     // Maximum diameter in direction of second principal component.
				titles.push_back("bounding scale [1]");
				return titles;
			}

			virtual std::vector<double> analyze(const std::vector<POINT>& points) const override
			{
				std::vector<double> results;


				Vec3d centroid3 = mean<POINT, Vec3d, double>(points);
				centroid3.z = 0;
				Vec2d centroid(centroid3.x, centroid3.y);

				Matrix2x2d CI;
				for (const auto& p : points)
				{
					Vec2d d = Vec2d(p.x, p.y) - centroid;
					CI += Matrix2x2d::outer(d, d);
				}
				CI /= (double)points.size();

				Vec2d t1, t2;
				double lambda1, lambda2;
				CI.eigsym(t1, t2, lambda1, lambda2);


				if (lambda1 < lambda2)
					throw ITLException("Eigenvalues are unsorted.");

				double r1, phi1;

				toPolar(t1.x, t1.y, r1, phi1);

				double e = 0.0;
				if (abs(lambda1) > 1e-6)
					e = sqrt(1.0 - (lambda2 * lambda2) / (lambda1 * lambda1));

				double l1 = sqrt(lambda1);
				double l2 = sqrt(lambda2);

				results.push_back(centroid.x);
				results.push_back(centroid.y);

				results.push_back(e);

				results.push_back(l1);
				results.push_back(l2);

				results.push_back(phi1);


				// Calculate maximum radius, projections to principal axes, etc.
				double maxr = 0;
				t1.normalize();
				t2.normalize();
				double maxd1 = 0;
				double mind1 = 0;
				double maxd2 = 0;
				double mind2 = 0;
				double boundingScale = 0;
				for (size_t n = 0; n < points.size(); n++)
				{
					Vec2d p = Vec2d(points[n].x, points[n].y) - centroid;
					double r = p.norm();
					maxr = std::max(maxr, r);

					r = t1.dot(p);
					if (r > maxd1)
						maxd1 = r;
					if (r < mind1)
						mind1 = r;

					r = t2.dot(p);
					if (r > maxd2)
						maxd2 = r;
					if (r < mind2)
						mind2 = r;


					double f = getEllipseFunctionValue(Vec2d(points[n].x, points[n].y),
						centroid,
						l1, l2,
						phi1);
					double scale = sqrt(f);
					if (!std::isnan(scale) && !std::isinf(scale))
						boundingScale = std::max(boundingScale, scale);
				}
				results.push_back(maxr);
				results.push_back(maxd1 - mind1);
				results.push_back(maxd2 - mind2);
				results.push_back(boundingScale);

				return results;
			}
		};

		///**
		//Returns measurements of histogram (from another image) at the points of the particle.
		//*/
  //      template<class POINT, typename pixel_t, typename PTYPE2> class Histogram : public Analyzer<POINT, pixel_t>
  //      {
  //      private:
  //          Image<PTYPE2>& anotherImage;

		//	std::vector<std::vector<float> > pixelValues;

  //      public:

		//	const std::vector<std::vector<float> >& getExtraResults() const
		//	{
		//		return pixelValues;
		//	}

  //          Histogram(Image<PTYPE2>& anotherImage) : anotherImage(anotherImage)
  //          {
  //          }

  //          virtual std::vector<string> getTitles() const
  //          {
		//		std::vector<string> labels;
		//		labels.push_back("Mean");
		//		labels.push_back("Stdev");
		//		labels.push_back("Min");
		//		labels.push_back("Max");
  //              return labels;
  //          }

  //          virtual std::vector<double> analyze(const std::vector<POINT >& points)
  //          {
		//		pixelValues.push_back(std::vector<float>());
		//		std::vector<float>& vec = pixelValues[pixelValues.size()-1];
		//		vec.reserve(points.size());

		//		// Calculate statistics
		//		double mean = 0;
		//		double stdev = 0;
		//		double min = numeric_limits<double>::infinity();
		//		double max = -numeric_limits<double>::infinity();

		//		
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			double pix = anotherImage.getPixel(points[n]);
		//			vec.push_back((float)pix);
		//			mean += pix;
		//			if(pix < min)
		//				min = pix;
		//			if(pix > max)
		//				max = pix;
		//		}

		//		mean /= (double)points.size();

		//		// Calculate standard deviation
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			double pix = anotherImage.getPixel(points[n]);
		//			stdev += (pix - mean) * (pix - mean);
		//		}

		//		stdev = sqrt(stdev/(double)points.size());

		//		std::vector<double> results;
		//		results.push_back(mean);
		//		results.push_back(stdev);
		//		results.push_back(min);
		//		results.push_back(max);
  //              return results;
  //          }
  //      };


		/**
		Calculates axis-aligned bounding box of the particle.
		*/
		template<class POINT, typename pixel_t> class BoundingBox3D : public Analyzer<POINT, pixel_t>
		{
		private:


		public:
			virtual std::vector<std::string> getTitles() const override
			{
				std::vector<std::string> labels;
				labels.push_back("minx [pix]");
				labels.push_back("maxx [pix]");
				labels.push_back("miny [pix]");
				labels.push_back("maxy [pix]");
				labels.push_back("minz [pix]");
				labels.push_back("maxz [pix]");
				return labels;
			}

			virtual std::vector<double> analyze(const std::vector<POINT >& points) const override
			{
				// Calculate bounds
				POINT min = points[0];
				POINT max = points[0];
				for (size_t n = 1; n < points.size(); n++)
				{
					if (points[n].x < min.x)
						min.x = points[n].x;
					if (points[n].y < min.y)
						min.y = points[n].y;
					if (points[n].z < min.z)
						min.z = points[n].z;
					if (points[n].x > max.x)
						max.x = points[n].x;
					if (points[n].y > max.y)
						max.y = points[n].y;
					if (points[n].z > max.z)
						max.z = points[n].z;
				}

				std::vector<double> results;
				results.push_back((double)min.x);
				results.push_back((double)max.x);
				results.push_back((double)min.y);
				results.push_back((double)max.y);
				results.push_back((double)min.z);
				results.push_back((double)max.z);
				return results;
			}

			virtual std::string name() const override
			{
				return "bounds";
			}

			virtual std::string description() const override
			{
				return "Calculates axis-aligned bounding box of the particle. Returns minimum and maximum coordinates of particle points in all coordinate dimensions. Outputs columns 'minX' and 'maxX' for each dimension, where X is either x, y, or z.";
			}
		};


		/**
		Calculates axis-aligned bounding box of the particle.
		*/
		template<class POINT, typename pixel_t> class BoundingBox2D : public Analyzer<POINT, pixel_t>
		{
		private:


		public:
			virtual std::vector<std::string> getTitles() const override
			{
				std::vector<std::string> labels;
				labels.push_back("minx [pix]");
				labels.push_back("maxx [pix]");
				labels.push_back("miny [pix]");
				labels.push_back("maxy [pix]");
				return labels;
			}

			virtual std::vector<double> analyze(const std::vector<POINT >& points) const override
			{
				// Calculate bounds
				POINT min = points[0];
				POINT max = points[0];
				for (size_t n = 1; n < points.size(); n++)
				{
					if (points[n].x < min.x)
						min.x = points[n].x;
					if (points[n].y < min.y)
						min.y = points[n].y;
					if (points[n].x > max.x)
						max.x = points[n].x;
					if (points[n].y > max.y)
						max.y = points[n].y;
				}

				std::vector<double> results;
				results.push_back((double)min.x);
				results.push_back((double)max.x);
				results.push_back((double)min.y);
				results.push_back((double)max.y);
				return results;
			}

			virtual std::string name() const override
			{
				return "bounds2d";
			}

			virtual std::string description() const override
			{
				return "Calculates axis-aligned bounding box of the particle. Returns minimum and maximum coordinates of particle points in all coordinate dimensions. Outputs columns 'minX' and 'maxX' for each dimension, where X is either x, or y.";
			}
		};


		/**
		Calculates minimum bounding sphere of the particle
		*/
		template<class POINT, typename pixel_t> class BoundingSphere : public Analyzer<POINT, pixel_t>
		{
		private:


		public:
			virtual std::vector<std::string> getTitles() const override
			{
				std::vector<std::string> labels;
				labels.push_back("bounding sphere X [pixel]");
				labels.push_back("bounding sphere Y [pixel]");
				labels.push_back("bounding sphere Z [pixel]");
				labels.push_back("bounding sphere radius [pixel]");
				return labels;
			}

			virtual std::vector<double> analyze(const std::vector<POINT >& points) const override
			{
				Sphere<double> bs = Sphere<double>::miniball(points);

				std::vector<double> results;
				results.push_back(bs.center.x);
				results.push_back(bs.center.y);
				results.push_back(bs.center.z);
				results.push_back(bs.radius);
				return results;
			}

			virtual std::string name() const override
			{
				return "boundingsphere";
			}

			virtual std::string description() const override
			{
				return "Shows position and radius of the smallest possible sphere that contains all the points in the particle. Outputs columns 'bounding sphere X', 'bounding sphere Y', 'bounding sphere Z', and 'bounding sphere radius'. Uses the MiniBall algorithm from Welzl, E. - Smallest enclosing disks (balls and ellipsoids).";
			}
		};
	}

	/**
	Creates list of all particle analyzers for 3D particles.
	*/
	template<typename pixel_t> AnalyzerSet<Vec3sc, pixel_t> allAnalyzers(const Vec3c& dimensions)
	{
		AnalyzerSet<Vec3sc, pixel_t> analyzers;
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::Coordinates3D<Vec3sc, pixel_t>()));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::IsOnEdge<Vec3sc, pixel_t>(dimensions)));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::Volume<Vec3sc, pixel_t>()));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::PCA3D<Vec3sc, pixel_t>()));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::SurfaceArea3D<Vec3sc, pixel_t>()));
		//analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::Color<Vec3sc, pixel_t>()));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::BoundingBox3D<Vec3sc, pixel_t>()));
		analyzers.push_back(std::shared_ptr<Analyzer<Vec3sc, pixel_t> >(new analyzers::BoundingSphere<Vec3sc, pixel_t>()));

		return analyzers;
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

	/**
	Creates list of all particle analyzers.
	*/
	template<typename pixel_t> AnalyzerSet<Vec3sc, pixel_t> allAnalyzers(const Image<pixel_t>& img)
	{
		return allAnalyzers<pixel_t>(img.dimensions());
	}

	namespace internals
	{
		inline int isnotalnum(int ch)
		{
			return !isalnum(ch);
		}

		/**
		Get analyzer from list by its name.
		*/
		template<typename pixel_t> std::shared_ptr<Analyzer<Vec3sc, pixel_t> > getAnalyzer(AnalyzerSet<Vec3sc, pixel_t>& all, const string& name, bool throwIfNotFound)
		{
			string loname = name;
			toLower(loname);

			for (size_t n = 0; n < all.size(); n++)
			{
				string cname = all[n]->name();
				toLower(cname);
				if (cname == loname)
					return all[n];
			}

			if(throwIfNotFound)
				throw ITLException(string("Invalid analyzer name: ") + name);
			return nullptr;
		}
	}

	/**
	Converts list of analyzer names to analyzer set.
	*/
	template<typename pixel_t> AnalyzerSet<Vec3sc, pixel_t> createAnalyzers(string names, const Vec3c& dimensions)
	{
		// Replace all possible delimiters by space
		replace_if(names.begin(), names.end(), internals::isnotalnum, ' ');

		std::vector<string> parts = split(names, false, ' ');

		AnalyzerSet<Vec3sc, pixel_t> all = allAnalyzers<pixel_t>(dimensions);
		AnalyzerSet<Vec3sc, pixel_t> all2d = allCrossSectionAnalyzers<pixel_t>();
		AnalyzerSet<Vec3sc, pixel_t> result;

		for (size_t n = 0; n < parts.size(); n++)
		{
			string name = parts[n];
			std::shared_ptr<Analyzer<Vec3sc, pixel_t> > analyzer = internals::getAnalyzer<pixel_t>(all, name, false);
			if (analyzer)
			{
				// 3D analyzer found
				result.push_back(analyzer);
			}
			else
			{
				// 3D analyzer not found, try 2D
				analyzer = internals::getAnalyzer<pixel_t>(all2d, name, false);
				
				if (analyzer)
				{
					// 2D analyzer found.
					if (dimensions.z <= 1)
						result.push_back(analyzer);
					else
						throw ITLException(string("Analyzer ") + name + " is for 2-dimensional images, but the image to be analyzed is not 2-dimensional.");
				}
				else
				{
					// 2D analyzer not found
					ITLException(string("Invalid analyzer name: ") + name);
				}
				
			}
		}

		return result;
	}


	/**
	Converts list of analyzer names to analyzer set.
	*/
	template<typename pixel_t> AnalyzerSet<Vec3sc, pixel_t> createCrossSectionAnalyzers(string names)
	{
		// Replace all possible delimiters by space
		replace_if(names.begin(), names.end(), internals::isnotalnum, ' ');

		std::vector<string> parts = split(names, false, ' ');

		AnalyzerSet<Vec3sc, pixel_t> all = allCrossSectionAnalyzers<pixel_t>();
		AnalyzerSet<Vec3sc, pixel_t> result;

		for (size_t n = 0; n < parts.size(); n++)
		{
			std::shared_ptr<Analyzer<Vec3sc, pixel_t> > tmp = internals::getAnalyzer<pixel_t>(all, parts[n], true);
			result.push_back(tmp);
		}

		return result;
	}

}
