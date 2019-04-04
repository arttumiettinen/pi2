#pragma once

#include <string>
#include <vector>
#include <memory>
#include <algorithm>

#include "math/vec3.h"
#include "image.h"
#include "resultstable.h"
#include "sphere.h"

//#include "math/matrix3x3.h"
//#include "convexhull.h"
//#include "marchingcubes.h"
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
		//	virtual vector<string> getTitles() const
		//	{
		//		vector<string> labels;
		//		labels.push_back("X [pixel]");
		//		labels.push_back("Y [pixel]");
		//		return labels;
		//	}

		//	virtual vector<double> analyze(const vector<POINT >& points) const
		//	{
		//		vector<double> results;
		//		results.push_back(points[0][0]);
		//		results.push_back(points[0][1]);
		//		return results;
		//	}
		//};

		//      /**
			  //Returns the color of the point where the particle was found.
			  //*/
		//      template<class POINT, typename pixel_t> class Color : public Analyzer<POINT, pixel_t>
		//      {
		//      public:
		//          virtual vector<string> getTitles() const
		//          {
			  //		vector<string> labels;
			  //		labels.push_back("Color");
		//              return labels;
		//          }

		//          virtual vector<double> analyze(const vector<POINT >& points) const
		//          {
			  //		vector<double> results;
			  //		results.push_back((double)pixelValue);
		//              return results;
		//          }
		//      };

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

			virtual vector<string> getTitles() const
			{
				vector<string> labels;
				labels.push_back(title);
				return labels;
			}

			virtual vector<double> analyze(const vector<POINT >& points) const
			{
				vector<double> results;
				results.push_back(value);
				return results;
			}

			void setValue(double newValue)
			{
				value = newValue;
			}

			virtual std::string name() const
			{
				return "value";
			}

			virtual std::string description() const
			{
				return "Shows user-specified value that is the same for all the particles. Output column name is 'value'.";
			}
		};

		/**
		Returns the coordinates of the first point on the particle.
		*/
		template<class POINT, typename pixel_t> class Coordinates : public Analyzer<POINT, pixel_t>
		{
		public:
			virtual vector<string> getTitles() const
			{
				vector<string> labels;
				labels.push_back("X [pixel]");
				labels.push_back("Y [pixel]");
				labels.push_back("Z [pixel]");
				return labels;
			}

			virtual vector<double> analyze(const vector<POINT >& points) const
			{
				// This is needed so that threaded and non-threaded particle analysis versions
				// give exactly the same results.
				// TODO: The comparer template is not compatible with all POINT arguments.
				auto pi = std::min_element(points.begin(), points.end(), math::vecComparer<int32_t>);
				POINT p = *pi;

				vector<double> results;
				results.push_back((double)p[0]);
				results.push_back((double)p[1]);
				results.push_back((double)p[2]);
				return results;
			}

			virtual std::string name() const
			{
				return "coordinates";
			}

			virtual std::string description() const
			{
				return "Shows (zero-based) coordinates of a pixel that is guaranteed to be inside the particle. Output column names are X, Y, and Z.";
			}
		};

		/**
		Calculates the volume of the particle.
		*/
		template<class POINT, typename pixel_t> class Volume : public Analyzer<POINT, pixel_t>
		{
		public:
			virtual vector<string> getTitles() const
			{
				vector<string> labels;
				labels.push_back("Volume [pixel]");
				return labels;
			}

			virtual vector<double> analyze(const vector<POINT >& points) const
			{
				vector<double> results;
				results.push_back((double)points.size());
				return results;
			}

			virtual std::string name() const
			{
				return "volume";
			}

			virtual std::string description() const
			{
				return "Shows total volume of the particle in pixels. Outputs one column with name 'volume'.";
			}
		};

		///**
		//Calculates the surface area of a 3D object using Marching Cubes.
		//*/
		//template<class POINT, typename pixel_t> class SurfaceArea3D : public Analyzer<POINT, pixel_t>
		//{
		//public:
  //          virtual vector<string> getTitles() const
  //          {
		//		vector<string> labels;
		//		labels.push_back("Surface area [pixel^2]");
  //              return labels;
  //          }

  //          virtual vector<double> analyze(const vector<POINT >& points)
  //          {
		//		// Calculate bounding box for the points.
		//		POINT min;
		//		POINT max;
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			POINT p = points[n];
		//			if(p.x < min.x)
		//				min.x = p.x;
		//			if(p.y < min.y)
		//				min.y = p.y;
		//			if(p.z < min.z)
		//				min.z = p.z;

		//			if(p.x > max.x)
		//				max.x = p.x;
		//			if(p.y > max.y)
		//				max.y = p.y;
		//			if(p.z > max.z)
		//				max.z = p.z;
		//		}

		//		coord_t w = (coord_t)round(max.x) - (coord_t)round(min.x) + 2;
		//		coord_t h = (coord_t)round(max.y) - (coord_t)round(min.y) + 2;
		//		coord_t d = (coord_t)round(max.z) - (coord_t)round(min.z) + 2;

		//		// Plot the points into temporary image
		//		Image<uint8_t> block(w, h, d);

		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			POINT p = points[n];
		//			p.x++;
		//			p.y++;
		//			p.z++;
		//			block.setPixel(p, 255);
		//		}

		//		// Calculate area of the particle using Marching Cubes
		//		double area = getMarchingCubesArea<uint8_t>(block, 128);

		//		vector<double> results;
		//		results.push_back(area);
		//		return results;
  //          }
		//};

		///**
		//Calculates measures of the convex hull of the particle.
		//*/
  //      template<class POINT, typename pixel_t> class ConvexHull2D : public Analyzer<POINT, pixel_t>
  //      {
		//private:
		//	/**
		//	Calculates area of convex 2D polygon.
		//	*/
		//	double convexPolyArea(const vector<POINT>& v)
		//	{
		//		double area = 0;
		//		for(size_t n = 1; n < v.size()-1; n++)
		//		{
		//			POINT A = v[0];
		//			POINT B = v[n];
		//			POINT C = v[n + 1];
		//			area += abs((A.x-C.x)*(B.y-A.y)-(A.x-B.x)*(C.y-A.y));
		//		}
		//		area *= 0.5;

		//		return area;
		//	}

		//	/**
		//	Calculates edge length of polygon.
		//	*/
		//	double polyPerimeter(const vector<POINT>& v)
		//	{
		//		double p = 0;
		//		for(size_t n = 0; n < v.size() - 1; n++)
		//			p += (v[n + 1] - v[n]).norm();
		//		p += (v[0] - v[v.size() - 1]).norm();
		//		return p;
		//	}

  //      public:
  //          virtual vector<string> getTitles() const
  //          {
		//		vector<string> labels;
		//		labels.push_back("Convex area [pixel^2]");
		//		labels.push_back("Convex perimeter [pixel]");
  //              return labels;
  //          }

  //          virtual vector<double> analyze(const vector<POINT >& points)
  //          {
		//		vector<double> results;

		//		// Calculate convex hull from edge points so that area of convex hull is always >= nonconvex area.
		//		vector<POINT > edgePoints;
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			POINT p = points[n];
		//			int ix = (int)p.x;
		//			int iy = (int)p.y;

		//			edgePoints.push_back(POINT(ix, iy));
		//			edgePoints.push_back(POINT(ix+1, iy));
		//			edgePoints.push_back(POINT(ix, iy+1));
		//			edgePoints.push_back(POINT(ix+1, iy+1));
		//		}

		//		vector<POINT> hull;
		//		convexHull2D(edgePoints, hull);
		//		double area = convexPolyArea(hull);
		//		double perimeter = polyPerimeter(hull);

		//		results.push_back(area);
		//		results.push_back(perimeter);
		//		return results;
  //          }
  //      };


		/**
		Tests whether the particle touches image edge.
		*/
		template<class POINT, typename pixel_t> class IsOnEdge : public Analyzer<POINT, pixel_t>
		{
		private:
			math::Vec3c dimensions;
		public:
			/**
			Constructor
			@param dimensions Dimensions of the image.
			*/
			IsOnEdge(const math::Vec3c& dimensions)
			{
				this->dimensions = dimensions;
			}

			virtual vector<string> getTitles() const
			{
				vector<string> titles;
				titles.push_back("Is on edge [0/1]");
				return titles;
			}

			/**
			Tests if any of the given points is on the edge of image whose dimensions are also given.
			*/
			static bool isOnEdge(const vector<POINT>& points, const math::Vec3c& dimensions)
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

			virtual vector<double> analyze(const vector<POINT>& points) const
			{
				vector<double> results;
				results.push_back(0.0);

				if (isOnEdge(points, dimensions))
					results[0] = 1.0;

				return results;
			}

			virtual std::string name() const
			{
				return "isonedge";
			}

			virtual std::string description() const
			{
				return "Tests whether the particle touches image edge. If it does, returns 1.0 in a column 'isonedge'; if it doesn't, returns zero.";
			}
		};

		///**
		//Calculates principal components.
		//*/
		//template<class POINT, typename pixel_t> class PCA3D : public Analyzer<POINT, pixel_t>
		//{
		//	virtual vector<string> getTitles() const
		//	{
		//		vector<string> titles;
		//		titles.push_back("CX [pixel]");
		//		titles.push_back("CY [pixel]");
		//		titles.push_back("CZ [pixel]");
		//		titles.push_back("e (meridional)");
		//		titles.push_back("l1 [pixel]");
		//		titles.push_back("l2 [pixel]");
		//		titles.push_back("l3 [pixel]");
		//		titles.push_back("phi1 [rad]");
		//		titles.push_back("theta1 [rad]");
		//		titles.push_back("phi2 [rad]");
		//		titles.push_back("theta2 [rad]");
		//		titles.push_back("phi3 [rad]");
		//		titles.push_back("theta3 [rad]");
		//		titles.push_back("rmax [pixel]");	// Overall maximum radius from center point
		//		titles.push_back("d1 [pixel]");     // Maximum diameter in direction of first principal component.
		//		titles.push_back("d2 [pixel]");     // Maximum diameter in direction of second principal component.
		//		titles.push_back("d3 [pixel]");     // Maximum diameter in direction of third principal component.
		//		titles.push_back("bounding scale"); // Semi-axis of bounding ellipsoid are calculated as (bounding scale)*l1, (bounding scale)*l2, etc.
		//		return titles;
		//	}

		//	virtual vector<double> analyze(const vector<POINT >& points)
		//	{
		//		vector<double> results;

		//		// Calculate centroid
		//		double mx = 0.0;
		//		double my = 0.0;
		//		double mz = 0.0;

		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			mx += points[n][0];
		//			my += points[n][1];
		//			mz += points[n][2];
		//		}
		//		mx /= points.size();
		//		my /= points.size();
		//		mz /= points.size();

		//		results.push_back(mx);
		//		results.push_back(my);
		//		results.push_back(mz);

		//		// Calculate covariance matrix
		//		double sumxx = 0.0;
		//		double sumyy = 0.0;
		//		double sumzz = 0.0;
		//		double sumxy = 0.0;
		//		double sumxz = 0.0;
		//		double sumyz = 0.0;

		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			double currX = points[n][0] - mx;
		//			double currY = points[n][1] - my;
		//			double currZ = points[n][2] - mz;

		//			sumxx += currX * currX;
		//			sumyy += currY * currY;
		//			sumzz += currZ * currZ;

		//			sumxy += currX * currY;
		//			sumxz += currX * currZ;
		//			sumyz += currY * currZ;
		//		}

		//		double K = (double)points.size();
		//		sumxx /= K;
		//		sumyy /= K;
		//		sumzz /= K;
		//		sumxy /= K;
		//		sumxz /= K;
		//		sumyz /= K;

		//		//Matrix3d S;
		//		//SelfAdjointEigenSolver<Matrix3d> solver;

		//		//S(0, 0) = sumxx;
		//		//S(0, 1) = sumxy;
		//		//S(0, 2) = sumxz;

		//		//S(1, 0) = sumxy;
		//		//S(1, 1) = sumyy;
		//		//S(1, 2) = sumyz;

		//		//S(2, 0) = sumxz;
		//		//S(2, 1) = sumyz;
		//		//S(2, 2) = sumzz;

		//		//solver.compute(S);

		//		//double lambda1 = 0;
		//		//double lambda2 = 0;
		//		//double lambda3 = 0;

		//		////Matrix3d V;

		//		//if(solver.info() == Success)
		//		//{
		//		//	lambda1 = solver.eigenvalues()[2];
		//		//	lambda2 = solver.eigenvalues()[1];
		//		//	lambda3 = solver.eigenvalues()[0];

		//		//	//V = solver.eigenvectors();
		//		//}
		//		//else
		//		//{
		//		//	//V.setZero();
		//		//}

		//		double S[3][3];

		//		S[0][0] = sumxx;
		//		S[0][1] = sumxy;
		//		S[0][2] = sumxz;

		//		S[1][0] = sumxy;
		//		S[1][1] = sumyy;
		//		S[1][2] = sumyz;

		//		S[2][0] = sumxz;
		//		S[2][1] = sumyz;
		//		S[2][2] = sumzz;

		//		// Symmetric matrix A => eigenvectors in columns of V, corresponding eigenvalues in d.
		//		// See also analysis of structure tensor - approx the same code is there, too.
		//		double V[3][3];
		//		double d[3];
		//		eigen_decomposition(S, V, d);
		//		double lambda1 = d[2];
		//		double lambda2 = d[1];
		//		double lambda3 = d[0];

		//		if(lambda1 < lambda2 || lambda2 < lambda3 || lambda1 < lambda3)
		//			throw ITLException("Eigenvalues are unsorted.");

		//		double r1, phi1, theta1;
		//		double r2, phi2, theta2;
		//		double r3, phi3, theta3;

		//		double V1mult = V[0][0] < 0 ? -1 : 1;
		//		double V2mult = V[0][1] < 0 ? -1 : 1;
		//		double V3mult = V[0][2] < 0 ? -1 : 1;
		//		math::tospherical(V1mult * V[0][0], V1mult * V[1][0], V1mult * V[2][0], r3, phi3, theta3);
		//		math::tospherical(V2mult * V[0][1], V2mult * V[1][1], V2mult * V[2][1], r2, phi2, theta2);
		//		math::tospherical(V3mult * V[0][2], V3mult * V[1][2], V3mult * V[2][2], r1, phi1, theta1);

		//		double e = 0.0;
		//		if(abs(lambda1) > 1e-6)
		//			e = sqrt(1.0 - (lambda3 * lambda3) / (lambda1 * lambda1));

		//		results.push_back(e);

		//		results.push_back(lambda1);
		//		results.push_back(lambda2);
		//		results.push_back(lambda3);

		//		results.push_back(phi1);
		//		results.push_back(theta1);
		//		
		//		results.push_back(phi2);
		//		results.push_back(theta2);
		//		
		//		results.push_back(phi3);
		//		results.push_back(theta3);


		//		// Calculate maximum radius
		//		double maxr = 0;
		//		math::Vec3d t3dir(V1mult * V[0][0], V1mult * V[1][0], V1mult * V[2][0]);
		//		math::Vec3d t2dir(V2mult * V[0][1], V2mult * V[1][1], V2mult * V[2][1]);
		//		math::Vec3d t1dir(V3mult * V[0][2], V3mult * V[1][2], V3mult * V[2][2]);
		//		double l1 = sqrt(lambda1);
		//		double l2 = sqrt(lambda2);
		//		double l3 = sqrt(lambda3);
		//		math::Vec3d c(mx, my, mz);
		//		double maxd1 = 0;
		//		double mind1 = 0;
		//		double maxd2 = 0;
		//		double mind2 = 0;
		//		double maxd3 = 0;
		//		double mind3 = 0;
		//		double boundingScale = 0;
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			math::Vec3d p(points[n][0] - mx, points[n][1] - my, points[n][2] - mz);
		//			double r = p.norm();
		//			maxr = math::max(maxr, r);

		//			r = t1dir.dot(p);
		//			if(r > maxd1)
		//				maxd1 = r;
  //                  if(r < mind1)
  //                      mind1 = r;

		//			r = t2dir.dot(p);
		//			if(r > maxd2)
		//				maxd2 = r;
  //                  if(r < mind2)
  //                      mind2 = r;

		//			r = t3dir.dot(p);
		//			if(r > maxd3)
		//				maxd3 = r;
  //                  if(r < mind3)
  //                      mind3 = r;


		//			double f = getEllipsoidFunctionValue(math::Vec3d(points[n][0], points[n][1], points[n][2]),
		//						c,
		//						l1, l2, l3,
		//						phi1, theta1,
		//						phi2, theta2,
		//						phi3, theta3);
		//			double scale = sqrt(f);
		//			if(!isnan(scale) && !isinf(scale))
		//				boundingScale = math::max(boundingScale, scale);
		//		}
		//		results.push_back(maxr);
		//		results.push_back(maxd1 - mind1);
		//		results.push_back(maxd2 - mind2);
		//		results.push_back(maxd3 - mind3);
		//		results.push_back(boundingScale);

		//		return results;
		//	}
		//};

		///**
		//Returns measurements of histogram (from another image) at the points of the particle.
		//*/
  //      template<class POINT, typename pixel_t, typename PTYPE2> class Histogram : public Analyzer<POINT, pixel_t>
  //      {
  //      private:
  //          Image<PTYPE2>& anotherImage;

		//	vector<vector<float> > pixelValues;

  //      public:

		//	const vector<vector<float> >& getExtraResults() const
		//	{
		//		return pixelValues;
		//	}

  //          Histogram(Image<PTYPE2>& anotherImage) : anotherImage(anotherImage)
  //          {
  //          }

  //          virtual vector<string> getTitles() const
  //          {
		//		vector<string> labels;
		//		labels.push_back("Mean");
		//		labels.push_back("Stdev");
		//		labels.push_back("Min");
		//		labels.push_back("Max");
  //              return labels;
  //          }

  //          virtual vector<double> analyze(const vector<POINT >& points)
  //          {
		//		pixelValues.push_back(vector<float>());
		//		vector<float>& vec = pixelValues[pixelValues.size()-1];
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

		//		vector<double> results;
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
			virtual vector<string> getTitles() const
			{
				vector<string> labels;
				labels.push_back("minx");
				labels.push_back("maxx");
				labels.push_back("miny");
				labels.push_back("maxy");
				labels.push_back("minz");
				labels.push_back("maxz");
				return labels;
			}

			virtual vector<double> analyze(const vector<POINT >& points) const
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

				vector<double> results;
				results.push_back((double)min.x);
				results.push_back((double)max.x);
				results.push_back((double)min.y);
				results.push_back((double)max.y);
				results.push_back((double)min.z);
				results.push_back((double)max.z);
				return results;
			}

			virtual std::string name() const
			{
				return "bounds";
			}

			virtual std::string description() const
			{
				return "Calculates axis-aligned bounding box of the particle. Returns minimum and maximum coordinates of particle points in all coordinate dimensions. Outputs columns 'minX' and 'maxX' for each dimension, where X is either x, y, or z.";
			}
		};

		/**
		Calculates minimum bounding sphere of the particle
		*/
		template<class POINT, typename pixel_t> class BoundingSphere : public Analyzer<POINT, pixel_t>
		{
		private:


		public:
			virtual vector<string> getTitles() const
			{
				vector<string> labels;
				labels.push_back("bounding sphere X [pixel]");
				labels.push_back("bounding sphere Y [pixel]");
				labels.push_back("bounding sphere Z [pixel]");
				labels.push_back("bounding sphere radius [pixel]");
				return labels;
			}

			virtual vector<double> analyze(const vector<POINT >& points) const
			{
				Sphere<double> bs = Sphere<double>::miniball(points);

				vector<double> results;
				results.push_back(bs.center.x);
				results.push_back(bs.center.y);
				results.push_back(bs.center.z);
				results.push_back(bs.radius);
				return results;
			}

			virtual std::string name() const
			{
				return "boundingsphere";
			}

			virtual std::string description() const
			{
				return "Shows position and radius of the smallest possible sphere that contains all the points in the particle. Outputs columns 'bounding sphere X', 'bounding sphere Y', 'bounding sphere Z', and 'bounding sphere radius'.";
			}
		};
	}

	/**
	Creates list of all particle analyzers.
	*/
	template<typename pixel_t> AnalyzerSet<math::Vec3sc, pixel_t> allAnalyzers(const math::Vec3c& dimensions)
	{
		AnalyzerSet<math::Vec3sc, pixel_t> analyzers;
		analyzers.push_back(shared_ptr<Analyzer<math::Vec3sc, pixel_t> >(new analyzers::Coordinates<math::Vec3sc, pixel_t>()));
		analyzers.push_back(shared_ptr<Analyzer<math::Vec3sc, pixel_t> >(new analyzers::IsOnEdge<math::Vec3sc, pixel_t>(dimensions)));
		analyzers.push_back(shared_ptr<Analyzer<math::Vec3sc, pixel_t> >(new analyzers::Volume<math::Vec3sc, pixel_t>()));
		//analyzers.push_back(shared_ptr<Analyzer<math::Vec3sc, pixel_t> >(new analyzers::PCA3D<math::Vec3sc, pixel_t>()));
		//analyzers.push_back(shared_ptr<Analyzer<math::Vec3sc, pixel_t> >(new analyzers::SurfaceArea3D<math::Vec3sc, pixel_t>()));
		//analyzers.push_back(shared_ptr<Analyzer<math::Vec3sc, pixel_t> >(new analyzers::Color<math::Vec3sc, pixel_t>()));
		analyzers.push_back(shared_ptr<Analyzer<math::Vec3sc, pixel_t> >(new analyzers::BoundingBox3D<math::Vec3sc, pixel_t>()));
		analyzers.push_back(shared_ptr<Analyzer<math::Vec3sc, pixel_t> >(new analyzers::BoundingSphere<math::Vec3sc, pixel_t>()));

		return analyzers;
	}

	/**
	Creates list of all particle analyzers.
	*/
	template<typename pixel_t> AnalyzerSet<math::Vec3sc, pixel_t> allAnalyzers(const Image<pixel_t>& img)
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
		template<typename pixel_t> shared_ptr<Analyzer<math::Vec3sc, pixel_t> > getAnalyzer(AnalyzerSet<math::Vec3sc, pixel_t>& all, const string& name)
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

			throw ITLException(string("Invalid analyzer name: ") + name);
		}
	}

	/**
	Converts list of analyzer names to analyzer set.
	*/
	template<typename pixel_t> AnalyzerSet<math::Vec3sc, pixel_t> createAnalyzers(string names, const math::Vec3c& dimensions)
	{
		// Replace all possible delimiters by space
		replace_if(names.begin(), names.end(), internals::isnotalnum, ' ');

		vector<string> parts = split(names, false, ' ');

		AnalyzerSet<math::Vec3sc, pixel_t> all = allAnalyzers<pixel_t>(dimensions);
		AnalyzerSet<math::Vec3sc, pixel_t> result;

		for (size_t n = 0; n < parts.size(); n++)
		{
			shared_ptr<Analyzer<math::Vec3sc, pixel_t> > tmp = internals::getAnalyzer<pixel_t>(all, parts[n]);
			result.push_back(tmp);
		}

		return result;
	}

}
