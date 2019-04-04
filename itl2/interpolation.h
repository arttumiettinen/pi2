#pragma once

#include <memory>

#include "image.h"
#include "boundarycondition.h"
#include "math/mathutils.h"
#include "math/vec3.h"
#include "math/numberutils.h"
#include "interpolationmode.h"

using math::Vec3;
using math::pixelRound;
using math::NumberUtils;

namespace itl2
{

	/*
	Gets a pixel value from the image, obeys boundary conditions.
	*/
	template<typename input_t> input_t getPixelSafe(const Image<input_t>& img, coord_t x, coord_t y, coord_t z, BoundaryCondition bc)
	{
		if (img.isInImage(x, y, z))
			return img(x, y, z);

		if (bc == BoundaryCondition::Zero)
			return input_t();

		// bc == Nearest
		math::clamp<coord_t>(x, 0, img.width() - 1);
		math::clamp<coord_t>(y, 0, img.height() - 1);
		math::clamp<coord_t>(z, 0, img.depth() - 1);
		return img(x, y, z);
	}

	/**
	Base class for interpolator functors.
	*/
	template<typename output_t, typename input_t, typename real_t = typename NumberUtils<output_t>::RealFloatType> class Interpolator
	{
	private:
		BoundaryCondition bc;

	protected:
		/*
		Gets pixel from the image, obeys boundary conditions.
		*/
		input_t getPixelSafe(const Image<input_t>& img, coord_t x, coord_t y, coord_t z) const
		{
			return itl2::getPixelSafe(img, x, y, z, bc);
		}

	public:
		Interpolator(BoundaryCondition bc) : bc(bc)
		{

		}

		virtual ~Interpolator()
		{

		}

		/**
		Interpolates given image at the given location.
		*/
		virtual output_t operator()(const Image<input_t>& img, real_t x, real_t y, real_t z) const = 0;

		/**
		Interpolates given image at the given location.
		*/
		output_t operator()(const Image<input_t>& img, const Vec3<real_t>& x) const
		{
			return operator()(img, x.x, x.y, x.z);
		}
	};

	/**
	Nearest neighbour interpolation functor.
	*/
	template<typename output_t, typename input_t, typename real_t = typename NumberUtils<output_t>::RealFloatType> class NearestNeighbourInterpolator : public Interpolator<output_t, input_t, real_t>
	{
	public:
		NearestNeighbourInterpolator(BoundaryCondition bc) : Interpolator<output_t, input_t, real_t>(bc)
		{

		}

		virtual output_t operator()(const Image<input_t>& img, real_t x, real_t y, real_t z) const
		{
			coord_t ix = (coord_t)::round(x);
			coord_t iy = (coord_t)::round(y);
			coord_t iz = (coord_t)::round(z);
			return pixelRound<output_t>(this->getPixelSafe(img, ix, iy, iz));
		}
	};

	/**
	Linear interpolation functor.
	@param output_t Type of output values.
	@param input_t Pixel type of input image.
	@param real_t Scalar real number type, typically double.
	@param intermediate_t Type of intermediate values, typically double or Vec3d etc.
	*/
	template<typename output_t, typename input_t, typename real_t = typename NumberUtils<output_t>::RealFloatType, typename intermediate_t = typename NumberUtils<output_t>::FloatType > class LinearInterpolator : public Interpolator<output_t, input_t, real_t>
	{
	private:
		inline real_t w_lin(real_t dx) const
		{
			return 1 - fabs(dx);
		}

	public:
		LinearInterpolator(BoundaryCondition bc) : Interpolator<output_t, input_t, real_t>(bc)
		{

		}

		virtual output_t operator()(const Image<input_t>& img, real_t x, real_t y, real_t z) const
		{
			coord_t u0 = (coord_t)floor(x);
			coord_t v0 = (coord_t)floor(y);
			coord_t w0 = (coord_t)floor(z);

			intermediate_t r = intermediate_t();
			for (int k = 0; k <= 1; k++)
			{
				intermediate_t q = intermediate_t();
				coord_t w = w0 + k;

				for (int j = 0; j <= 1; j++)
				{
					intermediate_t p = intermediate_t();
					coord_t v = v0 + j;

					for (int i = 0; i <= 1; i++)
					{
						coord_t u = u0 + i;
						input_t pixval = this->getPixelSafe(img, u, v, w);
						p = p + pixval * w_lin(x - u);
					}

					q = q + p * w_lin(y - v);
				}

				r = r + q * w_lin(z - w);
			}

			return pixelRound<output_t>(r);
		}
	};


	/**
	Linear interpolation functor with invalid value support.
	Pixels having a specific value are not used in the interpolation process.
	@param output_t Type of output values.
	@param input_t Pixel type of input image.
	@param real_t Scalar real number type, typically double.
	@param intermediate_t Type of intermediate values, typically double or Vec3d etc.
	*/
	template<typename output_t, typename input_t, typename real_t = typename NumberUtils<output_t>::RealFloatType, typename intermediate_t = typename NumberUtils<output_t>::FloatType> class LinearInvalidValueInterpolator : public Interpolator<output_t, input_t, real_t>
	{
	private:
		inline real_t w_lin(real_t dx) const
		{
			return 1 - fabs(dx);
		}

		input_t invalidInputValue;
		output_t invalidOutputValue;

	public:
		LinearInvalidValueInterpolator(BoundaryCondition bc, input_t invalidInputValue, output_t invalidOutputValue) : Interpolator<output_t, input_t, real_t>(bc), invalidInputValue(invalidInputValue), invalidOutputValue(invalidOutputValue)
		{

		}

		virtual output_t operator()(const Image<input_t>& img, real_t x, real_t y, real_t z) const
		{
			coord_t u0 = (coord_t)floor(x);
			coord_t v0 = (coord_t)floor(y);
			coord_t w0 = (coord_t)floor(z);

			intermediate_t r = intermediate_t();
			real_t wTotR = 0;
			for (int k = 0; k <= 1; k++)
			{
				intermediate_t q = intermediate_t();
				coord_t w = w0 + k;
				real_t wTotQ = 0;

				for (int j = 0; j <= 1; j++)
				{
					intermediate_t p = intermediate_t();
					coord_t v = v0 + j;
					real_t wTotP = 0;

					for (int i = 0; i <= 1; i++)
					{
						coord_t u = u0 + i;
						input_t pixval = this->getPixelSafe(img, u, v, w);
						if (pixval != invalidInputValue)
						{
							real_t ww = w_lin(x - u);
							p = p + pixval * ww;
							wTotP += ww;
						}
					}

					if (wTotP > 0)
					{
						real_t ww = w_lin(y - v);
						q = q + p / wTotP * ww;
						wTotQ += ww;
					}
				}

				if (wTotQ > 0)
				{
					real_t ww = w_lin(z - w);
					r = r + q / wTotQ * ww;
					wTotR += ww;
				}
			}

			if (wTotR > 0)
			{
				r /= wTotR;
				return pixelRound<output_t>(r);
			}

			return invalidOutputValue;
		}
	};

	/**
	Cubic interpolation functor.
	This code is from
	https://github.com/imagingbook/imagingbook-common/blob/master/src/main/java/imagingbook/lib/interpolation/BicubicInterpolator.java
	See http://imagingbook.com
	*/
	template<typename output_t, typename input_t, typename real_t = typename NumberUtils<output_t>::RealFloatType, typename intermediate_t = typename NumberUtils<output_t>::FloatType> class CubicInterpolator : public Interpolator<output_t, input_t, real_t>
	{
	private:
		real_t a;

		/*
		Cubic weight funtion in 1D.
		*/
		real_t w_cub(real_t x) const
		{
			if (x < 0)
				x = -x;
			
			if (x < 1)
				return (-a + 2) * x * x * x + (a - 3) * x * x + 1;
			else if (x < 2)
				return -a * x * x * x + 5 * a * x * x - 8 * a * x + 4 * a;

			return 0;
		}

	public:
		CubicInterpolator(BoundaryCondition bc, real_t sharpness = (real_t)0.5) : Interpolator<output_t, input_t, real_t>(bc), a(sharpness)
		{

		}

		virtual output_t operator()(const Image<input_t>& img, real_t x, real_t y, real_t z) const
		{
			coord_t u0 = (coord_t)floor(x);
			coord_t v0 = (coord_t)floor(y);
			coord_t w0 = (coord_t)floor(z);

			intermediate_t r = intermediate_t();
			for (int k = 0; k <= 3; k++)
			{
				coord_t w = w0 - 1 + k;
				intermediate_t q = intermediate_t();

				for (int j = 0; j <= 3; j++)
				{
					coord_t v = v0 - 1 + j;
					intermediate_t p = intermediate_t();

					for (int i = 0; i <= 3; i++)
					{
						coord_t u = u0 - 1 + i;
						input_t pixval = this->getPixelSafe(img, u, v, w);
						p = p + pixval * w_cub(x - u);
					}

					q = q + p * w_cub(y - v);
				}

				r = r + q * w_cub(z - w);
			}

			return pixelRound<output_t>(r);
		}
	};


	/**
	Cubic interpolation functor with invalid value support.
	Pixels having a specific value are not used in the interpolation process.
	This code is partially from
	https://github.com/imagingbook/imagingbook-common/blob/master/src/main/java/imagingbook/lib/interpolation/BicubicInterpolator.java
	See http://imagingbook.com
	*/
	template<typename output_t, typename input_t, typename real_t = typename NumberUtils<output_t>::RealFloatType, typename intermediate_t = typename NumberUtils<output_t>::FloatType> class CubicInvalidValueInterpolator : public Interpolator<output_t, input_t, real_t>
	{
	private:
		real_t a;

		/*
		Cubic weight funtion in 1D.
		*/
		real_t w_cub(real_t x) const
		{
			if (x < 0)
				x = -x;

			if (x < 1)
				return (-a + 2) * x * x * x + (a - 3) * x * x + 1;
			else if (x < 2)
				return -a * x * x * x + 5 * a * x * x - 8 * a * x + 4 * a;

			return 0;
		}

		input_t invalidInputValue;
		output_t invalidOutputValue;

	public:
		CubicInvalidValueInterpolator(BoundaryCondition bc, input_t invalidInputValue, output_t invalidOutputValue, real_t sharpness = (real_t)0.5) : Interpolator<output_t, input_t, real_t>(bc), a(sharpness), invalidInputValue(invalidInputValue), invalidOutputValue(invalidOutputValue)
		{

		}

		virtual output_t operator()(const Image<input_t>& img, real_t x, real_t y, real_t z) const
		{
			coord_t u0 = (coord_t)floor(x);
			coord_t v0 = (coord_t)floor(y);
			coord_t w0 = (coord_t)floor(z);

			intermediate_t r = intermediate_t();
			real_t wTotR = 0;
			for (int k = 0; k <= 3; k++)
			{
				intermediate_t q = intermediate_t();
				coord_t w = w0 - 1 + k;
				real_t wTotQ = 0;

				for (int j = 0; j <= 3; j++)
				{
					intermediate_t p = intermediate_t();
					coord_t v = v0 - 1 + j;
					real_t wTotP = 0;

					for (int i = 0; i <= 3; i++)
					{
						coord_t u = u0 - 1 + i;
						input_t pixval = this->getPixelSafe(img, u, v, w);
						if (pixval != invalidInputValue)
						{
							real_t ww = w_cub(x - u);
							p = p + pixval * ww;
							wTotP += ww;
						}
					}

					if (wTotP > 0)
					{
						real_t ww = w_cub(y - v);
						q = q + p / wTotP * ww;
						wTotQ += ww;
					}
				}

				if (wTotQ > 0)
				{
					real_t ww = w_cub(z - w);
					r = r + q / wTotQ * ww;
					wTotR += ww;
				}
			}

			if (wTotR > 0)
			{
				r /= wTotR;
				return pixelRound<output_t>(r);
			}

			return invalidOutputValue;
		}
	};

	/**
	Create Interpolator object from interpolation mode and boundary condition.
	*/
	template<typename output_t, typename input_t, typename real_t = typename NumberUtils<output_t>::RealFloatType> std::shared_ptr<Interpolator<output_t, input_t, real_t> > createInterpolator(InterpolationMode mode, BoundaryCondition bc)
	{
		switch (mode)
		{
			case InterpolationMode::Nearest: return std::make_shared<NearestNeighbourInterpolator<output_t, input_t, real_t> >(bc);
			case InterpolationMode::Linear: return std::make_shared<LinearInterpolator<output_t, input_t, real_t> >(bc);
			case InterpolationMode::Cubic: return std::make_shared<CubicInterpolator<output_t, input_t, real_t>>(bc);
		}
		throw ITLException("Unsupported interpolation mode.");
	}
}
