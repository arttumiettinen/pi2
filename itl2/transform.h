#pragma once

#include <omp.h>
#include <vector>
#include <array>
#include "image.h"
#include "math/vec3.h"
#include "math/matrix3x3.h"
#include "interpolation.h"
#include "utilities.h"
#include "math/numberutils.h"
#include "math/vectoroperations.h"
#include "median.h"
#include "filters.h"
#include "conversions.h"
#include "iteration.h"

namespace itl2
{
	/**
	Possible directions for volume image re-slicing.
	*/
	enum class ResliceDirection
	{
		Top,
		Bottom,
		Left,
		Right
	};

	template<>
	inline std::string toString(const ResliceDirection& dir)
	{
		switch (dir)
		{
		case ResliceDirection::Top: return "Top";
		case ResliceDirection::Bottom: return "Bottom";
		case ResliceDirection::Left: return "Left";
		case ResliceDirection::Right: return "Right";
		default: throw ITLException("Unknown reslice direction.");
		}
	}

	template<>
	inline ResliceDirection fromString(const string& str)
	{
		string str2 = str;
		trim(str2);
		toLower(str2);
		if (str2 == "top")
			return ResliceDirection::Top;
		if (str2 == "bottom")
			return ResliceDirection::Bottom;
		if (str2 == "left")
			return ResliceDirection::Left;
		if (str2 == "right")
			return ResliceDirection::Right;
		
		throw ITLException("Invalid reslice direction.");
	}

	/**
	Re-slices the input image into the output image.
	Rotates image like a cube, and makes the face of the cube defined by ResliceDirection to be the first slice in the output image.
	*/
	template<typename pixel_t, typename out_t> void reslice(const Image<pixel_t>& in, Image<out_t>& out, ResliceDirection dir)
	{
		out.mustNotBe(in);
		
		if (dir == ResliceDirection::Top)
		{
			// Top:
			// x' = x
			// y' = -z
			// z' = y

			out.ensureSize(in.width(), in.depth(), in.height());

			forAllPixels(in, [&](coord_t x, coord_t y, coord_t z)
			{
				out(x, in.depth() - 1 - z, y) = pixelRound<out_t>(in(x, y, z));
			});
		}
		else if (dir == ResliceDirection::Bottom)
		{
			// Bottom:
			// x' = x
			// y' = z
			// z' = -y

			out.ensureSize(in.width(), in.depth(), in.height());

			forAllPixels(in, [&](coord_t x, coord_t y, coord_t z)
			{
				out(x, z, in.height() - 1 - y) = pixelRound<out_t>(in(x, y, z));
			});
		}
		else if (dir == ResliceDirection::Left)
		{
			// Left:
			// x' = -z
			// y' = y
			// z' = x

			out.ensureSize(in.depth(), in.height(), in.width());

			forAllPixels(in, [&](coord_t x, coord_t y, coord_t z)
			{
				out(in.depth() - 1 - z, y, x) = pixelRound<out_t>(in(x, y, z));
			});
		}
		else if (dir == ResliceDirection::Right)
		{
			// Right:
			// x' = z
			// y' = y
			// z' = -x

			out.ensureSize(in.depth(), in.height(), in.width());

			forAllPixels(in, [&](coord_t x, coord_t y, coord_t z)
			{
				out(z, y, in.width() - 1 - x) = pixelRound<out_t>(in(x, y, z));
			});
		}
		else
		{
			throw ITLException("Unsupported reslice direction.");
		}
		
	}

	/**
	Flips image 'in' in one or more dimensions and stores the output to 'out'.
	*/
	template<typename pixel_t, typename out_t> void flip(const Image<pixel_t>& in, Image<out_t>& out, bool flipX, bool flipY, bool flipZ)
	{
		out.mustNotBe(in);
		out.ensureSize(in);

		#pragma omp parallel for if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			coord_t oz = z;
			if (flipZ)
				oz = out.depth() - 1 - z;

			for (coord_t y = 0; y < out.height(); y++)
			{
				coord_t oy = y;
				if (flipY)
					oy = out.height() - 1 - y;

				for (coord_t x = 0; x < out.width(); x++)
				{
					coord_t ox = x;
					if (flipX)
						ox = out.width() - 1 - x;

					out(ox, oy, oz) = pixelRound<out_t>(in(x, y, z));
				}
			}
		}

	}

	/**
	Flips image in single dimension.
	*/
	template<typename pixel_t> void flip(Image<pixel_t>& img, size_t dimension = 0)
	{
		if (dimension == 0)
		{
			forAllInBox(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(img.width() / 2, img.height(), img.depth())), [&](coord_t x, coord_t y, coord_t z)
			{
				std::swap(img(x, y, z), img(img.width() - 1 - x, y, z));
			});
		}
		else if (dimension == 1)
		{
			forAllInBox(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(img.width(), img.height() / 2, img.depth())), [&](coord_t x, coord_t y, coord_t z)
			{
				std::swap(img(x, y, z), img(x, img.height() - 1 - y, z));
			});
		}
		else if (dimension == 2)
		{
			forAllInBox(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(img.width(), img.height(), img.depth() / 2)), [&](coord_t x, coord_t y, coord_t z)
			{
				std::swap(img(x, y, z), img(x, y, img.depth() - 1 - z));
			});
		}
		else
		{
			throw ITLException(string("Unsupported dimensionality: ") + toString(dimension));
		}
	}

	/**
	Rotates input image 90 degrees clockwise around the z-axis.
	Sets the size of the output image automatically.
	*/
	template<typename pixel_t, typename out_t> void rot90cw(const Image<pixel_t>& in, Image<out_t>& out)
	{
		out.mustNotBe(in);
		out.ensureSize(in.height(), in.width(), in.depth());

		forAllPixels(in, [&](coord_t x, coord_t y, coord_t z)
		{
			out(out.width() - 1 - y, x, z) = pixelRound<out_t>(in(x, y, z));
		});
	}

	/**
	Rotates input image 90 degrees counterclockwise around the z-axis.
	Sets the size of the output image automatically.
	*/
	template<typename pixel_t, typename out_t> void rot90ccw(const Image<pixel_t>& in, Image<out_t>& out)
	{
		out.mustNotBe(in);
		out.ensureSize(in.height(), in.width(), in.depth());

		forAllPixels(in, [&](coord_t x, coord_t y, coord_t z)
		{
			out(y, out.height() - 1 - x, z) = pixelRound<out_t>(in(x, y, z));
		});
	}

	/**
	Rotates image around the given axis by given angle (in radians).
	Does not change the size of the output image.
	This function tranforms pixel positions according to
	outPos = R * (inPos - inCenter) + outCenter,
	where inPos is position in the input image, outPos is position in the output image,
	R is rotation matrix that rotates around given axis by given amount, and
	inCenter and outCenter define center point for the rotation.
	Pixel at inCenter will map to position outCenter in the output image.
	*/
	template<typename pixel_t, typename out_t> void rotate(const Image<pixel_t>& in, Image<out_t>& out, double angle, const Vec3d& axis, const Vec3d& inCenter, const Vec3d& outCenter, const Interpolator<out_t, pixel_t>& interpolate = LinearInterpolator<out_t, pixel_t>(BoundaryCondition::Zero))
	{
		out.mustNotBe(in);

		using real_t = double;
		Matrix3x3<real_t> R = Matrix3x3<real_t>::rotationMatrix(angle, axis);
		R.transpose();

		forAllPixels(out, [&](coord_t x, coord_t y, coord_t z)
		{
			Vec3d outPos((real_t)x, (real_t)y, (real_t)z);
			Vec3d inPos = R * (outPos - outCenter) + inCenter;
			out(x, y, z) = interpolate(in, Vec3<typename NumberUtils<out_t>::RealFloatType>(inPos));
		});
	}

	/**
	Rotates image around the given axis by given angle (in radians).
	Does not change the size of the output image.
	This function tranforms pixel positions according to
	outPos = R * (inPos - inCenter) + outCenter,
	where inPos is position in the input image, outPos is position in the output image,
	R is rotation matrix that rotates around given axis by given amount, and
	inCenter and outCenter are the center points of the input and output images, respectively.
	Pixel at the center of the input image will map to the center of the output image.
	*/
	template<typename pixel_t, typename out_t> void rotate(const Image<pixel_t>& in, Image<out_t>& out, double angle, const Vec3d& axis = Vec3d(0, 0, 1), const Interpolator<out_t, pixel_t>& interpolate = LinearInterpolator<out_t, pixel_t>(BoundaryCondition::Zero))
	{
		rotate(in, out, angle, axis, Vec3d(in.dimensions()) / 2.0, Vec3d(out.dimensions()) / 2.0, interpolate);
	}


	/**
	Shifts the input image and stores the result in the output image.
	@param in Image to be shifted.
	@param out Output image.
	@param shift Shift vector.
	@param interpolate Interpolation type.
	*/
	template<typename pixel_t, typename out_t> void translate(const Image<pixel_t>& in, Image<out_t>& out, const Vec3d& shift, const Interpolator<out_t, pixel_t>& interpolate = LinearInterpolator<out_t, pixel_t>(BoundaryCondition::Zero))
	{
		using real_t = typename NumberUtils<out_t>::RealFloatType;
		out.mustNotBe(in);
		if(out.dimensions().max() <= 1)
			out.ensureSize(in);

		forAllPixels(out, [&](coord_t x, coord_t y, coord_t z)
		{
			real_t xs = (real_t)(x - shift.x);
			real_t ys = (real_t)(y - shift.y);
			real_t zs = (real_t)(z - shift.z);
			out(x, y, z) = interpolate(in, xs, ys, zs);
		});
	}
	/**
	Transposes the input image and stores the result in the output image.
	@param in Image to be transposed.
	@param out Output image.
	@param order Permutation indicating the order for transposing.
	*/
	template<typename pixel_t>
	void transpose(Image<pixel_t>& in, const Vec3c& order, int fillValue)
	{
		Vec3c transposedShape = in.dimensions().transposed(order);
		Image<pixel_t> temp(transposedShape, fillValue);
		forAllPixels(in, [&](coord_t x, coord_t y, coord_t z)
		{
		  Vec3c cords(x, y, z);
		  Vec3c transposedCords = cords.transposed(order);
		  temp(transposedCords) = in(cords);
		});

		in.ensureSize(transposedShape);
		forAllPixels(temp, [&](coord_t x, coord_t y, coord_t z)
		{
		  Vec3c cords(x, y, z);
		  in(cords) = temp(cords);
		});
	}

	/**
	Crops the input image to the size of the output image and places the result to the output image.
	Left-top corner of the output image is placed at the given position in the input image.
	*/
	template<typename pixel_t, typename out_t> void crop(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& outPos)
	{
		out.mustNotBe(in);

		setValue(out, out_t());

		AABoxc inBox = AABoxc::fromPosSize(Vec3c(0, 0, 0), in.dimensions());
		AABoxc outBox = AABoxc::fromPosSize(Vec3c(0, 0, 0), out.dimensions());
		AABoxc clippedBox = outBox.translate(outPos).intersection(inBox).translate(-outPos);

		forAllInBox(clippedBox, [&](coord_t x, coord_t y, coord_t z)
		{
			Vec3c xi = Vec3c(x, y, z) + outPos;
			out(x, y, z) = pixelRound<out_t>(in(xi));
		});
	}

	/**
	Copies 'block' to 'target', to given position.
	In other words, makes the 'inverse' of crop operation:
	crop(img, block, pos)
	assigns pixels in img to pixels of block but
	copyValues(img, block, pos)
	assigns pixels in block to pixels of img.
	Does not resize target image.
	@param target Target image where the pixels are written to.
	@param block Source image where the pixels are copied from.
	@param pos Position of source image data in the target image.
	*/
	template<typename pixel_t, typename out_t> void copyValues(Image<pixel_t>& target, const Image<out_t>& block, const Vec3c& pos = Vec3c(0,0,0))
	{
		target.mustNotBe(block);
		std::cout << "copyValues pos: " << pos << std::endl;
		AABox<coord_t> sourceBox = AABox<coord_t>::fromPosSize(Vec3c(0, 0, 0), block.dimensions());
		AABox<coord_t> targetBox = AABox<coord_t>::fromPosSize(Vec3c(0, 0, 0), target.dimensions());
		AABox<coord_t> clippedBox = sourceBox.translate(pos).intersection(targetBox).translate(-pos);

		// TODO: This optimization is not correct at the moment. It does not account for clipping above.
		//			I will remove it for now, and debug it later if it is really necessary.
		//if (pos.x == 0 && pos.y == 0 && block.dimensions().x == target.dimensions().x && block.dimensions().y == target.dimensions().y)
		//{
		//	// Shift in z only. This can be done very fast as memory copy-style operation.

		//	size_t targetStartIndex = target.getLinearIndex(0, 0, pos.z);

		//	#pragma omp parallel for if(block.pixelCount() > PARALLELIZATION_THRESHOLD)
		//	for (coord_t n = 0; n < block.pixelCount(); n++)
		//	{
		//		target(targetStartIndex + n) = pixelRound<pixel_t, out_t>(block(n));
		//	}
		//}
		//else
		//{
			// General shift
			forAllInBox(clippedBox, [&](coord_t x, coord_t y, coord_t z)
			{
				Vec3c xi = Vec3c(x, y, z) + pos;
				pixel_t result = pixelRound<pixel_t>(block(x, y, z));
				std::cout << "setting target(= "<<xi << ") to block("<<x<<", "<<y<<", "<<z<<")=?"<<std::endl;
				target(xi) = result;
			});
		//}
	}


    /**
	Copies some pixels from 'source' to 'target', to given position.
	@param target Target image where the pixels are written to.
	@param block Source image where the pixels are copied from.
	@param targetPos Position of source image data in the target image.
	@param blockPos Position of the first pixel of the block to copy.
	@param copySize Size of the block to copy.
	*/
	template<typename pixel_t, typename out_t> void copyValues(Image<pixel_t>& target, const Image<out_t>& block, const Vec3c& targetPos, const Vec3c& sourcePos, const Vec3c& copySize)
	{
		target.mustNotBe(block);

		// The region where we are going to copy from.
		AABox<coord_t> sourceBox = AABox<coord_t>::fromPosSize(sourcePos, copySize);
		
		// Clip it to the available region in the source block
		AABox<coord_t> fullSourceBox = AABox<coord_t>::fromPosSize(Vec3c(0, 0, 0), block.dimensions());
		sourceBox = sourceBox.intersection(fullSourceBox);

		// Clip it to the target box so that we don't copy values out of target image.
		AABox<coord_t> targetBox = AABox<coord_t>::fromPosSize(Vec3c(0, 0, 0), target.dimensions());
		
		AABox<coord_t> clippedSourceBox = sourceBox.translate(-sourcePos).translate(targetPos).intersection(targetBox).translate(-targetPos).translate(sourcePos);

		// General shift
		forAllInBox(clippedSourceBox, [&](coord_t x, coord_t y, coord_t z)
			{
				Vec3c xi = Vec3c(x, y, z) - sourcePos + targetPos;
				target(xi) = pixelRound<pixel_t>(block(x, y, z));
			});
	}

	namespace internals
	{
		/**
		Gets a block of pixel values from image.
		*/
		template<typename pixel_t> void getBlock(const Image<pixel_t>& img, std::vector<pixel_t>& block, const Vec3c& pos, const Vec3c& size)
		{
			Vec3c end = min(pos + size, img.dimensions());

			coord_t n = 0;
			for (coord_t z = pos.z; z < end.z; z++)
			{
				for (coord_t y = pos.y; y < end.y; y++)
				{
					for (coord_t x = pos.x; x < end.x; x++)
					{
						block.push_back(img(x, y, z));
					}
				}
			}
		}
	}

	namespace binningop
	{
		template<typename pixel_t, typename out_t> out_t mean(const std::vector<pixel_t>& block)
		{
			using real_t = typename NumberUtils<pixel_t>::RealFloatType;
			using float_t = typename NumberUtils<pixel_t>::FloatType;

			float_t val = 0;
			for (pixel_t p : block)
				val += (float_t)p;
			val /= (real_t)block.size();
			return pixelRound<out_t>(val);
		}

		template<typename pixel_t, typename out_t> pixel_t max(const std::vector<pixel_t>& block)
		{
			return pixelRound<out_t>(itl2::max(block));
		}

		template<typename pixel_t, typename out_t> pixel_t min(const std::vector<pixel_t>& block)
		{
			return pixelRound<out_t>(itl2::min(block));
		}
	}

	/**
	Converts input image to smaller scale by transforming each binSize block of input to one pixel in output.
	The supported transforms are defined in binningop namespace.
	@param in Input image.
	@param out Output image. The image is automatically initialized to correct size.
	@param binSize Bin size. 2 makes the output image dimensions half of the input image dimensions, 3 makes them one third etc.
	@param indicateProgress Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t, out_t operation(const std::vector<pixel_t>&)> void binning(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& binSize, bool indicateProgress = true)
	{
		if (binSize.min() <= 0)
			throw ITLException("Bins size must be positive.");

		out.mustNotBe(in);
		out.ensureSize(in.dimensions().componentwiseDivide(binSize));

		// TODO: This is separable operation for most operations.

		size_t counter = 0;
		#pragma omp parallel if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{
			std::vector<pixel_t> block;
			block.reserve(binSize.x * binSize.y * binSize.z);

			#pragma omp for
			for (coord_t z = 0; z < out.depth(); z++)
			{
				coord_t inz = z * binSize.z;
				coord_t iny = 0;
				for (coord_t y = 0; y < out.height(); y++, iny += binSize.y)
				{
					coord_t inx = 0;
					for (coord_t x = 0; x < out.width(); x++, inx += binSize.x)
					{
						block.clear();
						internals::getBlock(in, block, Vec3c(inx, iny, inz), binSize);
						out(x, y, z) = operation(block);
					}
				}

				showThreadProgress(counter, out.depth(), indicateProgress);
			}
		}
	}

	/**
	Converts input image to smaller scale by averaging each binSize block of input and placing the average to the corresponding pixel of the output image.
	@param in Input image.
	@param out Output image. The image is automatically initialized to correct size.
	@param binSize Bin size. 2 makes the output image dimensions half of the input image dimensions, 3 makes them one third etc.
	@param indicateProgress Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t> void binning(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& binSize, bool indicateProgress = true)
	{
		binning<pixel_t, out_t, binningop::mean<pixel_t> >(in, out, binSize, indicateProgress);
	}

	/**
	Converts input image to smaller scale by averaging binSize^dimensionality blocks.
	@param in Input image.
	@param out Output image. The image is automatically initialized to correct size.
	@param binSize Bin size. 2 makes the output image dimensions half of the input image dimensions, 3 makes them one third etc.
	@param indicateProgress Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t, out_t operation(const std::vector<pixel_t>&)> void binning(const Image<pixel_t>& in, Image<out_t>& out, size_t binSize, bool indicateProgress = true)
	{
		binning<pixel_t, out_t, operation>(in, out, Vec3c(binSize, binSize, binSize), indicateProgress);
	}

	/**
	Converts input image to smaller scale by averaging each binSize^dimensionality block of input and placing the average to the corresponding pixel of the output image.
	@param in Input image.
	@param out Output image. The image is automatically initialized to correct size.
	@param binSize Bin size. 2 makes the output image dimensions half of the input image dimensions, 3 makes them one third etc.
	@param indicateProgress Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t> void binning(const Image<pixel_t>& in, Image<out_t>& out, size_t binSize, bool indicateProgress = true)
	{
		binning<pixel_t, out_t, binningop::mean<pixel_t> >(in, out, Vec3c(binSize, binSize, binSize), indicateProgress);
	}

	/**
	Converts input image to smaller scale by averaging binSize^dimensionality blocks.
	Does not average value specified as an argument.
	@param in Input image.
	@param out Output image. The image is automatically initialized to correct size.
	@param amount Bin size. 2 makes the output image dimensions half of the input image dimensions, 3 makes them one third etc.
	@param badValue Value that should not be considered in the averaging calculations.
	@param undefinedValue Value that is placed to those pixels of the output image that do not correspond to any valid pixels in the input image.
	*/
	template<typename pixel_t, typename out_t> void maskedBinning(const Image<pixel_t>& in, Image<out_t>& out, size_t amount, pixel_t badValue, out_t undefinedValue, bool indicateProgress = true)
	{
		out.mustNotBe(in);
		coord_t binSize = (coord_t)amount;
		out.ensureSize(round(Vec3d(in.dimensions()) / (double)binSize));

		size_t counter = 0;
#pragma omp parallel for if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			coord_t inz = z * binSize;
			coord_t iny = 0;
			for (coord_t y = 0; y < out.height(); y++, iny += binSize)
			{
				coord_t inx = 0;
				for (coord_t x = 0; x < out.width(); x++, inx += binSize)
				{

					coord_t inzEnd = std::min(inz + binSize, in.depth());
					coord_t inyEnd = std::min(iny + binSize, in.height());
					coord_t inxEnd = std::min(inx + binSize, in.width());
					typename NumberUtils<out_t>::FloatType M = 0;
					typename NumberUtils<out_t>::RealFloatType count = 0;
					for (coord_t zz = inz; zz < inzEnd; zz++)
					{
						for (coord_t yy = iny; yy < inyEnd; yy++)
						{
							for (coord_t xx = inx; xx < inxEnd; xx++)
							{
								pixel_t pix = in(xx, yy, zz);
								if (pix != badValue)
								{
									M += pix;
									count++;
								}
							}
						}
					}

					if (count > 0)
					{
						M /= count;
						out(x, y, z) = pixelRound<out_t>(M);
					}
					else
					{
						out(x, y, z) = undefinedValue;
					}

				}
			}

			showThreadProgress(counter, out.depth(), indicateProgress);
		}
	}

	namespace internals
	{
		template<typename in_t, typename out_t, typename real_t> void scaleHelper(const Image<in_t>& in, Image<out_t>& out, const Interpolator<out_t, in_t, real_t>& interpolate = LinearInterpolator<out_t, in_t>(BoundaryCondition::Zero), bool indicateProgress = true, const Vec3d& factor = Vec3d(0, 0, 0), const Vec3d& delta = Vec3d(0, 0, 0))
		{
			size_t counter = 0;
			#pragma omp parallel for if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t z = 0; z < out.depth(); z++)
			{
				real_t sz = (real_t)(z / factor.z) + (real_t)delta.z;
				for (coord_t y = 0; y < out.height(); y++)
				{
					real_t sy = (real_t)(y / factor.y) + (real_t)delta.y;
					for (coord_t x = 0; x < out.width(); x++)
					{
						real_t sx = (real_t)(x / factor.x) + (real_t)delta.x;
						out(x, y, z) = interpolate(in, sx, sy, sz);
					}
				}

				showThreadProgress(counter, out.depth(), indicateProgress);
			}

		}
	}

	/**
	Scales input image to the size of the output image and replaces output image by the scaled image.
	@param in Input image.
	@param out Output image.
	@param interp Interpolation type.
	@param indicateProgress Set to true to show a progress bar.
	@param factor Scaling factor. If zero or negative, determined from dimensions of input and output image.
	@param delta Shift that is added to coordinates of each input point. Used in distributed processing.
	*/
	template<typename in_t, typename out_t> void scale(const Image<in_t>& in, Image<out_t>& out, bool averageWhenDownSizing = true, InterpolationMode interp = InterpolationMode::Linear, BoundaryCondition bc = BoundaryCondition::Zero, bool indicateProgress = true, Vec3d factor = Vec3d(0, 0, 0), const Vec3d& delta = Vec3d(0, 0, 0))
	{
		in.mustNotBe(out);

		if (factor.max() <= 0)
		{
			factor = Vec3d(out.dimensions()).componentwiseDivide(Vec3d(in.dimensions()));
		}


		// This could work but it feels a bit stupid and requires interpolation between two binning values.
		//const Image<in_t>* pSrc = &in;
		//Image<in_t> temp;
		//if (averageWhenDownSizing)
		//{
		//	if (factor.min() < 0.5)
		//	{
		//		Vec3c binSize(1, 1, 1);
		//		while (factor.min() < 0.5)
		//		{
		//			for (size_t i = 0; i < factor.size(); i++)
		//			{
		//				if (factor[i] < 0.5)
		//				{
		//					binSize[i] *= 2;
		//					factor[i] *= 2;
		//				}
		//			}
		//		}

		//		binning(in, temp, binSize, false);
		//		pSrc = &temp;
		//	}
		//}

		// TODO: This is a bit simplistic approach. Both scaling and filtering operations
		// are separable so we could be a lot faster by taking advantage of that. Additionally we
		// calculate too many filtered samples here if factor is small.
		// Additionally this requires a lot of memory (sizeof(float32_t) * size of input in pixels)
		// Good side is that this should be fine with distributed processing and 'delta' parameter.
		if (averageWhenDownSizing && factor.min() < 1)
		{
			Vec3d sigma = 0.5 * (Vec3d(1 / factor.x, 1 / factor.y, 1 / factor.z) - Vec3d(1, 1, 1));
			Image<typename NumberUtils<in_t>::FloatType> temp;
			convert(in, temp);
			internals::sepgauss(temp, sigma, -1, -1, bc);
			internals::scaleHelper<typename NumberUtils<in_t>::FloatType, out_t, typename NumberUtils<in_t>::RealFloatType>(temp, out, *createInterpolator<out_t, typename NumberUtils<in_t>::FloatType>(interp, bc), indicateProgress, factor, delta);
		}
		else
		{
			internals::scaleHelper<in_t, out_t>(in, out, *createInterpolator<out_t, in_t>(interp, bc), indicateProgress, factor, delta);
		}
	}

	

	/**
	Scale binary or label image so that region borders remain smooth even when significant upscaling is made.
	Equals to scaling with nearest neighbour interpolation followed by median filtering.
	Supports only positive integer scaling factors, i.e. only upscaling.
	@param in Input image.
	@param out Output image.
	@param indicateProgress Set to true to show a progress bar.
	@param factor Scaling factor. If zero or negative, determined from dimensions of input and output image.
	@param outputOrigin If processing a block of a full image, value of this argument is the origin of the output image in the coordinates of the full image. Used in distributed processing.
	*/
	template<typename pixel_t> void scaleLabels(const Image<pixel_t>& in, Image<pixel_t>& out, bool indicateProgress = true, Vec3c factor = Vec3c(0, 0, 0), const Vec3c& outputOrigin = Vec3c(0, 0, 0))
	{
		// TODO: This might work better if the radius of the median filtering was 0.5 * scale and not 1.0 * scale.

		in.mustNotBe(out);

		if (factor.max() <= 0)
		{
			factor = out.dimensions().componentwiseDivide(in.dimensions());
		}

		if((out.dimensions() - in.dimensions().componentwiseMultiply(factor)).abs().max() > 1)
			out.ensureSize(in.dimensions().componentwiseMultiply(factor));

		size_t counter = 0;
		#pragma omp parallel if(!omp_in_parallel() && out.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			std::array<coord_t, 27> counts;
			std::array<pixel_t, 27> values;
			std::vector<size_t> tmpSpace;

			#pragma omp for
			for (coord_t zz = 0; zz < out.depth(); zz++) // zz is in output block coordinates
			{
				// Convert to full image coordinates
				coord_t z = zz + outputOrigin.z;

				coord_t zi = z / factor.z;
				coord_t Nz1 = zi * factor.z - z + factor.z;
				coord_t Nz2 = factor.z;
				coord_t Nz3 = z - zi * factor.z + 1;

				// Convert to input image coordinates
				zi = zi - outputOrigin.z / factor.z;

				coord_t zim = zi - 1;
				coord_t zip = zi + 1;

				clamp(zip, (coord_t)0, in.depth() - 1);
				clamp(zim, (coord_t)0, in.depth() - 1);

				clamp(zi, (coord_t)0, in.depth() - 1);

				for (coord_t yy = 0; yy < out.height(); yy++)
				{
					coord_t y = yy;

					coord_t yi = y / factor.y;
					coord_t Ny1 = yi * factor.y - y + factor.y;
					coord_t Ny2 = factor.y;
					coord_t Ny3 = y - yi * factor.y + 1;

					yi = yi - outputOrigin.y / factor.y;

					coord_t yim = yi - 1;
					coord_t yip = yi + 1;

					clamp(yip, (coord_t)0, in.height() - 1);
					clamp(yim, (coord_t)0, in.height() - 1);

					clamp(yi, (coord_t)0, in.height() - 1);

					for (coord_t xx = 0; xx < out.width(); xx++)
					{
						coord_t x = xx;

						coord_t xi = x / factor.x;
						coord_t Nx1 = xi * factor.x - x + factor.x;
						coord_t Nx2 = factor.x;
						coord_t Nx3 = x - xi * factor.x + 1;

						xi = xi - outputOrigin.x / factor.x;

						coord_t xim = xi - 1;
						coord_t xip = xi + 1;

						clamp(xip, (coord_t)0, in.width() - 1);
						clamp(xim, (coord_t)0, in.width() - 1);

						clamp(xi, (coord_t)0, in.width() - 1);

						
						counts[0] = Nx1 * Ny1 * Nz1;
						counts[1] = Nx2 * Ny1 * Nz1;
						counts[2] = Nx3 * Ny1 * Nz1;
						counts[3] = Nx1 * Ny2 * Nz1;
						counts[4] = Nx2 * Ny2 * Nz1;
						counts[5] = Nx3 * Ny2 * Nz1;
						counts[6] = Nx1 * Ny3 * Nz1;
						counts[7] = Nx2 * Ny3 * Nz1;
						counts[8] = Nx3 * Ny3 * Nz1;

						counts[9] = Nx1 * Ny1 * Nz2;
						counts[10] = Nx2 * Ny1 * Nz2;
						counts[11] = Nx3 * Ny1 * Nz2;
						counts[12] = Nx1 * Ny2 * Nz2;
						counts[13] = Nx2 * Ny2 * Nz2;
						counts[14] = Nx3 * Ny2 * Nz2;
						counts[15] = Nx1 * Ny3 * Nz2;
						counts[16] = Nx2 * Ny3 * Nz2;
						counts[17] = Nx3 * Ny3 * Nz2;

						counts[18] = Nx1 * Ny1 * Nz3;
						counts[19] = Nx2 * Ny1 * Nz3;
						counts[20] = Nx3 * Ny1 * Nz3;
						counts[21] = Nx1 * Ny2 * Nz3;
						counts[22] = Nx2 * Ny2 * Nz3;
						counts[23] = Nx3 * Ny2 * Nz3;
						counts[24] = Nx1 * Ny3 * Nz3;
						counts[25] = Nx2 * Ny3 * Nz3;
						counts[26] = Nx3 * Ny3 * Nz3;

						values[0] = in(xim, yim, zim);
						values[1] = in(xi , yim, zim);
						values[2] = in(xip, yim, zim);
						values[3] = in(xim, yi , zim);
						values[4] = in(xi , yi , zim);
						values[5] = in(xip, yi , zim);
						values[6] = in(xim, yip, zim);
						values[7] = in(xi , yip, zim);
						values[8] = in(xip, yip, zim);

						values[9] = in(xim, yim, zi);
						values[10] = in(xi , yim, zi);
						values[11] = in(xip, yim, zi);
						values[12] = in(xim, yi , zi);
						values[13] = in(xi , yi , zi);
						values[14] = in(xip, yi , zi);
						values[15] = in(xim, yip, zi);
						values[16] = in(xi , yip, zi);
						values[17] = in(xip, yip, zi);

						values[18] = in(xim, yim, zip);
						values[19] = in(xi , yim, zip);
						values[20] = in(xip, yim, zip);
						values[21] = in(xim, yi , zip);
						values[22] = in(xi , yi , zip);
						values[23] = in(xip, yi , zip);
						values[24] = in(xim, yip, zip);
						values[25] = in(xi , yip, zip);
						values[26] = in(xip, yip, zip);


						pixel_t outVal = weightedMedian(values, counts, tmpSpace);

						out(xx, yy, zz) = outVal;
					}
				}

				showThreadProgress(counter, out.depth(), indicateProgress);
			}
		}
	}


	namespace internals
	{
		/**
		Inverse distance interpolation.
		@param p Smoothing exponent.
		*/
		inline Vec3f inverseDistanceInterpolate(const std::vector<Vec3f>& refPoints, const std::vector<Vec3f>& values, const Vec3f& x, float p = 2.5)
		{
			if (refPoints.size() != values.size())
				throw ITLException("refPoints and values lists must have the same size.");

			if (refPoints.size() <= 0)
				return x;

			Vec3f sum = Vec3f(0, 0, 0);
			float wsum = 0;
			for (size_t n = 0; n < refPoints.size(); n++)
			{
				float dist = (x - refPoints[n]).norm();

				if (dist < 1e-7)
				{
					return values[n];
				}

				float w = 1 / ::pow(dist, p);

				sum += w * values[n];
				wsum += w;
			}

			return sum / wsum;
		}
	}

	/**
	Transforms img using a free-form point-to-point transformation.
	@param img Input image.
	@param out Output image. Size of this image must be set to the size of desired output.
	@param outPos Position of the output image relative to the origin of the input image.
	@param refPoints, defPoints List of point pairs. Element of refPoints defines a position in the input image and the corresponding element in defPoints gives the corresponding deformed position.
	*/
	template<typename pixel_t, typename out_t> void genericTransform(const Image<pixel_t>& img,
		Image<out_t>& out,
		const Vec3c& outPos,
		const std::vector<Vec3f>& refPoints,
		const std::vector<Vec3f>& defPoints,
		float exponent = 2.5f,
		const Interpolator<out_t, pixel_t>& interpolate = LinearInterpolator<out_t, pixel_t>(BoundaryCondition::Zero))
	{
		out.mustNotBe(img);

		std::vector<Vec3f> shifts = defPoints - refPoints;

		size_t counter = 0;
		#pragma omp parallel for if (!omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			for (coord_t y = 0; y < out.height(); y++)
			{
				for (coord_t x = 0; x < out.width(); x++)
				{
					Vec3c oix(x, y, z);

					Vec3f ix = Vec3f(oix + outPos);
					Vec3f transformed = ix + internals::inverseDistanceInterpolate(refPoints, shifts, ix, exponent);

					out(oix) = interpolate(img, transformed.x, transformed.y, transformed.z);
				}
			}

			showThreadProgress(counter, out.pixelCount());
		}

	}

	namespace tests
	{
		void scale();
		void translate();
		void binning();
		void genericTransform();
		void scaleLabels();
		void rot90();
		void rotate();
		void reslice();
		void crop();
	}

}
