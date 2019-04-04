#pragma once

#include "transform.h"

#include "command.h"
#include "commandsbase.h"
#include "math/vec3.h"
#include "overlapdistributable.h"
#include "pointprocesscommands.h"

using math::Vec3c;

namespace pilib
{

	template<typename pixel_t> class BinCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	public:
		BinCommand() : TwoImageInputOutputCommand<pixel_t>("bin", "Reduces size of input image by given integer factor. Each output pixel corresponds to average of factor^dimensionality block of pixels in the input image.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "factor", "Binning factor in each coordinate direction. Value 2 makes the output image dimension half of the input image dimension, 3 makes it one third etc.")
			})
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const
		{
			Vec3c binSize = pop<Vec3c>(args);
			
			if (binSize.min() <= 0)
				throw ITLException("Bins size must be positive.");

			binning(in, out, binSize);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
			Vec3c binSize = get<Vec3c>(args[2]);

			if(binSize.min() <= 0)
				throw ITLException("Bins size must be positive.");

			Vec3c outDimensions = in.dimensions().componentwiseDivide(binSize);
			out.ensureSize(outDimensions);

			// distribute in z, use overlap
			return distributor.distribute(this, args, 2, binSize);
		}

		virtual void getCorrespondingBlock(vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 0)
			{
				// Calculate block of input image that corresponds to given block of output image.

				//size_t binSize = get<size_t>(args[2]);

				//readStart *= binSize;
				//readSize *= binSize;
				//writeFilePos *= binSize;
				//writeImPos *= binSize;
				//writeSize *= binSize;

				Vec3c binSize = get<Vec3c>(args[2]);

				readStart = readStart.componentwiseMultiply(binSize);
				readSize = readSize.componentwiseMultiply(binSize);
				writeFilePos = writeFilePos.componentwiseMultiply(binSize);
				writeImPos = writeImPos.componentwiseMultiply(binSize);
				writeSize = writeSize.componentwiseMultiply(binSize);
			}
		}

	};

	template<typename input_t, typename output_t = input_t> class ScaleCommand : public TwoImageInputOutputCommand<input_t, output_t>, public Distributable
	{
	private:
		void adjustOutSize(const Vec3d& factor, const Vec3c& inSize, Vec3c& outSize) const
		{
			for (size_t n = 0; n < outSize.size(); n++)
				outSize[n] = factor[n] > 0 ? math::round(factor[n] * inSize[n]) : outSize[n];
		}
	public:
		ScaleCommand() : TwoImageInputOutputCommand<input_t, output_t>("scale", "Scales input image and places the result into the output image. Set size of output image before calling this command or specify scaling factor as argument. Does not suppress aliasing artifacts when downscaling.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "scaling factor", "Scaling factor in each coordinate direction. If zero, the current size of the output image in that dimension is used.", Vec3d(0, 0, 0)),
				CommandArgument<InterpolationMode>(ParameterDirection::In, "interpolation mode", "Interpolation mode.", InterpolationMode::Linear),
				CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", "Boundary condition.", BoundaryCondition::Nearest),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Specifies origin of current calculation block. This parameter is used internally in distributed processing and should be set to (-1, -1, -1) in normal use.", Distributor::BLOCK_ORIGIN_ARG_TYPE(-1, -1, -1))
			})
		{
		}

		virtual void run(Image<input_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const
		{
			Vec3d factor = pop<Vec3d>(args);
			InterpolationMode ip = pop<InterpolationMode>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE outputOrigin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

			Vec3c inSize = in.dimensions();
			Vec3c outSize = out.dimensions();
			
			Vec3d delta(0, 0, 0);
			if (outputOrigin.x < 0)
			{
				// This is non-distributed call of this command.
				// Adjust output size according to scale factor.
				adjustOutSize(factor, inSize, outSize);
			}
			else
			{
				// This is distributed call of this command (as outputOrigin is set)
				// Do not adjust output size as that has been set already,
				// but calculate delta

				// Shift required in block-wise distributed processing to correct for shift between block origin and image origin in
				// scaled coordinates.
				delta = Vec3d(outputOrigin).componentwiseDivide(factor) - Vec3d(round(Vec3d(outputOrigin).componentwiseDivide(factor)));
			}

			out.ensureSize(outSize);

			shared_ptr<Interpolator<output_t, input_t>> pInterp = createInterpolator<output_t, input_t>(ip, bc);
			scale(in, out, *pInterp, true, factor, delta);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<input_t>& in = *get<DistributedImage<input_t>* >(args[0]);
			DistributedImage<output_t>& out = *get<DistributedImage<output_t>* >(args[1]);
			Vec3d factor = get<Vec3d>(args[2]);

			Vec3c inSize = in.dimensions();
			Vec3c outSize = out.dimensions();

			adjustOutSize(factor, inSize, outSize);

			out.ensureSize(outSize);

			return distributor.distribute(this, args, 2, Vec3c(0, 0, 5), 1, &args);
		}

		virtual void getCorrespondingBlock(vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 0)
			{
				// Calculate block of input image that corresponds to given block of output image.

				Vec3d factor = get<Vec3d>(args[2]);

				readStart = math::round(Vec3d(readStart).componentwiseDivide(factor));
				readSize = math::round(Vec3d(readSize).componentwiseDivide(factor));
			}
		}
	};

	template<typename pixel_t> class CropCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	public:
		CropCommand() : TwoImageInputOutputCommand<pixel_t>("crop", "Crops the image into size of output. Set size of output before calling this command.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position in input image where the top-left corner of the cropped image is placed.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "size", "Size of output image. Specify zeroes or nothing to crop to current size of output image.", Vec3c(0, 0, 0))
			})
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const
		{
			Vec3c pos = pop<Vec3c>(args);
			Vec3c size = pop<Vec3c>(args);

			if (size.x > 0 && size.y > 0 && size.z > 0)
				out.ensureSize(size);

			crop(in, out, pos);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
			Vec3c pos = get<Vec3c>(args[2]);
			Vec3c size = get<Vec3c>(args[3]);

			if (size.x > 0 && size.y > 0 && size.z > 0)
				out.ensureSize(size);

			// The block coordinates are calculated in getCorrespondingBlock so each block is only copied from input to output.
			vector<ParamVariant> newArgs;
			ParamVariant p1, p2, p3, p4;
			p1 = &in;
			p2 = &out;
			p3 = Vec3c(0, 0, 0);
			p4 = Vec3c(0, 0, 0);
			newArgs.push_back(p1);
			newArgs.push_back(p2);
			newArgs.push_back(p3);
			newArgs.push_back(p4);

			return distributor.distribute(this, newArgs, 2, Vec3c(0, 0, 0), 1, &args);
		}

		virtual void getCorrespondingBlock(vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 0)
			{
				// Calculate block of input image that corresponds to given block of output image.

				Vec3c pos = get<Vec3c>(args[2]);
				readStart += pos;
			}
		}
	};


	template<typename input_t> class TranslateCommand : public TwoImageInputOutputCommand<input_t>
	{
	public:
		TranslateCommand() : TwoImageInputOutputCommand<input_t>("translate", "Translates input image by specified amount.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "shift", "Translation that will be applied to the input image."),
			})
		{
		}

		virtual void run(Image<input_t>& in, Image<input_t>& out, vector<ParamVariant>& args) const
		{
			Vec3d shift = pop<Vec3d>(args);
			translate(in, out, shift);
		}
	};


	template<typename pixel_t> class GenericTransformCommand : public Command
	{
	public:
		GenericTransformCommand() : Command("generictransform", "Transforms image based on point-to-point transformation data.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "image", "Image that will be transformed."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "transformed image", "The result of the transformation is set to this image. Size of this image must be set before calling this command."),
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position of the transformed image in coordinates of the original."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "reference points", "Points in the original image as 3xN image where each row contains (x, y, z)-coordinates of a single point, and there are N points in total."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "deformed points", "Locations of points in reference points image after the deformation has been applied. Encoded similarly to reference points image."),
				CommandArgument<double>(ParameterDirection::In, "exponent", "Smoothing exponent. Smaller values smooth more.", 2.5)
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& out = *pop<Image<pixel_t>* >(args);
			Vec3c outPos = pop<Vec3c>(args);
			Image<float32_t>& refPointImg = *pop<Image<float32_t>* >(args);
			Image<float32_t>& defPointImg = *pop<Image<float32_t>* >(args);
			double p = pop<double>(args);

			if (refPointImg.width() != 3 || defPointImg.width() != 3)
				throw ITLException("Reference and deformed point matrices must have the same size.");

			refPointImg.checkSize(defPointImg);

			coord_t N = refPointImg.height();
			vector<Vec3f> refPoints(N), defPoints(N);
			for (coord_t n = 0; n < N; n++)
			{
				refPoints[n] = Vec3f(refPointImg(0, n), refPointImg(1, n), refPointImg(2, n));
				defPoints[n] = Vec3f(defPointImg(0, n), defPointImg(1, n), defPointImg(2, n));
			}

			genericTransform(in, out, outPos, refPoints, defPoints, (float)p);
		}
	};

}
