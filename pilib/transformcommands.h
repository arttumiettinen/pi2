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
				CommandArgument<size_t>(In, "factor", "Binning factor. Value 2 makes the output image dimensions half of the input image dimensions, 3 makes them one third etc.")
			})
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const
		{
			size_t amount = pop<size_t>(args);
			binning(in, out, amount);
		}

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
			size_t binSize = get<size_t>(args[2]);

			if (binSize <= 0)
				throw ITLException("Bins size must be positive.");

			out.ensureSize(in.dimensions() / binSize);

			// distribute in z, use overlap
			distributor.distribute(this, args, 2, Vec3c(binSize, binSize, binSize));
		}

		virtual void getCorrespondingBlock(vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 0)
			{
				// Calculate block of input image that corresponds to given block of output image.

				size_t binSize = get<size_t>(args[2]);

				readStart *= binSize;
				readSize *= binSize;
				writeFilePos *= binSize;
				writeImPos *= binSize;
				writeSize *= binSize;
			}
		}

	};

	template<typename input_t, typename output_t> class ScaleCommand : public TwoImageInputOutputCommand<input_t, output_t>
	{
	public:
		ScaleCommand() : TwoImageInputOutputCommand<input_t, output_t>("scale", "Scales input image into the output image. Set size of output image before calling this command. Does not suppress aliasing artifacts when downscaling.",
			{
			})
		{
		}

		virtual void run(Image<input_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const
		{
			scale(in, out);
		}
	};

	template<typename pixel_t> class CropCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	public:
		CropCommand() : TwoImageInputOutputCommand<pixel_t>("crop", "Crops the image into size of output. Set size of output before calling this command.",
			{
				CommandArgument<Vec3c>(In, "position", "Position in input image where the top-left corner of the cropped image is placed.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(In, "size", "Size of output image. Specify zeroes or nothing to crop to current size of output image.", Vec3c(0, 0, 0))
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

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
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
			p1.dimgval = &in;
			p2.dimgval = &out;
			p3.vix = 0;
			p3.viy = 0;
			p3.viz = 0;
			p4.vix = 0;
			p4.viy = 0;
			p4.viz = 0;
			newArgs.push_back(p1);
			newArgs.push_back(p2);
			newArgs.push_back(p3);
			newArgs.push_back(p4);

			// Distribute in z, use overlap
			distributor.distribute(this, newArgs, 2, Vec3c(0, 0, 0), 0, 1, &args);
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

	template<typename pixel_t> class GenericTransformCommand : public Command
	{
	public:
		GenericTransformCommand() : Command("generictransform", "Transforms image based on point-to-point transformation data.",
			{
				CommandArgument<Image<pixel_t> >(In, "image", "Image that will be transformed."),
				CommandArgument<Image<pixel_t> >(Out, "transformed image", "The result of the transformation is set to this image. Size of this image must be set before calling this command."),
				CommandArgument<Vec3c>(In, "position", "Position of the transformed image in coordinates of the original."),
				CommandArgument<Image<float32_t> >(In, "reference points", "Points in the original image as 3xN image where each row contains (x, y, z)-coordinates of a single point, and there are N points in total."),
				CommandArgument<Image<float32_t> >(In, "deformed points", "Locations of points in reference points image after the deformation has been applied. Encoded similarly to reference points image."),
				CommandArgument<double>(In, "exponent", "Smoothing exponent. Smaller values smooth more.", 2.5)
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
