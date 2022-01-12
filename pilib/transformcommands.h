#pragma once

#include "transform.h"

#include "command.h"
#include "commandsbase.h"
#include "commandlist.h"
#include "math/vec3.h"
#include "overlapdistributable.h"
#include "pointprocesscommands.h"
#include "standardhelp.h"


namespace pilib
{

	inline std::string transformSeeAlso()
	{
		return "rot90cw, rot90ccw, rotate, flip, reslice, crop, copy, scalelabels";
	}

	template<typename pixel_t> class Rotate90CWCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		Rotate90CWCommand() : TwoImageInputOutputCommand<pixel_t>("rot90cw", "Rotates input image clockwise 90 degrees around the $z$-axis.",
			{
			},
			transformSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			rot90cw(in, out);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>*>(args[1]);
			out.ensureSize(in.height(), in.width(), in.depth());
			return distributor.distribute(this, args);
		}

		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 0)
			{
				// readStart and readSize are those for the output image.
				// Convert them to the input coordinates

				DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>*>(args[1]);

				readStart = Vec3c(readStart.y, out.width() - (readStart.x + readSize.x), readStart.z);
				readSize = Vec3c(readSize.y, readSize.x, readSize.z);
			}
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const
		{
			return JobType::Fast;
		}
	};



	template<typename pixel_t> class Rotate90CCWCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		Rotate90CCWCommand() : TwoImageInputOutputCommand<pixel_t>("rot90ccw", "Rotates input image counterclockwise 90 degrees around the $z$-axis.",
			{
			},
			transformSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			rot90ccw(in, out);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>*>(args[1]);
			out.ensureSize(in.height(), in.width(), in.depth());
			return distributor.distribute(this, args);
		}

		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 0)
			{
				// readStart and readSize are those for the output image.
				// Convert them to the input coordinates

				DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>*>(args[1]);

				// x' = -y
				// y' = x
				// z' = z
				readStart = Vec3c(out.height() - (readStart.y + readSize.y), readStart.x, readStart.z);
				readSize = Vec3c(readSize.y, readSize.x, readSize.z);
			}
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const
		{
			return JobType::Fast;
		}
	};




	template<typename pixel_t> class FlipCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		FlipCommand() : OneImageInPlaceCommand<pixel_t>("flip", "Flips image in the given dimension, e.g. if dimension is zero, the left edge of the image becomes the right edge and vice versa.",
			{
				CommandArgument<size_t>(ParameterDirection::In, "dimension", "The dimension to flip. Zero corresponds to $x$, one to $y$, etc."),
			},
			transformSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const override
		{
			size_t dim = pop<size_t>(args);
			flip(img, dim);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const
		{
			return distributor.distribute(this, args);
		}

		virtual size_t getDistributionDirection1(const std::vector<ParamVariant>& args) const
		{
			// Return any other direction than the one we are flipping in.
			// This is to make sure that the scanlines are full in the flipping dimension.
			size_t dim = std::get<size_t>(args[1]);
			switch (dim)
			{
			case 0: return 2;
			case 1: return 2;
			case 2: return 1;
			default: throw ITLException(string("Invalid flipping dimension: ") + itl2::toString(dim));
			}
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const
		{
			return JobType::Fast;
		}
	};





	template<typename pixel_t> class ResliceCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		ResliceCommand() : TwoImageInputOutputCommand<pixel_t>("reslice", "Rotates the input image like a 3D cube such that the face of the cube defined by 'direction' argument will be the first slice in the output image.",
			{
				CommandArgument<string>(ParameterDirection::In, "direction", "The reslice direction. Can be Top, Bottom, Left, or Right."),
			},
			transformSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			string dirs = pop<string>(args);
			ResliceDirection dir = fromString<ResliceDirection>(dirs);

			reslice(in, out, dir);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>*>(args[1]);
			string dirs = std::get<string>(args[2]);
			ResliceDirection dir = fromString<ResliceDirection>(dirs);

			if (dir == ResliceDirection::Top)
				out.ensureSize(in.width(), in.depth(), in.height());
			else if (dir == ResliceDirection::Bottom)
				out.ensureSize(in.width(), in.depth(), in.height());
			else if (dir == ResliceDirection::Left)
				out.ensureSize(in.depth(), in.height(), in.width());
			else if (dir == ResliceDirection::Right)
				out.ensureSize(in.depth(), in.height(), in.width());
			else
				throw ITLException(string("Invalid reslice direction: ") + itl2::toString(dir));

			return distributor.distribute(this, args);
		}

		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 0)
			{
				// readStart and readSize are those for the output image.
				// Convert them to the input coordinates

				DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>*>(args[1]);
				string dirs = std::get<string>(args[2]);
				ResliceDirection dir = fromString<ResliceDirection>(dirs);

				if (dir == ResliceDirection::Top)
				{
					// Top:
					// x' = x
					// y' = -z
					// z' = y

					readStart = Vec3c(readStart.x, readStart.z, out.height() - (readStart.y + readSize.y));
					readSize = Vec3c(readSize.x, readSize.z, readSize.y);
				}
				else if (dir == ResliceDirection::Bottom)
				{
					// Bottom:
					// x' = x
					// y' = z
					// z' = -y

					readStart = Vec3c(readStart.x, out.depth() - (readStart.z + readSize.z), readStart.y);
					readSize = Vec3c(readSize.x, readSize.z, readSize.y);
				}
				else if (dir == ResliceDirection::Left)
				{
					// Left:
					// x' = -z
					// y' = y
					// z' = x

					readStart = Vec3c(readStart.z, readStart.y, out.width() - (readStart.x + readSize.x));
					readSize = Vec3c(readSize.z, readSize.y, readSize.x);
				}
				else if (dir == ResliceDirection::Right)
				{
					// Right:
					// x' = z
					// y' = y
					// z' = -x

					readStart = Vec3c(out.depth() - (readStart.z + readSize.z), readStart.y, readStart.x);
					readSize = Vec3c(readSize.z, readSize.y, readSize.x);
				}
				else
				{
					throw ITLException(string("Invalid reslice direction: ") + itl2::toString(dir));
				}

			}
		}

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const
		{
			return 1;
		}
	};








	template<typename pixel_t> class RotateCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		RotateCommand() : TwoImageInputOutputCommand<pixel_t>("rotate", "Rotates input image around given axis. NOTE: This command does not set the size of the output image automatically. Please set the size of the output to the desired value before calling this command.",
			{
				CommandArgument<double>(ParameterDirection::In, "angle", "Rotation angle in radians."),
				CommandArgument<Vec3d>(ParameterDirection::In, "axis", "Rotation axis. Does not need to be unit vector.", Vec3d(0, 0, 1)),
				CommandArgument<Vec3d>(ParameterDirection::In, "input center", "Rotation center in the input image."),
				CommandArgument<Vec3d>(ParameterDirection::In, "output center", "The rotation center in the input image is mapped to this point in the output image."),
				CommandArgument<InterpolationMode>(ParameterDirection::In, "interpolation mode", string("Interpolation mode. ") + interpolationHelp(), InterpolationMode::Linear),
				CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Zero),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "This argument is used internally in distributed processing. It is assigned the origin of the current calculation block. In normal operation it should be assigned to zero vector.", Distributor::BLOCK_ORIGIN_ARG_TYPE()),
				CommandArgument<Vec3c>(ParameterDirection::In, "full input dimensions", "This argument is used internally in distributed processing. It is assigned the full dimensions of the input image. In normal operation it should be assigned to zero vector.", Distributor::BLOCK_ORIGIN_ARG_TYPE()),
			},
			transformSeeAlso())
		{
		}

	private:

		static Vec3d transformPoint(const Vec3d& outPos, const Matrix3x3<double>& R, const Vec3d& inCenter, const Vec3d& outCenter)
		{
			// Calculate inPos
			return R * (outPos - outCenter) + inCenter;
		}

		static void inputPosAndSize(double angle, const Vec3d& axis, const Vec3c& outputBlockPos, const Vec3c& outputBlockSize, const Vec3c& inputDimensions, const Vec3d& inCenter, const Vec3d& outCenter, Vec3c& start, Vec3c& end)
		{
			// First find locations of all 8 corners of the rotated block, and convert them to the input coordinates.
			Matrix3x3<double> R = Matrix3x3<double>::rotationMatrix(angle, axis);
			R.transpose();

			Vec3d ps[] = {
				transformPoint(Vec3d(Vec3c(outputBlockPos.x,                     outputBlockPos.y,                     outputBlockPos.z)), R, inCenter, outCenter),
				transformPoint(Vec3d(Vec3c(outputBlockPos.x + outputBlockSize.x, outputBlockPos.y,                     outputBlockPos.z)), R, inCenter, outCenter),
				transformPoint(Vec3d(Vec3c(outputBlockPos.x,                     outputBlockPos.y + outputBlockSize.y, outputBlockPos.z)), R, inCenter, outCenter),
				transformPoint(Vec3d(Vec3c(outputBlockPos.x + outputBlockSize.x, outputBlockPos.y + outputBlockSize.y, outputBlockPos.z)), R, inCenter, outCenter),
				transformPoint(Vec3d(Vec3c(outputBlockPos.x,                     outputBlockPos.y,                     outputBlockPos.z + outputBlockSize.z)), R, inCenter, outCenter),
				transformPoint(Vec3d(Vec3c(outputBlockPos.x + outputBlockSize.x, outputBlockPos.y,                     outputBlockPos.z + outputBlockSize.z)), R, inCenter, outCenter),
				transformPoint(Vec3d(Vec3c(outputBlockPos.x,                     outputBlockPos.y + outputBlockSize.y, outputBlockPos.z + outputBlockSize.z)), R, inCenter, outCenter),
				transformPoint(Vec3d(Vec3c(outputBlockPos.x + outputBlockSize.x, outputBlockPos.y + outputBlockSize.y, outputBlockPos.z + outputBlockSize.z)), R, inCenter, outCenter),
			};

			// Find minimum and maximum elementwise coordinates
			Vec3d minp = ps[0];
			Vec3d maxp = ps[0];
			for (size_t n = 0; n < 8; n++)
			{
				minp = min(minp, ps[n]);
				maxp = max(maxp, ps[n]);
			}

			// Add margin for interpolation
			minp -= Vec3d(2, 2, 2);
			maxp += Vec3d(2, 2, 2);

			start = floor(minp);
			end = ceil(maxp);

			Vec3c M = inputDimensions - Vec3c(1, 1, 1);
			clamp(start, Vec3c(0, 0, 0), M);
			clamp(end, Vec3c(0, 0, 0), M);
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			double angle = pop<double>(args);
			Vec3d axis = pop<Vec3d>(args);
			Vec3d inCenter = pop<Vec3d>(args);
			Vec3d outCenter = pop<Vec3d>(args);
			InterpolationMode imode = pop<InterpolationMode>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);
			Vec3c outputOrigin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			Vec3c inputDimensions = pop<Vec3c>(args);

			// blockOrigin is the origin of the calculation block in the output image
			// Get origin of the block in the input image
			Vec3c inputOrigin(0, 0, 0);
			Vec3c inputEnd;
			if(inputDimensions.min() > 0)
				inputPosAndSize(angle, axis, outputOrigin, out.dimensions(), inputDimensions, inCenter, outCenter, inputOrigin, inputEnd);

			itl2::rotate<pixel_t, pixel_t>(in, out, angle, axis, inCenter - Vec3d(inputOrigin), outCenter - Vec3d(outputOrigin), *createInterpolator<pixel_t, pixel_t>(imode, bc));
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);
			args[9] = in.dimensions();
			return distributor.distribute(this, args);
		}

		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const
		{
			if (argIndex == 0)
			{
				// readStart and readSize are those for the output image.
				// Convert them to the input coordinates

				DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);

				double angle = std::get<double>(args[2]);
				Vec3d axis = std::get<Vec3d>(args[3]);
				Vec3d inCenter = std::get<Vec3d>(args[4]);
				Vec3d outCenter = std::get<Vec3d>(args[5]);

				Vec3c start, end;
				inputPosAndSize(angle, axis, readStart, readSize, in.dimensions(), inCenter, outCenter, start, end);

				// Assign read start and size
				readStart = start;
				readSize = end - start + Vec3c(1, 1, 1);
			}
		}

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const
		{
			return 1;
		}

	};


	template<typename pixel_t> class Rotate2Command : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		Rotate2Command() : TwoImageInputOutputCommand<pixel_t>("rotate", "Rotates input image around given axis. Rotation center is in the center of the input image, and it is mapped to the center of the output image. NOTE: This command does not set the size of the output image automatically. Please set the size of the output to the desired value before calling this command.",
			{
				CommandArgument<double>(ParameterDirection::In, "angle", "Rotation angle in radians."),
				CommandArgument<Vec3d>(ParameterDirection::In, "axis", "Rotation axis. Does not need to be unit vector.", Vec3d(0, 0, 1)),
				CommandArgument<InterpolationMode>(ParameterDirection::In, "interpolation mode", string("Interpolation mode. ") + interpolationHelp(), InterpolationMode::Linear),
				CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Zero),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "This argument is used internally in distributed processing. It is assigned the origin of the current calculation block. In normal operation it should be assigned to zero vector.", Distributor::BLOCK_ORIGIN_ARG_TYPE()),
				CommandArgument<Vec3c>(ParameterDirection::In, "full input dimensions", "This argument is used internally in distributed processing. It is assigned the full dimensions of the input image. In normal operation it should be assigned to zero vector.", Distributor::BLOCK_ORIGIN_ARG_TYPE()),
			},
			transformSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			args.insert(args.begin() + 2, Vec3d(in.dimensions() / 2));
			args.insert(args.begin() + 3, Vec3d(out.dimensions() / 2));

			CommandList::get<RotateCommand<pixel_t>>().run(in, out, args);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>*>(args[1]);

			args.insert(args.begin() + 2 + 2, Vec3d(in.dimensions() / 2));
			args.insert(args.begin() + 2 + 3, Vec3d(out.dimensions() / 2));

			return CommandList::get<RotateCommand<pixel_t>>().runDistributed(distributor, args);
		}
	};



	inline std::string scaleSeeAlso()
	{
		return "scale, bin, maskedbin, scalelabels";
	}


	template<typename pixel_t> class MaskedBinCommand : public TwoImageInputOutputCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		MaskedBinCommand() : TwoImageInputOutputCommand<pixel_t>("maskedbin", "Reduces size of input image by given integer factor. Each output pixel corresponds to factor^dimensionality block of pixels in the input image. Supports treating one value in the input image as 'bad' such that pixels having that value do not contribute to the output at all.",
			{
				CommandArgument<size_t>(ParameterDirection::In, "factor", "Binning factor in each coordinate direction. Value 2 makes the output image dimension half of the input image dimension, 3 makes them one third etc."),
				CommandArgument<double>(ParameterDirection::In, "bad value", "Value that should not be considered in the averaging calculations.", 0.0),
				CommandArgument<double>(ParameterDirection::In, "undefined value", " Value that is placed to those pixels of the output image that do not correspond to any valid pixels in the input image.", 0.0),
			},
			scaleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			size_t binSize = pop<size_t>(args);
			pixel_t badValue = pixelRound<pixel_t>(pop<double>(args));
			pixel_t undefinedValue = pixelRound<pixel_t>(pop<double>(args));

			if (binSize <= 0)
				throw ITLException("Bins size must be positive.");

			maskedBinning<pixel_t>(in, out, binSize, badValue, undefinedValue);
		}
	};

	template<typename pixel_t> class BinCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		BinCommand() : TwoImageInputOutputCommand<pixel_t>("bin", "Reduces size of input image by given integer factor. Each output pixel corresponds to factor^dimensionality block of pixels in the input image.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "factor", "Binning factor in each coordinate direction. Value 2 makes the output image dimension half of the input image dimension, 3 makes them one third etc."),
				CommandArgument<string>(ParameterDirection::In, "binning type", "Name of binning type to be performed. Currently 'mean', 'min' and 'max' are supported.", "mean")
			},
			scaleSeeAlso())
		{
		}

	private:

		enum class BinType
		{
			Mean,
			Min,
			Max
		};

		static BinType fromString(string s)
		{
			toLower(s);
			if (s == "mean" || s == "average" || s == "avg" || s == "normal")
				return BinType::Mean;
			if (s == "max" || s == "maximum")
				return BinType::Max;
			if (s == "min" || s == "minimum")
				return BinType::Min;

			throw ITLException(string("Invalid binning type: ") + s);
		}

		static string toString(BinType b)
		{
			switch (b)
			{
			case BinType::Mean: return "mean";
			case BinType::Max: return "max";
			case BinType::Min: return "min";
			default: throw ITLException("Unsupported binning type.");
			}
		}


	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			Vec3c binSize = pop<Vec3c>(args);
			string binType = pop<string>(args);
			
			if (binSize.min() <= 0)
				throw ITLException("Bins size must be positive.");

			BinType type = fromString(binType);
			switch (type)
			{
			case BinType::Mean:
				binning<pixel_t, pixel_t, binningop::mean<pixel_t, pixel_t> >(in, out, binSize);
				break;
			case BinType::Max:
				binning<pixel_t, pixel_t, binningop::max<pixel_t, pixel_t> >(in, out, binSize);
				break;
			case BinType::Min:
				binning<pixel_t, pixel_t, binningop::min<pixel_t, pixel_t> >(in, out, binSize);
				break;
			default:
				throw ITLException(string("Unsupported binning type: ") + toString(type));
			}
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = * std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = * std::get<DistributedImage<pixel_t>* >(args[1]);
			Vec3c binSize =  std::get<Vec3c>(args[2]);
			string binType =  std::get<string>(args[3]);

			in.mustNotBe(out);

			if(binSize.min() <= 0)
				throw ITLException("Bin size must be positive.");

			// Test that binning type is ok.
			fromString(binType);

			Vec3c outDimensions = in.dimensions().componentwiseDivide(binSize);
			out.ensureSize(outDimensions);

			return distributor.distribute(this, args);
		}

		virtual Vec3c getMargin(const vector<ParamVariant>& args) const override
		{
			// NOTE: Margin is given in output image coordinates
			return Vec3c(1, 1, 1);
			//Vec3c binSize =  std::get<Vec3c>(args[2]);
			//return binSize;
		}

		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 0)
			{
				// Calculate block of input image that corresponds to given block of output image.

				Vec3c binSize =  std::get<Vec3c>(args[2]);

				readStart = readStart.componentwiseMultiply(binSize);
				readSize = readSize.componentwiseMultiply(binSize);
				writeFilePos = writeFilePos.componentwiseMultiply(binSize);
				writeImPos = writeImPos.componentwiseMultiply(binSize);
				writeSize = writeSize.componentwiseMultiply(binSize);
			}
		}

		virtual double calculateExtraMemory(const vector<ParamVariant>& args) const override
		{
			// Calculate expected output size and if current output image size is not correct, calculate amount of extra memory needed to create it.
			DistributedImage<pixel_t>& in = * std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = * std::get<DistributedImage<pixel_t>* >(args[1]);
			Vec3c binSize =  std::get<Vec3c>(args[2]);

			Vec3c expectedOutSize = in.dimensions().componentwiseDivide(binSize);
			if (out.dimensions() != expectedOutSize)
			{
				double currentTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + out.pixelCount() * out.pixelSize());
				double trueTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + expectedOutSize.x * expectedOutSize.y * expectedOutSize.z * out.pixelSize());
				return std::max(0.0, trueTotalDataSize / currentTotalDataSize - 1.0);
			}

			return 0.0;
		}

	};

	template<typename input_t, typename output_t = input_t> class ScaleCommand : public TwoImageInputOutputCommand<input_t, output_t>, public Distributable
	{
	private:
		void adjustOutSize(const Vec3d& factor, const Vec3c& inSize, Vec3c& outSize) const
		{
			for (size_t n = 0; n < outSize.size(); n++)
				outSize[n] = factor[n] > 0 ? itl2::round(factor[n] * inSize[n]) : outSize[n];
		}


		void adjustFactor(Vec3d& factor, const Vec3c& inSize, const Vec3c& outSize) const
		{
			for (size_t n = 0; n < outSize.size(); n++)
				if (factor[n] <= 0)
					factor[n] = (double)outSize[n] / (double)inSize[n];
		}

	protected:
		friend class CommandList;

		ScaleCommand() : TwoImageInputOutputCommand<input_t, output_t>("scale", "Scales input image and places the result into the output image. Set size of output image before calling this command or specify scaling factor as an argument. Does not suppress aliasing artifacts when downscaling unless average when downsizing-parameter is set to true.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "scaling factor", "Scaling factor in each coordinate direction. If zero in some dimension, the factor is calculated from current size of the output image and the input image in that dimension.", Vec3d(0, 0, 0)),
				CommandArgument<bool>(ParameterDirection::In, "average when downsizing", "Set to true to average when downsizing.", false),
				CommandArgument<InterpolationMode>(ParameterDirection::In, "interpolation mode", string("Interpolation mode. ") + interpolationHelp(), InterpolationMode::Linear),
				CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Specifies origin of current calculation block. This parameter is used internally in distributed processing and should be set to (-1, -1, -1) in normal use.", Distributor::BLOCK_ORIGIN_ARG_TYPE(-1, -1, -1))
			},
			scaleSeeAlso())
		{
		}

	public:
		virtual void run(Image<input_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const override
		{
			Vec3d factor = pop<Vec3d>(args);
			bool averageWhenDownSizing = pop<bool>(args);
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
				// but calculate delta.

				// Shift required in block-wise distributed processing to correct for shift between block origin and image origin in
				// scaled coordinates.
				adjustFactor(factor, inSize, outSize);
				delta = Vec3d(outputOrigin).componentwiseDivide(factor) - Vec3d(floor(Vec3d(outputOrigin).componentwiseDivide(factor)));
			}

			out.ensureSize(outSize);

			scale(in, out, averageWhenDownSizing, ip, bc, true, factor, delta);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<input_t>& in = * std::get<DistributedImage<input_t>* >(args[0]);
			DistributedImage<output_t>& out = * std::get<DistributedImage<output_t>* >(args[1]);
			Vec3d factor =  std::get<Vec3d>(args[2]);

			in.mustNotBe(out);

			Vec3c inSize = in.dimensions();
			Vec3c outSize = out.dimensions();

			adjustOutSize(factor, inSize, outSize);

			out.ensureSize(outSize);

			return distributor.distribute(this, args);
		}

		virtual Vec3c getMargin(const vector<ParamVariant>& args) const override
		{
			return Vec3c(0, 0, 5);
		}

		virtual size_t getRefIndex(const vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 0)
			{
				// Calculate block of input image that corresponds to given block of output image.

				DistributedImage<input_t>& in = * std::get<DistributedImage<input_t>* >(args[0]);
				DistributedImage<output_t>& out = * std::get<DistributedImage<output_t>* >(args[1]);
				Vec3d factor =  std::get<Vec3d>(args[2]);

				adjustFactor(factor, in.dimensions(), out.dimensions());

				readStart = floor(Vec3d(readStart).componentwiseDivide(factor));
				readSize = ceil(Vec3d(readSize).componentwiseDivide(factor));
				Vec3c dims = in.dimensions();
				for (size_t n = 0; n < 3; n++)
				{
					if (readStart[n] + readSize[n] > dims[n])
						readSize[n] = dims[n] - readStart[n];
				}
			}
		}

		virtual double calculateExtraMemory(const vector<ParamVariant>& args) const override
		{
			// Calculate expected output size and if current output image size is not correct, calculate amount of extra memory needed to create it.
			DistributedImage<input_t>& in = * std::get<DistributedImage<input_t>* >(args[0]);
			DistributedImage<output_t>& out = * std::get<DistributedImage<output_t>* >(args[1]);
			Vec3d factor =  std::get<Vec3d>(args[2]);
			bool averageWhenDownSizing =  std::get<bool>(args[3]);

			Vec3c inSize = in.dimensions();
			Vec3c outSize = out.dimensions();

			adjustOutSize(factor, inSize, outSize);

			double currentTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + out.pixelCount() * out.pixelSize());
			double trueTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + outSize.x * outSize.y * outSize.z * out.pixelSize());
			if (averageWhenDownSizing &&
				(factor.x < 1 || factor.y < 1 || factor.z < 1))
				trueTotalDataSize += (double)(in.pixelCount() * sizeof(float32_t));
			return std::max(0.0, trueTotalDataSize / currentTotalDataSize - 1.0);

			//if (out.dimensions() != outSize)
			//{
			//	double currentTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + out.pixelCount() * out.pixelSize());
			//	double trueTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + outSize.x * outSize.y * outSize.z * out.pixelSize());
			//	return std::max(0.0, trueTotalDataSize / currentTotalDataSize - 1.0);
			//}

			//return 0.0;
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Fast;
		}
	};



	template<typename pixel_t> class ScaleLabelsCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	private:
		void adjustOutSize(const Vec3c& factor, const Vec3c& inSize, Vec3c& outSize) const
		{
			for (size_t n = 0; n < outSize.size(); n++)
				outSize[n] = factor[n] > 0 ? factor[n] * inSize[n] : outSize[n];
		}

	protected:
		friend class CommandList;

		ScaleLabelsCommand() : TwoImageInputOutputCommand<pixel_t>("scalelabels", "Scales input image and places the result into the output image. Assumes that input image is a binary image or contains label values, and performs scaling such that edges of regions do not become jagged. Supports only upscaling by integer scaling factor. Set size of output image before calling this command or specify scaling factor as an argument.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "scaling factor", "Scaling factor in each coordinate direction. If zero in some dimension, the factor is calculated from current size of the output image and the input image in that dimension.", Vec3c(0, 0, 0)),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Specifies origin of current calculation block. This parameter is used internally in distributed processing and should be set to (-1, -1, -1) in normal use.", Distributor::BLOCK_ORIGIN_ARG_TYPE(-1, -1, -1))
			},
			scaleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			Vec3c factor = pop<Vec3c>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE outputOrigin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

			Vec3c inSize = in.dimensions();
			Vec3c outSize = out.dimensions();

			Vec3d delta(0, 0, 0);
			if (outputOrigin.x < 0)
			{
				// This is non-distributed call of this command.
				// Adjust output size according to scale factor.
				adjustOutSize(factor, inSize, outSize);
				outputOrigin = Vec3c(0, 0, 0);
			}
			else
			{
				// This is distributed call of this command (as outputOrigin is set)
				// Do not adjust output size as that has been set already.
			}

			out.ensureSize(outSize);

			scaleLabels(in, out, true, factor, outputOrigin);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = * std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = * std::get<DistributedImage<pixel_t>* >(args[1]);
			Vec3c factor = std::get<Vec3c>(args[2]);

			in.mustNotBe(out);

			Vec3c inSize = in.dimensions();
			Vec3c outSize = out.dimensions();

			adjustOutSize(factor, inSize, outSize);

			out.ensureSize(outSize);

			return distributor.distribute(this, args);
		}

		virtual Vec3c getMargin(const vector<ParamVariant>& args) const override
		{
			return Vec3c(0, 0, 6);
		}

		virtual size_t getRefIndex(const vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 0)
			{
				// Calculate block of input image that corresponds to given block of output image.

				DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
				Vec3c factor = std::get<Vec3c>(args[2]);

				readStart = readStart.componentwiseDivide(factor);
				readSize = readSize.componentwiseDivide(factor) + Vec3c(1, 1, 1); // Add 1 in all directions to account for possible rounding inaccuracies
				Vec3c readEnd = readStart + readSize;
				clamp(readEnd, Vec3c(0, 0, 0), in.dimensions());
				readSize = readEnd - readStart;

			}
		}

		virtual double calculateExtraMemory(const vector<ParamVariant>& args) const override
		{
			// Calculate expected output size and if current output image size is not correct, calculate amount of extra memory needed to create it.
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			Vec3c factor = std::get<Vec3c>(args[2]);

			Vec3c inSize = in.dimensions();
			Vec3c outSize = out.dimensions();

			adjustOutSize(factor, inSize, outSize);

			if (out.dimensions() != outSize)
			{
				double currentTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + out.pixelCount() * out.pixelSize());
				double trueTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + outSize.x * outSize.y * outSize.z * out.pixelSize());
				return std::max(0.0, trueTotalDataSize / currentTotalDataSize - 1.0);
			}

			return 0.0;
		}
	};




	template<typename pixel_t> class CropCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		CropCommand() : TwoImageInputOutputCommand<pixel_t>("crop", "Crops the input image. If the size parameter is set to zero, crops to current size of the output image. NOTE: Input and output images must have the same data type. If not, output image is converted to the correct type.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position in input image where the top-left corner of the cropped image is placed.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "size", "Size of output image. Specify zeroes or nothing to crop to the current size of the output image.", Vec3c(0, 0, 0))
			},
			scaleSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			Vec3c pos = pop<Vec3c>(args);
			Vec3c size = pop<Vec3c>(args);

			if (size.x > 0 && size.y > 0 && size.z > 0)
				out.ensureSize(size);

			crop(in, out, pos);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			Vec3c pos = std::get<Vec3c>(args[2]);
			Vec3c size = std::get<Vec3c>(args[3]);

			in.mustNotBe(out);

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

			return distributor.distribute(this, newArgs, &args);
		}

		virtual size_t getRefIndex(const vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 0)
			{
				// Calculate block of input image that corresponds to given block of output image.

				Vec3c pos = std::get<Vec3c>(args[2]);
				readStart += pos;
			}
		}

		virtual double calculateExtraMemory(const vector<ParamVariant>& args) const override
		{
			// Calculate expected output size and if current output image size is not correct, calculate amount of extra memory needed to create it.
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);

			Vec3c size = std::get<Vec3c>(args[3]);

			if (size.x > 0 && size.y > 0 && size.z > 0)
			{
				if (out.dimensions() != size)
				{
					double currentTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + out.pixelCount() * out.pixelSize());
					double trueTotalDataSize = (double)(in.pixelCount() * in.pixelSize() + size.x * size.y * size.z * out.pixelSize());
					return std::max(0.0, trueTotalDataSize / currentTotalDataSize - 1.0);
				}
			}

			return 0.0;
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Fast;
		}
	};


	template<typename pixel_t> class Copy2Command : public Command
	{
	protected:
		friend class CommandList;

		Copy2Command() : Command("copy", "Copies pixel values from source image to target image to specified location. The size of the target image is not changed and out-of-bounds pixels are not copied. See also `set` command.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "source image", "Image that is copied to the target image."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "target image", "Image whose values are set."),
				CommandArgument<Vec3c>(ParameterDirection::In, "location", "Location where the target image is placed in the source image.")
			},
			scaleSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& source = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& target = *pop<Image<pixel_t>* >(args);
			Vec3c pos = pop<Vec3c>(args);

			itl2::copyValues(target, source, pos);
		}
	};




	template<typename input_t> class TranslateCommand : public TwoImageInputOutputCommand<input_t>
	{
	protected:
		friend class CommandList;

		TranslateCommand() : TwoImageInputOutputCommand<input_t>("translate", "Translates input image by specified amount.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "shift", "Translation that will be applied to the input image."),
			},
			transformSeeAlso())
		{
		}

	public:
		virtual void run(Image<input_t>& in, Image<input_t>& out, vector<ParamVariant>& args) const override
		{
			Vec3d shift = pop<Vec3d>(args);
			translate(in, out, shift);
		}
	};


	template<typename pixel_t> class GenericTransformCommand : public Command
	{
	protected:
		friend class CommandList;

		GenericTransformCommand() : Command("generictransform", "Transforms image based on point-to-point correspondence data. Transformation between the points is interpolated from the point data using inverse distance interpolation.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "image", "Image that will be transformed."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "transformed image", "The result of the transformation is set to this image. Size of this image must be set before calling this command."),
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position of the transformed image in coordinates of the original."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "reference points", "Points in the original image as 3xN image where each row contains (x, y, z)-coordinates of a single point, and there are N points in total."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "deformed points", "Locations of points in reference points image after the deformation has been applied. Encoded similarly to reference points image."),
				CommandArgument<double>(ParameterDirection::In, "exponent", "Smoothing exponent in the inverse distance interpolation. Smaller values smooth more.", 2.5)
			})
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& out = *pop<Image<pixel_t>* >(args);
			Vec3c outPos = pop<Vec3c>(args);
			Image<float32_t>& refPointImg = *pop<Image<float32_t>* >(args);
			Image<float32_t>& defPointImg = *pop<Image<float32_t>* >(args);
			double p = pop<double>(args);

			in.mustNotBe(out);

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
