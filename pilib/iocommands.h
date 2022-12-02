#pragma once

#include "command.h"
#include "trivialdistributable.h"
#include "parseexception.h"
#include "io/io.h"
#include "commandlist.h"
#include "standardhelp.h"
#include "timing.h"

#include <vector>

using namespace itl2;

namespace pilib
{
	class IsImageFileCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		IsImageFileCommand() : Command("isimagefile", "Checks if a file with given name is an image file.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name of image file."),
				CommandArgument<Image<uint8_t> >(ParameterDirection::In, "result", "Set to 1 if the given file name is a readable image file, and to 0 otherwise.")
			},
			"",
			"In Python/pi2py2, the result parameter is not specified, but the test result is returned as a boolean return value.")
		{

		}

	public:

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			std::string fname = pop<std::string>(args);
			Image<uint8_t>& result = *pop<Image<uint8_t>*>(args);

			Vec3c dims;
			ImageDataType dt;
			std::string reason;
			result(0) = io::getInfo(fname, dims, dt, reason);
		}
	};

	class FileInfoCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		FileInfoCommand() : Command("fileinfo", "Reads metadata of given image file.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name of image file."),
				CommandArgument<Image<uint32_t>>(ParameterDirection::In, "info", "At output, this image will contain four pixels corresponding to width, height, depth, and data type index of the image. If the image has only zero pixels, the given file is not an image that can be read into the pi system. The data type index one correspond to uint8, two to uint16, etc. in this sequence: uint8, uint16, uint32, uint64, float32, complex32, int8, int16, int32, int64.")
			})
		{

		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			std::string fname = pop<std::string>(args);
			Image<uint32_t>& results = *pop<Image<uint32_t>*>(args);

			results.ensureSize(4);
			setValue(results, 0);

			Vec3c dims;
			ImageDataType dt;
			std::string reason;
			if (io::getInfo(fname, dims, dt, reason))
			{
				results(0) = (uint32_t)dims[0];
				results(1) = (uint32_t)dims[1];
				results(2) = (uint32_t)dims[2];
				results(3) = (uint32_t)dt;
			}
		}
	};


	class ShowFileInfoCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		ShowFileInfoCommand() : Command("fileinfo", "Displays metadata of given image file.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name of image file.")
			})
		{

		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			std::string fname = pop<std::string>(args);

			Vec3c dims;
			ImageDataType dt;
			std::string reason;
			if (io::getInfo(fname, dims, dt, reason))
			{
				std::cout << dims << std::endl;
				std::cout << itl2::toString(dt) << std::endl;
			}
			else
			{
				std::cout << "Not an image file: " << fname << std::endl;
				std::cout << reason << std::endl;
			}
		}
	};

	class ShowRawInfoCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		ShowRawInfoCommand() : Command("rawinfo", "Displays metadata of .raw file.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name of .raw file.")
			})
		{

		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			std::string fname = pop<std::string>(args);

			Vec3c dims;
			ImageDataType dt;
			std::string reason;
			if (raw::getInfo(fname, dims, dt, reason))
			{
				std::cout << dims << std::endl;
				std::cout << itl2::toString(dt) << std::endl;
			}
			else
			{
				std::cout << "Not a .raw file: " << fname << std::endl;
				std::cout << reason << std::endl;
			}
		}
	};

	class ShowSequenceInfoCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		ShowSequenceInfoCommand() : Command("sequenceinfo", "Displays metadata of image sequence.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "filename template", "Filename template corresponding to the sequence, as described in readsequence command documentation.")
			})
		{

		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			std::string fname = pop<std::string>(args);

			Vec3c dims;
			ImageDataType dt;
			std::string reason;
			if (sequence::getInfo(fname, dims, dt, reason))
			{
				std::cout << dims << std::endl;
				std::cout << itl2::toString(dt) << std::endl;
			}
			else
			{
				std::cout << "Unable to parse sequence: " << fname << std::endl;
				std::cout << reason << std::endl;
			}
		}
	};

	template<typename pixel_t> class NopSingleImageCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		NopSingleImageCommand() : Command("nop", "Does nothing. This command is used internally as a helper that tags image that must be saved but not processed.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::InOut, "image", "Name of image that is not processed.")
			})
		{

		}

	public:
		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(std::vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Fast;
		}
	};

	template<typename pixel_t> class WriteTiffCommand : public Command
	{
	protected:
		friend class CommandList;

		WriteTiffCommand() : Command("writetif", "Write an image to a .tif file.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the file to write. If the file exists, its current contents are erased. Extension .tif is automatically appended to the name of the file.")
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);

			itl2::tiff::writed(in, fname);
		}
	};


	template<typename pixel_t> class WritePngCommand : public Command
	{
	protected:
		friend class CommandList;

		WritePngCommand() : Command("writepng", "Write an image to a .png file. This command supports 1- and 2-dimensional images only.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the file to write. If the file exists, its current contents are erased. Extension .png is automatically appended to the name of the file.")
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);

			if (in.depth() > 1)
				throw ITLException("3-dimensional images cannot be saved to a .png file. Consider saving to an image sequence instead.");

			itl2::png::writed(in, fname);
		}
	};


	template<typename pixel_t> class WriteNRRDCommand : public Command
	{
	protected:
		friend class CommandList;

		WriteNRRDCommand() : Command("writenrrd", "Write an image to an .nrrd file.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the file to write. If the file exists, its current contents are erased. Extension .nrrd is automatically appended to the name of the file.")
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);

			itl2::nrrd::writed(in, fname);
		}
	};

	template<typename pixel_t> class WriteLZ4Command : public Command
	{
	protected:
		friend class CommandList;

		WriteLZ4Command() : Command("writelz4", "Write an image to an .lz4raw file.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the file to write. If the file exists, its current contents are erased. Extension .lz4raw is automatically appended to the name of the file.")
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);

			itl2::lz4::writed(in, fname);
		}
	};

	template<typename pixel_t> class WriteNN5Command : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		WriteNN5Command() : Command("writenn5", "Write an image to an .nn5 dataset.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "path", "Name (and path) of the dataset to write. If the dataset exists, its current contents are erased."),
				CommandArgument<Vec3c>(ParameterDirection::In, "chunk size", "Chunk size for the NN5 dataset to be written.", nn5::DEFAULT_CHUNK_SIZE)
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);
			Vec3c chunkSize = pop<Vec3c>(args);

			itl2::nn5::write(in, fname, chunkSize);
		}

		using Distributable::runDistributed;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			distributor.flush();

			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);
			Vec3c chunkSize = pop<Vec3c>(args);
			
			//if (in.isSavedToDisk() && in.isRaw())
			if (in.isSavedToDisk() && in.currentReadSourceType() == DistributedImageStorageType::NN5)
			{
				// Effectively copy data from input image to output image.
				// First move input image current data source file to output file.
				// Then point input image current data source file to the copied file.
				// This way no data needs to be copied (if all the data is stored on single partition)
				// and input image can still be used as further changes are anyway saved to a temp file.
				// NOTE: This does not work if the input and output image are stored in "cloud", but
				// that is not supported at the moment anyway.

				if (in.isSavedToTemp())
				{
					if (fs::exists(in.currentReadSource()))
					{
						// The image has been saved to a temporary file
						// Just move the temporary file to new location (and name) and sets read source to that file.
						moveFile(in.currentReadSource(), fname);
						in.setReadSource(fname, false);
					}
					else
					{
						// The input image is empty, so just create an empty file.
						setFileSize(fname, in.dimensions().x * in.dimensions().y * in.dimensions().z * sizeof(pixel_t));
					}
				}
				else
				{
					if (fs::exists(in.currentReadSource()))
					{
						// The image has been saved to a non-temporary file, so we cannot just move it.
						// The file must be copied.
						copyFile(in.currentReadSource(), fname, true);
					}
					else
					{
						// The input image has been saved to a non-temporary file but that file does not exist.
						// This is impossible situation: the source data is not available anymore.
						throw ITLException("The image " + in.varName() + " is loaded from " + in.currentReadSource() + " but that file does not exist anymore.");
					}
				}
			}
			else
			{
				// The input is not stored as .nn5 dataset so actual copying must be made.

				auto oldStorageType = in.currentWriteTargetType();
				in.setWriteTarget(fname, DistributedImageStorageType::NN5);

				std::vector<ParamVariant> args2;
				ParamVariant p;
				p = &in;
				args2.push_back(p);
				NopSingleImageCommand<pixel_t>& cmd = CommandList::get<NopSingleImageCommand<pixel_t> >();
				std::vector<std::string> result = cmd.runDistributed(distributor, args2);

				// Ensure that following operations do not write to the same file.
				// Retain old write target storage type.
				in.newWriteTarget(oldStorageType);
				return result;
			}

			return std::vector<std::string>();
		}
	};


	template<typename pixel_t> class WriteRawCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		WriteRawCommand() : Command("writeraw", "Write an image to .raw file. The dimensions of the image are automatically appended to the file name. If distributed processing is enabled, uses optimized implementation that does not need to copy or write any image data if the image is saved to a temporary image storage location and that location is on the same partition than the file being written.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the file to write. If the file exists, its current contents are erased (unless append parameter is set to true)."),
				CommandArgument<bool>(ParameterDirection::In, "append", "Set to true to append to existing .raw file. This parameter must be set to false in distributed mode.", false)
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);
			bool append = pop<bool>(args);

			raw::writed(in, fname, !append);
		}

		using Distributable::runDistributed;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			distributor.flush();

			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);
			bool append = pop<bool>(args);

			if (append)
				throw ITLException("Unable to append in distributed mode.");


			fname = concatDimensions(fname, in.dimensions());

			//if (in.isSavedToDisk() && in.isRaw())
			if (in.isSavedToDisk() && in.currentReadSourceType() == DistributedImageStorageType::Raw)
			{
				// Effectively copy data from input image to output image.
				// First move input image current data source file to output file.
				// Then point input image current data source file to the copied file.
				// This way no data needs to be copied (if all the data is stored on single partition)
				// and input image can still be used as further changes are anyway saved to a temp file.
				// NOTE: This does not work if the input and output image are stored in "cloud", but
				// that is not supported at the moment anyway.

				if (in.isSavedToTemp())
				{
					if (fs::exists(in.currentReadSource()))
					{
						// The image has been saved to a temporary file
						// Just move the temporary file to new location (and name) and sets read source to that file.
						moveFile(in.currentReadSource(), fname);
						in.setReadSource(fname, false);
					}
					else
					{
						// The input image is empty, so just create an empty file.
						setFileSize(fname, in.dimensions().x * in.dimensions().y * in.dimensions().z * sizeof(pixel_t));
					}
				}
				else
				{
					if (fs::exists(in.currentReadSource()))
					{
						// The image has been saved to a non-temporary file, so we cannot just move it.
						// The file must be copied.
						copyFile(in.currentReadSource(), fname, true);
					}
					else
					{
						// The input image has been saved to a non-temporary file but that file does not exist.
						// This is impossible situation: the source data is not available anymore.
						throw ITLException("The image " + in.varName() + " is loaded from " + in.currentReadSource() + " but that file does not exist anymore.");
					}
				}
			}
			else
			{
				// The input is not stored as a .raw file so actual copying must be made.
				
				auto oldStorageType = in.currentWriteTargetType();
				in.setWriteTarget(fname, DistributedImageStorageType::Raw);
				
				std::vector<ParamVariant> args2;
				ParamVariant p;
				p = &in;
				args2.push_back(p);
				NopSingleImageCommand<pixel_t>& cmd = CommandList::get<NopSingleImageCommand<pixel_t> >();
				std::vector<std::string> result = cmd.runDistributed(distributor, args2);

				// Ensure that following operations do not write to the same file.
				// Retain old write target storage type.
				in.newWriteTarget(oldStorageType);
				return result;
			}

			return std::vector<std::string>();
		}
	};


	class WriteRGBRawCommand : virtual public Command
	{
	protected:
		friend class CommandList;

		WriteRGBRawCommand() : Command("writeraw", "Writes an RGB image to a .raw file. The dimensions of the image are automatically appended to the file name.",
			{
				CommandArgument<Image<uint8_t> >(ParameterDirection::In, "r", "Red component image."),
				CommandArgument<Image<uint8_t> >(ParameterDirection::In, "g", "Green component image."),
				CommandArgument<Image<uint8_t> >(ParameterDirection::In, "b", "Blue component image."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the file to write. If the file exists, its current contents are erased (unless append parameter is set to true)."),
				CommandArgument<bool>(ParameterDirection::In, "append", "Set to true to append to existing .raw file.", false)
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<uint8_t>& r = *pop<Image<uint8_t>* >(args);
			Image<uint8_t>& g = *pop<Image<uint8_t>* >(args);
			Image<uint8_t>& b = *pop<Image<uint8_t>* >(args);
			std::string fname = pop<std::string>(args);
			bool append = pop<bool>(args);

			raw::writed(r, g, b, fname, !append);
		}
	};


	template<typename pixel_t> class WriteSequenceCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		WriteSequenceCommand() : Command("writesequence", "Write an image sequence to disk.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of file to write. Existing files are overwritten. Specify file name as a template where @ signifies the location where the slice number should be added. Use @(n) to specify field width, e.g., if n = 4, slice number 7 would be labeled 0007. Specify @(-) to select suitable field width automatically. Specify empty file name to save using default template '@.png' that corresponds to file names 0.png, 1.png, 2.png etc.")
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);
			sequence::write(in, fname);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			distributor.flush();

			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);

			//if (in.isSavedToDisk() && in.isSequence())
			if (in.isSavedToDisk() && in.currentReadSourceType() == DistributedImageStorageType::Sequence)
			{
				if (in.isSavedToTemp())
				{
					if (sequence::isSequence(in.currentReadSource()))
					{
						// The image has been saved to temporary files
						// Just move the temporary files to new location (and name) and set read source to that file.
						sequence::moveSequence(in.currentReadSource(), fname);
						in.setReadSource(fname, false);
						return std::vector<std::string>();
					}
					else
					{
						// The input image is empty, use default processing to create a new sequence.
					}
				}
				else
				{
					if (sequence::isSequence(in.currentReadSource()))
					{
						// The image has been saved to non-temporary files, so we cannot just move them.
						// The files must be copied.
						sequence::copySequence(in.currentReadSource(), fname);
						return std::vector<std::string>();
					}
					else
					{
						// The input image has been saved to a non-temporary files but those do not form a sequence.
						// This is impossible situation: the source data is not available anymore.
						throw ITLException("The image " + in.varName() + " is loaded from " + in.currentReadSource() + " but a valid image sequence is not found given that template.");
					}
				}
			}
			
			auto oldStorageType = in.currentWriteTargetType();
			in.setWriteTarget(fname, DistributedImageStorageType::Sequence);

			std::vector<ParamVariant> args2;
			ParamVariant p;
			p = &in;
			args2.push_back(p);
			NopSingleImageCommand<pixel_t>& cmd = CommandList::get<NopSingleImageCommand<pixel_t> >();
			std::vector<std::string> result = cmd.runDistributed(distributor, args2);

			// Ensure that following operations do not write to the same file.
			// Retain old write target storage type.
			in.newWriteTarget(oldStorageType);
			return result;
		}
	};

	//template<typename pixel_t> class WriteRawBlockCommand : public Command
	//{
	//protected:
	//	friend class CommandList;

	//	WriteRawBlockCommand() : Command("writerawblock", "Write an image to a specified position in a .raw file. Optionally can write only a block of the source image. " + rawFilenameFormatHelp(),
	//		{
	//			CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
	//			CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of file to write."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "x", "X-position of the image in the target file."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-position of the image in the target file."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-position of the image in the target file."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the output file. Specify zero to parse dimensions from the file name.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the output file. Specify zero to parse dimensions from the file name.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the output file. Specify zero to parse dimensions from the file name.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "source x", "X-position of the block of the source image to write.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "source y", "Y-position of the block of the source image to write.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "source z", "Z-position of the block of the source image to write.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "source width", "Width of the block to write. Specify a negative value to write the whole source image.", -1),
	//			CommandArgument<coord_t>(ParameterDirection::In, "source height", "Height of the block to write. Specify a negative value to write the whole source image.", -1),
	//			CommandArgument<coord_t>(ParameterDirection::In, "source depth", "Depth of the block to write. Specify a negative value to write the whole source image.", -1),
	//		})
	//	{
	//	}

	//public:
	//	virtual void run(std::vector<ParamVariant>& args) const override
	//	{
	//		Image<pixel_t>& img = *pop<Image<pixel_t>* >(args);
	//		std::string fname = pop<std::string>(args);

	//		coord_t x = pop<coord_t>(args);
	//		coord_t y = pop<coord_t>(args);
	//		coord_t z = pop<coord_t>(args);
	//		
	//		coord_t w = pop<coord_t>(args);
	//		coord_t h = pop<coord_t>(args);
	//		coord_t d = pop<coord_t>(args);

	//		coord_t ix = pop<coord_t>(args);
	//		coord_t iy = pop<coord_t>(args);
	//		coord_t iz = pop<coord_t>(args);
	//		coord_t iw = pop<coord_t>(args);
	//		coord_t ih = pop<coord_t>(args);
	//		coord_t id = pop<coord_t>(args);


	//		if (iw <= 0 || ih <= 0 || id <= 0)
	//		{
	//			iw = img.width();
	//			ih = img.height();
	//			id = img.depth();
	//		}

	//		// Parse dimensions from file name if no dimensions are provided
	//		if (w <= 0 || h <= 0 || d <= 0)
	//		{
	//			Vec3c dims;
	//			if (!raw::internals::parseDimensions(fname, dims))
	//				throw ParseException(std::string("Unable to find dimensions from file name: ") + fname);
	//			w = dims.x;
	//			h = dims.y;
	//			d = dims.z;
	//		}

	//		raw::writeBlock(img, fname, Vec3c(x, y, z), Vec3c(w, h, d), Vec3c(ix, iy, iz), Vec3c(iw, ih, id), true);
	//	}
	//};

	template<typename pixel_t> class WriteRawBlock2Command : public Command
	{
	protected:
		friend class CommandList;

		WriteRawBlock2Command() : Command("writerawblock", "Write an image to a specified position in a .raw file. Optionally can write only a block of the source image. " + rawFilenameFormatHelp(),
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of file to write."),
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position of the image in the target file."),
				CommandArgument<Vec3c>(ParameterDirection::In, "file dimensions", "Dimensions of the output file. Specify zero to parse dimensions from the file name.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "source position", "Position of the block of the source image to write.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "source block size", "Size of the block to write. Specify zero to write the whole source image.", Vec3c(0, 0, 0)),
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& img = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);

			Vec3c position = pop<Vec3c>(args);
			Vec3c fileSize = pop<Vec3c>(args);
			Vec3c blockPosition = pop<Vec3c>(args);
			Vec3c blockSize = pop<Vec3c>(args);


			if (blockSize.x <= 0 || blockSize.y <= 0 || blockSize.z <= 0)
				blockSize = img.dimensions();

			// Parse dimensions from file name if no dimensions are provided
			if (fileSize.x <= 0 || fileSize.y <= 0 || fileSize.z <= 0)
			{
				Vec3c dims;
				if (!raw::internals::parseDimensions(fname, dims))
					throw ParseException(std::string("Unable to find dimensions from file name: ") + fname);
				fileSize = dims;
			}

			raw::writeBlock(img, fname, position, fileSize, blockPosition, blockSize, true);
		}
	};


	//template<typename pixel_t> class WriteSequenceBlockCommand : public Command
	//{
	//protected:
	//	friend class CommandList;

	//	WriteSequenceBlockCommand() : Command("writesequenceblock", "Write an image to a specified position in an image sequence. Optionally can write only a block of the source image.",
	//		{
	//			CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
	//			CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of file to write."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "x", "X-position of the image in the target file."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-position of the image in the target file."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-position of the image in the target file."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the output file. Specify zero to parse dimensions from the files.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the output file. Specify zero to parse dimensions from the files.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the output file. Specify zero to parse dimensions from the files.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "ix", "X-position of the block of the source image to write.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "iy", "Y-position of the block of the source image to write.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "iz", "Z-position of the block of the source image to write.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "iwidth", "Width of the block to write. Specify a negative value to write the whole source image.", -1),
	//			CommandArgument<coord_t>(ParameterDirection::In, "iheight", "Height of the block to write. Specify a negative value to write the whole source image.", -1),
	//			CommandArgument<coord_t>(ParameterDirection::In, "idepth", "Depth of the block to write. Specify a negative value to write the whole source image.", -1),
	//		})
	//	{
	//	}

	//public:
	//	virtual void run(std::vector<ParamVariant>& args) const override
	//	{
	//		Image<pixel_t>& img = *pop<Image<pixel_t>* >(args);
	//		std::string fname = pop<std::string>(args);

	//		coord_t x = pop<coord_t>(args);
	//		coord_t y = pop<coord_t>(args);
	//		coord_t z = pop<coord_t>(args);

	//		coord_t w = pop<coord_t>(args);
	//		coord_t h = pop<coord_t>(args);
	//		coord_t d = pop<coord_t>(args);

	//		coord_t ix = pop<coord_t>(args);
	//		coord_t iy = pop<coord_t>(args);
	//		coord_t iz = pop<coord_t>(args);
	//		coord_t iw = pop<coord_t>(args);
	//		coord_t ih = pop<coord_t>(args);
	//		coord_t id = pop<coord_t>(args);

	//		if (iw <= 0 || ih <= 0 || id <= 0)
	//		{
	//			iw = img.width();
	//			ih = img.height();
	//			id = img.depth();
	//		}

	//		// Parse dimensions from files if no dimensions are provided
	//		if (w <= 0 || h <= 0 || d <= 0)
	//		{
	//			Vec3c dims;
	//			itl2::ImageDataType dt2;
	//			std::string reason;
	//			if (!sequence::getInfo(fname, dims, dt2, reason))
	//				throw ParseException(std::string("Unable to find metadata from sequence with template: ") + fname + ". " + reason);
	//			w = dims.x;
	//			h = dims.y;
	//			d = dims.z;
	//		}

	//		sequence::writeBlock(img, fname, Vec3c(x, y, z), Vec3c(w, h, d), Vec3c(ix, iy, iz), Vec3c(iw, ih, id), true);
	//	}
	//};

	template<typename pixel_t> class WriteSequenceBlock2Command : public Command
	{
	protected:
		friend class CommandList;

		WriteSequenceBlock2Command() : Command("writesequenceblock", "Write an image to a specified position in an image sequence. Optionally can write only a block of the source image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the file to write."),
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position of the image in the target file."),
				CommandArgument<Vec3c>(ParameterDirection::In, "file dimensions", "Dimensions of the output file. Specify zero to parse dimensions from the file. In this case it must exist.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "source position", "Position of the block of the source image to write.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "source block size", "Size of the block to write. Specify zero to write the whole source image.", Vec3c(0, 0, 0)),
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& img = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);

			Vec3c position = pop<Vec3c>(args);
			Vec3c fileSize = pop<Vec3c>(args);
			Vec3c blockPosition = pop<Vec3c>(args);
			Vec3c blockSize = pop<Vec3c>(args);

			if (blockSize.x <= 0 || blockSize.y <= 0 || blockSize.z <= 0)
				blockSize = img.dimensions();

			// Parse dimensions from files if no dimensions are provided
			if (fileSize.x <= 0 || fileSize.y <= 0 || fileSize.z <= 0)
			{
				Vec3c dims;
				itl2::ImageDataType dt2;
				std::string reason;
				if (!sequence::getInfo(fname, dims, dt2, reason))
					throw ParseException(std::string("Unable to find metadata from sequence with template: ") + fname + ". " + reason);
				fileSize = dims;
			}

			sequence::writeBlock(img, fname, position, fileSize, blockPosition, blockSize, true);
		}
	};




	//template<typename pixel_t> class WriteNN5BlockCommand : public Command
	//{
	//protected:
	//	friend class CommandList;

	//	WriteNN5BlockCommand() : Command("writenn5block", "Write an image to a specified position in an NN5 dataset. Optionally can write only a block of the source image.",
	//		{
	//			CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
	//			CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the dataset to write."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "x", "X-position of the image in the target file."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-position of the image in the target file."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-position of the image in the target file."),
	//			CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the output file. Specify zero to parse dimensions from the files.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the output file. Specify zero to parse dimensions from the files.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the output file. Specify zero to parse dimensions from the files.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "ix", "X-position of the block of the source image to write.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "iy", "Y-position of the block of the source image to write.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "iz", "Z-position of the block of the source image to write.", 0),
	//			CommandArgument<coord_t>(ParameterDirection::In, "iwidth", "Width of the block to write. Specify a negative value to write the whole source image.", -1),
	//			CommandArgument<coord_t>(ParameterDirection::In, "iheight", "Height of the block to write. Specify a negative value to write the whole source image.", -1),
	//			CommandArgument<coord_t>(ParameterDirection::In, "idepth", "Depth of the block to write. Specify a negative value to write the whole source image.", -1),
	//		})
	//	{
	//	}

	//public:
	//	virtual void run(std::vector<ParamVariant>& args) const override
	//	{
	//		Image<pixel_t>& img = *pop<Image<pixel_t>* >(args);
	//		std::string fname = pop<std::string>(args);

	//		coord_t x = pop<coord_t>(args);
	//		coord_t y = pop<coord_t>(args);
	//		coord_t z = pop<coord_t>(args);

	//		coord_t w = pop<coord_t>(args);
	//		coord_t h = pop<coord_t>(args);
	//		coord_t d = pop<coord_t>(args);

	//		coord_t ix = pop<coord_t>(args);
	//		coord_t iy = pop<coord_t>(args);
	//		coord_t iz = pop<coord_t>(args);
	//		coord_t iw = pop<coord_t>(args);
	//		coord_t ih = pop<coord_t>(args);
	//		coord_t id = pop<coord_t>(args);

	//		if (iw <= 0 || ih <= 0 || id <= 0)
	//		{
	//			iw = img.width();
	//			ih = img.height();
	//			id = img.depth();
	//		}

	//		// Parse dimensions from files if no dimensions are provided
	//		if (w <= 0 || h <= 0 || d <= 0)
	//		{
	//			Vec3c dims;
	//			itl2::ImageDataType dt2;
	//			std::string reason;
	//			if (!nn5::getInfo(fname, dims, dt2, reason))
	//				throw ParseException(std::string("Unable to find metadata from NN5 dataset: ") + fname + ". " + reason);
	//			w = dims.x;
	//			h = dims.y;
	//			d = dims.z;
	//		}

	//		nn5::writeBlock(img, fname, nn5::DEFAULT_CHUNK_SIZE, nn5::NN5Compression::LZ4, Vec3c(x, y, z), Vec3c(w, h, d), Vec3c(ix, iy, iz), Vec3c(iw, ih, id));
	//	}
	//};

	template<typename pixel_t> class WriteNN5Block2Command : public Command
	{
	protected:
		friend class CommandList;

		WriteNN5Block2Command() : Command("writenn5block", "Write an image to a specified position in an NN5 dataset. Optionally can write only a block of the source image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image to save."),
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the dataset to write."),
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position of the image in the target file."),
				CommandArgument<Vec3c>(ParameterDirection::In, "file dimensions", "Dimensions of the output file. Specify zero to parse dimensions from the file. In this case it must exist.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "source position", "Position of the block of the source image to write.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "source block size", "Size of the block to write. Specify zero to write the whole source image.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "chunk size", "Size of chunks in the NN5 dataset.", nn5::DEFAULT_CHUNK_SIZE)
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			Image<pixel_t>& img = *pop<Image<pixel_t>* >(args);
			std::string fname = pop<std::string>(args);

			Vec3c position = pop<Vec3c>(args);
			Vec3c fileSize = pop<Vec3c>(args);
			Vec3c blockPosition = pop<Vec3c>(args);
			Vec3c blockSize = pop<Vec3c>(args);
			Vec3c chunkSize = pop<Vec3c>(args);

			if (blockSize.x <= 0 || blockSize.y <= 0 || blockSize.z <= 0)
				blockSize = img.dimensions();

			// Parse dimensions from files if no dimensions are provided
			if (fileSize.x <= 0 || fileSize.y <= 0 || fileSize.z <= 0)
			{
				Vec3c dims;
				itl2::ImageDataType dt2;
				std::string reason;
				if (!nn5::getInfo(fname, dims, dt2, reason))
					throw ParseException(std::string("Unable to find metadata from NN5 dataset: ") + fname + ". " + reason);
				fileSize = dims;
			}

			nn5::writeBlock(img, fname, chunkSize, nn5::NN5Compression::LZ4, position, fileSize, blockPosition, blockSize, true);
		}
	};


	class EndConcurrentWriteCommand : public Command
	{
	protected:
		friend class CommandList;

		EndConcurrentWriteCommand() : Command("endconcurrentwrite", "End concurrent writing to an NN5 dataset. This command is used internally in distributed processing.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "filename", "Name (and path) of the dataset."),
				CommandArgument<Vec3c>(ParameterDirection::In, "chunk index", "Index of the chunk to process.")
			})
		{
		}

	public:

		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			TimingFlag flag(TimeClass::IO);

			std::string fname = pop<std::string>(args);
			Vec3c chunkIndex = pop<Vec3c>(args);

			nn5::endConcurrentWrite(fname, chunkIndex);
		}
	};
}
