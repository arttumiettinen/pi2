#pragma once

#include "command.h"
#include "trivialdistributable.h"
#include "parseexception.h"

#include <vector>

#include "itl2.h"

using namespace std;
using namespace itl2;

namespace pilib
{
	class RawInfoCommand : virtual public Command, public TrivialDistributable
	{
	public:
		RawInfoCommand() : Command("rawinfo", "Displays metadata of .raw file.",
			{
				CommandArgument<string>(In, "filename", "Name of .raw file.")
			})
		{

		}

		virtual void run(vector<ParamVariant>& args) const
		{
			string fname = pop<string>(args);

			Vec3c dims;
			ImageDataType dt;
			if (!raw::internals::parseDimensions(fname, dims, dt))
				throw ITLException(string("Unable to parse file name ") + fname);

			cout << itl2::toString(dt) << endl;
			cout << dims.x << endl;
			cout << dims.y << endl;
			cout << dims.z << endl;
		}
	};

	class SequenceInfoCommand : virtual public Command, public TrivialDistributable
	{
	public:
		SequenceInfoCommand() : Command("sequenceinfo", "Displays metadata of image sequence.",
			{
				CommandArgument<string>(In, "filename template", "Filename template corresponding to the sequence, as described in readsequence command documentation.")
			})
		{

		}

		virtual void run(vector<ParamVariant>& args) const
		{
			string fname = pop<string>(args);

			Vec3c dims;
			ImageDataType dt;

			if (!sequence::getInfo(fname, dims, dt))
				throw ITLException(string("Unable to parse sequence ") + fname);

			cout << itl2::toString(dt) << endl;
			cout << dims.x << endl;
			cout << dims.y << endl;
			cout << dims.z << endl;
		}
	};

	template<typename pixel_t> class NopSingleImageCommand : virtual public Command, public TrivialDistributable
	{
	public:
		NopSingleImageCommand() : Command("nop", "Does nothing. This command is used internally as a helper that tags image that must be saved but not processed.",
			{
				CommandArgument<Image<pixel_t> >(InOut, "image", "Name of image that is not processed.")
			})
		{

		}

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};

	template<typename pixel_t> class WriteRawCommand : virtual public Command, public Distributable
	{
	public:
		WriteRawCommand() : Command("writeraw", "Write an image to .raw file. The dimensions of the image are automatically appended to the file name. If distributed processing is enabled, uses optimized implementation that does not need to copy or write any image data if the image is saved to a temporary image storage location and that location is on the same partition than the file being written.",
			{
				CommandArgument<Image<pixel_t> >(In, "input image", "Name of image to save."),
				CommandArgument<string>(In, "filename", "Name (and path) of file to write. If the file exists, its current contents are erased.")
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			string fname = pop<string>(args);

			raw::writed(in, fname);
		}

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			string fname = pop<string>(args);

			fname = raw::internals::concatDimensions(fname, in.dimensions());

			if (in.isRaw())
			{
				// Effectively copy data from input image to output image.
				// First move input image current data source file to output file.
				// Then point input image current data source file to the copied file.
				// This way no data needs to be copied (if all the data is stored on single partition)
				// and input image can still be used as further changes are anyway saved to a temp file.
				// NOTE: This does not work if the input and output image are stored in "cloud", but
				// that is not supported at the moment anyway.

				if (in.savedToTemp())
				{
					if (fileExists(in.currentReadSource()))
					{
						// The image has been saved to a temporary file
						// Just move the temporary file to new location (and name) and sets read source to that file.
						moveFile(in.currentReadSource(), fname);
						in.setReadSource(fname);
					}
					else
					{
						// The input image is empty, so just create an empty file.
						setFileSize(fname, in.dimensions().x * in.dimensions().y * in.dimensions().z * sizeof(pixel_t));
					}
				}
				else
				{
					if (fileExists(in.currentReadSource()))
					{
						// The image has been saved to a non-temporary file, so we cannot just move it.
						// The file must be copied.
						copyFile(in.currentReadSource(), fname);
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
				// The input is not stored as .raw file so actual copying must be made.
				//distributor.distributeReadWrite(in, &WriteRawBlockCommand(), fname);
				vector<ParamVariant> args2;
				ParamVariant p;
				p.dimgval = &in;
				args2.push_back(p);
				NopSingleImageCommand<pixel_t> cmd;
				distributor.distribute(&cmd, args2, 2, Vec3c(0, 0, 0), &fname);
			}

		}
	};

	template<typename pixel_t> class WriteSequenceCommand : public Command, public Distributable
	{
	public:
		WriteSequenceCommand() : Command("writesequence", "Write an image sequence to disk.",
			{
				CommandArgument<Image<pixel_t> >(In, "input image", "Name of image to save."),
				CommandArgument<string>(In, "filename", "Name (and path) of file to write. Existing files are overwritten. Specify file name as a template where @ signifies the location where the slice number should be added. Use @(n) to specify field width, e.g., if n = 4, slice number 7 would be labeled 0007. Specify @(-) to select suitable field width automatically. Specify empty file name to save using default template '@.png' that corresponds to file names 0.png, 1.png, 2.png etc.")
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			string fname = pop<string>(args);
			sequence::write(in, fname);
		}

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			string fname = pop<string>(args);

			if (in.isSequence())
			{
				if (in.savedToTemp())
				{
					if (sequence::isSequence(in.currentReadSource()))
					{
						// The image has been saved to temporary files
						// Just move the temporary files to new location (and name) and set read source to that file.
						sequence::moveSequence(in.currentReadSource(), fname);
						in.setReadSource(fname);
						return;
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
						return;
					}
					else
					{
						// The input image has been saved to a non-temporary files but those do not form a sequence.
						// This is impossible situation: the source data is not available anymore.
						throw ITLException("The image " + in.varName() + " is loaded from " + in.currentReadSource() + " but a valid image sequence is not found given that template.");
					}
				}
			}

			vector<ParamVariant> args2;
			ParamVariant p;
			p.dimgval = &in;
			args2.push_back(p);
			NopSingleImageCommand<pixel_t> cmd;
			distributor.distribute(&cmd, args2, 2, Vec3c(0, 0, 0), &fname);
		}
	};

	template<typename pixel_t> class WriteRawBlockCommand : public Command
	{
	public:
		WriteRawBlockCommand() : Command("writerawblock", "Write an image to a specified position in a .raw file. Optionally can write only a block of the source image.",
			{
				CommandArgument<Image<pixel_t> >(In, "input image", "Name of image file. If the file does not exist, it is created. If the size of the file is incorrect, its size is changed but the contents are not erased."),
				CommandArgument<string>(In, "filename", "Name (and path) of file to write."),
				CommandArgument<coord_t>(In, "x", "X-position of the image in the target file."),
				CommandArgument<coord_t>(In, "y", "Y-position of the image in the target file."),
				CommandArgument<coord_t>(In, "z", "Z-position of the image in the target file."),
				CommandArgument<coord_t>(In, "width", "Width of the output file. Specify zero to parse dimensions from the file name.", 0),
				CommandArgument<coord_t>(In, "height", "Height of the output file. Specify zero to parse dimensions from the file name.", 0),
				CommandArgument<coord_t>(In, "depth", "Depth of the output file. Specify zero to parse dimensions from the file name.", 0),
				CommandArgument<coord_t>(In, "ix", "X-position of the block of the source image to write.", 0),
				CommandArgument<coord_t>(In, "iy", "Y-position of the block of the source image to write.", 0),
				CommandArgument<coord_t>(In, "iz", "Z-position of the block of the source image to write.", 0),
				CommandArgument<coord_t>(In, "iwidth", "Width of the block to write. Specify a negative value to write the whole source image.", -1),
				CommandArgument<coord_t>(In, "iheight", "Height of the block to write. Specify a negative value to write the whole source image.", -1),
				CommandArgument<coord_t>(In, "idepth", "Depth of the block to write. Specify a negative value to write the whole source image.", -1),
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& img = *pop<Image<pixel_t>* >(args);
			string fname = pop<string>(args);

			coord_t x = pop<coord_t>(args);
			coord_t y = pop<coord_t>(args);
			coord_t z = pop<coord_t>(args);
			
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);

			coord_t ix = pop<coord_t>(args);
			coord_t iy = pop<coord_t>(args);
			coord_t iz = pop<coord_t>(args);
			coord_t iw = pop<coord_t>(args);
			coord_t ih = pop<coord_t>(args);
			coord_t id = pop<coord_t>(args);


			if (iw <= 0 || ih <= 0 || id <= 0)
			{
				iw = img.width();
				ih = img.height();
				id = img.depth();
			}

			// Parse dimensions from file name if no dimensions are provided
			if (w <= 0 || h <= 0 || d <= 0)
			{
				Vec3c dims;
				itl2::ImageDataType dt2;
				if (!raw::internals::parseDimensions(fname, dims, dt2))
					throw ParseException(string("Unable to find dimensions from file name: ") + fname);
				w = dims.x;
				h = dims.y;
				d = dims.z;
			}

			raw::writeBlock(img, fname, Vec3c(x, y, z), Vec3c(w, h, d), Vec3c(ix, iy, iz), Vec3c(iw, ih, id), true);
		}
	};

	template<typename pixel_t> class WriteSequenceBlockCommand : public Command
	{
	public:
		WriteSequenceBlockCommand() : Command("writesequenceblock", "Write an image to a specified position in an image sequence. Optionally can write only a block of the source image.",
			{
				CommandArgument<Image<pixel_t> >(In, "input image", "Name of image file. If the file does not exist, it is created. If the size of the file is incorrect, its size is changed but the contents are not erased."),
				CommandArgument<string>(In, "filename", "Name (and path) of file to write."),
				CommandArgument<coord_t>(In, "x", "X-position of the image in the target file."),
				CommandArgument<coord_t>(In, "y", "Y-position of the image in the target file."),
				CommandArgument<coord_t>(In, "z", "Z-position of the image in the target file."),
				CommandArgument<coord_t>(In, "width", "Width of the output file. Specify zero to parse dimensions from the files.", 0),
				CommandArgument<coord_t>(In, "height", "Height of the output file. Specify zero to parse dimensions from the files.", 0),
				CommandArgument<coord_t>(In, "depth", "Depth of the output file. Specify zero to parse dimensions from the files.", 0),
				CommandArgument<coord_t>(In, "ix", "X-position of the block of the source image to write.", 0),
				CommandArgument<coord_t>(In, "iy", "Y-position of the block of the source image to write.", 0),
				CommandArgument<coord_t>(In, "iz", "Z-position of the block of the source image to write.", 0),
				CommandArgument<coord_t>(In, "iwidth", "Width of the block to write. Specify a negative value to write the whole source image.", -1),
				CommandArgument<coord_t>(In, "iheight", "Height of the block to write. Specify a negative value to write the whole source image.", -1),
				CommandArgument<coord_t>(In, "idepth", "Depth of the block to write. Specify a negative value to write the whole source image.", -1),
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& img = *pop<Image<pixel_t>* >(args);
			string fname = pop<string>(args);

			coord_t x = pop<coord_t>(args);
			coord_t y = pop<coord_t>(args);
			coord_t z = pop<coord_t>(args);

			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);

			coord_t ix = pop<coord_t>(args);
			coord_t iy = pop<coord_t>(args);
			coord_t iz = pop<coord_t>(args);
			coord_t iw = pop<coord_t>(args);
			coord_t ih = pop<coord_t>(args);
			coord_t id = pop<coord_t>(args);

			if (iw <= 0 || ih <= 0 || id <= 0)
			{
				iw = img.width();
				ih = img.height();
				id = img.depth();
			}

			// Parse dimensions from files if no dimensions are provided
			if (w <= 0 || h <= 0 || d <= 0)
			{
				Vec3c dims;
				itl2::ImageDataType dt2;
				if (!sequence::getInfo(fname, dims, dt2))
					throw ParseException(string("Unable to find metadata from sequence with template: ") + fname);
				w = dims.x;
				h = dims.y;
				d = dims.z;
			}

			sequence::writeBlock(img, fname, Vec3c(x, y, z), Vec3c(w, h, d), Vec3c(ix, iy, iz), Vec3c(iw, ih, id), true);
		}
	};
}
