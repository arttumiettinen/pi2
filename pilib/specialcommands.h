#pragma once

#include "command.h"
#include "distributable.h"
#include "trivialdistributable.h"
#include "commandsbase.h"
#include "standardhelp.h"

namespace pilib
{

	inline std::string helpSeeAlso()
	{
		return "help, info, license, echo, print, waitreturn, hello, timing, savetiming, resettiming";
	}

	inline std::string timeClassHelp()
	{
		return
			"The output includes the following time classes.\n"
			"\n"
			"**Overhead**\n"
			"\n"
			"General overhead, e.g. parsing inputs, finding correct commands to run etc.\n"
			"This includes the total overhead time spent in the main process, and in possible cluster job processes.\n"
			"\n"
			"**I/O**\n"
			"\n"
			"Time spent in I/O-bound processing. This is the time when the disk I/O is the bottleneck.\n"
			"This includes the total I/O time spent in the main process, and in possible cluster job processes.\n"
			"Time spent in output data compression is counted to this time class.\n"
			"\n"
			"**Computation**\n"
			"\n"
			"Time spent in CPU/GPU-bound processing. This is the time when the CPU/GPU is the bottleneck.\n"
			"This includes the total computation time spent in the main process, and in possible cluster job processes.\n"
			"This is the default mode for all commands.\n"
			"\n"
			"**Job execution**\n"
			"\n"
			"Total distributed job execution time.\n"
			"This value includes Overhead+IO+Computation of all jobs, plus workload manager node reservation, process starting, etc. overhead.\n"
			"This value does not include time spent in workload manager queue.\n"
			"\n"
			"**Job queuing**\n"
			"\n"
			"Total distributed job queuing time.\n"
			"This is the total time all jobs have spent in the workload manager queue, waiting to be executed.\n"
			"\n"
			"**Total job waiting**\n"
			"\n"
			"Total time from submitting the first distributed job until all of them are found to be finished.\n"
			"This is the total time spent in the job execution process, from submission to jobs until all of them are done.\n"
			"\n"
			"**Write preparation**\n"
			"\n"
			"Time spent in preparing for writing output images (e.g. NN5 write preparation).\n"
			"This is the total time spent in the submitting process while preparing writing of output images.\n"
			"\n"
			"**Total write finalization waiting**\n"
			"\n"
			"Time spent in write finalization jobs, including queuing (e.g. NN5 write finalization jobs).\n"
			"This is the total time spent in the submitting process, from submission of write finalization jobs until all of them are done.\n";
	}

	class TimingCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		TimingCommand() : Command("timing", string("Prints information about wall-clock time taken by various sub-processes. Running this command causes all delayed commands to be executed. ") + timeClassHelp(),
			{
			},
			helpSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override;

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return false;
		}
	};

	class SaveTimingCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		SaveTimingCommand() : Command("savetiming", string("Saves timing information to a file. Running this command causes all delayed commands to be executed. ") + timeClassHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "file name", "The name and path of the file where the information is to be saved.")
			},
			helpSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override;

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return false;
		}
	};

	class ResetTimingCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		ResetTimingCommand() : Command("resettiming", string("Zeroes all existing timing data."),
			{
			},
			helpSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override;

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return false;
		}
	};

	
	class HelloCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		HelloCommand() : Command("hello", "Shows a greetings message. This command is used for testing things.",
			{
				CommandArgument<string>(ParameterDirection::In, "name", "Name of caller. Omit to output a generic greeting.", "there")
			},
			helpSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			std::cout << "Hello " << pop<string>(args) << "!" << std::endl;
		}
	};


	class PrintCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		PrintCommand() : Command("print", "Prints a message to the output. This is mainly used internally in distributed processing.",
			{
				CommandArgument<string>(ParameterDirection::In, "string", "String to print")
			},
			helpSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			std::cout << pop<string>(args) << std::endl;
		}
	};


	class WaitReturnCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		WaitReturnCommand() : Command("waitreturn", "Waits until the user presses return key. Does not show any prompts on screen.",
			{},
			helpSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			//std::cout << "Press return to exit..." << std::endl;
			char dummy;
			std::cin.getline(&dummy, 1);
		}
	};

	class LicenseCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		LicenseCommand() : Command("license", "Displays license of this program and PI system.", {}, helpSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override;
	};

	class ListCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;


		ListCommand() : Command("list", "Displays information about images loaded to memory, or images prepared for distributed processing.", {}, helpSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{

		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};


	inline std::string sequenceDefinitionHelp()
	{
		return "If a directory name is given, all files in the directory will be read. If the a file name is given, it can contain wildcards $*$, $?$ and $@$. Wildcard $*$ corresponds to any sequence of characters, wildcard $?$ corresponds to any character, and wildcard $@$ corresponds to sequence of numerical digits. For example, sequence containing files xyz_000.png, xyz_001.png, xyz_002.png, etc. could be read with template xyz_@.png.";
	}

	class ReadCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		ReadCommand() : Command("read", "Reads an image or image sequence from disk. Determines type of file automatically.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read or a sequence definition. " + sequenceDefinitionHelp()),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ". Specify empty value to infer data type from file content.", "")
			},
			"",
			"In Python/pi2py2, the image name parameter is not specified, and the function returns the newly created image read from the disk.")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;

		
	};


	class ReadRawCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		ReadRawCommand() : Command("readraw", "Reads a .raw image from a file. If the image is not in the native byte order of the host computer, the byte order may be changed using `swapbyteorder` command. " + rawFilenameFormatHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image to create."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ". Specify empty value to infer data type from file size."),
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the image. Omit width, height and depth to infer dimensions from file name."),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the image. Omit width, height and depth to infer dimensions from file name."),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the image. Omit width, height and depth to infer dimensions from file name.")
			},
			"mapraw, readrawblock")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};


	class ReadRaw2Command : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		ReadRaw2Command() : Command("readraw", "Reads a .raw image from a file. If the image is not in the native byte order of the host computer, the byte order may be changed using `swapbyteorder` command. " + rawFilenameFormatHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image to create."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ". Specify empty value to infer data type from file size.", ""),
				CommandArgument<Vec3c>(ParameterDirection::In, "dimensions", "Size of the image. Set to zero to infer dimensions from file name.", Vec3c(0, 0, 0))
			},
			"mapraw, readrawblock")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};


	class ReadSequenceCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		ReadSequenceCommand() : Command("readsequence", "Reads an image sequence from disk. Supports any image formats supported by the back end (at least .png).",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename template", "Name (and path) template of the sequence to be read. " + sequenceDefinitionHelp()),
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};

	class ReadVolCommand : public Command
	{
	protected:
		friend class CommandList;

		ReadVolCommand() : Command("readvol", "Reads a .vol image from a file.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read.")
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};

	class ReadBlockCommand : public Command
	{
	protected:
		friend class CommandList;

		ReadBlockCommand() : Command("readblock", "Reads a block of an image from a file. If the file is a .raw file, special conditions apply: " + rawFilenameFormatHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "x", "X-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block width", "Width of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block height", "Height of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block depth", "Depth of block to read."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Used only if reading .raw images. Leave empty to guess data type based on file size.", "")
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};

	class ReadBlock2Command : public Command
	{
	protected:
		friend class CommandList;

		ReadBlock2Command() : Command("readblock", "Reads a block of an image from a file. If the file is a .raw file, special conditions apply: " + rawFilenameFormatHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Coordinates of the first pixel to read."),
				CommandArgument<Vec3c>(ParameterDirection::In, "block size", "Dimensions of the block to read."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Used only if reading .raw images. Leave empty to guess data type based on file size.", "")
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class ReadRawBlockCommand : public Command
	{
	protected:
		friend class CommandList;

		ReadRawBlockCommand() : Command("readrawblock", "Reads a block of a .raw image file. " + rawFilenameFormatHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "x", "X-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block width", "Width of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block height", "Height of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block depth", "Depth of block to read."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ". Specify empty value to infer data type from image dimensions", ""),
				CommandArgument<coord_t>(ParameterDirection::In, "total width", "Width of the image. Omit width, height and depth to infer dimensions from file name.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "total height", "Height of the image. Omit width, height and depth to infer dimensions from file name.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "total depth", "Depth of the image. Omit width, height and depth to infer dimensions from file name.", 0)
			},
			"mapraw, readraw")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};

	class ReadRawBlock2Command : public Command
	{
	protected:
		friend class CommandList;

		ReadRawBlock2Command() : Command("readrawblock", "Reads a block of a .raw image file. " + rawFilenameFormatHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Coordinates of the first pixel to read."),
				CommandArgument<Vec3c>(ParameterDirection::In, "block size", "The size of the block to read."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ". Specify empty value to infer data type from image dimensions", ""),
				CommandArgument<Vec3c>(ParameterDirection::In, "image size", "Dimensions of the full image. Set to zero to infer dimensions from file name.", Vec3c(0, 0, 0)),
			},
			"mapraw, readraw")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};



	class ReadSequenceBlockCommand : public Command
	{
	protected:
		friend class CommandList;

		ReadSequenceBlockCommand() : Command("readsequenceblock", "Reads a block of an image sequence.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "x", "X-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block width", "Width of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block height", "Height of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block depth", "Depth of block to read.")
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class ReadSequenceBlock2Command : public Command
	{
	protected:
		friend class CommandList;

		ReadSequenceBlock2Command() : Command("readsequenceblock", "Reads a block of an image sequence.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Coordinates of the first pixel to read."),
				CommandArgument<Vec3c>(ParameterDirection::In, "block size", "The size of the block to read."),
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class ReadNN5BlockCommand : public Command
	{
	protected:
		friend class CommandList;

		ReadNN5BlockCommand() : Command("readnn5block", "Reads a block of an NN5 dataset.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of the dataset to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "x", "X-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block width", "Width of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block height", "Height of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block depth", "Depth of block to read.")
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class ReadNN5Block2Command : public Command
	{
	protected:
		friend class CommandList;

		ReadNN5Block2Command() : Command("readnn5block", "Reads a block of an NN5 dataset.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of the dataset to read."),
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Coordinates of the first pixel to read."),
				CommandArgument<Vec3c>(ParameterDirection::In, "block size", "The size of the block to read."),
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class MapRawCommand : public Command
	{
	protected:
		friend class CommandList;

		MapRawCommand() : Command("mapraw", "Maps a .raw image file to an image. If the image file does not exist, it is created. Changes made to the file are IMMEDIATELY reflected to disk, so there is no need to save or load the image. The image must be in the native byte order of the host computer. " + rawFilenameFormatHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to map."),
				CommandArgument<string>(ParameterDirection::In, "data type", string("Data type of the image. Can be ") + listSupportedImageDataTypes() + ". Specify empty value to infer data type from image dimensions."),
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the image. Omit width, height and depth to infer dimensions from file name."),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the image. Omit width, height and depth to infer dimensions from file name."),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the image. Omit width, height and depth to infer dimensions from file name."),
				CommandArgument<bool>(ParameterDirection::In, "read only", "Set to true to do read-only mapping. This might be beneficial if the image file is accessed through a network share. WARNING: If set to true, writes to the image result in undefined behaviour, probably program crash to access violation or equivalent error.", false)
			},
			"readrawblock, readraw, getmapfile",
			"In Python/pi2py2, the image name parameter is not specified, and the function returns a newly created image mapped to the .raw image file.")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};

	class MapRaw2Command : public Command
	{
	protected:
		friend class CommandList;

		MapRaw2Command() : Command("mapraw", "Maps a .raw image file to an image. If the image file does not exist, it is created. Changes made to the file are IMMEDIATELY reflected to disk, so there is no need to save or load the image. The image must be in the native byte order of the host computer. " + rawFilenameFormatHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to map."),
				CommandArgument<string>(ParameterDirection::In, "data type", string("Data type of the image. Can be ") + listSupportedImageDataTypes() + ". Specify empty value to infer data type from image dimensions.", ""),
				CommandArgument<Vec3c>(ParameterDirection::In, "dimensions", "Dimensions of the image. Set to zero to infer dimensions from file name.", Vec3c(0, 0, 0)),
				CommandArgument<bool>(ParameterDirection::In, "read only", "Set to true to do read-only mapping. This might be beneficial if the image file is accessed through a network share. WARNING: If set to true, writes to the image result in undefined behaviour, probably program crash to access violation or equivalent error.", false)
			},
			"readrawblock, readraw, getmapfile",
			"In Python/pi2py2, the image name parameter is not specified, and the function returns a newly created image mapped to the .raw image file.")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	template<typename pixel_t> class GetMapFileCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		GetMapFileCommand() : OneImageInPlaceCommand<pixel_t>("getmapfile", "Get a path to the file where the argument image has been memory-mapped to. Returns empty string if no mapping has been made for the argument image.",
			{
				CommandArgument<string>(ParameterDirection::Out, "mapfile", "The name and path to the memory-mapped file."),
			},
			"mapraw")
		{
		}
	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			std::string* out = std::get<std::string*>(args[0]);
			*out = in.mappedFile();
		}
	};



	class SetStringCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		SetStringCommand() : Command("set", "Sets value of a string variable",
			{
				CommandArgument<string>(ParameterDirection::Out, "name", "Variable to set."),
				CommandArgument<string>(ParameterDirection::In, "value", "New value."),
			},
			"set, clear")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};




	class SetIntCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		SetIntCommand() : Command("set", "Sets value of an integer variable",
			{
				CommandArgument<coord_t>(ParameterDirection::Out, "name", "Variable to set."),
				CommandArgument<coord_t>(ParameterDirection::In, "value", "New value."),
			},
			"set, clear")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};



	class SetRealCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		SetRealCommand() : Command("set", "Sets value of an real number variable",
			{
				CommandArgument<double>(ParameterDirection::Out, "name", "Variable to set."),
				CommandArgument<double>(ParameterDirection::In, "value", "New value."),
			},
			"set, clear")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};


	class SetBoolCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		SetBoolCommand() : Command("set", "Sets value of an boolean variable",
			{
				CommandArgument<bool>(ParameterDirection::Out, "name", "Variable to set."),
				CommandArgument<bool>(ParameterDirection::In, "value", "New value."),
			},
			"set, clear")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};



	class NewValueCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		NewValueCommand() : Command("newvalue", "Creates a new variable.",
			{
				CommandArgument<string>(ParameterDirection::In, "name", "Name of the variable in the system."),
				CommandArgument<string>(ParameterDirection::In, "type", "Data type of the variable. Can be 'string', 'int', 'real', or 'bool'.", "string"),
				CommandArgument<string>(ParameterDirection::In, "value", "Initial value of the variable", ""),
			},
			"set, clear",
			"In Python/pi2py2, one should use the newstring command.")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};


	class NewImageCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		NewImageCommand() : Command("newimage", "Creates a new, empty image.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of the image in the system."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ".", "uint8"),
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the image.", 1),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the image.", 1),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the image.", 1)
			},
			"ensuresize, newlike",
			"In Python/pi2py2, the image name parameter is not specified, and the return value is a Pi2Image object that can be passed to any command expecting an image name as an argument.")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};

	class NewImage2Command : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		NewImage2Command() : Command("newimage", "Creates a new, empty image.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of the image in the system."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + "."),
				CommandArgument<Vec3c>(ParameterDirection::In, "dimensions", "Dimensions of the image."),
			},
			"ensuresize, newlike",
			"In Python/pi2py2, the image name parameter is not specified, and the return value is a Pi2Image object that can be passed to any command expecting an image name as an argument.")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};
	

	class ClearCommand : virtual public Command, public Distributable
	{
	protected:
		friend class CommandList;

		ClearCommand() : Command("clear", "Dispose image or value from the system (and free memory that the object consumes).",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image or value to erase. Specify empty name to clear everything.", "")
			},
			"list")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};

	class HelpCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		HelpCommand() : Command("help", "Shows usage information.",
			{
				CommandArgument<string>(ParameterDirection::In, "topic", "Name of command whose help is to be retrieved. Specify nothing to show general help.", ""),
				//CommandArgument<string>(ParameterDirection::In, "search term", "Search term. Specify something to show only those matches that contain this term. E.g. specify uint8 to show help only for functions that process uint8 images.", ""),
				CommandArgument<string>(ParameterDirection::In, "format", "Output format. Can be text for textual output or rst for ReStructuredText output. The latter is more complex format.", "text"),
			},
			helpSeeAlso())
		{
		}

	public:

		string run(const string& name, HelpFormat format) const;

		virtual void run(vector<ParamVariant>& args) const override;
	};

	class CommandReferenceCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		CommandReferenceCommand() : Command("commandreference", "Generates command reference in ReStructuredText format for building help.",
			{
				CommandArgument<string>(ParameterDirection::In, "folder", "Target folder where the RST files are placed."),
			})
		{
		}

	public:
		virtual bool isInternal() const override
		{
			return true;
		}

		virtual void run(vector<ParamVariant>& args) const override;
	};

	class EchoCommandsCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		EchoCommandsCommand() : Command("echo", "Enables or disables echoing of commands and timing information to output.",
			{
				CommandArgument<bool>(ParameterDirection::In, "commands", "Set to true to show commands that have been run.", true),
				CommandArgument<bool>(ParameterDirection::In, "timing", "Set to true to show timing information on commands that have been run.", false),
			},
			helpSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	inline std::string distributeSeeAlso()
	{
		return "distribute, delaying, maxmemory, maxjobs, chunksize, printscripts";
	}

	class DistributeCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		DistributeCommand() : Command("distribute", "Enables or disables distributed processing of commands. Run this command before commands that you would like to run using distributed processing. Images used during distributed processing are not available for local processing and vice versa, unless they are loaded again from disk. All commands do not support distributed processing.",
			{
				CommandArgument<string>(ParameterDirection::In, "workload manager system name", "Set to 'SLURM' to use SLURM workload manager; set to 'LOCAL' to process tasks sequentially using the local computer; set to empty string to disable distributed processing (default).", "")
			},
			distributeSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class MaxMemoryCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		MaxMemoryCommand() : Command("maxmemory", "Sets the maximum memory setting used in distributed processing. This command overrides the value read from the configuration file. The maximum memory is the amount of memory that can be used either on the local computer (Local distribution mode) or in a compute node (Slurm etc. distribution modes).",
			{
				CommandArgument<double>(ParameterDirection::In, "maximum memory", "Maximum amount of memory to use, in megabytes. Specify zero to calculate the value automatically.", 0.0)
			},
			distributeSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class GetMaxMemoryCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		GetMaxMemoryCommand() : Command("getmaxmemory", "Gets the maximum memory setting used in distributed processing. The maximum memory is the amount of memory that can be used either on the local computer (Local distribution mode) or in a compute node (Slurm etc. distribution modes).",
			{
				CommandArgument<double>(ParameterDirection::Out, "maximum memory", "Maximum amount of memory to use, in megabytes.")
			},
			distributeSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class MaxJobsCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		MaxJobsCommand() : Command("maxjobs", "Sets the maximum number of jobs that will be submitted in parallel in distributed mode. If there are more jobs to be submitted at once, the jobs are combined into larger jobs until the total number of jobs is below or equal to the speficied maximum. This command overrides max_parallel_submit_count value read from the distributor configuration file.",
			{
				CommandArgument<size_t>(ParameterDirection::In, "max number of jobs", "Maximum number of jobs to submit in parallel. If the analysis task requires more jobs than specified, the jobs are combined until only maximum number of jobs are left. Specify zero for unlimited number of jobs.", 0)
			},
			distributeSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class ChunkSizeCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		ChunkSizeCommand() : Command("chunksize", "Sets the NN5 dataset chunk size used in distributed processing. This command overrides the value read from the configuration file.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "chunk size", "New NN5 dataset chunk size.", nn5::DEFAULT_CHUNK_SIZE)
			},
			distributeSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	class DelayingCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		DelayingCommand() : Command("delaying", "Enables or disables possibility for combination of similar commands before sending them to compute processes/nodes in distributed processing (kind of lazy evaluation). Delaying decreases amount of disk space and I/O used for temporary images. When delaying is enabled, command run times (see echo command) are not accurate. Delaying is enabled by default.",
			{
				CommandArgument<bool>(ParameterDirection::In, "enable", "Set to true to enable delaying.", true)
			},
			distributeSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};

	class PrintTaskScriptsCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		PrintTaskScriptsCommand() : Command("printscripts", "Enables or disables printing of pi2 script run for each task in compute node/process. Has effect only if distributed processing is enabled.",
			{
				CommandArgument<bool>(ParameterDirection::In, "enable", "Set to true to enable printing.", true)
			},
			distributeSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
	};


	template<typename pixel_t> class EnsureSizeCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		EnsureSizeCommand() : OneImageInPlaceCommand<pixel_t>("ensuresize", "Makes sure that the size of the parameter image equals dimensions given as an arguments. The parameter image is re-allocated only if its size must be changed.",
			{
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Desired width of the image.", 1),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Desired height of the image.", 1),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Desired depth of the image.", 1)
			},
			"newimage")
		{
		}

	public:

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);
			in.ensureSize(w, h, d);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>*>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);
			img.ensureSize(w, h, d);
			return vector<string>();
		}
	};

	template<typename pixel_t> class EnsureSize2Command : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		EnsureSize2Command() : OneImageInPlaceCommand<pixel_t>("ensuresize", "Makes sure that the size of the parameter image equals dimensions given as an argument. The parameter image is re-allocated only if its size must be changed.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "size", "Desired image dimensions.")
			},
			"newimage")
		{
		}

	public:

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Vec3c d = pop<Vec3c>(args);
			in.ensureSize(d);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>*>(args);
			Vec3c d = pop<Vec3c>(args);
			img.ensureSize(d);
			return vector<string>();
		}
	};

}
