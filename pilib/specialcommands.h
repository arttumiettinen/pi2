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
		return "help, info, license, echo, print, waitreturn, hello";
	}
	
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

	class InfoCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		InfoCommand() : Command("info", "Displays information about the computer and the PI system.",
			{},
			helpSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override;
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

		ReadRawBlockCommand() : Command("readrawblock", "Reads a block of a .raw image from a file. " + rawFilenameFormatHelp(),
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

		ReadRawBlock2Command() : Command("readrawblock", "Reads a block of a .raw image from a file. " + rawFilenameFormatHelp(),
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

		ReadSequenceBlockCommand() : Command("readsequenceblock", "Reads a block of an image sequence from a file.",
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

		ReadSequenceBlock2Command() : Command("readsequenceblock", "Reads a block of an image sequence from a file.",
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
			"readrawblock, readraw")
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
			"readrawblock, readraw")
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override;

		virtual void run(vector<ParamVariant>& args) const override
		{
		}
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
			"ensuresize, newlike")
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
			"ensuresize, newlike")
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

		ClearCommand() : Command("clear", "Dispose image from the system (and free memory that the image consumes).",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image to erase. Specify empty name to clear everything.", "")
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
		return "distribute, delaying, maxmemory, printscripts";
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

		MaxMemoryCommand() : Command("maxmemory", "Sets the maximum memory setting used in distributed processing. This command overrides the value read from the configuration file. The maximum memory is the amount of memory that can be used either on the local computer (Local distribution mode) or in a compute node (Slurm distribution mode).",
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
