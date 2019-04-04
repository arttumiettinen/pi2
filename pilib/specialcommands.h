#pragma once

#include "command.h"
#include "distributable.h"
#include "trivialdistributable.h"

namespace pilib
{
	
	class HelloCommand : virtual public Command, public TrivialDistributable
	{
	public:

		HelloCommand() : Command("hello", "Shows greetings message. This command is used for testing things.",
		{
			CommandArgument<string>(ParameterDirection::In, "name", "Name of caller. Omit to output a generic greeting.", "there")
		})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			cout << "Hello " << pop<string>(args) << "!" << endl;
		}
	};

	class WaitReturnCommand : virtual public Command, public TrivialDistributable
	{
	public:

		WaitReturnCommand() : Command("waitreturn", "Waits until the user presses return key. Does not show any prompts on screen.")
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			//cout << "Press return to exit..." << endl;
			char dummy;
			cin.getline(&dummy, 1);
		}
	};

	class InfoCommand : virtual public Command, public TrivialDistributable
	{
	public:

		InfoCommand() : Command("info", "Displays information about the computer and the PI system.")
		{
		}

		virtual void run(vector<ParamVariant>& args) const;
	};

	class LicenseCommand : virtual public Command, public TrivialDistributable
	{
	public:

		LicenseCommand() : Command("license", "Displays license of this program and PI system.")
		{
		}

		virtual void run(vector<ParamVariant>& args) const;
	};

	class ListCommand : virtual public Command, public Distributable
	{
	public:

		ListCommand() : Command("list", "Displays information about images loaded to memory, or images prepared for distributed processing.")
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		

		virtual void run(vector<ParamVariant>& args) const
		{

		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const;
	};

	class ReadCommand : virtual public Command, public Distributable
	{
	public:

		ReadCommand() : Command("read", "Reads an image or image sequence from disk. Determines type of file automatically.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read or a sequence definition. If a directory name is given, all files in the directory will be read. If the a file name is given, it can contain wildcards *, ? and @. Wildcard * corresponds to any sequence of characters, wildcard ? corresponds to any character, and wildcard @ corresponds to sequence of numerical digits. For example, sequence containing files xyz_000.png, xyz_001.png, xyz_002.png, etc. could be read with template 'xyz_@.png'.")
			})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const;
	};

	class ReadRawCommand : virtual public Command, public Distributable
	{
	private:
		static void parseArgs(vector<ParamVariant>& args, string& name, string& fname, coord_t& w, coord_t& h, coord_t& d, ImageDataType& dt);

	public:

		ReadRawCommand() : Command("readraw", "Reads a .raw image from a file. The image must be in native byte order of the host computer.",
		{
			CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
			CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
			CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be uint8, uint16, uint32, uint64, float32, or complex32. Specify empty value to infer data type from file size.", ""),
			CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the image. Omit width, height and depth to infer dimensions from file name.", 0),
			CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the image. Omit width, height and depth to infer dimensions from file name.", 0),
			CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the image. Omit width, height and depth to infer dimensions from file name.", 0)
		})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const;
	};

	class ReadSequenceCommand : public Command, public Distributable
	{
	public:

		ReadSequenceCommand() : Command("readsequence", "Reads an image sequence from disk. Supports any image formats supported by the back end (at least .png).",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename template", "Name (and path) template of the sequence to be read. If template is directory name, all files in the directory will be read. If the template is file name, it can contain wildcards *, ? and @. Wildcard * corresponds to any sequence of characters, wildcard ? corresponds to any character, and wildcard @ corresponds to sequence of numerical digits. For example, sequence containing files xyz_000.png, xyz_001.png, xyz_002.png, etc. could be read with template 'xyz_@.png'."),
			})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const;
	};

	class ReadVolCommand : public Command
	{
	public:

		ReadVolCommand() : Command("readvol", "Reads a .vol image from a file.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read.")
			})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};

	class ReadBlockCommand : public Command
	{
	public:

		ReadBlockCommand() : Command("readblock", "Reads a block of an image sequence from a file.",
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

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};

	class ReadRawBlockCommand : public Command
	{
	public:

		ReadRawBlockCommand() : Command("readrawblock", "Reads a block of a .raw image from a file. The image must be in native byte order of the host computer.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "x", "X-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "y", "Y-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "z", "Z-coordinate of the first pixel to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block width", "Width of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block height", "Height of block to read."),
				CommandArgument<coord_t>(ParameterDirection::In, "block depth", "Depth of block to read."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be uint8, uint16, or float32. Specify empty value to infer data type from image dimensions", ""),
				CommandArgument<coord_t>(ParameterDirection::In, "total width", "Width of the image. Omit width, height and depth to infer dimensions from file name.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "total height", "Height of the image. Omit width, height and depth to infer dimensions from file name.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "total depth", "Depth of the image. Omit width, height and depth to infer dimensions from file name.", 0)
			})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};

	class ReadSequenceBlockCommand : public Command
	{
	public:

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

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};

	class MapRawCommand : public Command
	{
	public:

		MapRawCommand() : Command("mapraw", "Maps a .raw image file to an image. The image file does not need to exist. Changes made to the file are IMMEDIATELY reflected to disk, so there is no need to save or load the image. The image must be in native byte order of the host computer.",
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of image in the system."),
				CommandArgument<string>(ParameterDirection::In, "filename", "Name (and path) of file to map."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image (see command newimage). Specify empty value to infer data type from image dimensions.", ""),
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the image. Omit width, height and depth to infer dimensions from file name.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the image. Omit width, height and depth to infer dimensions from file name.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the image. Omit width, height and depth to infer dimensions from file name.", 0)
			})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}

		// TODO: This could be used to make non-temporary disk-stored images in distributed processing
	};


	class NewImageCommand : virtual public Command, public Distributable
	{
	public:
		NewImageCommand() : Command("newimage", "Creates a new, empty image.",
		{
			CommandArgument<string>(ParameterDirection::In, "image name", "Name of the image in the system."),
			CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be uint8, uint16, uint32, float32, or complex32."),
			CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the image.", 1),
			CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the image.", 1),
			CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the image.", 1)
		})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const;
	};


	

	class ClearCommand : virtual public Command, public Distributable
	{
	public:
		ClearCommand() : Command("clear", "Dispose image from the system (and free memory that the image consumes).",
		{
			CommandArgument<string>(ParameterDirection::In, "image name", "Name of image to erase. Specify empty name to clear everything.", "")
		})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const;
	};

	class HelpCommand : virtual public Command, public TrivialDistributable
	{
	public:
		HelpCommand() : Command("help", "Shows usage information.",
		{
			CommandArgument<string>(ParameterDirection::In, "topic", "Name of command whose help is to be retrieved. Specify nothing to show general help.", ""),
			CommandArgument<string>(ParameterDirection::In, "search term", "Search term. Specify something to show only those matches that contain this term. E.g. specify uint8 to show help only for functions that process uint8 images.", ""),
		})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};

	class EchoCommandsCommand : virtual public Command, public TrivialDistributable
	{
	public:
		EchoCommandsCommand() : Command("echo", "Enables or disables echoing of commands and timing information to output.",
		{
			CommandArgument<bool>(ParameterDirection::In, "commands", "Set to true to show commands that have been run.", true),
			CommandArgument<bool>(ParameterDirection::In, "timing", "Set to true to show timing information on commands that have been run.", false),
		})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};

	class DistributeCommand : virtual public Command, public TrivialDistributable
	{
	public:
		DistributeCommand() : Command("distribute", "Enables or disables distributed processing of commands. Run this command before commands that you would like to run using distributed processing. Images used during distributed processing are not available for local processing and vice versa, unless they are loaded again from disk. All commands do not support distributed processing.",
			{
				CommandArgument<string>(ParameterDirection::In, "workload manager system name", "Set to SLURM to use SLURM workload manager; set to Local to process tasks sequentially using the local computer; set to empty string to disable distributed processing (default).", "")
			})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};
}
