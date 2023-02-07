#pragma once

#include "commandsbase.h"
#include "trivialdistributable.h"

namespace pilib
{
	inline std::string metaSeeAlso()
	{
		return "setmeta, getmeta, writemeta, readmeta, clearmeta, listmeta, copymeta, metarowcount, metacolumncount";
	}

	template<typename input_t, typename output_t> class CopyMetadataCommand : public TwoImageInputOutputCommand<input_t, output_t>
	{
	protected:
		friend class CommandList;

		CopyMetadataCommand() : TwoImageInputOutputCommand<input_t, output_t>("copymeta", "Copies metadata from one image to another.",
			{
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(Image<input_t>& in, Image<output_t>& out, std::vector<ParamVariant>& args) const override
		{
			out.metadata = in.metadata;
		}
	};


	template<typename input_t> class MetaRowCountCommand : public OneImageCommand<input_t>
	{
	protected:
		friend class CommandList;

		MetaRowCountCommand() : OneImageCommand<input_t>("metarowcount", "Retrieves count of rows in a metadata item.",
			{
				CommandArgument<string>(ParameterDirection::In, "key", "Name of the metadata item."),
				CommandArgument<coord_t>(ParameterDirection::Out, "count", "Item count.")
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(const Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			const string& key = std::get<string>(args[0]);
			coord_t* out = std::get<coord_t*>(args[1]);
			*out = in.metadata.rowCount(key);
		}
	};


	template<typename input_t> class MetaColumnCountCommand : public OneImageCommand<input_t>
	{
	protected:
		friend class CommandList;

		MetaColumnCountCommand() : OneImageCommand<input_t>("metacolumncount", "Retrieves count of columns in a specific row of a metadata item.",
			{
				CommandArgument<string>(ParameterDirection::In, "key", "Name of the metadata item."),
				CommandArgument<coord_t>(ParameterDirection::In, "row index", "Index of the row whose column count is to be returned."),
				CommandArgument<coord_t>(ParameterDirection::Out, "count", "Item count.")
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(const Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			const string& key = std::get<string>(args[0]);
			coord_t row = std::get<coord_t>(args[1]);
			coord_t* out = std::get<coord_t*>(args[1]);
			*out = in.metadata.columnCount(key, row);
		}
	};


	template<typename input_t> class SetMetadataCommand : public OneImageInPlaceCommand<input_t>
	{
	protected:
		friend class CommandList;

		SetMetadataCommand() : OneImageInPlaceCommand<input_t>("setmeta", "Sets metadata item of an image.",
			{
				CommandArgument<string>(ParameterDirection::In, "key", "Name of the metadata item."),
				CommandArgument<string>(ParameterDirection::In, "value", "Value of the metadata item."),
				CommandArgument<size_t>(ParameterDirection::In, "i", "Row index of the item to set in the data matrix.", 0),
				CommandArgument<size_t>(ParameterDirection::In, "j", "Column index of the item to set in the data matrix.", 0),
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			const string& key = std::get<string>(args[0]);
			const string& value = std::get<string>(args[1]);
			size_t i = std::get<size_t>(args[2]);
			size_t j = std::get<size_t>(args[3]);
			in.metadata.set(key, value, i, j);
		}
	};


	template<typename input_t> class GetMetadataCommand : public OneImageCommand<input_t>
	{
	protected:
		friend class CommandList;

		GetMetadataCommand() : OneImageCommand<input_t>("getmeta", "Gets metadata item from an image.",
			{
				CommandArgument<string>(ParameterDirection::In, "key", "Name of the metadata item."),
				CommandArgument<string>(ParameterDirection::Out, "value", "Value of the metadata item is placed into this string."),
				CommandArgument<size_t>(ParameterDirection::In, "i", "Row index of the item to retrieve from the data matrix.", 0),
				CommandArgument<size_t>(ParameterDirection::In, "j", "Column index of the item to retrieve from the data matrix.", 0),
				CommandArgument<string>(ParameterDirection::In, "default", "This value is returned if the key is not found", ""),
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(const Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			const string& key = std::get<string>(args[0]);
			string* value = std::get<string*>(args[1]);
			size_t i = std::get<size_t>(args[2]);
			size_t j = std::get<size_t>(args[3]);
			const string& def = std::get<string>(args[4]);

			*value = in.metadata.get(key, def, i, j);
		}
	};


	template<typename pixel_t> class WriteMetadataCommand : public OneImageCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		WriteMetadataCommand() : OneImageCommand<pixel_t>("writemeta", "Writes image metadata to disk.",
			{
				CommandArgument<string>(ParameterDirection::In, "filename", "Name and path of file to write. The file will be replaced."),
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(const Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			const string& filename = std::get<string>(args[0]);
			in.metadata.writeToFile(filename);
		}
	};

	template<typename pixel_t> class ReadMetadataCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		ReadMetadataCommand() : OneImageInPlaceCommand<pixel_t>("readmeta", "Reads image metadata from disk.",
			{
				CommandArgument<string>(ParameterDirection::In, "filename", "Name and path of file to read from."),
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			const string& filename = std::get<string>(args[0]);
			in.metadata.readFromFile(filename);
		}
	};


	template<typename pixel_t> class ClearMetadataCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		ClearMetadataCommand() : OneImageInPlaceCommand<pixel_t>("clearmeta", "Clears metadata of an image",
			{
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			in.metadata.clear();
		}
	};


	template<typename pixel_t> class ListMetadataCommand : public OneImageCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		ListMetadataCommand() : OneImageCommand<pixel_t>("listmeta", "Builds a comma-separated list of the names of all the metadata items in an image.",
			{
				CommandArgument<string>(ParameterDirection::Out, "names", "Names of metadata items will be stored in this string."),
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(const Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			string& s = *std::get<string*>(args[0]);
			s = "";
			const auto& c = in.metadata.keys();
			for (size_t n = 0; n < c.size(); n++)
			{
				if (n > 0)
					s += ", ";
				s += c[n];
			}
		}
	};
}