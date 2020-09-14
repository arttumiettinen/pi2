#pragma once

#include "commandsbase.h"
#include "trivialdistributable.h"

namespace pilib
{
	inline std::string metaSeeAlso()
	{
		return "setmeta, getmeta";
	}

	template<typename input_t> class SetMetadataCommand : public OneImageInPlaceCommand<input_t>
	{
	protected:
		friend class CommandList;

		SetMetadataCommand() : OneImageInPlaceCommand<input_t>("setmeta", "Sets metadata item of an image.",
			{
				CommandArgument<string>(ParameterDirection::In, "key", "Name of the metadata item."),
				CommandArgument<string>(ParameterDirection::In, "value", "Value of the metadata item."),
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			const string& key = std::get<string>(args[0]);
			const string& value = std::get<string>(args[1]);
			in.metadata.set(key, value);
		}
	};


	template<typename input_t> class GetMetadataCommand : public OneImageInPlaceCommand<input_t>
	{
	protected:
		friend class CommandList;

		GetMetadataCommand() : OneImageInPlaceCommand<input_t>("getmeta", "Gets metadata item from an image.",
			{
				CommandArgument<string>(ParameterDirection::In, "key", "Name of the metadata item."),
				CommandArgument<string>(ParameterDirection::Out, "value", "Value of the metadata item is placed into this string."),
				CommandArgument<string>(ParameterDirection::In, "default", "This value is returned if the key is not found", ""),
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			const string& key = std::get<string>(args[0]);
			string* value = std::get<string*>(args[1]);
			const string& def = std::get<string>(args[2]);
			*value = in.metadata.get(key, def);
		}
	};

}