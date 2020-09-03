#pragma once

#include "commandsbase.h"
#include "trivialdistributable.h"

namespace pilib
{
	inline std::string metaSeeAlso()
	{
		return "line, capsule, sphere, ellipsoid, set, get, ramp";
	}

	template<typename input_t> class SetMetadataCommand : public OneImageInPlaceCommand<input_t>, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		SetMetadataCommand() : Command("setmeta", "Sets metadata item of an image.",
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
			const string& key = args[0];
			const string& value = args[1];
			in.metadata.set(key, value);
		}
	};


	template<typename input_t> class GetMetadataCommand : public OneImageInPlaceCommand<input_t>, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		GetMetadataCommand() : Command("getmeta", "Gets metadata item from an image.",
			{
				CommandArgument<string>(ParameterDirection::In, "key", "Name of the metadata item."),
			},
			metaSeeAlso())
		{
		}

	public:
		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			const string& key = args[0];
			const string& value = args[1];
			in.metadata.set(key, value);
		}
	};

}