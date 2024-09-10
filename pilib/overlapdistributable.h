#pragma once

#include "distributable.h"

namespace pilib
{

	/**
	Base class for processes that distribute using simple overlap.
	@param CMDBASE The parent class that is derived from Command.
	*/
	template<class CMDBASE>
	class OverlapDistributable : public CMDBASE, public Distributable
	{
	protected:
		friend class CommandList;

		OverlapDistributable(const string& name, const string& help, const vector<CommandArgumentBase>& extraArgs = {}, const string& seeAlso = "") : CMDBASE(name, help, extraArgs, seeAlso)
		{
		}

	public:
		/**
		Calculate and return required overlap margin.
		Ensure that possible output image and input image have correct size.
		*/
		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const = 0;


		virtual Vec3c getMargin(const vector<ParamVariant>& args) const override
		{
			return calculateOverlap(args);
		}

		virtual size_t getDistributionDirection1(const std::vector<ParamVariant>& args) const override
		{
			return 2;
		}

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual size_t getDistributionDirection3(const std::vector<ParamVariant>& args) const override
		{
			return 0;
		}

		virtual bool canDelay(const vector<ParamVariant>& args) const override
		{
			return true;
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}
	};

}
