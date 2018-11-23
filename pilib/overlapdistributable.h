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
	public:

		OverlapDistributable(const string& name, const string& help, const vector<CommandArgumentBase>& extraArgs = {}) : CMDBASE(name, help, extraArgs)
		{
		}

		/**
		Calculate and return required overlap margin.
		Ensure that possible output image and input image have correct size.
		*/
		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const = 0;

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			Vec3c margin = calculateOverlap(args);

			// distribute in z, use overlap
			distributor.distribute(this, args, 2, margin, 0);
		}
	};

}
