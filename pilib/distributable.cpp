
#include "distributable.h"

namespace pilib
{
	void Distributable::runDistributedInternal(PISystem* system, Distributor& distributor, vector<ParamVariant>& args) const
	{
		runDistributed(distributor, args);
	}
}
