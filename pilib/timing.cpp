
#include "timing.h"
#include <sstream>
#include <iomanip>

using namespace std;

namespace pilib
{
	std::map<TimeClass, double> Timing::times = {
		{TimeClass::JobsInclQueuing, 0.0},
		{TimeClass::WriteFinalizationInclQueuing, 0.0},
		{TimeClass::WritePreparation, 0.0}
	};

	void Timing::Add(TimeClass timeClass, double seconds)
	{
		times[timeClass] += seconds;
	}

	std::string Timing::toString()
	{
		stringstream s;

		s << "Jobs incl. queuing: " << setprecision(2) << times[TimeClass::JobsInclQueuing] << " s" << endl
			<< "Write preparation: " << setprecision(2) << times[TimeClass::WritePreparation] << " s" << endl
			<< "Write finalization incl. queuing: " << setprecision(2) << times[TimeClass::WriteFinalizationInclQueuing] << " s";

		return s.str();
	}
}
