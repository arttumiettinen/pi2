
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

		s << "Waiting for jobs to finish: " << setprecision(2) << times[TimeClass::JobsInclQueuing] << " s" << endl
			<< "Write preparation: " << setprecision(2) << times[TimeClass::WritePreparation] << " s" << endl
			<< "Waiting for write finalization: " << setprecision(2) << times[TimeClass::WriteFinalizationInclQueuing] << " s" << endl
			<< "Total job queuing time: " << setprecision(2) << times[TimeClass::JobQueueing] << " s" << endl
			<< "Total job execution time: " << setprecision(2) << times[TimeClass::JobExecution] << " s" << endl;

		return s.str();
	}
}
