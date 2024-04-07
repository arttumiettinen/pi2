
#include "timing.h"
#include <sstream>
#include <iomanip>
#include "stringutils.h"

using namespace std;

namespace pilib
{

	TimingResults GlobalTimer::tresults = TimingResults();
	itl2::Timer GlobalTimer::timer = itl2::Timer();
	TimeClass GlobalTimer::mode = TimeClass::Overhead;

	void GlobalTimer::setMode(TimeClass mode)
	{
		double elapsed = timer.lap();
		tresults.add(GlobalTimer::mode, elapsed);
		GlobalTimer::mode = mode;
	}

	void GlobalTimer::start()
	{
		mode = TimeClass::Overhead;
		timer.start();
	}

	void GlobalTimer::stop()
	{
		timer.stop();
		double elapsed = timer.getSeconds();
		tresults.add(GlobalTimer::mode, elapsed);
	}

	void GlobalTimer::reset()
	{
		tresults = TimingResults();
	}


	TimingFlag::TimingFlag(TimeClass mode) :
		oldMode(GlobalTimer::currentMode())
	{
		GlobalTimer::setMode(mode);
	}

	TimingFlag::~TimingFlag()
	{
		GlobalTimer::setMode(oldMode);
	}

	void TimingFlag::changeMode(TimeClass mode)
	{
		GlobalTimer::setMode(mode);
	}

	void TimingResults::add(TimeClass timeClass, double seconds)
	{
		times[timeClass] += seconds;
	}

	std::string TimingResults::toString() const
	{
		stringstream s;

		bool first = true;
		for (const auto& entry : times)
		{
			if (!first)
			{
				s << endl;
			}
			first = false;
			s << pilib::toString(entry.first) << ": " << entry.second << " s";
		}

		return s.str();
	}


	TimingResults TimingResults::parse(const std::string& data)
	{
		TimingResults result;
		std::vector<string> lines = itl2::split(data, false);
		for (const string& line : lines)
		{
			std::vector<string> parts = itl2::split(line, true, ':');
			if (parts.size() == 2)
			{
				TimeClass mode = itl2::fromString<TimeClass>(parts[0]);
				double time = itl2::fromString<double>(parts[1]);
				result.times[mode] = time;
			}
		}
		return result;
	}

	void TimingResults::accumulate(const TimingResults& r)
	{
		for (const auto& entry : r.times)
		{
			times[entry.first] += entry.second;
		}
	}

	void TimingResults::toFile(const std::string& filename) const
	{
		itl2::writeText(filename, toString());
	}

	TimingResults TimingResults::fromFile(const std::string& filename)
	{
		return TimingResults::parse(itl2::readText(filename, true));
	}
}
