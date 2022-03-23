#pragma once

#include <map>
#include <string>

namespace pilib
{
	enum class TimeClass
	{
		JobsInclQueuing,
		WritePreparation,
		WriteFinalizationInclQueuing
	};

	class Timing
	{
	private:
		static std::map<TimeClass, double> times;
	public:
		/**
		Add the given amount of seconds to timing of given operation class.
		*/
		static void Add(TimeClass timeClass, double seconds);

		/**
		Retrieve a timing report as a string.
		*/
		static std::string toString();
	};
}