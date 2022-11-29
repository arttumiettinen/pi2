#pragma once

#include <map>
#include <string>

namespace pilib
{
	enum class TimeClass
	{
		/**
		Total distributed job execution time.
		*/
		JobExecution,
		/**
		Total distributed job queuing time.
		*/
		JobQueueing,
		/**
		Total time from starting of all the distributed jobs until they are found to be finished.
		*/
		JobsInclQueuing,
		/**
		Time spent in preparing for writing output images (e.g. NN5 write preparation).
		*/
		WritePreparation,
		/**
		Time spent in write finalization jobs, including queuing (e.g. NN5 write finalization jobs).
		*/
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