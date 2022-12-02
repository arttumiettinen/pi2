#pragma once

#include <map>
#include <string>
#include <timer.h>
#include "utilities.h"

namespace pilib
{
	enum class TimeClass
	{
		/**
		General overhead, e.g. parsing inputs, finding correct commands to run etc.
		This includes the total overhead time spent in the main process, and in possible cluster job processes.
		*/
		Overhead,
		/**
		Time spent in I/O-bound processing. This is the time when the disk I/O is the bottleneck.
		This includes the total I/O time spent in the main process, and in possible cluster job processes.
		Time spent in output data compression is counted to this time class.
		*/
		IO,
		/**
		Time spent in CPU/GPU-bound processing. This is the time when the CPU/GPU is the bottleneck.
		This includes the total computation time spent in the main process, and in possible cluster job processes.
		This is the default mode for all commands.
		*/
		Computation,

		/**
		Total distributed job execution time.
		This value includes Overhead+IO+Computation of all jobs, plus workload manager job starting overhead, time spent in starting pi2 process etc.
		This value does not include time spent in workload manager queue.
		*/
		JobExecution,
		/**
		Total distributed job queuing time.
		This is the total time all jobs have spent in the workload manager queue, waiting to be executed.
		*/
		JobQueueing,
		/**
		Total time from starting of all the distributed jobs until they are found to be finished.
		This is the total time spent in the submitting process, from submission to jobs until all of them are done.
		*/
		JobsInclQueuing,
		/**
		Time spent in preparing for writing output images (e.g. NN5 write preparation).
		This is the total time spent in the submitting process while preparing writing of output images.
		*/
		WritePreparation,
		/**
		Time spent in write finalization jobs, including queuing (e.g. NN5 write finalization jobs).
		This is the total time spent in the submitting process, from submission of write finalization jobs until all of them are done.
		*/
		WriteFinalizationInclQueuing
	};

	/*
	Convert time class to string.
	*/
	static inline string toString(TimeClass t)
	{
		if (t == TimeClass::Overhead)
			return "Overhead";
		if (t == TimeClass::IO)
			return "I/O";
		if (t == TimeClass::Computation)
			return "Computation";
		if (t == TimeClass::JobExecution)
			return "Job execution";
		if (t == TimeClass::JobQueueing)
			return "Job queuing";
		if (t == TimeClass::JobsInclQueuing)
			return "Total job waiting";
		if (t == TimeClass::WritePreparation)
			return "Write preparation";
		if (t == TimeClass::WriteFinalizationInclQueuing)
			return "Total write finalization waiting";
		throw std::runtime_error("Not implemented");
	}
}

namespace itl2
{
	template<>
	inline pilib::TimeClass fromString(const std::string& dt)
	{
		std::string dt2 = dt;
		trim(dt2);
		toLower(dt2);
		if (dt2 == "overhead")
			return pilib::TimeClass::Overhead;
		if (dt2 == "i/o")
			return pilib::TimeClass::IO;
		if (dt2 == "computation")
			return pilib::TimeClass::Computation;
		if (dt2 == "job execution")
			return pilib::TimeClass::JobExecution;
		if (dt2 == "job queuing")
			return pilib::TimeClass::JobQueueing;
		if (dt2 == "total job waiting time")
			return pilib::TimeClass::JobsInclQueuing;
		if (dt2 == "write preparation")
			return pilib::TimeClass::WritePreparation;
		if (dt2 == "total write finalization waiting")
			return pilib::TimeClass::WriteFinalizationInclQueuing;

		throw ITLException("Invalid job type: " + dt);
	}
}

namespace pilib
{
	/**
	Class that stores timing results.
	*/
	class TimingResults
	{
	private:
		std::map<TimeClass, double> times;
	public:
		/**
		Add the given amount of seconds to timing of given operation class.
		*/
		void add(TimeClass timeClass, double seconds);

		/**
		Retrieve a timing report as a string.
		*/
		std::string toString() const;

		/**
		Parse TimingResults from a string generated with toString method.
		*/
		static TimingResults parse(const std::string& data);

		/**
		Add the given timing results to this instance.
		*/
		void accumulate(const TimingResults& r);

		/**
		Save timing results to a file.
		*/
		void toFile(const std::string& filename) const;

		/**
		Load timing results from a file.
		Throws an exception if an I/O error is encountered.
		*/
		static TimingResults fromFile(const std::string& filename);
	};

	/**
	Class used to measure I/O and processing time separately.
	*/
	class GlobalTimer
	{
	private:
		static TimingResults tresults;
		static itl2::Timer timer;
		static TimeClass mode;
	public:
		/**
		Sets current timing mode to a new value.
		When the timing mode changes, the time accumulated so far is added to the total time spent in the previous active timing mode.
		*/
		static void setMode(TimeClass mode);

		/**
		Starts the timing in Overhead mode.
		*/
		static void start();
		
		/**
		Stops timing.
		*/
		static void stop();

		/**
		Gets a reference to the timing results.
		*/
		static inline TimingResults& results()
		{
			return tresults;
		}

		/**
		Gets the active timing mode.
		*/
		static inline TimeClass currentMode()
		{
			return mode;
		}
	};

	/**
	RAII class that is used to change timing mode for a scope.
	*/
	class TimingFlag
	{
	private:
		TimeClass oldMode;
	public:
		TimingFlag(TimeClass mode);

		virtual ~TimingFlag();

		/**
		Change current timing mode for the scope.
		*/
		void changeMode(TimeClass mode);
	};
}