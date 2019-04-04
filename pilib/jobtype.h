#pragma once

namespace pilib
{
	/**
	Classifies calculation jobs according to expected runtime.
	*/
	enum class JobType
	{
		Fast,
		Normal,
		Slow
	};

	inline string toString(const JobType type)
	{
		switch (type)
		{
		case JobType::Fast: return "fast";
		case JobType::Normal : return "normal";
		case JobType::Slow: return "slow";
		default: throw ITLException("Invalid job type.");
		}
	}
}
