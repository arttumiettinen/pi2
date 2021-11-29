#pragma once

#include <string>
#include "utilities.h"

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
		case JobType::Normal: return "normal";
		case JobType::Slow: return "slow";
		default: throw ITLException("Invalid job type.");
		}
	}

}

namespace itl2
{
	template<>
	inline pilib::JobType fromString(const std::string& dt)
	{
		std::string dt2 = dt;
		trim(dt2);
		toLower(dt2);
		if (dt2 == "fast")
			return pilib::JobType::Fast;
		if (dt2 == "normal")
			return pilib::JobType::Normal;
		if (dt2 == "slow")
			return pilib::JobType::Slow;

		throw ITLException("Invalid job type: " + dt);
	}
}
