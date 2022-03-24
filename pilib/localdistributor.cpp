
#include "localdistributor.h"

#include "stringutils.h"
#include "exeutils.h"

#include <algorithm>
#include "filesystem.h"

using namespace itl2;
using namespace std;

namespace pilib
{
	LocalDistributor::LocalDistributor(PISystem* piSystem) : Distributor(piSystem), allowedMem(0)
	{
		fs::path configPath = getPiCommand();
		size_t mem = 0;
		if (configPath.has_filename())
		{
			configPath = configPath.replace_filename("local_config.txt");

			INIReader reader(configPath.string());

			mem = (size_t)(reader.get<double>("max_memory", 0) * 1024 * 1024);
			
			readSettings(reader);
		}

		allowedMemory(mem);
	}

	void LocalDistributor::allowedMemory(size_t maxMem)
	{
		allowedMem = maxMem;

		if (allowedMem <= 0)
			allowedMem = (size_t)(0.85 * itl2::memorySize());

		cout << "Using " << bytesToString((double)allowedMem) << " RAM per task." << endl;
	}

	void LocalDistributor::submitJob(const string& piCode, JobType jobType)
	{
		// Write the code to (temporary) file
		{
			ofstream f("pi2_local_job.txt");
			f << piCode << endl;
			f << "print(Everything done.)" << endl;
		}

		string output = execute(getPiCommand(), "pi2_local_job.txt", true);

		outputs.push_back(output);
	}

	vector<string> LocalDistributor::waitForJobs()
	{
		// Do not wait as the jobs have already finished (they are run sequentially on submit).

		ostringstream msg;
		for(size_t n = 0; n < outputs.size(); n++)
		{
			string line = lastLine(outputs[n]);

			if (startsWith(line, "Error"))
			{
				msg << "Job " << n << " failed with message '" << line << "'" << endl;
			}
			else if (line != "Everything done.")
			{
				msg << "Job " << n << " failed without error message." << endl;
			}
		}

		string s = msg.str();
		if (s.length() > 0)
			throw ITLException(s.substr(0, s.length() - 1));

		vector<string> result = outputs;
		outputs.clear();
		return result;
	}

}
