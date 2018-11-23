
#include "localdistributor.h"

#include "stringutils.h"
#include "exeutils.h"
#include "inireader.h"

#include <algorithm>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;


namespace pilib
{
	LocalDistributor::LocalDistributor(PISystem* piSystem) : Distributor(piSystem)
	{
		fs::path configPath = getPiCommand();
		if (configPath.has_filename())
			configPath = configPath.replace_filename("local_config.txt");

		INIReader reader(configPath.string());

		allowedMem = reader.get<size_t>("max_memory", 0) * 1024 * 1024;

		if(allowedMem <= 0)
			allowedMem = (size_t)(0.85 * itl2::memorySize());

		cout << "Using " << bytesToString((double)allowedMem) << " RAM per task." << endl;
	}

	void LocalDistributor::submitJob(const string& piCode)
	{
		// Write the code to (temporary) file
		{
			ofstream f("pi2_local_job.txt");
			f << piCode;
		}

		string output = execute(getPiCommand(), "pi2_local_job.txt");

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
		}

		string s = msg.str();
		if (s.length() > 0)
			throw ITLException(s.substr(0, s.length() - 1));

		vector<string> result = outputs;
		outputs.clear();
		return result;
	}

}
