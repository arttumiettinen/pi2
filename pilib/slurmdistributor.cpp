
#include "slurmdistributor.h"

#include "exeutils.h"
#include "inireader.h"
#include "utilities.h"
#include "stringutils.h"
#include "math/vectoroperations.h"

#if defined(__linux__)
    #include <sys/types.h>
       #include <dirent.h>
#endif

#include <algorithm>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

using namespace itl2;

namespace pilib
{
	

	SLURMDistributor::SLURMDistributor(PISystem* system) : Distributor(system), allowedMem(0)
	{
		// Read config file, try to read from the current folder and from the folder of pi2 executable.
		fs::path configPath = getPiCommand();
		if(configPath.has_filename())
			configPath = configPath.replace_filename("slurm_config.txt");

		INIReader reader(configPath.string());

		extraArgs = reader.get<string>("extra_args", "");
		jobInitCommands = reader.get<string>("job_init_commands", "");
		allowedMem = reader.get<size_t>("max_memory", 0) * 1024 * 1024;
		
		if(allowedMem <= 0)
		{
		    // Run cluster info command to get amount of memory in megabytes
		    string cmd = string("sinfo --Node --Format=memory ") + extraArgs;
		    string output = execute(cmd);

		    // Split lines
		    trim(output);
		    vector<string> lines = split(output, false);

		    // Remove header
		    lines.erase(lines.begin() + 0);

		    // Convert list to list of ints
		    vector<size_t> values = fromString<size_t>(lines);

		    // Calculate minimum in bytes
		    if(values.size() > 0)
    		    allowedMem = (size_t)round(math::min(values) * 1024 * 1024 * 0.9);
		    
		    if(allowedMem <= 0)
		        throw ITLException(string("Unable to determine maximum memory available in the SLURM cluster. Are extra parameters in configuration file correct. The file is read from ") + configPath.string());
		}
		
		cout << "Per node memory in the SLURM cluster: " << bytesToString((double)allowedMem) << endl;
	}

	string makeJobName(size_t jobIndex)
	{
		return "pi2-" + itl2::toString<size_t>(jobIndex);
	}
	
	string makeInputName(size_t jobIndex)
	{
	    return "./slurm-io-files/" + makeJobName(jobIndex) + "-in.txt";
	}

	string makeOutputName(size_t jobIndex)
	{
		return "./slurm-io-files/" + makeJobName(jobIndex) + "-out.txt";
	}
	
	string makeErrorName(size_t jobIndex)
	{
		return "./slurm-io-files/" + makeJobName(jobIndex) + "-err.txt";
	}

	void SLURMDistributor::submitJob(const string& piCode)
	{
		size_t jobIndex = submittedJobs.size();
		string jobName = makeJobName(jobIndex);
		
		string inputName = makeInputName(jobIndex);
        string outputName = makeOutputName(jobIndex);
        string errorName = makeErrorName(jobIndex);
        fs::remove(inputName);
        fs::remove(outputName);
        fs::remove(errorName);
        
        writeText(inputName, piCode);

		string jobCmdLine;
		if(jobInitCommands.length() > 0)
			jobCmdLine = jobInitCommands + "; ";
		jobCmdLine += getPiCommand() + " " + inputName;

		string sbatchArgs = "--job-name=" + jobName + " --output=" + outputName + " --error=" + errorName + " --wrap=\"" + jobCmdLine + "\" " + extraArgs;
		
		//cout << "sbatch input: " << sbatchArgs << endl;

		string result = execute("sbatch", sbatchArgs);
		
		//cout << "sbatch output: " << result << endl;
		
		vector<string> lines = split(result);
		if(lines.size() == 1)
		{
		    vector<string> parts = split(*lines.rbegin(), false, ' ');
		    if(parts.size() < 1)
		    {
		        cout << result << endl;
    		    throw ITLException("SLURM returned no batch job id.");
    		}
		    size_t slurmId = fromString<int>(parts[parts.size() - 1]);
		    submittedJobs.push_back(slurmId);
		}
		else
		{
		    cout << result << endl;
		    throw ITLException("Unexpected sbatch output. The output has been printed to standard output.");
		}
	}

	bool SLURMDistributor::isJobDone(size_t jobIndex) const
	{
		string slurmId = itl2::toString(submittedJobs[jobIndex]);

		string result = execute("squeue", string("--noheader --jobs=") + slurmId);
		trim(result);

		return !startsWith(result, slurmId);
	}
	
	void flushCache(const string& filename)
	{
#if defined(__linux__)
	    // This may flush NFS cache on the files of the folder where the log file lives.
	    int fd = open(filename.c_str(), O_RDONLY);
	    fsync(fd);
	    close(fd);
	    
	    fs::path dir(filename);
        dir = dir.parent_path();
        DIR* dr = opendir(dir.string().c_str());
        closedir(dr);
#endif
	}

	string SLURMDistributor::getLog(size_t jobIndex, bool flush) const
	{
	    string outputName = makeOutputName(jobIndex);
	    flushCache(outputName);
		return readText(outputName);
	}
	
	constexpr int JOB_FAILED = -2;
	constexpr int JOB_WAITING = -1;

	int SLURMDistributor::getJobProgressFromLog(size_t jobIndex) const
	{
		string log = getLog(jobIndex);
		if (log.length() <= 0)
			return JOB_WAITING;
		string line = lastLine(log);
		
		if (startsWith(line, "Error"))
			return JOB_FAILED;

        // NOTE: This counting method does not work if showProgress function is changed!
		return (int)std::count(line.begin(), line.end(), '=') * 10;
	}

	vector<int> SLURMDistributor::getJobProgress() const
	{
		vector<int> progress;
		progress.reserve(submittedJobs.size());
		for (size_t n = 0; n < submittedJobs.size(); n++)
		{
			int state;
			if (!isJobDone(n))
			{
				// Read progress from log file or -1 if it does not exist.
				state = getJobProgressFromLog(n);
				// Do not let the progress read from log ever reach 100 as that is interpreted
				// as job done.
				if(state >= 100)
    				state = 99;
			}
			else
			{
				// Job is ready (or failed)
				state = getJobProgressFromLog(n);
				if(state != JOB_FAILED)
					state = 100;
			}

			progress.push_back(state);
		}
		return progress;
	}

	char SLURMDistributor::getProgressChar(int progress) const
	{
		if (progress == JOB_FAILED)
			return '!';
		if (progress == JOB_WAITING)
			return 'W';


		// 20	40	60	80	100
		// _	.	o	O	*

		if (progress < 20)
			return '_';
		if (progress < 40)
			return '.';
		if (progress < 60)
			return 'o';
		if (progress < 80)
			return 'O';

		return '*';
	}

	vector<string> SLURMDistributor::waitForJobs()
	{
		vector<string> result;
		if (submittedJobs.size() <= 0)
			return result;

		// Wait until jobs are done
		bool done = false;
		do
		{
			vector<int> progress = getJobProgress();

			// Show progress report
			for (size_t n = 0; n < progress.size(); n++)
			{
				cout << getProgressChar(progress[n]);
			}
			cout << '\r' << flush;

			// Check if all jobs are done
			size_t doneCount = 0;
			for (size_t n = 0; n < progress.size(); n++)
			{
				if (progress[n] == JOB_FAILED || progress[n] >= 100)
					doneCount++;
			}

			done = doneCount == progress.size();
		} while (!done);

		// Remove progress bar
		for (size_t n = 0; n < submittedJobs.size(); n++)
			cout << ' ';
		cout << '\r' << flush;

		// Read logs and generate error if something went wrong.
		result.reserve(submittedJobs.size());
		ostringstream msg;
		for (size_t n = 0; n < submittedJobs.size(); n++)
		{
		    string log = getLog(n, true);
		    
		    if(log.length() <= 0)
		    {
		        msg << "Job " << n << " failed without running pi2. SLURM log " << makeErrorName(n) << " may contain further details. If it does not exist, please make sure that SLURM has read and write access to the current folder." << endl;
		    }
			else
			{
			    string line = lastLine(log);
			    if (startsWith(line, "Error"))
			    {
				    msg << "Job " << n << " failed with message '" << line << "'" << endl;
			    }
		    }
		    
			result.push_back(log);
		}
		submittedJobs.clear();
		
		string s = msg.str();
		if (s.length() > 0)
			throw ITLException(s.substr(0, s.length() - 1));

		return result;
	}


}
