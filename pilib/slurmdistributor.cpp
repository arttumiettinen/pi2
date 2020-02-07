
#include "slurmdistributor.h"

#include "exeutils.h"
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

#include <random>

using namespace itl2;
using namespace std;

namespace pilib
{
	

	SLURMDistributor::SLURMDistributor(PISystem* system) : Distributor(system), allowedMem(0)
	{
		// Create unique name for this distributor
		std::random_device dev;
		std::mt19937 rng(dev());
		std::uniform_int_distribution<std::mt19937::result_type> dist(0, 10000);
		myName = itl2::toString(dist(rng));

		// Read config file, try to read from the current folder and from the folder of pi2 executable.
		fs::path configPath = getPiCommand();
		size_t mem = 0;
		if (configPath.has_filename())
		{		
			configPath = configPath.replace_filename("slurm_config.txt");
			
//			cout << "Reading settings from " << configPath << endl;

			INIReader reader(configPath.string());

			extraArgsFastJobs = reader.get<string>("extra_args_fast_jobs", "");
			extraArgsNormalJobs = reader.get<string>("extra_args_normal_jobs", "");
			extraArgsSlowJobs = reader.get<string>("extra_args_slow_jobs", "");
			jobInitCommands = reader.get<string>("job_init_commands", "");
			mem = (size_t)(reader.get<double>("max_memory", 0) * 1024 * 1024);
			maxSubmissions = reader.get<size_t>("max_resubmit_count", 5) + 1;
			sbatchCommand = reader.get<string>("sbatch_command", "sbatch");
			squeueCommand = reader.get<string>("squeue_command", "squeue");
			scancelCommand = reader.get<string>("scancel_command", "scancel");
			sinfoCommand = reader.get<string>("sinfo_command", "sinfo");

			readSettings(reader);
		}
		
		allowedMemory(mem);
	}

	void SLURMDistributor::allowedMemory(size_t maxMem)
	{
		allowedMem = maxMem;

		if (allowedMem <= 0)
		{

		}

		if (allowedMem <= 0)
		{
			// Run cluster info command to get amount of memory in megabytes
			//string cmd = string("sinfo --Node --Format=freemem ") + extraArgs; // Memory available to start new programs. Jobs running on nodes may affect this.
			string cmd = sinfoCommand + string(" --Node --Format=memory ") + extraArgsNormalJobs; // Total installed memory. Not all of that may be available!
			string output = execute(cmd);

			// Split lines
			trim(output);
			vector<string> lines = split(output, false);

			// Remove header
			lines.erase(lines.begin() + 0);

			// Convert list to list of ints
			vector<size_t> values = fromString<size_t>(lines);

			// Calculate minimum in bytes
			if (values.size() > 0)
				allowedMem = (size_t)std::round(min(values) * 1024 * 1024 * 0.75);

			if (allowedMem <= 0)
				throw ITLException("Unable to determine maximum memory available in the SLURM cluster. Are the extra_args_* parameters in the slurm_config.txt file correct?");
					//The file is read from ") + configPath.string());
		}

		cout << "Per node memory in the SLURM cluster: " << bytesToString((double)allowedMem) << endl;
	}

	string SLURMDistributor::makeJobName(size_t jobIndex) const
	{
		return "pi2-" + itl2::toString<size_t>(jobIndex) + "-" + myName;
	}
	
	string SLURMDistributor::makeInputName(size_t jobIndex) const
	{
	    return "./slurm-io-files/" + makeJobName(jobIndex) + "-in.txt";
	}

	string SLURMDistributor::makeOutputName(size_t jobIndex) const
	{
		return "./slurm-io-files/" + makeJobName(jobIndex) + "-out.txt";
	}
	
	string SLURMDistributor::makeErrorName(size_t jobIndex) const
	{
		return "./slurm-io-files/" + makeJobName(jobIndex) + "-err.txt";
	}

	void SLURMDistributor::resubmit(size_t jobIndex)
	{
		string jobName = makeJobName(jobIndex);
		string inputName = makeInputName(jobIndex);
		string outputName = makeOutputName(jobIndex);
		string errorName = makeErrorName(jobIndex);
		JobType jobType = get<1>(submittedJobs[jobIndex]);

		fs::remove(outputName);
		fs::remove(errorName);

		string jobCmdLine;
		if (jobInitCommands.length() > 0)
			jobCmdLine = jobInitCommands + "; ";
		jobCmdLine += "'" + getPiCommand() + "' " + inputName;

		string sbatchArgs = string("--no-requeue") + " --job-name=" + jobName + " --output=" + outputName + " --error=" + errorName + " " + extraArgs(jobType) + " --wrap=\"" + jobCmdLine + "\"";

	    //cout << "sbatch input: " << sbatchArgs << endl;

		string result = execute("sbatch", sbatchArgs);

		//cout << "sbatch output: " << result << endl;

		vector<string> lines = split(result);
		if (lines.size() == 1)
		{
			vector<string> parts = split(*lines.rbegin(), false, ' ');
			if (parts.size() < 1)
			{
				cout << result << endl;
				throw ITLException("SLURM returned no batch job id.");
			}
			size_t slurmId = fromString<int>(parts[parts.size() - 1]);
			get<0>(submittedJobs[jobIndex]) = slurmId;
			get<2>(submittedJobs[jobIndex])++;
		}
		else
		{
			if (result.length() > 0)
			{
				cout << "Command" << endl;
				cout << "sbatch " << sbatchArgs << endl;
				cout << "returned" << endl;
				cout << result << endl;
				throw ITLException("Unexpected sbatch output. The output has been printed to standard output.");
			}
			else
			{
				throw ITLException("Empty sbatch output.");
			}
		}
	}

	void SLURMDistributor::submitJob(const string& piCode, JobType jobType)
	{
		// Add job completion marker
		string piCode2 = piCode + "\nprint(Everything done.);\n";

		// Create slot for the job
		size_t jobIndex = submittedJobs.size();
		submittedJobs.push_back(make_tuple(0, jobType, 0));

		// Write input file
		string inputName = makeInputName(jobIndex);
		fs::remove(inputName);
		writeText(inputName, piCode2);

		// Submit "again"
		resubmit(jobIndex);
	}

	bool SLURMDistributor::isJobDone(size_t jobIndex) const
	{
		string slurmId = itl2::toString(get<0>(submittedJobs[jobIndex]));

		string result = execute(squeueCommand, string("--noheader --jobs=") + slurmId);
		trim(result);
		
		// Empty result means the job is done
		if(result.length() <= 0)
		    return true;

        // If the message starts with job id, the job is running or in queue.
		bool running = startsWith(result, slurmId);
		
		// If the message contains "invalid job id...", the job data has been already erased.
		bool notRunning = contains(result, "Invalid job id specified");

        // Check that the job is running or its data has been erased.
        // Otherwise we don't know what is going on, and assume this is a temporary error message.
		if ((running && notRunning) ||
			(!running && !notRunning))
		{
			// Erroneous squeue output. Assume that the job is not done.
			cout << "Warning: Unexpected squeue output '" << result << "' for job " << jobIndex << " (SLURM id " << slurmId << "). Assuming the job is still running." << endl;
			return false;
		}

		return !running;
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

	string SLURMDistributor::getSlurmErrorLog(size_t jobIndex, bool flush) const
	{
		string outputName = makeErrorName(jobIndex);
		flushCache(outputName);
		return readText(outputName);
	}
	
	constexpr int JOB_FAILED = -2;
	constexpr int JOB_WAITING = -1;

	int SLURMDistributor::getJobProgressFromLog(size_t jobIndex) const
	{
	    // Check if the job has been cancelled?
        string slurmLog = getSlurmErrorLog(jobIndex);
        if(contains(slurmLog, "CANCELLED"))
            return JOB_FAILED;
	
	    // Get progress from output
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
       			// The job is not in squeue anymore
				// Job is ready or failed
				state = getJobProgressFromLog(n);
				
				// If the job is in waiting state something is wrong!
				if (state == JOB_WAITING)
				{
					state = JOB_FAILED;
				}
				else
				{
					// Check that the job log ends with success message, if not set job to failed state.
					string last = lastLine(getLog(n));
					if (last != "Everything done.")
						state = JOB_FAILED;
				}
				    
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

	string SLURMDistributor::createProgressBar(const vector<int>& progress)
	{
		int readyCount = 0;
		int waitingCount = 0;

		for (size_t n = 0; n < progress.size(); n++)
		{
			if (progress[n] == JOB_WAITING)
				waitingCount++;
			else if (progress[n] >= 100)
				readyCount++;
		}

		int runningCount = (int)progress.size() - waitingCount - readyCount;

		stringstream s;
		s << "waiting: " << waitingCount << ", ready: " << readyCount << ", running: " << runningCount << "; ";

		for (size_t n = 0; n < progress.size(); n++)
		{
			if (progress[n] != JOB_WAITING && progress[n] < 100)
				s << getProgressChar(progress[n]);
		}
		
		return s.str();
	}

	void showProgressBar(const string& bar, size_t& barLength)
	{
		// Remove old bar
		for (size_t n = 0; n < barLength; n++)
			cout << ' ';
		if(barLength > 0)
			cout << '\r';

		// Store length of new bar
		barLength = bar.length();

		// Show new bar
		cout << bar;
		if (barLength > 0)
			cout << '\r';

		cout << flush;
	}

	bool isCancelled(const string& errorMessage)
	{
		return startsWith(errorMessage, "Job ") && contains(errorMessage, "was cancelled");
	}

	bool isCancelledDueToTimeLimit(const string& errorMessage)
	{
		return startsWith(errorMessage, "Job ") && contains(errorMessage, "was cancelled due to time limit");
	}

	/**
	Extracts error message from pi2 log string.
	@param jobIndex Index of the job whose log we are parsing. The value is printed to some error messages.
	*/
	string SLURMDistributor::getErrorMessage(size_t jobIndex) const
	{
		stringstream msg;

		string log = getLog(jobIndex, true);
		string slurmErrors = getSlurmErrorLog(jobIndex, true);
		
		// First check if the job has been cancelled
		if (contains(slurmErrors, "CANCELLED"))
		{
			// The job has been cancelled by user or by system.
			if (contains(slurmErrors, "DUE TO TIME LIMIT"))
			{
				// Cancellation due to time limit
				msg << "Job " << jobIndex << " was cancelled due to time limit.";
			}
			else
			{
				// Cancellation due to unknown reason
				msg << "Job " << jobIndex << " was cancelled: " << lastLine(slurmErrors);
			}
			return msg.str();
		}

		if (log.length() <= 0)
		{
			msg << "SLURM did not run pi2. SLURM log " << makeErrorName(jobIndex) << " may contain further details. If it does not exist, please make sure that SLURM has read and write access to the current folder.";
		}
		else
		{
			// Check for internal pi2 errors in the output
			string line = lastLine(log);
			if (startsWith(line, "Error"))
			{
				msg << line;
			}
			else
			{
				// No internal error message.

				// Output info from slurm error log if it is not empty. If it is, we have no clue what the error could be, so output "no error message available".
				if (slurmErrors != "")
					msg << slurmErrors << endl;
				else
					msg << "No error message available.";
			}
		}

		return msg.str();
	}

	void SLURMDistributor::cancelJob(size_t slurmId) const
	{
		execute(scancelCommand, itl2::toString(slurmId));
	}

	void SLURMDistributor::cancelAll() const
	{
		for (size_t n = 0; n < submittedJobs.size(); n++)
		{
			cancelJob(get<0>(submittedJobs[n]));
		}
	}

	vector<string> SLURMDistributor::waitForJobs()
	{
		vector<string> result;
		if (submittedJobs.size() <= 0)
			return result;

		size_t barLength = 0;
		try
		{
			// Wait until jobs are done
			bool done = false;
			do
			{
				itl2::sleep(500);

				vector<int> progress = getJobProgress();

				// Re-submit failed jobs
				for (size_t n = 0; n < progress.size(); n++)
				{
					if (progress[n] == JOB_FAILED)
					{
						// Remove progress bar
						showProgressBar("", barLength);

						// Get error message
						string errorMessage = getErrorMessage(n);

						size_t submissionCount = get<2>(submittedJobs[n]);
						if (submissionCount < maxSubmissions)
						{
							// We can re-submit

							if (isCancelled(errorMessage))
							{
								if (isCancelledDueToTimeLimit(errorMessage))
								{
									// Try to move the job to slower queue
									JobType& type = get<1>(submittedJobs[n]);
									JobType oldType = type;
									if (type == JobType::Fast)
									{
										if (extraArgsFastJobs != extraArgsNormalJobs)
											type = JobType::Normal;
										else if (extraArgsFastJobs != extraArgsSlowJobs)
											type = JobType::Slow;
									}
									else if (type == JobType::Normal)
									{
										if (extraArgsNormalJobs != extraArgsSlowJobs)
											type = JobType::Slow;
									}

									if (type == oldType)
									{
										// The job was cancelled due to time limit and it is in the Slow queue, throw fatal error.
										// Fail and notify the user to change settings.
										throw ITLException(string("Job ") + itl2::toString(n) + " was cancelled due to time limit and no queue with higher time limits is available. Consider increasing time limits for jobs in the slurm_settings.txt file.");
									}

									cout << "Job " << n << " was cancelled due to time limit. Moving it from " << toString(oldType) << " to " << toString(type) << " queue." << endl;
								}
								else
								{
									// Job was cancelled, throw fatal error as the system is failing or somebody is cancelling the jobs.
									throw ITLException(errorMessage);
								}
							}
							else
							{
								cout << "Re-submitting failed job " << n << ". (" << errorMessage << ")" << endl;
							}

							resubmit(n);
							progress[n] = JOB_WAITING;
						}
						else
						{
							throw ITLException(string("Job ") + itl2::toString(n) + " has failed (" + errorMessage + ") Unable to re-submit as the job has been re-submitted too many times.");
						}
					}
				}

				// Show progress bar
				string bar = createProgressBar(progress);
				showProgressBar(bar, barLength);

				// Check if all jobs are done
				size_t doneCount = 0;
				for (size_t n = 0; n < progress.size(); n++)
				{
					//if (progress[n] == JOB_FAILED || progress[n] >= 100)
					if (progress[n] >= 100)
						doneCount++;
				}

				done = doneCount == progress.size();

			} while (!done);
		}
		catch (const ITLException&)
		{
			// Remove progress bar
			showProgressBar("", barLength);

			// Something went badly wrong, cancel remaining jobs.
			cout << "Cancelling remaining jobs..." << endl;
			cancelAll();

			throw;
		}

		// Remove progress bar
		showProgressBar("", barLength);

		// Read logs and generate error if something went wrong.
		vector<int> progress = getJobProgress();
		result.reserve(submittedJobs.size());
		ostringstream msg;
		for (size_t n = 0; n < submittedJobs.size(); n++)
		{
			if(progress[n] == JOB_FAILED)
				msg << "Job " << n << " failed: " << getErrorMessage(n) << endl;
		    
			string log = getLog(n, true);
			result.push_back(log);
		}
		submittedJobs.clear();
		
		string s = msg.str();
		if (s.length() > 0)
			throw ITLException(s.substr(0, s.length() - 1));

		return result;
	}


}
