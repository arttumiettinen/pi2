
#include "lsfdistributor.h"

#include "exeutils.h"
#include "utilities.h"
#include "stringutils.h"
#include "math/vectoroperations.h"
#include <algorithm>
#include "filesystem.h"

#include <random>

using namespace itl2;
using namespace std;

namespace pilib
{


	LSFDistributor::LSFDistributor(PISystem* system) : Distributor(system), allowedMem(0)
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
			configPath = configPath.replace_filename("lsf_config.txt");

			//cout << "Reading settings from " << configPath << endl;

			INIReader reader(configPath.string());

			extraArgsFastJobs = reader.get<string>("extra_args_fast_jobs", "");
			extraArgsNormalJobs = reader.get<string>("extra_args_normal_jobs", "");
			extraArgsSlowJobs = reader.get<string>("extra_args_slow_jobs", "");
			jobInitCommands = reader.get<string>("job_init_commands", "");
			mem = (size_t)(reader.get<double>("max_memory", 0) * 1024 * 1024);
			maxSubmissions = reader.get<size_t>("max_resubmit_count", 5) + 1;
			bsubCommand = reader.get<string>("bsub_command", "bsub");
			bjobsCommand = reader.get<string>("bjobs_command", "bjobs");
			bkillCommand = reader.get<string>("bkill_command", "bkill");

			readSettings(reader);
		}

		allowedMemory(mem);
	}

	void LSFDistributor::allowedMemory(size_t maxMem)
	{
		allowedMem = maxMem;

		if (allowedMem <= 0)
		{

		}

		// TODO: Convert this from Slurm to LSF
		//if (allowedMem <= 0)
		//{
		//	// Run cluster info command to get amount of memory in megabytes
		//	//string cmd = string("sinfo --Node --Format=freemem ") + extraArgs; // Memory available to start new programs. Jobs running on nodes may affect this.
		//	string cmd = sinfoCommand + string(" --Node --Format=memory ") + extraArgsNormalJobsSInfo; // Total installed memory. Not all of that may be available!
		//	string output = execute(cmd);

		//	// Split lines
		//	trim(output);
		//	vector<string> lines = split(output, false);

		//	// Remove header
		//	lines.erase(lines.begin() + 0);

		//	// Convert list to list of ints
		//	vector<size_t> values = fromString<size_t>(lines);

		//	// Calculate minimum in bytes
		//	if (values.size() > 0)
		//		allowedMem = (size_t)std::round(min(values) * 1024 * 1024 * 0.75);

		//	if (allowedMem <= 0)
		//		throw ITLException("Unable to determine maximum memory available in the SLURM cluster. Are the extra_args_* parameters in the slurm_config.txt file correct?");
		//	//The file is read from ") + configPath.string());
		//}

		cout << "Memory per node in the LSF cluster: " << bytesToString((double)allowedMem) << endl;
	}

	string LSFDistributor::makeJobName(size_t jobIndex) const
	{
		return "pi2-" + itl2::toString<size_t>(jobIndex) + "-" + myName;
	}

	string LSFDistributor::makeInputName(size_t jobIndex) const
	{
		return "./lsf-io-files/" + makeJobName(jobIndex) + "-in.txt";
	}

	string LSFDistributor::makeOutputName(size_t jobIndex) const
	{
		return "./lsf-io-files/" + makeJobName(jobIndex) + "-out.txt";
	}

	string LSFDistributor::makeErrorName(size_t jobIndex) const
	{
		return "./lsf-io-files/" + makeJobName(jobIndex) + "-err.txt";
	}

	void LSFDistributor::resubmit(size_t jobIndex)
	{
		string jobName = makeJobName(jobIndex);
		string inputName = makeInputName(jobIndex);
		string outputName = makeOutputName(jobIndex);
		string errorName = makeErrorName(jobIndex);
		JobType jobType = get<1>(submittedJobs[jobIndex]);

		fs::remove(outputName);
		fs::remove(errorName);

		string initStr = " ";
		if (jobInitCommands.length() > 0)
			initStr = string(" -E ") + jobInitCommands + " ";
		string jobCmdLine = "'" + getPiCommand() + "' " + inputName;

		string bsubArgs = string("") + "-J " + jobName + " -o " + outputName + " -e " + errorName + initStr + extraArgs(jobType) + " " + jobCmdLine;

	cout << "bsub arguments: " << bsubArgs << endl;

		string result = execute(bsubCommand, bsubArgs);

	cout << "bsub output: " << result << endl;

		// Here I assume that the output looks like this:
		// Job <930> is submitted to default queue <normal>.

		vector<string> lines = split(result);
		try
		{
			if (lines.size() == 1)
			{
				string line = lines[0];
				if (startsWith(line, "Job <"))
				{
					size_t idStart = string("Job <").length() + 1;
					size_t idEnd = line.find('>');
					if (idEnd != string::npos)
					{
						string idStr = line.substr(idStart, idEnd - idStart + 1);
						size_t id = fromString<int>(idStr);
						get<0>(submittedJobs[jobIndex]) = id;
						get<2>(submittedJobs[jobIndex])++;
					}
					else
					{
						throw ITLException(bsubCommand + " output does contain '>'.");
					}
				}
				else
				{
					throw ITLException(bsubCommand + " output does not start with 'Job <'.");
				}
			}
			else
			{
				throw ITLException(string("Unxpected ") + bsubCommand + " output.");
			}
		}
		catch (ITLException)
		{
			cout << "Command" << endl;
			cout << bsubCommand << " " << bsubArgs << endl;
			cout << "returned" << endl;
			cout << result << endl;
			throw;
		}
	}

	void LSFDistributor::submitJob(const string& piCode, JobType jobType)
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

	int LSFDistributor::getJobStatus(size_t jobIndex) const
	{
		string id = itl2::toString(get<0>(submittedJobs[jobIndex]));

		string bjobsArgs = "-X -noheader -o \"STAT\" " + id;

	cout << "bjobs arguments: " << bjobsArgs << endl;

		string result = execute(bjobsCommand, bjobsArgs);

	cout << "bjobs output: " << result << endl;

		trim(result);

		bool waiting = startsWith(result, "PEND") ||
			startsWith(result, "PROV") ||
			startsWith(result, "PSUSP") ||
			startsWith(result, "USUSP") ||
			startsWith(result, "SSUSP") ||
			startsWith(result, "WAIT");

		bool running = startsWith(result, "RUN");

		bool finished = startsWith(result, "DONE");

		bool error = startsWith(result, "EXIT") ||
			startsWith(result, "UNKWN") ||
			startsWith(result, "ZOMBI");

		if (waiting)
			return JOB_WAITING;
		if (running)
			return 0;
		if (finished)
			return 100;
		
		return JOB_FAILED;
	}

	

	string LSFDistributor::getLog(size_t jobIndex, bool flush) const
	{
		string outputName = makeOutputName(jobIndex);
		flushCache(outputName);
		return readText(outputName);
	}

	string LSFDistributor::getErrorLog(size_t jobIndex, bool flush) const
	{
		string outputName = makeErrorName(jobIndex);
		flushCache(outputName);
		return readText(outputName);
	}


	int LSFDistributor::getJobProgressFromLog(size_t jobIndex) const
	{
		// Check if the job has been cancelled?
		string errorLog = getErrorLog(jobIndex);

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

	vector<int> LSFDistributor::getJobProgress() const
	{
		vector<int> progress;
		progress.reserve(submittedJobs.size());
		for (size_t n = 0; n < submittedJobs.size(); n++)
		{
			int state = getJobStatus(n);

			if (state == 0)
			{
				// Job is running, read progress from the log file.
				state = getJobProgressFromLog(n);
				// Do not let the progress read from log ever reach 100 as that is interpreted
				// as job done.
				if (state >= 100)
					state = 99;
			}
			else if (state == 100)
			{
				// Job is done

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

				if (state != JOB_FAILED)
					state = 100;
			}

			progress.push_back(state);
		}
		return progress;
	}


	bool isCancelledLSF(const string& errorMessage)
	{
		// TODO
		return startsWith(errorMessage, "Job ") && contains(errorMessage, "was cancelled");
	}

	bool isCancelledDueToTimeLimitLSF(const string& errorMessage)
	{
		// TODO
		return startsWith(errorMessage, "Job ") && contains(errorMessage, "was cancelled due to time limit");
	}

	/**
	Extracts error message from pi2 log string.
	@param jobIndex Index of the job whose log we are parsing. The value is printed to some error messages.
	*/
	string LSFDistributor::getErrorMessage(size_t jobIndex) const
	{
		stringstream msg;

		// TODO: How LSF indicates that a job has been cancelled due to time limit?

		string log = getLog(jobIndex, true);
		string slurmErrors = getErrorLog(jobIndex, true);

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
			msg << "LSF did not run pi2. Log " << makeErrorName(jobIndex) << " may contain further details. If it does not exist, please make sure that LSF has read and write access to the current folder.";
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

	void LSFDistributor::cancelJob(size_t slurmId) const
	{
		execute(bkillCommand, itl2::toString(slurmId));
	}

	void LSFDistributor::cancelAll() const
	{
		for (size_t n = 0; n < submittedJobs.size(); n++)
		{
			cancelJob(get<0>(submittedJobs[n]));
		}
	}

	vector<string> LSFDistributor::waitForJobs()
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

							if (isCancelledLSF(errorMessage))
							{
								if (isCancelledDueToTimeLimitLSF(errorMessage))
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
			if (progress[n] == JOB_FAILED)
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
