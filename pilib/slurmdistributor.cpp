
#include "slurmdistributor.h"

#include "exeutils.h"
#include "utilities.h"
#include "stringutils.h"
#include "math/vectoroperations.h"
#include <algorithm>
#include "filesystem.h"

using namespace itl2;
using namespace std;

namespace pilib
{
	

	SLURMDistributor::SLURMDistributor(PISystem* system) : Distributor(system), allowedMem(0)
	{
		INIReader reader = Distributor::readConfig("slurm_config");

		extraArgsFastJobsSBatch = reader.get<string>("extra_args_fast_jobs_sbatch", "");
		extraArgsNormalJobsSBatch = reader.get<string>("extra_args_normal_jobs_sbatch", "");
		extraArgsSlowJobsSBatch = reader.get<string>("extra_args_slow_jobs_sbatch", "");
		extraArgsFastJobsSInfo = reader.get<string>("extra_args_fast_jobs_sinfo", "");
		extraArgsNormalJobsSInfo = reader.get<string>("extra_args_normal_jobs_sinfo", "");
		extraArgsSlowJobsSInfo = reader.get<string>("extra_args_slow_jobs_sinfo", "");
		jobInitCommands = reader.get<string>("job_init_commands", "");
		size_t mem = (size_t)(reader.get<double>("max_memory", 0) * 1024 * 1024);
		maxSubmissions = reader.get<size_t>("max_resubmit_count", 5) + 1;
		sbatchCommand = reader.get<string>("sbatch_command", "sbatch");
		squeueCommand = reader.get<string>("squeue_command", "squeue");
		scancelCommand = reader.get<string>("scancel_command", "scancel");
		sinfoCommand = reader.get<string>("sinfo_command", "sinfo");
		progressPollInterval = reader.get<int>("progress_poll_interval", 1000);

		if (progressPollInterval < 1)
			progressPollInterval = 1;
		
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
			string cmd = sinfoCommand + string(" --Node --Format=memory ") + extraArgsNormalJobsSInfo; // Total installed memory. Not all of that may be available!
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

		cout << "Memory per node in the SLURM cluster: " << bytesToString((double)allowedMem) << endl;
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

	string SLURMDistributor::makeSbatchName(size_t jobIndex) const
	{
		return "./slurm-io-files/" + makeJobName(jobIndex) + "-sbatch.sh";
	}

	void SLURMDistributor::resubmit(size_t jobIndex)
	{
		string jobName = makeJobName(jobIndex);
		string inputName = makeInputName(jobIndex);
		string outputName = makeOutputName(jobIndex);
		string errorName = makeErrorName(jobIndex);
		string sbatchName = makeSbatchName(jobIndex);
		JobType jobType = get<1>(submittedJobs[jobIndex]);

		fs::remove(outputName);
		fs::remove(errorName);

		//string jobCmdLine;
		//if (jobInitCommands.length() > 0)
		//	jobCmdLine = jobInitCommands + "; ";
		//jobCmdLine += "'" + getPiCommand() + "' " + inputName;
		//string sbatchArgs = string("--no-requeue") + " --job-name=" + jobName + " --output=" + outputName + " --error=" + errorName + " " + extraArgsSBatch(jobType) + " --wrap=\"" + jobCmdLine + "\"";

		string sbatchCode;
		sbatchCode += "#!/bin/bash\n";
		sbatchCode += "#SBATCH --no-requeue\n";
		sbatchCode += "#SBATCH --job-name=" + jobName + "\n";
		sbatchCode += "#SBATCH --output=" + outputName + "\n";
		sbatchCode += "#SBATCH --error=" + errorName + "\n";
		
		string sbatchExtra = extraArgsSBatch(jobType);
		if(!isWhitespace(sbatchExtra))
			sbatchCode += "#SBATCH " + sbatchExtra + "\n";
		sbatchCode += "\n";
		
		if(!isWhitespace(jobInitCommands))
			sbatchCode += jobInitCommands + "\n";
		
		sbatchCode += getJobPiCommand() + " \"" + inputName + "\"\n";
		sbatchCode += "\n";
		writeText(sbatchName, sbatchCode);

		string sbatchArgs = sbatchName;
		string result = execute(sbatchCommand, sbatchArgs);

		vector<string> lines = split(result);
		if (lines.size() == 1)
		{
			vector<string> parts = split(*lines.rbegin(), false, ' ');
			
			size_t slurmId;
			try
			{
			    if (parts.size() < 1)
				    throw ITLException("SLURM returned no batch job id.");
			    
			    slurmId = fromString<int>(parts[parts.size() - 1]);
			}
			catch(ITLException e)
			{
			    cout << "Command" << endl;
				cout << sbatchCommand << " " << sbatchArgs << endl;
				cout << "returned" << endl;
				cout << result << endl;
				throw ITLException(string("Command ") + sbatchCommand + " did not return a job id. The received output has been printed to standard output.");
			}
			
			get<0>(submittedJobs[jobIndex]) = slurmId;
			get<2>(submittedJobs[jobIndex])++;

			cout << "Submitted job " << jobName << ", SLURM id = " << slurmId << endl;
		}
		else
		{
			if (result.length() > 0)
			{
				cout << "Command" << endl;
				cout << sbatchCommand << " " << sbatchArgs << endl;
				cout << "returned" << endl;
				cout << result << endl;
				throw ITLException(string("Unexpected ") + sbatchCommand + " output. The output has been printed to standard output.");
			}
			else
			{
				throw ITLException(string("Empty ") + sbatchCommand + " output.");
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

		// TODO: Sometimes empty result means that slurm is somehow intermittently unavailable, but jobs are
		// still running. How to detect that? Perhaps wait a second and try again?
		
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
			cout << "Warning: Unexpected " << squeueCommand << " output '" << result << "' for job " << jobIndex << " (SLURM id " << slurmId << "). Assuming the job is still running." << endl;
			return false;
		}

		return !running;
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

	double parseSlurmTime(string str)
	{
		try
		{
			// Parse from format [days-]hours:minutes:seconds
			int days = 0;
			size_t dashPos = str.find('-');
			if (dashPos != string::npos)
			{
				string daysStr = str.substr(0, dashPos);
				str = str.substr(dashPos + 1);

				days = fromString<int32_t>(daysStr);
			}

			size_t colonPos = str.find(':');
			if (colonPos == string::npos)
				return -1;
			string hoursStr = str.substr(0, colonPos);
			str = str.substr(colonPos + 1);

			colonPos = str.find(':');
			if (colonPos == string::npos)
				return -1;
			string minsStr = str.substr(0, colonPos);
			string secsStr = str.substr(colonPos + 1);

			int hours = fromString<int32_t>(hoursStr);
			int mins = fromString<int32_t>(minsStr);
			int secs = fromString<int32_t>(secsStr);

			return days * 24 * 60 * 60 + hours * 60 * 60 + mins * 60 + secs;
		}
		catch (ITLException)
		{
			// Invalid conversion. The input is invalid. Return -1.
			return -1;
		}
	}

	double parseRawTime(const string& str)
	{
		try
		{
			return fromString<int32_t>(str);
		}
		catch (ITLException)
		{
			// Invalid conversion. The input is invalid. Return -1.
			return -1;
		}
	}

	tuple<double, double> SLURMDistributor::getJobTimes(size_t jobIndex) const
	{
		string slurmId = itl2::toString(get<0>(submittedJobs[jobIndex]));

		string result = execute("sacct", string("-X --noheader -o reserved,elapsedraw -j ") + slurmId);
		trim(result);

		vector<string> parts = split(result, false, ' ');
		
		if (parts.size() != 2)
			return make_tuple(-1.0, -1.0); // Output is invalid, unable to parse it.

		return make_tuple(parseSlurmTime(parts[0]), parseRawTime(parts[1]));
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
				itl2::sleep(progressPollInterval);

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
										if (extraArgsFastJobsSBatch != extraArgsNormalJobsSBatch)
											type = JobType::Normal;
										else if (extraArgsFastJobsSBatch != extraArgsSlowJobsSBatch)
											type = JobType::Slow;
									}
									else if (type == JobType::Normal)
									{
										if (extraArgsNormalJobsSBatch != extraArgsSlowJobsSBatch)
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
					if (progress[n] >= 100)
					{
						doneCount++;

						if (jobTiming.size() == progress.size())
						{
							if (get<0>(jobTiming[n]) < -1 || get<1>(jobTiming[n]) < -1)
							{
								jobTiming[n] = getJobTimes(n);
							}
						}
					}
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
