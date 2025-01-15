#pragma once

#if defined(__linux__) || defined(__APPLE__)
#include <sys/time.h>
#elif defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#else
#error timer.h not configured for this platform.
#endif

namespace itl2
{

	// TODO: This platform-dependent stuff could probaly be replaced by std::chrono (as such a thing is well available nowadays).

#if defined(__linux__) || defined(__APPLE__)

	/**
	 * Stopwatch class.
	 */
	class Timer
	{
	private:
		timeval startTime;
		timeval endTime;

		bool isRunning;

	public:
		Timer()
		{
			isRunning = false;
			gettimeofday(&startTime, 0);
			gettimeofday(&endTime, 0);
		}

		void start()
		{
			gettimeofday(&startTime, 0);
			isRunning = true;
		}

		void stop()
		{
			gettimeofday(&endTime, 0);
			isRunning = false;
		}

		/**
		Re-starts the timer and returns the elapsed time in seconds.
		*/
		double lap()
		{
			stop();
			double result = getSeconds();
			start();
			return result;
		}

		/**
		 * Get time between start and this call or start and stop if stop has been called.
		 * Returns time in milliseconds.
		 */
		double getTime()
		{
			double sec, usec;
			if(isRunning)
			{
				timeval c;
				gettimeofday(&c, 0);
				sec = c.tv_sec;
				usec = c.tv_usec;
			}
			else
			{
				sec = endTime.tv_sec;
				usec = endTime.tv_usec;
			}
			return (sec * 1000 + usec / 1000) - (startTime.tv_sec * 1000 + startTime.tv_usec / 1000);
		}

		/**
		 * Get elapsed time in seconds.
		 */
		double getSeconds()
		{
			return getTime() / 1000.0;
		}
	};

#elif defined(_WIN32)

	/**
	 * Stopwatch class.
	 */
	class Timer
	{
	private:
		/**
		Frequency
		*/
		LARGE_INTEGER m_frequency;

		/**
		Start time.
		*/
		LARGE_INTEGER m_startTime;

		/**
		Stop time.
		*/
		LARGE_INTEGER m_endTime;

		bool isRunning;

	public:
		Timer()
		{
			isRunning = false;
			QueryPerformanceFrequency(&m_frequency);
			QueryPerformanceCounter(&m_startTime);
			QueryPerformanceCounter(&m_endTime);
		}

		void start()
		{
			QueryPerformanceCounter(&m_startTime);
			isRunning = true;
		}

		void stop()
		{
			QueryPerformanceCounter(&m_endTime);
			isRunning = false;
		}

		/**
		Re-starts the timer and returns the elapsed time in seconds.
		*/
		double lap()
		{
			stop();
			double result = getSeconds();
			start();
			return result;
		}

		/**
		 * Get time between start and this call or start and stop if stop has been called.
		 * Returns time in milliseconds.
		 */
		double getTime()
		{
			double current;
			if(isRunning)
			{
				LARGE_INTEGER time;
				QueryPerformanceCounter(&time);
				current = (double)time.QuadPart;
			}
			else
			{
				current = (double)m_endTime.QuadPart;
			}

			return (current - (double)m_startTime.QuadPart) / (double)m_frequency.QuadPart * 1000;
		}

		/**
		 * Get elapsed time in seconds.
		 */
		double getSeconds()
		{
			return getTime() / 1000.0;
		}
	};

#endif
}
