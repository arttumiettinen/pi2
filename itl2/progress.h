#pragma once

#include <iostream>

#include "utilities.h"

namespace itl2
{
	/**
	Progress indicator that clears itself when it goes out of scope.
	*/
	class ProgressIndicator
	{
	private:
		float counter;

		float maxSteps;
		bool showIndicator;
		bool isTerm;

	public:
		template<typename T> ProgressIndicator(T maxSteps, bool showIndicator = true) : 
			counter(0),
			maxSteps((float)maxSteps),
			showIndicator(showIndicator)
		{
			isTerm = isTerminal();
			// Always show the indicator
			if(showIndicator && isTerm)
				std::cout << "0 %\r" << std::flush;
		}


		virtual ~ProgressIndicator()
		{
			if (showIndicator)
			{
				if (isTerm)
				{
					// Remove the indicator
					std::cout << "                   \r" << std::flush;
				}
				else
				{
					// Move to next line so that the progress bar does not conflict with future printing.
					std::cout << std::endl;
				}
			}
		}

		/**
		Step progress indicator to the next position
		*/
		void step()
		{
			if (showIndicator)
			{

				float localCounter;

#if defined(_MSC_VER)
#pragma omp atomic
				counter++;

				localCounter = counter;
#else
#pragma omp atomic capture
				localCounter = ++counter;
#endif

				if (localCounter > maxSteps)
					localCounter = maxSteps;

				if (isTerm)
				{
					coord_t prevProgress = round((localCounter - 1) / maxSteps * 100);
					coord_t currProgress = round(localCounter / maxSteps * 100);
					if (currProgress != prevProgress)
					{
#pragma omp critical(showProgress)
						std::cout << currProgress << " %              \r" << std::flush;
					}
				}
				else
				{
					coord_t prevProgress = round((localCounter - 1) / (maxSteps) * 10);
					coord_t currProgress = round(localCounter / (maxSteps) * 10);
#pragma omp critical(showProgress)
					{
						for (coord_t n = prevProgress; n < currProgress; n++)
						{
							std::cout << "=";
						}
						std::cout << std::flush;
					}
				}
			}
		}

	};

	namespace tests
	{
		void progress();
	}
}
