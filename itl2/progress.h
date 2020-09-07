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
		float currSteps;
		float maxSteps;
		bool showIndicator;
		bool isTerm;

	public:
		template<typename T> ProgressIndicator(T maxSteps, bool showIndicator = true) : 
			currSteps(-(float)maxSteps),
			maxSteps((float)maxSteps),
			showIndicator(showIndicator)
		{
			isTerm = isTerminal();
			// Always show the indicator
			show(0);
		}


		virtual ~ProgressIndicator()
		{
			if (showIndicator)
			{
				// Move to next line so that the progress bar does not conflict with future printing.
				std::cout << std::endl;
			}
		}

		/**
		Step progress indicator to the next position
		*/
		void step()
		{
			show(currSteps + 1);
		}

		/**
		Show progress indicator in given position.
		*/
		template<typename T> void show(T progress)
		{
			if (showIndicator)
			{
				if (isTerm)
				{
					coord_t prevProgress = round(currSteps / (maxSteps) * 100);
					coord_t currProgress = round((float)progress / (maxSteps) * 100);
					if (currProgress != prevProgress)
					{
						std::cout << currProgress << " %              \r" << std::flush;
					}
				}
				else
				{
					coord_t prevProgress = round(currSteps / (maxSteps) * 10);
					coord_t currProgress = round((float)progress / (maxSteps) * 10);
					for (coord_t n = prevProgress; n < currProgress; n++)
					{
						std::cout << "=";
					}
					std::cout << std::flush;
				}

				currSteps = (float)progress;
				if (currSteps > maxSteps)
					currSteps = maxSteps;
			}
		}

	};

	namespace tests
	{
		void progress();
	}
}
