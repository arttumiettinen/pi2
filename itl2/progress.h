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
		bool show;
		bool isTerm;

	public:
		//ProgressIndicator(float maxSteps, bool show = true) :
		template<typename T> ProgressIndicator(T maxSteps, bool show = true) : 
			currSteps(-(float)maxSteps),
			maxSteps((float)maxSteps),
			show(show)
		{
			isTerm = isTerminal();
			// Always show the indicator
			Show(0);
		}


		virtual ~ProgressIndicator()
		{
			if (show)
			{
				// Move to next line so that the progress bar does not conflict with future printing.
				std::cout << std::endl;
			}
		}

		void Step()
		{
			Show(currSteps + 1);
		}

		//void Show(float progress)
		template<typename T> void Show(T progress)
		{
			if (show)
			{
				if (isTerm)
				{
					coord_t prevProgress = round(currSteps / (maxSteps - 1) * 100);
					coord_t currProgress = round((float)progress / (maxSteps - 1) * 100);
					if (currProgress != prevProgress)
					{
						std::cout << currProgress << " %              \r" << std::flush;
					}
				}
				else
				{
					coord_t prevProgress = round(currSteps / (maxSteps - 1) * 10);
					coord_t currProgress = round((float)progress / (maxSteps - 1) * 10);
					for (coord_t n = prevProgress; n < currProgress; n++)
					{
						std::cout << "=";
					}
					std::cout << std::flush;
				}

				currSteps = (float)progress;
			}
		}

	};

	namespace tests
	{
		void progress();
	}
}
