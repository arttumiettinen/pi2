#pragma once

#include <iostream>

#include "utilities.h"

namespace itl2
{
	/**
	Progress indicator that clears itself when it goes out of scope.
	Progress indicators can be nested and only the top-most one will be shown.
	*/
	class ProgressIndicator
	{
	private:
		size_t counter;

		size_t maxSteps;
		bool showIndicator;
		bool isTerm;
		size_t messageLength;

		static int nestingCount;

		void eraseMessage()
		{
			for (size_t n = 0; n < messageLength; n++)
			{
				std::cout << " ";
			}
			std::cout << "\r" << std::flush;
		}

	public:
		
		/**
		Constructor
		Creates a progress indicator that shows a bar/percentage that has maxSteps steps.
		Proceed to next step by calling step() function.
		@param maxSteps The total number of times step() function will be called.
		@param showIndicator Set to false to not show the indicator at all.
		*/
		ProgressIndicator(size_t maxSteps, bool showIndicator = true);

		/**
		Constructor
		Creates a progress indicator that shows a message instead of counting steps until completion.
		Show a message by calling showMessage function.
		*/
		ProgressIndicator(const string& initialMessage, bool showIndicator = true);


		virtual ~ProgressIndicator();

		/**
		Show a progress message.
		Progress message and bar cannot be used interchangeably.
		*/
		void showMessage(const string& message)
		{
			if (showIndicator)
			{
				// TODO: Don't update if the message does not change.
#pragma omp critical(showProgress)
				{
					eraseMessage();
					messageLength = message.length();
					std::cout << message << "\r" << std::flush;
				}
			}
		}

		/**
		Step progress indicator to the next position.
		Progress message and bar cannot be used interchangeably.
		This function is thread-safe.
		*/
		void step()
		{
			if (showIndicator)
			{

				size_t localCounter;

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
					coord_t prevProgress = round((float)(localCounter - 1) / maxSteps * 100);
					coord_t currProgress = round((float)localCounter / maxSteps * 100);
					if (currProgress != prevProgress)
					{
#pragma omp critical(showProgress)
						std::cout << currProgress << " %              \r" << std::flush;
					}
				}
				else
				{
					coord_t prevProgress = round((float)(localCounter - 1) / (maxSteps) * 10);
					coord_t currProgress = round((float)localCounter / (maxSteps) * 10);
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
