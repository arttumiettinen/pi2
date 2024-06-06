
#include "progress.h"

namespace itl2
{
	int ProgressIndicator::nestingCount = 0;

	ProgressIndicator::ProgressIndicator(size_t maxSteps, bool showIndicator) :
		counter(0),
		maxSteps(maxSteps),
		showIndicator(showIndicator),
		messageLength(0)
	{
		nestingCount++;
		isTerm = isTerminal();

		if (nestingCount > 1)
			this->showIndicator = false;

		if (this->showIndicator && isTerm)
			std::cout << "0 %\r" << std::flush;
	}

	ProgressIndicator::ProgressIndicator(const string& initialMessage, bool showIndicator) :
		counter(0),
		maxSteps(0),
		showIndicator(showIndicator),
		messageLength(0)
	{
		nestingCount++;
		isTerm = isTerminal();

		if (nestingCount > 1)
			this->showIndicator = false;

		showMessage(initialMessage);
	}


	ProgressIndicator::~ProgressIndicator()
	{
		nestingCount--;
		if (nestingCount < 0)
			nestingCount = 0; // TODO: Deleting more objects than created; indicate error in code.

		if (showIndicator)
		{
			if (isTerm)
			{
				// Remove the indicator
				eraseMessage();
			}
			else
			{
				// Move to next line so that the progress bar does not conflict with future printing.
				std::cout << std::endl;
			}
		}
	}

	namespace tests
	{
		void progress()
		{
			std::cout << "Message test..." << std::endl;
			{
				ProgressIndicator ind("initial message");
				sleep(1000);
				ind.showMessage("Second message");
				sleep(1000);
				ind.showMessage("Third message, last until erase");
				sleep(1000);
			}

			std::cout << "Nesting test..." << std::endl;
			{
				ProgressIndicator prog(50000);
				for (size_t n = 0; n < 50000; n++)
				{
					// Nesting
					ProgressIndicator prog2(10000);
					for (size_t m = 0; m < 10000; m++)
					{
						prog2.step();
					}

					prog.step();
				}
			}

			std::cout << "Long progress indicator test..." << std::endl;
			{
				ProgressIndicator prog(450 * 450 * 450);
				for (size_t n = 0; n < 450 * 450 * 450; n++)
					prog.step();
			}

			std::cout << "Interrupted progress indicator test..." << std::endl;
			size_t MAX = 2235450000;
			{
				ProgressIndicator prog(MAX);
				for (float n = 0; n < 3 * MAX / 4; n += 173)
					prog.step();
			}
			std::cout << "this message interrupts the progress indicator..." << std::endl;
		}
	}
}