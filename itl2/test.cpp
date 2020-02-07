
#include "test.h"

#include "itlexception.h"

#include <iostream>

using namespace std;

namespace itl2
{
	size_t okTestCount = 0;
	size_t failedTestCount = 0;

	bool testOk = true;

	bool throwOnAssertion = false;

	void throwOnFailedAssertion(bool shouldThrow)
	{
		throwOnAssertion = shouldThrow;
	}

	void test(void (*testfunc)(), const string& testName)
	{
		cout << "Testing " << testName << endl;
		testOk = true;
		try
		{
			(*testfunc)();
		}
		catch (const ITLException& e)
		{
			testOk = false;
			cout << "Unexpected exception: " << e.message() << endl;
		}
		catch (const exception& e)
		{
			testOk = false;
			cout << "Unexpected exception: " << e.what() << endl;
		}
		catch(...)
		{
			testOk = false;
			cout << "Unexpected exception." << endl;
		}

		if(testOk)
		{
			cout << "Test OK" << endl;
			okTestCount++;
		}
		else
		{
			cout << "Test FAILED" << endl;
			failedTestCount++;
		}

		throwOnFailedAssertion(false);
	}

	bool testAssert(bool condition, const string& assertName)
	{
		if(!condition)
		{
			testOk = false;
			cout << "ASSERTION FAILED: " << assertName << endl;

			if (throwOnAssertion)
				throw ITLException("Exception requested on assertion failure.");
		}

		return condition;
	}

	void testReport()
	{
		cout << okTestCount << " tests passed." << endl;
		cout << failedTestCount << " tests failed." << endl;
	}
}
