#pragma once

#include <string>

using namespace std;

namespace itl2
{


	/**
	 * Run a given test function and report result.
	 */
	void test(void (*testfunc)(), const string& testName);

	/**
	 * If condition is false, report it to user.
	 */
	void testAssert(bool condition, const string& assertName);

	/**
	 * Shows count of failed and passed tests.
	 */
	void testReport();
}
