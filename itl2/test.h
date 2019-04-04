#pragma once

#include <string>

namespace itl2
{


	/**
	Run a given test function and report result.
	*/
	void test(void (*testfunc)(), const std::string& testName);

	/**
	If condition is false, report it to user.
	@return condition
	 */
	bool testAssert(bool condition, const std::string& assertName);

	/**
	Shows count of failed and passed tests.
	*/
	void testReport();

}
