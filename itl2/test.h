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
	Sets or clears a flag that indicates whether failed assertion should result in canceled test.
	The flag is automatically cleared after a test is finished.
	*/
	void throwOnFailedAssertion(bool shouldThrow);

	/**
	Shows count of failed and passed tests.
	*/
	void testReport();

}
