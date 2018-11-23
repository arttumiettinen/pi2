#pragma once

#include <string>

using std::string;

namespace itl2
{
	/**
	* Exception for itl system.
	*/
	class ITLException
	{
	private:
		string m_message;

	public:
		/**
		* Constructor
		*/
		ITLException(const string& msg) :
			m_message(msg)
		{
		}

		/**
		* Gets the error message.
		*/
		const string& message() const
		{
			return m_message;
		}
	};
}
