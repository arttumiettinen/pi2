#pragma once

#include <string>

namespace itl2
{
	/**
	* Exception for itl system.
	*/
	class ITLException
	{
	private:
		std::string m_message;

	public:
		/**
		* Constructor
		*/
		ITLException(const std::string& msg) :
			m_message(msg)
		{
		}

		/**
		* Gets the error message.
		*/
		const std::string& message() const
		{
			return m_message;
		}
	};
}
