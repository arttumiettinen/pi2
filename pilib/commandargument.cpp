
#include "commandargument.h"
#include "pilibutilities.h"
#include "stringutils.h"

using namespace std;

namespace pilib
{

	string CommandArgumentBase::defaultValueQuoted() const
	{
		if (!defaultAllowed())
			return "";

		if (defaultValue() == "" || itl2::contains(defaultValue(), " "))
			return "\"" + defaultValue() + "\"";

		return defaultValue();
	}

	string CommandArgumentBase::toString() const
	{
		stringstream msg;
		msg << "[" << pilib::toString(direction()) << ", " << pilib::toString(dataType()) << "] " << name();
		if (defaultAllowed())
		{
			msg << " = " << defaultValueQuoted();
		}
		return msg.str();
	}

	string CommandArgumentBase::toSimpleString() const
	{
		stringstream msg;
		msg << name();
		return msg.str();
	}

	string CommandArgumentBase::dataTypeHelpStringNoDefault(const string* dataTypeOverride) const
	{
		if (!dataTypeOverride)
		{
			return pilib::toString(dataType());
		}
		else if (dataTypeOverride->length() > 0)
		{
			return *dataTypeOverride;
		}
		return "";
	}

	string CommandArgumentBase::dataTypeHelpString(const string* dataTypeOverride) const
	{
		stringstream msg;

		msg << dataTypeHelpStringNoDefault(dataTypeOverride);

		if (direction() == ParameterDirection::In && defAllowed && !dataTypeOverride)
			msg << "(" << defaultValueQuoted() << ")";
		

		//bool printDefault = false;
		//if (!dataTypeOverride)
		//{
		//	msg << pilib::toString(dataType());
		//	printDefault = true;
		//}
		//else if (dataTypeOverride->length() > 0)
		//{
		//	msg << *dataTypeOverride;
		//}

		//if (direction() == ParameterDirection::In && defAllowed && printDefault)
		//	msg << "(\"" << defaultValue() << "\")";

		return msg.str();
	}

	string CommandArgumentBase::helpString(HelpFormat format, const string* dataTypeOverride) const
	{
		stringstream msg;

		string dtstr = dataTypeHelpString(dataTypeOverride);
		string dtstrNoDefault = dataTypeHelpStringNoDefault(dataTypeOverride);
		string title;

		switch (format)
		{
		case pilib::HelpFormat::Text:
			msg << name() << " [";
			if (direction() == ParameterDirection::In)
				msg << "input";
			else if (direction() == ParameterDirection::Out)
				msg << "output";
			else
				msg << "input & output";

			if (dtstr != "")
				msg << "; " << dtstr;

			msg << "]";
			msg << ": ";
			msg << endl;
			msg << help << endl;
			msg << endl;
			break;

		case pilib::HelpFormat::Rst:
			title = name() + " [";
			
			if (direction() == ParameterDirection::In)
				title += "input]";
			else if (direction() == ParameterDirection::Out)
				title += "output]";
			else
				title += "input & output]";
			
			printTitle(msg, title, 3);

			if (dtstrNoDefault != "")
				msg << "**Data type:** " << dtstrNoDefault << endl << endl;

			if (direction() == ParameterDirection::In && defAllowed)
			{
				if (!dataTypeOverride || *dataTypeOverride == pilib::toString(dataType()))
					msg << "**Default value:** " << defaultValueQuoted() << endl << endl;
				else
					msg << "**Default value:** Shown along data types." << endl << endl;
			}
			

			msg << help << endl;
			msg << endl;
			break;
		default:
			throw logic_error("Unsupported help format.");
		}
		
		
		return msg.str();
	}
}