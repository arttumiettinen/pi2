
#include "pilib.h"

#include "pisystem.h"

using namespace pilib;

void* createPI()
{
	return new PISystem();
}

void destroyPI(void* pi)
{
	delete (PISystem*)pi;
}

uint8_t run(void* pi, const char* commands)
{
	return ((PISystem*)pi)->run(commands) ? 1 : 0;
}

const char* lastErrorMessage(void* pi)
{
	return ((PISystem*)pi)->getLastErrorMessage();
}

int32_t lastErrorLine(void* pi)
{
	return (int32_t)((PISystem*)pi)->getLastErrorLine();
}

void clearLastError(void* pi)
{
	((PISystem*)pi)->clearLastError();
}

const char* commandList(void* pi)
{
	static string cmdList = ((PISystem*)pi)->commandList(false);
	return cmdList.c_str();
}

const char* help(void* pi, const char* commandName)
{
	vector<string> v = ((PISystem*)pi)->getHelp(commandName);
	static string hlp;
	hlp = "";
	for (size_t n = 0; n < v.size(); n++)
	{
		hlp += v[n] + "\n\n";
	}
	return hlp.c_str();
}

void* getImage(void* pi, const char* imgName, int64_t* width, int64_t* height, int64_t* depth, int32_t* dataType)
{
	PISystem* sys = (PISystem*)pi;
	ImageBase* img = sys->getImageNoThrow(imgName);

	if (!img)
	{
		*width = 0;
		*height = 0;
		*depth = 0;
		*dataType = Unknown;
		return 0;
	}

	*width = img->width();
	*height = img->height();
	*depth = img->depth();
	*dataType = (int)img->dataType();

	return img->getRawData();
}
