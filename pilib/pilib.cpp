
#include <mutex>

#include "pilib.h"

#include "pisystem.h"

using namespace pilib;

std::mutex mutex;

void* createPI()
{
	return new PISystem();
}

void destroyPI(void* pi)
{
	std::lock_guard<std::mutex> lock(mutex);
	delete (PISystem*)pi;
}

uint8_t run(void* pi, const char* commands)
{
	std::lock_guard<std::mutex> lock(mutex);
	return ((PISystem*)pi)->run(commands) ? 1 : 0;
}

const char* lastErrorMessage(void* pi)
{
	std::lock_guard<std::mutex> lock(mutex);
	return ((PISystem*)pi)->getLastErrorMessage();
}

int32_t lastErrorLine(void* pi)
{
	std::lock_guard<std::mutex> lock(mutex);
	return (int32_t)((PISystem*)pi)->getLastErrorLine();
}

void clearLastError(void* pi)
{
	std::lock_guard<std::mutex> lock(mutex);
	((PISystem*)pi)->clearLastError();
}

void getImageInfo(void* pi, const char* imgName, int64_t* width, int64_t* height, int64_t* depth, int32_t* dataType)
{
	std::lock_guard<std::mutex> lock(mutex);
	PISystem* sys = (PISystem*)pi;
	coord_t w = 0, h = 0, d = 0;
	ImageDataType dt = ImageDataType::Unknown;
	sys->getImageInfoNoThrow(imgName, w, h, d, dt);
	*width = w;
	*height = h;
	*depth = d;
	*dataType = (int32_t)dt;
}

void* getImage(void* pi, const char* imgName, int64_t* width, int64_t* height, int64_t* depth, int32_t* dataType)
{
	std::lock_guard<std::mutex> lock(mutex);
	PISystem* sys = (PISystem*)pi;
	ImageBase* img = sys->getImageNoThrow(imgName);

	if (!img)
	{
		*width = 0;
		*height = 0;
		*depth = 0;
		*dataType = (int32_t)ImageDataType::Unknown;
		return 0;
	}

	*width = img->width();
	*height = img->height();
	*depth = img->depth();
	*dataType = (int)img->dataType();

	return img->getRawData();
}

uint8_t finishUpdate(void* pi, const char* imgName)
{
	std::lock_guard<std::mutex> lock(mutex);
	return ((PISystem*)pi)->flushIfDistributedNoThrow(imgName) ? 1 : 0;
}

const char* commandList(void* pi)
{
	// Lock not required
	static string cmdList = CommandList::list(false);
	return cmdList.c_str();
}

const char* help(void* pi, const char* commandName)
{
	// Lock not required
	vector<string> v = CommandList::help(commandName, HelpFormat::Text);
	static string hlp;
	hlp = "";
	for (size_t n = 0; n < v.size(); n++)
	{
		hlp += v[n] + "\n\n";
	}
	return hlp.c_str();
}

const char* getString(void* pi, const char* name)
{
	std::lock_guard<std::mutex> lock(mutex);
	PISystem* sys = (PISystem*)pi;

	std::string* s = sys->getStringNoThrow(name);
	if (s)
		return s->c_str();
	return nullptr;
}