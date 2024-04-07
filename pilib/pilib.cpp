
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

	Value* val = sys->getValueNoThrow(name);

	if (val != nullptr && val->getType() == ValueType::String)
		return val->stringValue.c_str();
	return nullptr;
}

PILIB_API const int64_t getInt(void* pi, const char* name)
{
	std::lock_guard<std::mutex> lock(mutex);
	PISystem* sys = (PISystem*)pi;

	Value* val = sys->getValueNoThrow(name);

	if (val != nullptr && val->getType() == ValueType::Int)
		return val->intValue;
	return std::numeric_limits<int64_t>::lowest();
}

PILIB_API const float32_t getReal(void* pi, const char* name)
{
	std::lock_guard<std::mutex> lock(mutex);
	PISystem* sys = (PISystem*)pi;

	Value* val = sys->getValueNoThrow(name);

	if (val != nullptr && val->getType() == ValueType::Real)
		return (float)val->realValue;
	return std::numeric_limits<float32_t>::lowest();
}

PILIB_API const uint8_t getBool(void* pi, const char* name)
{
	std::lock_guard<std::mutex> lock(mutex);
	PISystem* sys = (PISystem*)pi;

	Value* val = sys->getValueNoThrow(name);

	if (val != nullptr && val->getType() == ValueType::Bool)
		return val->boolValue ? 1 : 0;
	return std::numeric_limits<uint8_t>::max();
}