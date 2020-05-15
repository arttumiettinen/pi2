#pragma once

#if defined(__linux__)

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#elif defined(_WIN32)

#include <filesystem>
namespace fs = std::filesystem;

#else

#error filesystem.h not configured for this platform.

#endif