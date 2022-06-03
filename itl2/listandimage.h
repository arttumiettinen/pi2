#pragma once

#include <vector>
#include "image.h"
#include <network.h>

namespace itl2
{
	void listToImage(const std::vector<Vec3f>& list, Image<float32_t>& image);

	void listToImage(const std::vector<Vec2f>& list, Image<float32_t>& image);

	void listToImage(const std::vector<Edge>& list, Image<uint64_t>& image);

	void imageToList(const Image<float32_t>& img, std::vector<Vec3f>& list);
}