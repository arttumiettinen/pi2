
#include "listandimage.h"

namespace itl2
{
	void listToImage(const std::vector<Vec3f>& list, Image<float32_t>& image)
	{
		image.ensureSize(3, list.size());
#pragma omp parallel for if(!omp_in_parallel() && list.size() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < (coord_t)list.size(); n++)
		{
			const Vec3f& p = list[n];
			image(0, n) = p.x;
			image(1, n) = p.y;
			image(2, n) = p.z;
		}
	}

	void listToImage(const std::vector<Vec2f>& list, Image<float32_t>& image)
	{
		image.ensureSize(2, list.size());
#pragma omp parallel for if(!omp_in_parallel() && list.size() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < (coord_t)list.size(); n++)
		{
			const Vec2f& p = list[n];
			image(0, n) = p.x;
			image(1, n) = p.y;
		}
	}

	void listToImage(const std::vector<Edge>& list, Image<uint64_t>& image)
	{
		image.ensureSize(2, list.size());
#pragma omp parallel for if(!omp_in_parallel() && list.size() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < (coord_t)list.size(); n++)
		{
			const Edge& p = list[n];
			image(0, n) = p.verts.x;
			image(1, n) = p.verts.y;
		}
	}

	void imageToList(const Image<float32_t>& img, std::vector<Vec3f>& list)
	{
		list.clear();
		list.resize(img.height());
#pragma omp parallel for if(!omp_in_parallel() && img.height() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < img.height(); n++)
		{
			list[n].x = img(0, n);
			list[n].y = img(1, n);
			list[n].z = img(2, n);
		}
	}
}