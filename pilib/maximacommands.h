#pragma once

#include "commandsbase.h"
#include "standardhelp.h"

#include "maxima.h"

using namespace itl2;

namespace pilib
{
	inline std::string maximaSeeAlso()
	{
		return "localmaxima, cleanmaxima, labelmaxima, growlabels, grow";
	}

	/**
	Packs region list into an image.
	Format: [count of lists][count of items in list 1][x1][y1][z1][zx2][y2][z2]...[count of items in list 2][x1][y1][z1][zx2][y2][z2]...
	*/
	inline void packToImage(const vector<vector<Vec3sc>>& lists, Image<int32_t>& img)
	{
		if (lists.size() > std::numeric_limits<int32_t>::max())
			throw ITLException("Too many regions.");

		size_t count = 0;
		count++; // Particle count
		for (const auto& list : lists)
		{
			count += 1 + list.size() * 3; // List size + elements
		}

		img.ensureSize(count);

		coord_t n = 0;

		// Count of lists
		img(n) = (int32_t)lists.size();
		n++;

		for (const auto& list : lists)
		{
			// Count of items
			img(n) = (int32_t)list.size();
			n++;

			// 3 values per item
			for (const Vec3sc& v : list)
			{
				img(n) = v.x;
				n++;
				img(n) = v.y;
				n++;
				img(n) = v.z;
				n++;
			}
			
		}
	}

	/**
	Unpacks region list from an image previously packed by packToImage.
	*/
	inline vector<vector<Vec3sc>> unpackFromImage(const Image<int32_t>& img)
	{
		vector<vector<Vec3sc>> lists;

		if (img.pixelCount() < 1)
			throw ITLException("Region list contains no count information.");

		size_t n = 0;
		int32_t count = img(n);
		n++;

		for (int32_t i = 0; i < count; i++)
		{
			vector<Vec3sc> list;

			int32_t lcount = img(n);
			n++;

			if ((size_t)img.pixelCount() < n + 3 * (size_t)lcount)
				throw ITLException("Region list is inconsistent. It does not contain enough values.");

			for (int32_t m = 0; m < lcount; m++)
			{
				Vec3sc v;
				v.x = img(n);
				n++;
				v.y = img(n);
				n++;
				v.z = img(n);
				n++;
				list.push_back(v);
			}

			lists.push_back(list);
		}

		return lists;
	}

	template<typename pixel_t> class LocalMaximaCommand : public TwoImageInputOutputCommand<pixel_t, int32_t>
	{
	protected:
		friend class CommandList;

		LocalMaximaCommand() : TwoImageInputOutputCommand<pixel_t, int32_t>("localmaxima",
"Finds local maxima in the input image. "
"Maxima migh be individual pixels or larger regions that have the same value and that are bordered by pixels of smaller value. "
"The output image will be in format [count of regions][count of pixels in region 1][x1][y1][z1][zx2][y2][z2]...[count of items in region 2][x1][y1][z1][zx2][y2][z2]...",
			{
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", "Connectivity of the maxima regions. " + connectivityHelp(), Connectivity::AllNeighbours)
			},
			maximaSeeAlso())
		{
		}
	public:
		virtual void run(Image<pixel_t>& in, Image<int32_t>& out, std::vector<ParamVariant>& args) const override
		{
			Connectivity connectivity = pop<Connectivity>(args);

			vector<vector<Vec3sc> > maxima = findLocalMaxima(in, connectivity);

			packToImage(maxima, out);
		}
	};


	template<typename pixel_t> class CleanMaximaCommand : public Command
	{
	protected:
		friend class CommandList;

		CleanMaximaCommand() : Command("cleanmaxima",
"Removes all maxima that are smaller in radius than neighbouring maximum. "
"Maximum is neighbour to another maximum if distance between them is less than radius of the larger maximum multiplied by radiusMultiplier. "
"Removes all maxima $m$ that satisfy $distance(m, n) < radiusMultiplier * radius(n)$ and $radius(n) > radius(m)$ for some $n$. "
"The distance is measured between centroids of the maxima. "
"The maxima are removed by combining them to the larger maxima.",
			{
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "image", "The image where the maxima have been extracted (by `localmaxima` command)."),
				CommandArgument<Image<int32_t>>(ParameterDirection::InOut, "maxima", "Image that contains the maxima. See output from `localmaxima` command."),
				CommandArgument<double>(ParameterDirection::In, "radius multiplier", "Maxima are enlarged by this multiplier.", 1.0),
			},
			maximaSeeAlso())
		{
		}
	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& img = *pop<Image<pixel_t>*>(args);
			Image<int32_t>& maxima = *pop<Image<int32_t>*>(args);
			double mul = pop<double>(args);

			vector<vector<Vec3sc>> arr = unpackFromImage(maxima);

			removeMaximaInsideLargerOnes(arr, img, mul);

			packToImage(arr, maxima);
		}
	};



	template<typename pixel_t> class DrawMaximaCommand : public Command
	{
	protected:
		friend class CommandList;

		DrawMaximaCommand() : Command("labelmaxima",
			"Draws all regions in the given list, each with different color. If there are not enough colors available, an error is shown.",
			{
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "image", "Image where the maxima are to be drawn."),
				CommandArgument<Image<int32_t>>(ParameterDirection::In, "regions", "Image that contains the regions. The format of this image is described in the `localmaxima`' command.")
			},
			maximaSeeAlso())
		{
		}
	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& img = *pop<Image<pixel_t>*>(args);
			Image<int32_t>& maxima = *pop<Image<int32_t>*>(args);
			
			vector<vector<Vec3sc>> arr = unpackFromImage(maxima);

			draw(img, arr);
		}
	};

}
