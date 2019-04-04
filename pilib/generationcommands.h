#pragma once

#include "commandsbase.h"
#include "generation.h"

namespace pilib
{

	template<typename pixel_t> class RampCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		RampCommand() : OneImageInPlaceCommand<pixel_t>("ramp", "Fills image with ramp in given dimension, i.e. performs img[r] = r[dimension].",
			{
				CommandArgument<coord_t>(ParameterDirection::In, "dimension", "Dimension of the ramp.", 0)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			coord_t dim = pop<coord_t>(args);

			ramp(in, dim);
		}
	};


	template<typename pixel_t> class BoxCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		BoxCommand() : OneImageInPlaceCommand<pixel_t>("box", "Draws a filled box into the image.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "position", "Position of the left-top corner of the box.", Vec3c(0, 0, 0)),
				CommandArgument<Vec3c>(ParameterDirection::In, "size", "Size of the box.", Vec3c(10, 10, 10)),
				CommandArgument<double>(ParameterDirection::In, "value", "Value for pixels inside the box.", 1)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			Vec3c pos = pop<Vec3c>(args);
			Vec3c size = pop<Vec3c>(args);
			double value = pop<double>(args);

			draw(in, Box<coord_t>(pos, pos + size), pixelRound<pixel_t>(value));
		}
	};


	template<typename pixel_t> class SphereCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		SphereCommand() : OneImageInPlaceCommand<pixel_t>("sphere", "Draws a filled sphere into the image.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "position", "Position of the center point of the sphere.", Vec3d(0, 0, 0)),
				CommandArgument<double>(ParameterDirection::In, "radius", "Radius of the sphere.", 10.0),
				CommandArgument<double>(ParameterDirection::In, "value", "Value for pixels inside the box.", 1)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			Vec3d pos = pop<Vec3d>(args);
			double radius = pop<double>(args);
			double value = pop<double>(args);

			draw(in, Sphere(pos, radius), pixelRound<pixel_t>(value));
		}
	};

	template<typename pixel_t> class LineCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		LineCommand() : OneImageInPlaceCommand<pixel_t>("line", "Draws a single-pixel wide line into the image.",
			{
				CommandArgument<Vec3d>(ParameterDirection::In, "start", "Start position of the line.", Vec3d(0, 0, 0)),
				CommandArgument<Vec3d>(ParameterDirection::In, "end", "End position of the line.", Vec3d(1, 1, 1)),
				CommandArgument<double>(ParameterDirection::In, "value", "Value for pixels of the line.", 1)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const
		{
			Vec3d start = pop<Vec3d>(args);
			Vec3d end = pop<Vec3d>(args);
			double value = pop<double>(args);

			draw(in, Line(start, end), pixelRound<pixel_t>(value));
		}
	};
}
