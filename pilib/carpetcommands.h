#pragma once

#include "carpet.h"

#include "command.h"
#include "standardhelp.h"

using namespace itl2;

namespace pilib
{
	inline std::string surfaceSeeAlso()
	{
		return "findsurface, drawheightmap, setbeforeheightmap, setafterheightmap, shiftz";
	}

	template<typename pixel_t> class FindSurfaceCommand : public Command
	{
	protected:
		friend class CommandList;

		FindSurfaceCommand() : Command("findsurface",
			"Surface recognition algorithm 'Carpet' according to Turpeinen - Interface Detection Using a Quenched-Noise Version of the Edwards-Wilkinson Equation. "
			"The algorithm places a surface above (alternatively below) the image, and moves it towards larger (alternatively smaller) $z$ values while controlling its dynamics according to the Edwards-Wilkinson equation. "
			"The movement of the surface stops when it encounters enough pixels with value above specific stopping gray level. "
			"The surface does not move through small holes in the object as it has controllable amount of surface tension. ",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "geometry", "Geometry image. This image does not need to be binary image (but it can be)."),
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "height map", R"(Height map that defines the initial position of the surface, or an empty image. The height map gives $z$-position of the surface for each $(x,y)$-position of the geometry image. The final position of the surface is also saved into this height map. If used as an input, the size of the height map must be $w \times h$ where $w$ and $h$ are the width and the height of the geometry image. If the size is not correct, the height map will be zeroed and set to the correct size.)"),
				CommandArgument<double>(ParameterDirection::In, "stopping gray value", "The movement of the surface stops when it encounters pixels whose value is above this gray level."),
				CommandArgument<string>(ParameterDirection::In, "direction", "Direction where the surface moves. 'Down' corresponds to direction towards larger $z$ values, and 'Up' corresponds to direction towards smaller $z$ values.", "Down"),
				CommandArgument<double>(ParameterDirection::In, "surface tension", "Value that indicates how smooth the surface will be. Larger value results in smoother surface.", 1.0),
				CommandArgument<size_t>(ParameterDirection::In, "iterations", "Count of iterations to perform.", 150),

				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "visualization", R"(An image where a visualization of the evolution of the surface will be saved. The dimensions of the visualization will be set to $w \times d \times N$, where $w$ and $d$ are width and depth of the geometry image, and $N$ is the count of iterations to perform.)"),
				CommandArgument<size_t>(ParameterDirection::In, "visualization y", "Indicates the $y$-coordinate of the $xz$-slice that will be visualized."),
				CommandArgument<double>(ParameterDirection::In, "visualization color", "Color of the surface in the visualization. If set to zero, the color will be set to one above the maximum in geometry image."),
			},
			surfaceSeeAlso())
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& geometry = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& hmap = *pop<Image<float32_t>* >(args);
			double stoppingValue = pop<double>(args);
			string directions = pop<string>(args);
			double surfaceTension = pop<double>(args);
			size_t iterations = pop<size_t>(args);

			Image<pixel_t>* visualization = pop<Image<pixel_t>*>(args);
			size_t visY = pop<size_t>(args);
			double visColor = pop<double>(args);
			
			Direction dir = fromString<Direction>(directions);

			findSurface(geometry, hmap, stoppingValue, dir, surfaceTension, iterations, visualization, visY, pixelRound<pixel_t>(visColor));
		}
	};


	template<typename pixel_t> class FindSurface2Command : public Command
	{
	protected:
		friend class CommandList;

		FindSurface2Command() : Command("findsurface",
			"Surface recognition algorithm 'Carpet' according to Turpeinen - Interface Detection Using a Quenched-Noise Version of the Edwards-Wilkinson Equation. "
			"The algorithm places a surface above (alternatively below) the image, and moves it towards larger (alternatively smaller) $z$ values while controlling its dynamics according to the Edwards-Wilkinson equation. "
			"The movement of the surface stops when it encounters enough pixels with value above specific stopping gray level. "
			"The surface does not move through small holes in the object as it has controllable amount of surface tension. ",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "geometry", "Geometry image. This image does not need to be binary image (but it can be)."),
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "height map", R"(Height map that defines the initial position of the surface, or an empty image. The height map gives $z$-position of the surface for each $(x,y)$-position of the geometry image. The final position of the surface is also saved into this height map. The size of the height map must be $w \times h$ where $w$ and $h$ are the width and the height of the geometry image.)"),
				CommandArgument<double>(ParameterDirection::In, "stopping gray value", "The movement of the surface stops when it encounters pixels whose value is above this gray level."),
				CommandArgument<string>(ParameterDirection::In, "direction", "Direction where the surface moves. 'Down' corresponds to direction towards larger $z$ values, and 'Up' corresponds to direction towards smaller $z$ values.", "Down"),
				CommandArgument<double>(ParameterDirection::In, "surface tension", "Value that indicates how smooth the surface will be. Larger value results in smoother surface.", 1.0),
				CommandArgument<size_t>(ParameterDirection::In, "iterations", "Count of iterations to perform.", 150),
			},
			surfaceSeeAlso())
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& geometry = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& hmap = *pop<Image<float32_t>* >(args);
			double stoppingValue = pop<double>(args);
			string directions = pop<string>(args);
			double surfaceTension = pop<double>(args);
			size_t iterations = pop<size_t>(args);

			Direction dir = fromString<Direction>(directions);

			findSurface(geometry, hmap, stoppingValue, dir, surfaceTension, iterations);
		}
	};



	template<typename pixel_t, void (*OP)(Image<pixel_t>&, const Image<float32_t>&, pixel_t)> class HeightMapCommandBase : public Command
	{
	protected:
		friend class CommandList;

		HeightMapCommandBase(const string& name, const string& help) : Command(name,
			help,
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::InOut, "geometry", "Image where whose pixels are to be set."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "height map", R"(Height map. The size of the height map must be $w \times h$ where $w$ and $h$ are the width and the height of the geometry image.)"),
				CommandArgument<double>(ParameterDirection::In, "visualization color", "Color of the surface in the visualization."),
			},
			surfaceSeeAlso())
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& geometry = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& hmap = *pop<Image<float32_t>* >(args);
			double color = pop<double>(args);
			
			OP(geometry, hmap, itl2::pixelRound<pixel_t>(color));
		}
	};

	template<typename pixel_t> class DrawHeightMapCommand : public HeightMapCommandBase<pixel_t, itl2::drawHeightMap>
	{
	protected:
		friend class CommandList;

		DrawHeightMapCommand() : HeightMapCommandBase<pixel_t, itl2::drawHeightMap>("drawheightmap", "Draws a height map to an image.")
		{
		}
	};

	template<typename pixel_t> class SetBeforeHeightMapCommand : public HeightMapCommandBase<pixel_t, itl2::setBeforeHeightMap>
	{
	protected:
		friend class CommandList;

		SetBeforeHeightMapCommand() : HeightMapCommandBase<pixel_t, itl2::setBeforeHeightMap>("setbeforeheightmap", "Sets values of all pixel located before (above) the given height map.")
		{
		}
	};
	
	template<typename pixel_t> class SetAfterHeightMapCommand : public HeightMapCommandBase<pixel_t, itl2::setAfterHeightMap>
	{
	protected:
		friend class CommandList;

		SetAfterHeightMapCommand() : HeightMapCommandBase<pixel_t, itl2::setAfterHeightMap>("setafterheightmap", "Sets values of all pixels located after (below) the given height map.")
		{
		}
	};

	template<typename pixel_t> class ShiftZCommand : public Command
	{
	protected:
		friend class CommandList;

		ShiftZCommand() : Command("shiftz",
			"Shift each $z$-directional column of input image by amount given in the shift map.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::InOut, "image", "Image whose $z$-directional columns are to be shifted."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "shift map", R"(Shift map that gives the amount of shift to apply in each $z$-directional column. The size of the shift map must be $w \times h$ where $w$ and $h$ are the width and the height of the geometry image.)"),
				CommandArgument<bool>(ParameterDirection::In, "subtract mean", R"(Set to true to automatically subtract the average value of the shift map from each shift. This is useful if the shift map is negation of a surface map found using the findsurface command, and the intent is to make the surface straight.)"),
				CommandArgument<InterpolationMode>(ParameterDirection::In, "interpolation mode", string("Interpolation mode. ") + interpolationHelp(), InterpolationMode::Linear),
				CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Zero),
			},
			surfaceSeeAlso())
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& geometry = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& hmap = *pop<Image<float32_t>* >(args);
			bool subtractMean = pop<bool>(args);
			InterpolationMode imode = pop<InterpolationMode>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);

			itl2::shiftZ(geometry, hmap, subtractMean, *createInterpolator<pixel_t, pixel_t>(imode, bc));
		}
	};
}
