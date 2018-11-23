#pragma once

#include "command.h"

#include <vector>

#include "itl2.h"

using namespace std;
using namespace itl2;

namespace pilib
{
	/**
	Function used to concatenate command argument lists.
	*/
	template<typename T> vector<T> concat(const vector<T>& a, const vector<T>& b)
	{
		vector<T> res;
		res.insert(res.end(), a.begin(), a.end());
		res.insert(res.end(), b.begin(), b.end());
		return res;
	}

	/**
	Base class for commands that do processing in-place.
	*/
	template<typename input_t> class OneImageInPlaceCommand : public Command
	{
	public:
		OneImageInPlaceCommand(const string& name, const string& help, const vector<CommandArgumentBase>& extraArgs = {}) :
			Command(name, help,
				concat({
					CommandArgument<Image<input_t> >(InOut, "image", "Image to process.")
					}, extraArgs)
			)
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<input_t>& in = *pop<Image<input_t>* >(args);
			run(in, args);
		}

		virtual void run(Image<input_t>& in, vector<ParamVariant>& args) const = 0;
	};

	/**
	Base class for commands that have input and output image.
	*/
	template<typename input_t, typename output_t = input_t> class TwoImageInputOutputCommand : public Command
	{
	public:
		TwoImageInputOutputCommand(const string& name, const string& help, const vector<CommandArgumentBase>& extraArgs = {}) :
			Command(name, help,
				concat({
					CommandArgument<Image<input_t> >(In, "input image", "Input image."),
					CommandArgument<Image<output_t> >(Out, "output image", "Output image.")
					}, extraArgs)
				)
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<input_t>& in = *pop<Image<input_t>* >(args);
			Image<output_t>& out = *pop<Image<output_t>* >(args);
			run(in, out, args);
		}

		virtual void run(Image<input_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const = 0;
	};

	/**
	Base class for commands that have input and param image.
	*/
	template<typename input_t, typename param_t = input_t> class TwoImageInputParamCommand : public OneImageInPlaceCommand<input_t>
	{
	public:
		TwoImageInputParamCommand(const string& name, const string& help, const vector<CommandArgumentBase>& extraArgs = {}) :
			OneImageInPlaceCommand<input_t>(name, help,
				concat({
					CommandArgument<Image<param_t> >(In, "parameter image", "Parameter image.")
					}, extraArgs)
				)
		{
		}

		virtual void run(Image<input_t>& in, vector<ParamVariant>& args) const
		{
			Image<param_t>& out = *pop<Image<param_t>* >(args);
			run(in, out, args);
		}

		virtual void run(Image<input_t>& in, Image<param_t>& param, vector<ParamVariant>& args) const = 0;
	};


	/**
	Base class for commands that process image in-place, and have neighbourhood radius and neighbourhood type parameters.
	*/
	template<typename input_t> class BasicOneImageNeighbourhoodCommand : public OneImageInPlaceCommand<input_t>
	{
	public:
		BasicOneImageNeighbourhoodCommand(const string& name, const string& help, vector<CommandArgumentBase> extraArgs = {}) :
			OneImageInPlaceCommand<input_t>(name, help,
				concat(concat({ CommandArgument<coord_t>(In, "radius", "Radius of neighbourhood. Diameter will be 2*r+1.", 1) },
					extraArgs),
					{ CommandArgument<NeighbourhoodType>(In, "neighbourhood type", "Type of neighbourhood, either Ellipsoidal or Rectangular.", Ellipsoidal),
					  CommandArgument<BoundaryCondition>(In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", Zero) })
				)
		{
		}

		virtual void run(Image<input_t>& in, vector<ParamVariant>& args) const
		{
			coord_t r = pop<coord_t>(args);
			NeighbourhoodType nbtype = pop<NeighbourhoodType>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);

			run(in, r, nbtype, bc, args);
		}

		virtual void run(Image<input_t>& in, coord_t r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const = 0;
	};
}
