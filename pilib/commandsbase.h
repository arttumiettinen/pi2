#pragma once

#include "command.h"

#include <vector>


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
					CommandArgument<Image<input_t> >(ParameterDirection::InOut, "image", "Image to process.")
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
					CommandArgument<Image<input_t> >(ParameterDirection::In, "input image", "Input image."),
					CommandArgument<Image<output_t> >(ParameterDirection::Out, "output image", "Output image.")
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
					CommandArgument<Image<param_t> >(ParameterDirection::In, "parameter image", "Parameter image.")
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
				concat(concat({ CommandArgument<coord_t>(ParameterDirection::In, "radius", "Radius of neighbourhood. Diameter will be 2*r+1.", 1) },
					extraArgs),
					{ CommandArgument<NeighbourhoodType>(ParameterDirection::In, "neighbourhood type", "Type of neighbourhood, either Ellipsoidal or Rectangular.", NeighbourhoodType::Ellipsoidal),
					  CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", BoundaryCondition::Nearest) })
				)
		{
		}

		virtual void run(Image<input_t>& in, vector<ParamVariant>& args) const
		{
			coord_t r = pop<coord_t>(args);

			NeighbourhoodType nbtype = get<NeighbourhoodType>(args[args.size() - 2]);
			BoundaryCondition bc = get<BoundaryCondition>(args[args.size() - 1]);
			args.erase(args.begin() + args.size() - 1);
			args.erase(args.begin() + args.size() - 1);

			run(in, r, nbtype, bc, args);
		}

		virtual void run(Image<input_t>& in, coord_t r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const = 0;
	};
}
