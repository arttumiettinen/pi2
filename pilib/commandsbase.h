#pragma once

#include "command.h"
#include "standardhelp.h"

#include <vector>

using namespace itl2;

namespace pilib
{
	/**
	Function used to concatenate command argument lists.
	*/
	template<typename T> std::vector<T> concat(const std::vector<T>& a, const std::vector<T>& b)
	{
		std::vector<T> res;
		res.insert(res.end(), a.begin(), a.end());
		res.insert(res.end(), b.begin(), b.end());
		return res;
	}

	/**
	Base class for commands that do processing in-place.
	*/
	template<typename input_t> class OneImageInPlaceCommand : public Command
	{
	protected:
		friend class CommandList;

		OneImageInPlaceCommand(const string& name, const string& help, const std::vector<CommandArgumentBase>& extraArgs = {}, const string& seeAlso = "") :
			Command(name, help,
				concat({
					CommandArgument<Image<input_t> >(ParameterDirection::InOut, "image", "Image to process.")
					}, extraArgs),
				seeAlso
			)
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<input_t>& in = *pop<Image<input_t>* >(args);
			run(in, args);
		}

		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const = 0;
	};

	/**
	Base class for commands that have input and output image.
	*/
	template<typename input_t, typename output_t = input_t> class TwoImageInputOutputCommand : public Command
	{
	protected:
		friend class CommandList;

		TwoImageInputOutputCommand(const string& name, const string& help, const std::vector<CommandArgumentBase>& extraArgs = {}, const string& seeAlso = "") :
			Command(name, help,
				concat({
					CommandArgument<Image<input_t> >(ParameterDirection::In, "input image", "Input image."),
					CommandArgument<Image<output_t> >(ParameterDirection::Out, "output image", "Output image.")
					}, extraArgs),
				seeAlso
				)
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<input_t>& in = *pop<Image<input_t>* >(args);
			Image<output_t>& out = *pop<Image<output_t>* >(args);
			run(in, out, args);
		}

		virtual void run(Image<input_t>& in, Image<output_t>& out, std::vector<ParamVariant>& args) const = 0;
	};

	/**
	Base class for commands that have input and param image.
	*/
	template<typename input_t, typename param_t = input_t> class TwoImageInputParamCommand : public OneImageInPlaceCommand<input_t>
	{
	protected:
		friend class CommandList;

		TwoImageInputParamCommand(const string& name, const string& help, const std::vector<CommandArgumentBase>& extraArgs = {}, const std::string& seeAlso = "") :
			OneImageInPlaceCommand<input_t>(name, help,
				concat({
					CommandArgument<Image<param_t> >(ParameterDirection::In, "parameter image", "Parameter image.")
					}, extraArgs),
				seeAlso
				)
		{
		}

	public:
		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			Image<param_t>& out = *pop<Image<param_t>* >(args);
			run(in, out, args);
		}

		virtual void run(Image<input_t>& in, Image<param_t>& param, std::vector<ParamVariant>& args) const = 0;
	};


	/**
	Base class for commands that process image in-place, and have neighbourhood radius and neighbourhood type parameters.
	*/
	template<typename input_t> class BasicOneImageNeighbourhoodCommand : public OneImageInPlaceCommand<input_t>
	{
	protected:
		friend class CommandList;

		BasicOneImageNeighbourhoodCommand(const string& name, const string& help, std::vector<CommandArgumentBase> extraArgs = {}) :
			OneImageInPlaceCommand<input_t>(name, help,
				concat(concat({ CommandArgument<coord_t>(ParameterDirection::In, "radius", "Radius of neighbourhood. Diameter will be 2*r+1.", 1) },
					extraArgs),
					{ CommandArgument<NeighbourhoodType>(ParameterDirection::In, "neighbourhood type", string("Type of neighbourhood. ") + neighbourhoodTypeHelp(), NeighbourhoodType::Ellipsoidal),
					  CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest) })
				)
		{
		}

	public:
		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			coord_t r = pop<coord_t>(args);

			NeighbourhoodType nbtype = std::get<NeighbourhoodType>(args[args.size() - 2]);
			BoundaryCondition bc = std::get<BoundaryCondition>(args[args.size() - 1]);
			args.erase(args.begin() + args.size() - 1);
			args.erase(args.begin() + args.size() - 1);

			run(in, r, nbtype, bc, args);
		}

		virtual void run(Image<input_t>& in, coord_t r, NeighbourhoodType nbtype, BoundaryCondition bc, std::vector<ParamVariant>& args) const = 0;
	};
}
