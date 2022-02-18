#pragma once

#include <string>
#include <vector>
#include <tuple>

#include "image.h"
#include "exprtk/exprtk.hpp"
#include "iteration.h"

namespace itl2
{
	namespace internals
	{
		using symbol_table_t = exprtk::symbol_table<double>;
		using expression_t = exprtk::expression<double>;
		using parser_t = exprtk::parser<double>;

		/**
		Parses string expression into symbol table and expression object.
		Throws ITLException if parse fails.
		*/
		inline std::tuple<symbol_table_t, expression_t> parse(const std::string& expression, const std::vector<ImageBase*>& params, std::vector<double>& varValues)
		{
			symbol_table_t symbol_table;

			varValues.reserve(varValues.size() + params.size());
			for (size_t n = 0; n < params.size(); n++)
			{
				std::string varName = std::string("x") + toString(n);
				varValues.push_back(0);
				symbol_table.add_variable(varName, varValues[n]);
			}
			symbol_table.add_constants();
			expression_t expr;
			expr.register_symbol_table(symbol_table);
			parser_t parser;
			if (!parser.compile(expression, expr))
				throw ITLException(string("Unable to parse expression: ") + parser.error());

			return std::make_tuple(symbol_table, expr);
		}
	}

	/**
	Evaluates mathematical expression, given as a string, on each pixel of the target and the argument images.
	The result is assigned into the corresponding pixel of the target image.
	The parameter images can be referenced in the expression using terms x0, x1, x2, etc.
	*/
	template<typename target_t> void eval(const std::string& expression, Image<target_t>& target, const std::vector<ImageBase*>& params)
	{
		// Parse so that syntax can be checked before entering processing threads
		std::vector<double> varValuesDummy;
		internals::parse(expression, params, varValuesDummy);

		// Check image sizes
		for (size_t n = 0; n < params.size(); n++)
		{
			target.checkSize(*params[n]);
		}

		#pragma omp parallel if(target.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			// Parse again to make separate expression object for each thread.
			std::vector<double> varValues;
			internals::symbol_table_t symbols;
			internals::expression_t expr;
			// Note: This does not work in clang
			//auto [symbols, expr] = internals::parse(expression, params, varValues);
			std::tie(symbols, expr) = internals::parse(expression, params, varValues);

			#pragma omp for
			for (coord_t n = 0; n < target.pixelCount(); n++)
			{
				for (size_t m = 0; m < params.size(); m++)
					varValues[m] = params[m]->getf(n);

				double result = expr.value();
				target(n) = pixelRound<target_t>(result);
			}
		}
	}

	namespace tests
	{
		void eval();
	}
}