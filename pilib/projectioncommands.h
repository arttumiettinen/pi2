#pragma once

#include "command.h"
#include "commandsbase.h"
#include "overlapdistributable.h"
#include "projections.h"
#include "math/vectoroperations.h"

#include <string>

using std::string;


namespace pilib
{
	namespace internals
	{
		/*
		Get value at the first line that reads "key = value".
		*/
		inline double getValue(const string& s, const string& key)
		{
			size_t startPos = s.find(key);
			if (startPos != string::npos)
			{
				startPos += key.length();
				size_t lineEnd = s.find('\n', startPos);
				string number = s.substr(startPos, lineEnd - startPos);
				return fromString<double>(number);
			}
			else
			{
				throw ITLException("Key not found.");
			}
		}

		inline double sumReducer(const vector<double>& vals, const vector<double>& counts)
		{
			return math::sum(vals);
		}

		inline double minReducer(const vector<double>& vals, const vector<double>& counts)
		{
			return math::min(vals);
		}

		inline double maxReducer(const vector<double>& vals, const vector<double>& counts)
		{
			return math::max(vals);
		}

		inline double meanReducer(const vector<double>& vals, const vector<double>& counts)
		{
			return math::mean(vals, counts);
		}
	}

#define CONCAT(str1, str2) #str1 #str2
#define DEF_PROJECT(classname, funcname, cmdname) \
	template<typename in_t> class classname##AllPixelsCommand : public TwoImageInputOutputCommand<in_t, float32_t>, public Distributable	\
	{																										\
	public:																									\
		classname##AllPixelsCommand() : TwoImageInputOutputCommand<in_t, float32_t>(#cmdname, "Calculates " #funcname " of all pixels in the input image. The output is a 1x1x1 image.",	\
			{																								\
				CommandArgument<bool>(ParameterDirection::In, "print to log", "Set to true to print the results to the log.", false)	\
			}) {}																							\
																											\
		virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const			\
		{																									\
			bool print = pop<bool>(args);																	\
			double res = (double)itl2:: funcname (in);														\
			out.ensureSize(1, 1, 1);																		\
			out(0) = (float32_t)res;																		\
			if(print)																						\
			{																								\
				cout << #funcname << " = " << out(0) << endl;												\
				cout << "count = " << in.pixelCount() << endl;												\
			}																								\
		}																									\
																											\
		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const	\
		{																									\
			DistributedImage<float32_t>& out = *get<DistributedImage<float32_t>* >(args[1]);				\
			bool print = get<bool>(args[2]);																\
			out.ensureSize(1, 1, 1);																		\
																											\
			vector<ParamVariant> args2;																		\
			args2.push_back(args[0]);																		\
			args2.push_back(args[1]);																		\
			ParamVariant p;																					\
			p = true;																						\
			args2.push_back(p);																				\
																											\
			vector<string> results = distributor.distribute(this, args2, 2, Vec3c(0, 0, 0), 0);				\
			vector<double> vals, counts;																	\
			for (const string& s : results)																	\
			{																								\
				vals.push_back(internals::getValue(s, #funcname " = "));									\
				counts.push_back(internals::getValue(s, "count = "));										\
			}																								\
																											\
			Image<float32_t> tmp;																			\
			out.readTo(tmp);																				\
			tmp(0) = (float32_t)internals::funcname##Reducer(vals, counts);									\
			out.setData(tmp);																				\
																											\
			if (print)																						\
			{																								\
				cout << #funcname << " = " << tmp(0) << endl;												\
				cout << "count = " << sum(counts) << endl;													\
			}																								\
																											\
			return vector<string>();																		\
		}																									\
																											\
		virtual void getCorrespondingBlock(vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const	\
		{																									\
			if (argIndex == 1)																				\
			{																								\
				readStart = Vec3c(0, 0, 0);																	\
				readSize = Vec3c(1, 1, 1);																	\
				writeFilePos = Vec3c(0, 0, 0);																\
				writeImPos = Vec3c(0, 0, 0);																\
				writeSize = Vec3c(1, 1, 1);																	\
			}																								\
		}																									\
	};																										\
																											\
	template<typename in_t> class classname##ProjectCommand : public TwoImageInputOutputCommand<in_t, float32_t>, public Distributable	\
	{																										\
	public:																									\
		classname##ProjectCommand() : TwoImageInputOutputCommand<in_t, float32_t>(CONCAT(funcname, project), "Calculates projection of the input image. The dimensionality of the output image is the dimensionality of the input image subtracted by one.",	\
		{ CommandArgument<size_t>(ParameterDirection::In, "dimension", "Dimension to project over, zero corresponding to x, one corresponding to y, and 2 corresponding to z.", 2) }) {}	\
																											\
		virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const			\
		{																									\
			size_t dim = pop<size_t>(args);																    \
			if(dim > 2)																			            \
				throw ITLException("Invalid dimension specification.");										\
																											\
			itl2:: funcname (in, dim, out);																	\
		}																									\
																											\
		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const	\
		{																									\
			DistributedImage<in_t>& in = *get<DistributedImage<in_t>* >(args[0]);							\
			DistributedImage<float32_t>& out = *get<DistributedImage<float32_t>* >(args[1]);				\
			size_t dim = get<size_t>(args[2]);															    \
																											\
			size_t distrDim;																				\
			if (dim == 2)																					\
			{																								\
				/* Z project */																				\
				out.ensureSize(in.width(), in.height());													\
				distrDim = 1;																				\
			}																								\
			else if (dim == 1)																				\
			{																								\
				/* Y project */																				\
				out.ensureSize(in.width(), in.depth());														\
				distrDim = 1;																				\
			}																								\
			else if (dim == 0)																				\
			{																								\
				/* X project. Swap z and y in output to make the image more logical (in some sense...) */	\
				out.ensureSize(in.depth(), in.height());													\
				distrDim = 1;																				\
			}																								\
			else																							\
			{																								\
				throw ITLException("Invalid dimension specification.");										\
			}																								\
																											\
			return distributor.distribute(this, args, distrDim, Vec3c(0, 0, 0));							\
		}																									\
																											\
		virtual void getCorrespondingBlock(vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const	\
		{																									\
			if (argIndex == 0)																				\
			{																								\
				DistributedImage<in_t>& in = *get<DistributedImage<in_t>* >(args[0]);						\
				size_t dim = get<size_t>(args[2]);															\
				if(dim == 2)																				\
				{																							\
					readSize = Vec3c(readSize.x, readSize.y, in.depth());									\
					readStart = Vec3c(readStart.x, readStart.y, 0);											\
				}																							\
				else if (dim == 1)																			\
				{																							\
					readSize = Vec3c(readSize.x, in.height(), readSize.y);									\
					readStart = Vec3c(readStart.x, 0, readStart.y);											\
				}																							\
				else /* dim == 0 */																			\
				{																							\
					readSize = Vec3c(in.width(), readSize.y, readSize.x);									\
					readStart = Vec3c(0, readStart.y, readStart.x);											\
				}																							\
			}																								\
		}																									\
	};

#undef min
#undef max
	DEF_PROJECT(Sum, sum, sum)
	DEF_PROJECT(Min, min, minval)
	DEF_PROJECT(Max, max, maxval)
	DEF_PROJECT(Mean, mean, mean)
	




	// TODO: Variance is not available at the moment
	//template<typename in_t> class VarianceAllPixelsCommand : public TwoImageInputOutputCommand<in_t, float32_t>	
	//{																										
	//public:																									
	//	VarianceAllPixelsCommand() : TwoImageInputOutputCommand<in_t, float32_t>("variance", "Calculates variance of all pixels in the input image. The output is a 1x1x1 image.",
	//		{																								
	//			CommandArgument<bool>(ParameterDirection::In, "print to log", "Set to true to print the results to the log.", false)	
	//		}) {}																							
	//																										
	//	virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const			
	//	{																									
	//		bool print = pop<bool>(args);																	
	//		double res = (double)variance(in);																
	//		out.ensureSize(1, 1, 1);																		
	//		out(0) = (float32_t)res;																		
	//		if(print)																						
	//		{																								
	//			cout << "variance = " << out(0) << endl;												
	//			cout << "count = " << in.pixelCount() << endl;												
	//		}																								
	//	}																									
	//																										
	//};

	// TODO: Standard deviation projection is not available at the moment
	// TODO: Standard deviation is not distributed at the moment
	template<typename in_t> class StdDevAllPixelsCommand : public TwoImageInputOutputCommand<in_t, float32_t>
	{
	public:
		StdDevAllPixelsCommand() : TwoImageInputOutputCommand<in_t, float32_t>("stddev", "Calculates standard deviation of all pixels in the input image. The output is a 1x1x1 image.",
			{
				CommandArgument<bool>(ParameterDirection::In, "print to log", "Set to true to print the results to the log.", false)
			}) {}

		virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const
		{
			bool print = pop<bool>(args);
			double res = (double)itl2::stdDev(in);
			out.ensureSize(1, 1, 1);
			out(0) = (float32_t)res;
			if (print)
			{
				cout << "stddev = " << out(0) << endl;
				cout << "count = " << in.pixelCount() << endl;
			}
		}

	};
}