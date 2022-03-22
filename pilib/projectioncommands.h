#pragma once

#include "projections.h"
#include "math/vectoroperations.h"

#include "command.h"
#include "commandsbase.h"
#include "overlapdistributable.h"
#include "commandlist.h"
#include "distributedtempimage.h"

#include <string>

using std::string;


namespace pilib
{
	namespace internals
	{
		/*
		Get value at the first line that reads "key value".
		*/
		template<typename result_t> result_t getValue(const string& s, const string& key)
		{
			size_t startPos = s.find(key);
			if (startPos != string::npos)
			{
				startPos += key.length();
				size_t lineEnd = s.find('\n', startPos);
				string number = s.substr(startPos, lineEnd - startPos);
				return fromString<result_t>(number);
			}
			else
			{
				throw ITLException("Key not found.");
			}
		}

		template<typename result_t> result_t sumReducer(const vector<result_t>& vals, const vector<size_t>& counts)
		{
			//return sum<result_t, result_t>(vals);
			result_t res = 0;
			for (result_t val : vals)
				res = NumberUtils<result_t>::saturatingAdd(res, val);
			return res;
		}

		template<typename result_t> result_t squareSumReducer(const vector<result_t>& vals, const vector<size_t>& counts)
		{
			return sum<result_t, result_t>(vals);
		}

		template<typename result_t> result_t minReducer(const vector<result_t>& vals, const vector<size_t>& counts)
		{
			return min(vals);
		}

		template<typename result_t> result_t maxReducer(const vector<result_t>& vals, const vector<size_t>& counts)
		{
			return max(vals);
		}

		template<typename result_t> result_t meanReducer(const vector<result_t>& vals, const vector<size_t>& counts)
		{
			return mean<result_t, size_t>(vals, counts);
		}

		template<typename result_t> result_t maskedMeanReducer(const vector<result_t>& vals, const vector<result_t>& counts)
		{
			return mean<result_t, result_t>(vals, counts);
		}


		/**
		This is used to select suitable distributed accumulator type for sum projections.
		The difference to itl2::sum_intermediate_type is that we must use float32_t instead of double
		in the final accumulation for floating-point images as pi2 does not currently support double/float64_t pixels.
		*/
		template<class pixel_t> struct sum_distributed_intermediate_type {
			using type = typename std::conditional <
				std::is_floating_point_v<pixel_t>,
				float32_t,  // floating-point pixels -> float32_t accumulator
				typename std::conditional<std::is_signed_v<pixel_t>,
				int64_t, // signed integer pixels -> int64 accumulator
				uint64_t // unsigned integer pixels -> uint64 accumulator
				>::type
			>::type;
		};
	}

#define CONCAT(str1, str2) #str1 #str2
#define DEF_PROJECT_ALL(classname, funcname, cmdname, result_type) \
	template<typename in_t, typename result_t = result_type> class classname##AllPixelsCommand : public TwoImageInputOutputCommand<in_t, result_t>, public Distributable	\
	{																										\
	using temp_t = typename internals::sum_distributed_intermediate_type<in_t>::type;						\
	protected:																								\
		friend class CommandList;																			\
																											\
		classname##AllPixelsCommand() : TwoImageInputOutputCommand<in_t, result_t>(#cmdname, "Calculates " #funcname " of all pixels in the input image. The output is a 1x1x1 image.",	\
			{																								\
				CommandArgument<bool>(ParameterDirection::In, "print to log", "Set to true to print the results to the log.", false)	\
			}) {}																							\
	public:																									\
		using Distributable::runDistributed;																\
																											\
		virtual void run(Image<in_t>& in, Image<result_t>& out, vector<ParamVariant>& args) const override	\
		{																									\
			bool print = pop<bool>(args);																	\
			result_t res = itl2:: funcname <in_t, result_t>(in);											\
			out.ensureSize(1, 1, 1);																		\
			out(0) = res;																					\
			if(print)																						\
			{																								\
				std::cout << #funcname << " = " << out(0) << std::endl;										\
				std::cout << "count = " << in.pixelCount() << std::endl;									\
			}																								\
		}																									\
																											\
		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override	\
		{																									\
			DistributedImage<result_t>& out = *std::get<DistributedImage<result_t>* >(args[1]);				\
			bool print = std::get<bool>(args[2]);															\
			out.ensureSize(1, 1, 1);																		\
																											\
			DistributedTempImage<temp_t> dummy(distributor, "project_all_pixels_dummy", Vec3c(1, 1, 1), DistributedImageStorageType::Raw);	\
																											\
			vector<ParamVariant> args2;																		\
			args2.push_back(args[0]);																		\
			args2.push_back(&dummy.get());																	\
			ParamVariant p;																					\
			p = true;																						\
			args2.push_back(p);																				\
																											\
			auto& cmd = CommandList::get<classname##AllPixelsCommand<in_t, temp_t> >();						\
			vector<string> results = distributor.distribute(&cmd, args2);									\
			vector<temp_t> vals;																			\
			vector<size_t> counts;																			\
			for (const string& s : results)																	\
			{																								\
				vals.push_back(internals::getValue<temp_t>(s, #funcname " = "));							\
				counts.push_back(internals::getValue<size_t>(s, "count = "));								\
			}																								\
																											\
			Image<result_t> tmp;																			\
			tmp(0) = pixelRound<result_t>(internals::funcname##Reducer<temp_t>(vals, counts));				\
			out.setData(tmp);																				\
																											\
			if (print)																						\
			{																								\
				std::cout << #funcname << " = " << tmp(0) << std::endl;										\
				std::cout << "count = " << sum<size_t, size_t>(counts) << std::endl;						\
			}																								\
																											\
			return vector<string>();																		\
		}																									\
																											\
		virtual size_t getRefIndex(const vector<ParamVariant>& args) const override							\
		{																									\
			return 0;																						\
		}																									\
																											\
		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override	\
		{																									\
			if (argIndex == 1)																				\
			{																								\
				readStart = Vec3c(0, 0, 0);																	\
				readSize = Vec3c(1, 1, 1);																	\
				writeFilePos = Vec3c(0, 0, 0);																\
				writeImPos = Vec3c(0, 0, 0);																\
				/* Zero write size as the result is not needed in the distributed mode. */						\
				/* We even don't know where to write it, and no processes should write to the same location */	\
				writeSize = Vec3c(0, 0, 0);																	\
			}																								\
		}																									\
																											\
		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override					\
		{																									\
			return JobType::Fast;																			\
		}																									\
	};																										\

#define DEF_PROJECT(classname, funcname, cmdname, result_type)												\
	template<typename in_t> class classname##ProjectCommand : public TwoImageInputOutputCommand<in_t, result_type>, public Distributable	\
	{																										\
	protected:																								\
		friend class CommandList;																			\
																											\
		classname##ProjectCommand() : TwoImageInputOutputCommand<in_t, result_type>(CONCAT(funcname, project), "Calculates projection of the input image. The dimensionality of the output image is the dimensionality of the input image subtracted by one.",	\
		{ CommandArgument<size_t>(ParameterDirection::In, "dimension", "Dimension to project over, zero corresponding to $x$, one corresponding to $y$, and 2 corresponding to $z$.", 2) }) {}	\
	public:																									\
		using Distributable::runDistributed;																\
																											\
		virtual void run(Image<in_t>& in, Image<result_type>& out, vector<ParamVariant>& args) const override	\
		{																									\
			size_t dim = pop<size_t>(args);																    \
			if(dim > 2)																			            \
				throw ITLException("Invalid dimension specification.");										\
																											\
			itl2:: funcname (in, dim, out);																	\
		}																									\
																											\
		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override	\
		{																									\
			DistributedImage<in_t>& in = *std::get<DistributedImage<in_t>* >(args[0]);						\
			DistributedImage<result_type>& out = *std::get<DistributedImage<result_type>* >(args[1]);				\
			size_t dim = std::get<size_t>(args[2]);															\
																											\
			if (dim == 2)																					\
			{																								\
				/* Z project */																				\
				out.ensureSize(in.width(), in.height());													\
			}																								\
			else if (dim == 1)																				\
			{																								\
				/* Y project */																				\
				out.ensureSize(in.width(), in.depth());														\
			}																								\
			else if (dim == 0)																				\
			{																								\
				/* X project. Swap z and y in output to make the image more logical (in some sense...) */	\
				out.ensureSize(in.depth(), in.height());													\
			}																								\
			else																							\
			{																								\
				throw ITLException("Invalid dimension specification.");										\
			}																								\
																											\
			return distributor.distribute(this, args);														\
		}																									\
																											\
		virtual size_t getDistributionDirection1(const vector<ParamVariant>& args) const override			\
		{																									\
			return 1; /* Always distribute along second dimension of output. */								\
		}																									\
																											\
		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override	\
		{																									\
			if (argIndex == 0)																				\
			{																								\
				DistributedImage<in_t>& in = *std::get<DistributedImage<in_t>* >(args[0]);					\
				size_t dim = std::get<size_t>(args[2]);														\
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
																											\
		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override					\
		{																									\
			return JobType::Fast;																			\
		}																									\
	};



#define DEF_PROJECT_2IMAGE(classname, funcname, cmdname, result_t) \
	template<typename in_t, typename val_t> class classname##Project2ImageCommand : public TwoImageInputOutputCommand<in_t, result_t>, public Distributable	\
	{																										\
	protected:																								\
		friend class CommandList;																			\
																											\
		classname##Project2ImageCommand() : TwoImageInputOutputCommand<in_t, result_t>(CONCAT(funcname, project), "Calculates projection of the input image. The dimensionality of the output image is the dimensionality of the input image subtracted by one. Constructs another image from values of second input image taken at location of the extrema.",	\
		{																									\
			CommandArgument<size_t>(ParameterDirection::In, "dimension", "Dimension to project over, zero corresponding to x, one corresponding to y, and 2 corresponding to z.", 2),\
			CommandArgument<Image<val_t> >(ParameterDirection::In, "input value image", "Image where values of the extra output are taken."), \
			CommandArgument<Image<val_t> >(ParameterDirection::Out, "output value image", "Value of this image will equal value of 'input value image' at location of the extrema (whose value is stored in the corresponding output image).") \
		}) {}																								\
	public:																									\
		using Distributable::runDistributed;																\
																											\
		virtual void run(Image<in_t>& in, Image<result_t>& out, vector<ParamVariant>& args) const override	\
		{																									\
			size_t dim = pop<size_t>(args);																    \
			if(dim > 2)																			            \
				throw ITLException("Invalid dimension specification.");										\
			Image<val_t>& valImg = *pop<Image<val_t>*>(args);												\
			Image<val_t>& valOut = *pop<Image<val_t>*>(args);												\
																											\
			itl2:: funcname (in, valImg, dim, out, valOut);													\
		}																									\
																											\
		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override	\
		{																									\
			DistributedImage<in_t>& in = *std::get<DistributedImage<in_t>* >(args[0]);						\
			DistributedImage<result_t>& out = *std::get<DistributedImage<result_t>* >(args[1]);				\
			size_t dim = std::get<size_t>(args[2]);															\
			DistributedImage<val_t>& inVal = *std::get<DistributedImage<val_t>* >(args[3]);					\
			DistributedImage<val_t>& outVal = *std::get<DistributedImage<val_t>* >(args[4]);				\
																											\
			in.checkSize(inVal.dimensions());																\
																											\
			if (dim == 2)																					\
			{																								\
				/* Z project */																				\
				out.ensureSize(in.width(), in.height());													\
				outVal.ensureSize(in.width(), in.height());													\
			}																								\
			else if (dim == 1)																				\
			{																								\
				/* Y project */																				\
				out.ensureSize(in.width(), in.depth());														\
				outVal.ensureSize(in.width(), in.depth());													\
			}																								\
			else if (dim == 0)																				\
			{																								\
				/* X project. Swap z and y in output to make the image more logical (in some sense...) */	\
				out.ensureSize(in.depth(), in.height());													\
				outVal.ensureSize(in.depth(), in.height());													\
			}																								\
			else																							\
			{																								\
				throw ITLException("Invalid dimension specification.");										\
			}																								\
																											\
			return distributor.distribute(this, args)			;											\
		}																									\
																											\
		virtual size_t getDistributionDirection1(const vector<ParamVariant>& args) const override			\
		{																									\
			return 1; /* Always distribute along second dimension of output. */								\
		}																									\
																											\
		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override	\
		{																									\
			if (argIndex == 0 || argIndex == 3)																\
			{																								\
				DistributedImage<in_t>& in = *std::get<DistributedImage<in_t>* >(args[0]);					\
				size_t dim = std::get<size_t>(args[2]);														\
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
	DEF_PROJECT_ALL(Sum, sum, sum, float32_t)
	DEF_PROJECT(Sum, sum, sum, float32_t)
	DEF_PROJECT_ALL(SquareSum, squareSum, squaresum, float32_t)
	DEF_PROJECT(SquareSum, squareSum, squaresum, float32_t)

	DEF_PROJECT_ALL(Min, min, minval, in_t)
	DEF_PROJECT(Min, min, minval, in_t)
	DEF_PROJECT_ALL(Max, max, maxval, in_t)
	DEF_PROJECT(Max, max, maxval, in_t)
	
	DEF_PROJECT_2IMAGE(Min, min, minval, in_t)
	DEF_PROJECT_2IMAGE(Max, max, maxval, in_t)

	// Mean project all pixels is defined separately using sum project. This cuts down the
	// required number of commands.
	DEF_PROJECT(Mean, mean, mean, float32_t)


	template<typename in_t> class MaskedMeanAllPixelsCommand : public TwoImageInputOutputCommand<in_t, float32_t>, public Distributable	
	{																										
	protected:																								
		friend class CommandList;																			
																											
		MaskedMeanAllPixelsCommand() : TwoImageInputOutputCommand<in_t, float32_t>("maskedmean", "Calculates mean of all pixels in the input image, but skips specific value in the calculation. The output is a 1x1x1 image.",	
			{																								
				CommandArgument<double>(ParameterDirection::In, "ignored value", "This value is ignored in the calculation."),
				CommandArgument<bool>(ParameterDirection::In, "print to log", "Set to true to print the results to the log. This is a hack required because there is currently no good way to return single/simple values in distributed processing.", false),
				CommandArgument<bool>(ParameterDirection::In, "square root", "Set to true to take square root of the pixel values before calculating mean. This is a hack currently needed in distributed thickness map calculation. If set to true, the input image will contain its square root at output, except in distributed processing where it will retain its old value.", false)
			}) {}																							
	public:																									
		using Distributable::runDistributed;																
																											
		virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const override
		{	
			in_t ignoreValue = pixelRound<in_t>(pop<double>(args));
			bool print = pop<bool>(args);	
			bool doSqrt = pop<bool>(args);

			double count;
			if (doSqrt)
				squareRoot(in);

			double res = itl2::maskedMean<in_t, double>(in, ignoreValue, count);

			out.ensureSize(1, 1, 1);																		
			out(0) = (float32_t)res;																					
			
			if(print)																						
			{																								
				std::cout << "maskedmean = " << res << std::endl;												
				std::cout << "count = " << count << std::endl;												
			}																								
		}																									
																											
		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override	
		{																									
			DistributedImage<float32_t>& out = *std::get<DistributedImage<float32_t>* >(args[1]);
			in_t ignoreValue = pixelRound<in_t>(std::get<double>(args[2]));
			bool print = std::get<bool>(args[3]);																
			bool doSqrt = std::get<bool>(args[4]);
			out.ensureSize(1, 1, 1);																		
																											
			vector<ParamVariant> args2;																		
			args2.push_back(args[0]);																		
			args2.push_back(args[1]);																		
			args2.push_back(args[2]);
			ParamVariant p;																					
			p = true;																						
			args2.push_back(p);																				
			ParamVariant p2;
			p2 = doSqrt;
			args2.push_back(p2);
																											
			vector<string> results = distributor.distribute(this, args2);									
			vector<float32_t> vals;																			
			vector<float32_t> counts;
			for (const string& s : results)																	
			{
				float32_t val = internals::getValue<float32_t>(s, "maskedmean = ");
				float32_t count = internals::getValue<float32_t>(s, "count = ");
				if (count <= 0) // If count is zero, mean for the block is nan. Avoid that as other blocks might have non-nan means and single nan makes the final result nan.
					val = 0;
				vals.push_back(val);
				counts.push_back(count);
			}																								
																											
			Image<float32_t> tmp;
			out.readTo(tmp);																				
			tmp(0) = internals::maskedMeanReducer<float32_t>(vals, counts);
			out.setData(tmp);																				
																											
			if (print)																						
			{																								
				std::cout << "maskedmean = " << tmp(0) << std::endl;												
				std::cout << "count = " << sum<float32_t, float32_t>(counts) << std::endl;
			}																								
																											
			return vector<string>();																		
		}																									
																											
		virtual size_t getRefIndex(const vector<ParamVariant>& args) const override							
		{																									
			return 0;																						
		}																									
																											
		virtual void getCorrespondingBlock(const vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override	
		{																									
			if (argIndex == 1)																				
			{																								
				readStart = Vec3c(0, 0, 0);																	
				readSize = Vec3c(1, 1, 1);																	
				writeFilePos = Vec3c(0, 0, 0);																
				writeImPos = Vec3c(0, 0, 0);																
				writeSize = Vec3c(1, 1, 1);																	
			}																								
		}																									
	};																										



	// TODO: Variance is not available at the moment
	//template<typename in_t> class VarianceAllPixelsCommand : public TwoImageInputOutputCommand<in_t, float32_t>	
	//{	
	//protected:
	//	friend class CommandList;																						
	//	VarianceAllPixelsCommand() : TwoImageInputOutputCommand<in_t, float32_t>("variance", "Calculates variance of all pixels in the input image. The output is a 1x1x1 image.",
	//		{																								
	//			CommandArgument<bool>(ParameterDirection::In, "print to log", "Set to true to print the results to the log.", false)	
	//		}) {}																							
	//				
	//public:
	//	virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const			
	//	{																									
	//		bool print = pop<bool>(args);																	
	//		double res = (double)variance(in);																
	//		out.ensureSize(1, 1, 1);																		
	//		out(0) = (float32_t)res;																		
	//		if(print)																						
	//		{																								
	//			std::cout << "variance = " << out(0) << std::endl;												
	//			std::cout << "count = " << in.pixelCount() << std::endl;												
	//		}																								
	//	}																									
	//																										
	//};

	// TODO: Standard deviation projection is not available at the moment
	template<typename in_t> class StdDevAllPixelsCommand : public TwoImageInputOutputCommand<in_t, float32_t>, public Distributable
	{
	protected:
		friend class CommandList;

		StdDevAllPixelsCommand() : TwoImageInputOutputCommand<in_t, float32_t>("stddev", "Calculates standard deviation of all pixels in the input image. The output is a 1x1x1 image.",
			{
				CommandArgument<bool>(ParameterDirection::In, "print to log", "Set to true to print the results to the log.", false)
			}) {}

		static void print(bool doPrint, float32_t stddev, coord_t pixelCount)
		{
			if (doPrint)
			{
				std::cout << "stddev = " << stddev << std::endl;
				std::cout << "count = " << pixelCount << std::endl;
			}
		}

	public:
		virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const override
		{
			bool doPrint = pop<bool>(args);
			double res = (double)itl2::stdDev(in);
			out.ensureSize(1, 1, 1);
			out(0) = (float32_t)res;
			print(doPrint, out(0), in.pixelCount());
		}

		using Distributable::runDistributed;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<in_t>& in = *pop<DistributedImage<in_t>*>(args);
			DistributedImage<float32_t>& out = *pop<DistributedImage<float32_t>*>(args);
			double doPrint = pop<bool>(args);
			
			CommandList::get<SumAllPixelsCommand<in_t, float32_t>>().runDistributed(distributor, { &in, &out, false });
			double total = out.getValue();

			CommandList::get<SquareSumAllPixelsCommand<in_t, float32_t>>().runDistributed(distributor, {&in, &out, false});
			double total2 = out.getValue();

			double pixelCount = (double)in.pixelCount();
			double std = sqrt((total2 - (total * total / pixelCount)) / (pixelCount - 1));

			Image<float32_t> temp;
			temp(0) = (float32_t)std;
			out.setData(temp);

			print(doPrint, (float32_t)std, in.pixelCount());

			return std::vector<string>();
		}
	};


	template<typename in_t> class MeanAllPixelsCommand : public TwoImageInputOutputCommand<in_t, float32_t>, public Distributable
	{
	protected:
		friend class CommandList;

		MeanAllPixelsCommand() : TwoImageInputOutputCommand<in_t, float32_t>("mean", "Calculates mean of all pixels in the input image. The output is a 1x1x1 image.",
			{
				CommandArgument<bool>(ParameterDirection::In, "print to log", "Set to true to print the results to the log.", false)
			}) {}

		static void print(bool doPrint, float32_t stddev, coord_t pixelCount)
		{
			if (doPrint)
			{
				std::cout << "mean = " << stddev << std::endl;
				std::cout << "count = " << pixelCount << std::endl;
			}
		}

	public:
		virtual void run(Image<in_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const override
		{
			bool doPrint = pop<bool>(args);
			double res = (double)itl2::mean(in);
			out.ensureSize(1, 1, 1);
			out(0) = (float32_t)res;
			print(doPrint, out(0), in.pixelCount());
		}

		using Distributable::runDistributed;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<in_t>& in = *pop<DistributedImage<in_t>*>(args);
			DistributedImage<float32_t>& out = *pop<DistributedImage<float32_t>*>(args);
			double doPrint = pop<bool>(args);

			CommandList::get<SumAllPixelsCommand<in_t, float32_t>>().runDistributed(distributor, { &in, &out, false });
			double total = out.getValue();
			double pixelCount = (double)in.pixelCount();

			double mean = total / pixelCount;

			Image<float32_t> temp;
			temp(0) = (float32_t)mean;
			out.setData(temp);

			print(doPrint, (float32_t)mean, in.pixelCount());

			return std::vector<string>();
		}
	};
}
