#pragma once

#include "command.h"
#include "commandsbase.h"
#include "sfcm.h"


namespace pilib
{

	inline std::string sFCMSeeAlso()
	{
		return "-, -";
	}


	inline std::string sFCMMethodsHelp()
	{
		return
			"**sFCM**\n"
			"\n"
			"x's spatial version of Fuzzy C-Means. Spatial neighbourhood size is determined by nbApothem and _singular pixel_ to _neighbourhood_ weighting by p and q respectively.\n"
			"\n"
			"\n"
			"**FCM**\n"
			"\n"
			"y's Fuzzy C-Means. When q = 0.\n"
			"\n"
			"\n"
			"**sHCM**\n"
			"\n"
			"Same as sFCM, but with fuzzyParam = 1. 'spatial Hard C-means'.\n"
			"\n"
			"\n"
			"**K-Means**\n"
			"\n"
			"Traditional K-means when fuzzyParam = 1 and q = 0. 'Hard C-Means'.\n"
			"\n"
			"\n"
			"\n";
	}


	template<typename pixel_t> class sFCMCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		sFCMCommand() : OneImageInPlaceCommand<pixel_t>("sfcm",
			"Segments given image using sFCm, FCM, sHCM, or K-means.\n\n" +
			sFCMMethodsHelp(),
			{
				CommandArgument<coord_t>(ParameterDirection::In, "clusterCount", "How many classes to segment the image to."),
				CommandArgument<double>(ParameterDirection::In, "fuzzyParam", "Fuzziness parameter of FCM."),
				CommandArgument<coord_t>(ParameterDirection::In, "nbApothem", "Apothem/radius of how many pixels to include in spatial neighbourhood calculations."),
				CommandArgument<coord_t>(ParameterDirection::In, "iterMax", "Maxmimum number of iterations to run the clustering algorithm."),
				CommandArgument<double>(ParameterDirection::In, "stopParam", "Convergence tolerance threshold value."),
				CommandArgument<double>(ParameterDirection::In, "p", "Weight parameter of pixel value."),
				CommandArgument<double>(ParameterDirection::In, "q", "Weight parameter of spatial part.")
			})
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			int clusterCount = static_cast<int>(pop<coord_t>(args));
			float fuzzyParam = static_cast<float>(pop<double>(args));
			coord_t nbApothem = pop<coord_t>(args);
			int iterMax = static_cast<int>(pop<coord_t>(args));
			float stopParam = static_cast<float>(pop<double>(args));
			float p = static_cast<float>(pop<double>(args));
			float q = static_cast<float>(pop<double>(args));

			doSFCM(img, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q);

		}


	};


	template<typename pixel_t> class sFCMCommand2 : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		sFCMCommand2() : OneImageInPlaceCommand<pixel_t>("sfcm",
			"Segments given images using sFCm, FCM, sHCM, or K-means.\n\n" +
			sFCMMethodsHelp(),
			{
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "imageData2", "Image with feature 2."),
				CommandArgument<coord_t>(ParameterDirection::In, "clusterCount", "How many classes to segment the image to."),
				CommandArgument<double>(ParameterDirection::In, "fuzzyParam", "Fuzziness parameter of FCM."),
				CommandArgument<coord_t>(ParameterDirection::In, "nbApothem", "Apothem/radius of how many pixels to include in spatial neighbourhood calculations."),
				CommandArgument<coord_t>(ParameterDirection::In, "iterMax", "Maxmimum number of iterations to run the clustering algorithm."),
				CommandArgument<double>(ParameterDirection::In, "stopParam", "Convergence tolerance threshold value."),
				CommandArgument<double>(ParameterDirection::In, "p", "Weight parameter of pixel value."),
				CommandArgument<double>(ParameterDirection::In, "q", "Weight parameter of spatial part.")
			})
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& imageData2 = *pop<Image<pixel_t>* >(args);
			int clusterCount = static_cast<int>(pop<coord_t>(args));
			float fuzzyParam = static_cast<float>(pop<double>(args));
			coord_t nbApothem = pop<coord_t>(args);
			int iterMax = static_cast<int>(pop<coord_t>(args));
			float stopParam = static_cast<float>(pop<double>(args));
			float p = static_cast<float>(pop<double>(args));
			float q = static_cast<float>(pop<double>(args));

			doSFCM(img, imageData2, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q);
			//doSFCM(img, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q, &externalMembershipMatrix, &externalSpatialMembershipMatrix);

		}


	};


	template<typename pixel_t> class sFCMCommand3 : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		sFCMCommand3() : OneImageInPlaceCommand<pixel_t>("sfcm",
			"Segments given images using sFCm, FCM, sHCM, or K-means.\n\n" +
			sFCMMethodsHelp(),
			{
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "imageData2", "Image with feature 2."),
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "imageData3", "Image with feature 3."),
				CommandArgument<coord_t>(ParameterDirection::In, "clusterCount", "How many classes to segment the image to."),
				CommandArgument<double>(ParameterDirection::In, "fuzzyParam", "Fuzziness parameter of FCM."),
				CommandArgument<coord_t>(ParameterDirection::In, "nbApothem", "Apothem/radius of how many pixels to include in spatial neighbourhood calculations."),
				CommandArgument<coord_t>(ParameterDirection::In, "iterMax", "Maxmimum number of iterations to run the clustering algorithm."),
				CommandArgument<double>(ParameterDirection::In, "stopParam", "Convergence tolerance threshold value."),
				CommandArgument<double>(ParameterDirection::In, "p", "Weight parameter of pixel value."),
				CommandArgument<double>(ParameterDirection::In, "q", "Weight parameter of spatial part.")
				//CommandArgument<std::vector<std::unique_ptr<Image<float>>>*>(ParameterDirection::In, "externalMembershipMatrix", "membershipMatrix is saved to this if given."),
				//CommandArgument<std::vector<std::unique_ptr<Image<float>>>*>(ParameterDirection::In, "externalSpatialMembershipMatrix", "spatialMembershipMatrix is saved to this if given."),
			})
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& imageData2 = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& imageData3 = *pop<Image<pixel_t>* >(args);
			int clusterCount = static_cast<int>(pop<coord_t>(args));
			float fuzzyParam = static_cast<float>(pop<double>(args));
			coord_t nbApothem = pop<coord_t>(args);
			int iterMax = static_cast<int>(pop<coord_t>(args));
			float stopParam = static_cast<float>(pop<double>(args));
			float p = static_cast<float>(pop<double>(args));
			float q = static_cast<float>(pop<double>(args));
			//std::vector<std::unique_ptr<Image<float>>>* externalMembershipMatrix = pop<std::vector<std::unique_ptr<Image<float>>>*>(args);
			//std::vector<std::unique_ptr<Image<float>>>* externalSpatialMembershipMatrix = pop<std::vector<std::unique_ptr<Image<float>>>*>(args);

			doSFCM(img, imageData2, imageData3, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q);
			//doSFCM(img, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q, &externalMembershipMatrix, &externalSpatialMembershipMatrix);

		}


	};


	template<typename pixel_t> class sFCMCommand4 : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		sFCMCommand4() : OneImageInPlaceCommand<pixel_t>("sfcm",
			"Segments given images using sFCm, FCM, sHCM, or K-means.\n\n" +
			sFCMMethodsHelp(),
			{
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "imageData2", "Image with feature 2."),
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "imageData3", "Image with feature 3."),
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "imageData4", "Image with feature 4."),
				CommandArgument<coord_t>(ParameterDirection::In, "clusterCount", "How many classes to segment the image to."),
				CommandArgument<double>(ParameterDirection::In, "fuzzyParam", "Fuzziness parameter of FCM."),
				CommandArgument<coord_t>(ParameterDirection::In, "nbApothem", "Apothem/radius of how many pixels to include in spatial neighbourhood calculations."),
				CommandArgument<coord_t>(ParameterDirection::In, "iterMax", "Maxmimum number of iterations to run the clustering algorithm."),
				CommandArgument<double>(ParameterDirection::In, "stopParam", "Convergence tolerance threshold value."),
				CommandArgument<double>(ParameterDirection::In, "p", "Weight parameter of pixel value."),
				CommandArgument<double>(ParameterDirection::In, "q", "Weight parameter of spatial part.")
				//CommandArgument<std::vector<std::unique_ptr<Image<float>>>*>(ParameterDirection::In, "externalMembershipMatrix", "membershipMatrix is saved to this if given."),
				//CommandArgument<std::vector<std::unique_ptr<Image<float>>>*>(ParameterDirection::In, "externalSpatialMembershipMatrix", "spatialMembershipMatrix is saved to this if given."),
			})
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& imageData2 = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& imageData3 = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& imageData4 = *pop<Image<pixel_t>* >(args);
			int clusterCount = static_cast<int>(pop<coord_t>(args));
			float fuzzyParam = static_cast<float>(pop<double>(args));
			coord_t nbApothem = pop<coord_t>(args);
			int iterMax = static_cast<int>(pop<coord_t>(args));
			float stopParam = static_cast<float>(pop<double>(args));
			float p = static_cast<float>(pop<double>(args));
			float q = static_cast<float>(pop<double>(args));
			//std::vector<std::unique_ptr<Image<float>>>* externalMembershipMatrix = pop<std::vector<std::unique_ptr<Image<float>>>*>(args);
			//std::vector<std::unique_ptr<Image<float>>>* externalSpatialMembershipMatrix = pop<std::vector<std::unique_ptr<Image<float>>>*>(args);

			doSFCM(img, imageData2, imageData3, imageData4, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q);
			//doSFCM(img, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q, &externalMembershipMatrix, &externalSpatialMembershipMatrix);

		}


	};


}