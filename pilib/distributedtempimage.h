#pragma once

#include "math/vec3.h"
#include "io/imagedatatype.h"

#include "distributedimage.h"
#include "commandlist.h"
#include "specialcommands.h"
#include "pisystem.h"

namespace pilib
{
	/**
	Temporary distributed image.
	The image is deleted when the object goes out of scope.
	This class is used if runDistributed method of a command does something complicated
	that requires temporary image(s).
	*/
	template<typename pixel_t> class DistributedTempImage
	{
	private:

		std::string name;

		DistributedImage<pixel_t>* dimg;

		Distributor& distributor;

	public:
		/**
		Constructor
		@param distributor Distributor object.
		@param purpose Free-form purpose string. This string will be shown in the name of the temporary file(s) corresponding to this image.
		@param dimensions Dimensions of the temp image.
		*/
		DistributedTempImage(Distributor& distributor, const std::string& purpose, const Vec3c& dimensions) :
			distributor(distributor),
			name(purpose + "_" + itl2::toString(randc(10000)))
		{
			std::string dts = itl2::toString(imageDataType<pixel_t>());

			CommandList::get<NewImageCommand>().runDistributed(distributor, { name, dts, dimensions.x, dimensions.y, dimensions.z });
			
			dimg = (DistributedImage<pixel_t>*)distributor.getSystem()->getDistributedImage(name);
		}

		DistributedTempImage(Distributor& distributor, const std::string& purpose, coord_t width = 1, coord_t height = 1, coord_t depth = 1) :
			DistributedTempImage(distributor, purpose, Vec3c(width, height, depth))
		{

		}

		/**
		Clears this temporary image from the pi2 system
		*/
		virtual ~DistributedTempImage()
		{
			CommandList::get<ClearCommand>().runDistributed(distributor, { name });
		}

		/**
		Gets the actual distributed image object.
		*/
		DistributedImage<pixel_t>& get()
		{
			return *dimg;
		}
	};

}