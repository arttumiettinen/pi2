
#include "distributedimage.h"
#include "io/raw.h"
#include "io/sequence.h"
#include "io/io.h"
#include "math/mathutils.h"
#include "distributor.h"
#include "utilities.h"

#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

using namespace itl2;
using namespace std;

namespace pilib
{
	
	DistributedImageBase::DistributedImageBase(Distributor& distributor, const string& name, const Vec3c& dimensions, ImageDataType dataType, const string& sourceFilename) :
		dims(dimensions),
		name(name),
		pixelDataType(dataType),
		distributor(&distributor)
	{
		setReadSource(sourceFilename, false);
		createTempFilenames();
	}

	DistributedImageBase::~DistributedImageBase()
	{
		// Remove temporary files if they were created.
		fs::remove_all(tempFilename1);
		fs::remove_all(tempFilename2);
	}

	void DistributedImageBase::flush() const
	{
		distributor->flush();
	}

	void DistributedImageBase::newWriteTarget()
	{
		// If input is raw, output should be raw, too.
		// If input is sequence, output should be sequence, too.
		// If temp file names are not .raw, convert them to .raw
		// NOTE that if temps are not the same type than input file, the temps cannot be
		// the input file and we can get rid of them.
		if ((isRaw() && !endsWith(tempFilename1, ".raw")) ||
			(!isRaw() && endsWith(tempFilename1, ".raw")))
		{
			fs::remove_all(tempFilename1);
			fs::remove_all(tempFilename2);

			// This assigns writeTarget, too.
			createTempFilenames();
		}
		else
		{
			// Just select new temp file.
			if (writeTarget == tempFilename1)
				writeTarget = tempFilename2;
			else
				writeTarget = tempFilename1;
		}
	}

	void DistributedImageBase::createTempFilenames()
	{
		stringstream s1, s2;
		string path = "./tmp_images/";
		// Add randomness to the name so that two images with the same name do not get saved to same files.
		// This may happen if image is cleared and re-created in PISystem but Distributor still keeps references to the cleared image.
		// TODO: Name generated like this is not necessarily 100 % unique.
		uniqName = name + "_" + itl2::toString(randc(10000));

		//if (readSource == "" || isRaw())
		if(!fileExists(readSource) || isRaw())
		{
			s1 << path << uniqueName() << "-1_" << dims.x << "x" << dims.y << "x" << dims.z << ".raw";
			s2 << path << uniqueName() << "-2_" << dims.x << "x" << dims.y << "x" << dims.z << ".raw";
		}
		else
		{
			s1 << path << uniqueName() << "-1/";
			s2 << path << uniqueName() << "-2/";
		}
		this->tempFilename1 = s1.str();
		this->tempFilename2 = s2.str();
		this->writeTarget = this->tempFilename1;
	}

	void DistributedImageBase::ensureSize(const Vec3c& newDimensions)
	{
		if (dims != newDimensions)
		{
			// Remove old temporary file so that a new one is created.
			fs::remove_all(tempFilename1);
			fs::remove_all(tempFilename2);
			dims = newDimensions;
			createTempFilenames();
		}
	}

	void DistributedImageBase::setReadSource(const string& filename, bool check)
	{
		readSource = filename;

        if(filename != "")
        {
			ImageDataType dt;
			Vec3c newDims;
			string reason;
			isNewImage = io::getInfo(filename, newDims, dt, reason) == false;
			if (!isNewImage)
			{
				dims = newDims;
				// NOTE: This test is not good as getInfo will recognize uint32 files as float32 files etc.
				//if (dt != dataType())
				//	throw ITLException("Invalid distributed image source file. Data type does not match data type of distributed image object.");
				if (check)
				{
					if (itl2::pixelSize(dt) != pixelSize())
						throw ITLException("Invalid distributed image source file. Data type does not match data type of distributed image object.");
				}
			}
		}
		else
		{
		    isNewImage = true;
		}
	}
	
	string DistributedImageBase::emitReadBlock(const Vec3c& filePos, const Vec3c& blockSize, bool dataNeeded) const
	{
		stringstream s;
		if (!isNewImage && dataNeeded)
		{
			s << "readblock(\"" << uniqueName() << "\", \"" << currentReadSource() << "\", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ", " << toString(pixelDataType) << ");" << endl;
		}
		else
		{
			s << "newimage(\"" << uniqueName() << "\", \"" << toString(pixelDataType) << "\", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		}
		return s.str();
	}

	string DistributedImageBase::emitWriteBlock(const Vec3c& filePos, const Vec3c& imagePos, const Vec3c& blockSize)
	{
		stringstream s;
		if(isOutputRaw())
			s << "writerawblock(\"" << uniqueName() << "\", \"" << currentWriteTarget() << "\", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << dims.x << ", " << dims.y << ", " << dims.z << ", " << imagePos.x << ", " << imagePos.y << ", " << imagePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		else
			s << "writesequenceblock(\"" << uniqueName() << "\", \"" << currentWriteTarget() << "\", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << dims.x << ", " << dims.y << ", " << dims.z << ", " << imagePos.x << ", " << imagePos.y << ", " << imagePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
			
		return s.str();
	}

    void DistributedImageBase::writeComplete()
    {
        // Temporary image corresponding to old read source is not needed anymore as it
        // is not up to date (unless read source and write target are the same).
        if(currentReadSource() != currentWriteTarget())
        {
            if(currentReadSource() == this->tempFilename1)
                fs::remove_all(this->tempFilename1);
            if(currentReadSource() == this->tempFilename2)
                fs::remove_all(this->tempFilename2);
        }
        
	    setReadSource(currentWriteTarget());
    }
}
