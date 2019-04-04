
#include "distributedimage.h"
#include "io/raw.h"
#include "io/sequence.h"
#include "io/io.h"

#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

namespace pilib
{
	DistributedImageBase::~DistributedImageBase()
	{
		// Remove temporary files if they were created.
		fs::remove_all(tempFilename1);
		fs::remove_all(tempFilename2);
	}

	void DistributedImageBase::createTempFilenames()
	{
		stringstream s1, s2;
		string path = "./tmp_images/";

		if (readSource == "" || endsWith(readSource, ".raw"))
		{
			s1 << path << name << "-1_" << dims.x << "x" << dims.y << "x" << dims.z << ".raw";
			s2 << path << name << "-2_" << dims.x << "x" << dims.y << "x" << dims.z << ".raw";
		}
		else
		{
			s1 << path << name << "-1/";
			s2 << path << name << "-2/";
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

	void DistributedImageBase::setReadSource(const string& filename)
	{
		readSource = filename;

        if(filename != "")
        {
			ImageDataType dt;
			isNewImage = io::getInfo(filename, dims, dt) == false;
			if (!isNewImage)
			{
				if (dt != dataType())
					throw ITLException("Invalid distributed image source file. Data type does not match data type of distributed image object.");
			}
		    //if (endsWith(filename, ".raw"))
		    //{
			   // isNewImage = !fileExists(filename);
		    //}
		    //else
		    //{
			   // isNewImage = !sequence::isSequence(filename);
		    //}
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
			s << "readblock(" << name << ", " << currentReadSource() << ", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
			//if(endsWith(srcfile, ".raw"))
			//	s << "readrawblock(" << name << ", " << srcfile << ", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ", " << dataTypeStr << ", " << dims.x << ", " << dims.y << ", " << dims.z << ");" << endl;
			//else
			//	s << "readsequenceblock(" << name << ", " << srcfile << ", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		}
		else
		{
			s << "newimage(" << name << ", " << dataTypeStr << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		}
		return s.str();
	}

	string DistributedImageBase::emitWriteBlock(const Vec3c& filePos, const Vec3c& imagePos, const Vec3c& blockSize)
	{
		stringstream s;
		if(isOutputRaw())
			s << "writerawblock(" << name << ", " << currentWriteTarget() << ", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << dims.x << ", " << dims.y << ", " << dims.z << ", " << imagePos.x << ", " << imagePos.y << ", " << imagePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		else
			s << "writesequenceblock(" << name << ", " << currentWriteTarget() << ", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << dims.x << ", " << dims.y << ", " << dims.z << ", " << imagePos.x << ", " << imagePos.y << ", " << imagePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
			
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
