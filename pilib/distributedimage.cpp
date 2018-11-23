
#include "distributedimage.h"
#include "io/raw.h"
#include "io/sequence.h"

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
	    //cout << "Set read source to " << filename << endl;
	    
		readSource = filename;

        if(filename != "")
        {
		    if (endsWith(filename, ".raw"))
		    {
			    isNewImage = !fileExists(filename);
		    }
		    else
		    {
			    isNewImage = !sequence::isSequence(filename);
		    }
		}
		else
		{
		    isNewImage = true;
		}
		
		//cout << "Is new image = " << isNewImage << endl;
	}
	
	string DistributedImageBase::emitReadBlock(const Vec3c& filePos, const Vec3c& blockSize, bool dataNeeded) const
	{
		string srcfile = currentReadSource();
		stringstream s;
		//cout << "is new = " << isNewImage << endl;
		if (!isNewImage && dataNeeded)
		{
			if(endsWith(srcfile, ".raw"))
				s << "readrawblock(" << name << ", " << srcfile << ", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ", " << dataTypeStr << ", " << dims.x << ", " << dims.y << ", " << dims.z << ");" << endl;
			else
				s << "readsequenceblock(" << name << ", " << srcfile << ", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		}
		else
		{
			//s << "newimage(" << name << ", " << dataTypeStr << ", " << dims.x << ", " << dims.y << ", " << dims.z << ");" << endl;
			s << "newimage(" << name << ", " << dataTypeStr << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		}
		//cout << "Read command: " << s.str() << endl;
		return s.str();
	}

	string DistributedImageBase::emitWriteBlock(const Vec3c& filePos, const Vec3c& imagePos, const Vec3c& blockSize, const string* outputFile)
	{
		string outFile = currentWriteTarget();

		// Override target file?
		if (outputFile)
			outFile = *outputFile;

		stringstream s;
		if(endsWith(outFile, ".raw"))
			s << "writerawblock(" << name << ", " << outFile << ", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << dims.x << ", " << dims.y << ", " << dims.z << ", " << imagePos.x << ", " << imagePos.y << ", " << imagePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		else
			s << "writesequenceblock(" << name << ", " << outFile << ", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << dims.x << ", " << dims.y << ", " << dims.z << ", " << imagePos.x << ", " << imagePos.y << ", " << imagePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
			
		//cout << "Write command: " << s.str() << endl;
		return s.str();
	}

}
