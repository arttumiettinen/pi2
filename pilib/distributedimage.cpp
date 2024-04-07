
#include "distributedimage.h"
#include "io/raw.h"
#include "io/sequence.h"
#include "io/io.h"
#include "math/mathutils.h"
#include "distributor.h"
#include "utilities.h"

#include "filesystem.h"

using namespace itl2;
using namespace std;

namespace pilib
{
	template class DistributedImage<uint8_t>;
	template class DistributedImage<uint16_t>;
	template class DistributedImage<uint32_t>;
	template class DistributedImage<uint64_t>;
	template class DistributedImage<int8_t>;
	template class DistributedImage<int16_t>;
	template class DistributedImage<int32_t>;
	template class DistributedImage<int64_t>;
	template class DistributedImage<float32_t>;
	template class DistributedImage<complex32_t>;

	DistributedImageStorageType DistributedImageBase::suggestStorageType(const Distributor& distributor, const Vec3c& dimensions)
	{
		if (distributor.getUseNN5())
		{
			// This condition is not good. Images of sizes like 3x100000000x1 will be stored in NN5 and that
			// leads to bad performance.
			//if (dimensions.product() >= distributor.getChunkSize().product())
			// Let's try with a condition that requires all dimensions be 'chunkable'.
			Vec3c cs = distributor.getChunkSize();
			if(dimensions.x >= cs.x &&
				dimensions.y >= cs.y &&
				dimensions.z >= cs.z)
				return DistributedImageStorageType::NN5;
			else
				return DistributedImageStorageType::Raw;
		}
		else
		{
			return DistributedImageStorageType::Raw;
		}
	}
	
	DistributedImageBase::DistributedImageBase(Distributor& distributor, const string& name, const Vec3c& dimensions, ImageDataType dataType, const string& sourceFilename, DistributedImageStorageType storageType) :
		dims(dimensions),
		name(name),
		pixelDataType(dataType),
		distributor(&distributor),
		nn5ChunkSize(distributor.getChunkSize())
	{
		setReadSourceInternal(sourceFilename, false);
		createTempFilenames(storageType);
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

	void DistributedImageBase::newWriteTarget(DistributedImageStorageType storageType)
	{
		
		if (writeTargetType == storageType)
		{
			// No change in write target storage type.
			if (writeTarget == tempFilename1)
				writeTarget = tempFilename2;
			else
				writeTarget = tempFilename1;
		}
		else
		{
			// Write target type changed.
			// Create new temp files.
			fs::remove_all(tempFilename1);
			fs::remove_all(tempFilename2);
			// TODO: Should we keep or change the storage type here?
			writeTargetType = storageType;
			// This assigns writeTarget, too.
			createTempFilenames(storageType);
		}

		//// Temporary files are always NN5, so here we just select new temp file.
		//writeTargetType = DistributedImageStorageType::NN5;
		//if (writeTarget == tempFilename1)
		//	writeTarget = tempFilename2;
		//else
		//	writeTarget = tempFilename1;

		// If input is raw, output should be raw, too.
		// If input is sequence, output should be sequence, too.
		// If temp file names are not .raw, convert them to .raw
		// NOTE that if temps are not the same type than input file, the temps cannot be
		// the input file and we can get rid of them.
		//if ((isRaw() && !endsWith(tempFilename1, ".raw")) ||
		//	(!isRaw() && endsWith(tempFilename1, ".raw")))
		//{
		//	fs::remove_all(tempFilename1);
		//	fs::remove_all(tempFilename2);

		//	// This assigns writeTarget, too.
		//	createTempFilenames();
		//}
		//else
		//{
		//	// Just select new temp file.
		//	if (writeTarget == tempFilename1)
		//		writeTarget = tempFilename2;
		//	else
		//		writeTarget = tempFilename1;
		//}
	}

	void DistributedImageBase::createTempFilenames(DistributedImageStorageType storageType)
	{
		string path = "./tmp_images/";
		// Add randomness to the name so that two images with the same name do not get saved to same files.
		// This may happen if image is cleared and re-created in PISystem but Distributor still keeps references to the cleared image.
		// TODO: Name generated like this is not necessarily 100 % unique.
		uniqName = name + "_" + itl2::toString(randc(10000));

		//stringstream s1, s2;
		//if(!fs::exists(readSource) || isRaw())
		//{
		//	s1 << path << uniqueName() << "-1_" << dims.x << "x" << dims.y << "x" << dims.z << ".raw";
		//	s2 << path << uniqueName() << "-2_" << dims.x << "x" << dims.y << "x" << dims.z << ".raw";
		//}
		//else
		//{
		//	s1 << path << uniqueName() << "-1/";
		//	s2 << path << uniqueName() << "-2/";
		//}
		
		//// We will always use NN5 datasets as temporary files
		//s1 << path << uniqueName() << "-1";
		//s2 << path << uniqueName() << "-2";

		//this->tempFilename1 = s1.str();
		//this->tempFilename2 = s2.str();
		//this->writeTargetType = DistributedImageStorageType::NN5;

		if (storageType == DistributedImageStorageType::Raw)
		{
			this->tempFilename1 = concatDimensions(path + uniqueName() + "-1", dims);
			this->tempFilename2 = concatDimensions(path + uniqueName() + "-2", dims);
		}
		else
		{
			this->tempFilename1 = path + uniqueName() + "-1";
			this->tempFilename2 = path + uniqueName() + "-2";
		}

		this->writeTarget = this->tempFilename1;
		this->writeTargetType = storageType;
	}

	void DistributedImageBase::ensureSize(const Vec3c& newDimensions)
	{
		if (dims != newDimensions)
		{
			// Remove old temporary file so that a new one is created.
			fs::remove_all(tempFilename1);
			fs::remove_all(tempFilename2);
			dims = newDimensions;
			createTempFilenames(currentWriteTargetType());
		}
	}

	void DistributedImageBase::setReadSourceInternal(const string& filename, bool check)
	{
		readSource = filename;
		// Reset read source type to some default even if input file does not exist.
		readSourceType = DistributedImageStorageType::NN5;

        if(filename != "")
        {
			ImageDataType dt;
			Vec3c newDims;
			string reason;

			// Try NN5
			isNewImage = nn5::getInfo(filename, newDims, dt, reason) == false;
			if (!isNewImage)
			{
				if (check)
				{
					if (dt != dataType())
						throw ITLException("Invalid distributed image source file. Data type does not match data type of distributed image object.");
				}
				dims = newDims;
				readSourceType = DistributedImageStorageType::NN5;
			}
			else
			{
				// Try Raw
				isNewImage = raw::getInfo(filename, newDims, dt, reason) == false;
				if (!isNewImage)
				{
					if (check)
					{
						if (itl2::pixelSize(dt) != pixelSize())
							throw ITLException("Invalid distributed image source file. Data type does not match data type of distributed image object.");
					}
					dims = newDims;
					readSourceType = DistributedImageStorageType::Raw;
					// This ensures that the file name is the full one, not the one without dimensions.
					raw::internals::expandRawFilename(readSource);
				}
				else
				{
					// Try sequence
					isNewImage = sequence::getInfo(filename, newDims, dt, reason) == false;
					if (!isNewImage)
					{
						if (check)
						{
							if (dt != dataType())
								throw ITLException("Invalid distributed image source file. Data type does not match data type of distributed image object.");
						}
						dims = newDims;
						readSourceType = DistributedImageStorageType::Sequence;
					}
					else
					{
						// The input file is nothing supported. We'll assume it is a new file.
					}
				}
			}

			//ImageDataType dt;
			//Vec3c newDims;
			//string reason;
			//isNewImage = io::getInfo(filename, newDims, dt, reason) == false;
			//if (!isNewImage)
			//{
			//	dims = newDims;
			//	// NOTE: This test is not good as getInfo will recognize uint32 files as float32 files etc.
			//	//if (dt != dataType())
			//	//	throw ITLException("Invalid distributed image source file. Data type does not match data type of distributed image object.");
			//	if (check)
			//	{
			//		if (itl2::pixelSize(dt) != pixelSize())
			//			throw ITLException("Invalid distributed image source file. Data type does not match data type of distributed image object.");
			//	}
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
			s << "readblock(\"" << uniqueName() << "\", \"" << currentReadSource() << "\", " << filePos << ", " << blockSize << ", " << toString(pixelDataType) << ");" << endl;
		}
		else
		{
			s << "newimage(\"" << uniqueName() << "\", \"" << toString(pixelDataType) << "\", " << blockSize << ");" << endl;
		}
		return s.str();
	}

	string DistributedImageBase::emitWriteBlock(const Vec3c& filePos, const Vec3c& imagePos, const Vec3c& blockSize) const
	{
		stringstream s;
		//if(isOutputRaw())
		//	s << "writerawblock(\"" << uniqueName() << "\", \"" << currentWriteTarget() << "\", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << dims.x << ", " << dims.y << ", " << dims.z << ", " << imagePos.x << ", " << imagePos.y << ", " << imagePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		//else
		//	s << "writesequenceblock(\"" << uniqueName() << "\", \"" << currentWriteTarget() << "\", " << filePos.x << ", " << filePos.y << ", " << filePos.z << ", " << dims.x << ", " << dims.y << ", " << dims.z << ", " << imagePos.x << ", " << imagePos.y << ", " << imagePos.z << ", " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << ");" << endl;
		
		switch (currentWriteTargetType())
		{
		case DistributedImageStorageType::NN5:
		{
			// TODO: No NN5 compression specified.
			s << "writenn5block(\"" << uniqueName() << "\", \"" << currentWriteTarget() << "\", " << filePos << ", " << dims << ", " << imagePos << ", " << blockSize << ", " << nn5ChunkSize << ");" << endl;
			break;
		}
		case DistributedImageStorageType::Raw:
		{
			s << "writerawblock(\"" << uniqueName() << "\", \"" << currentWriteTarget() << "\", " << filePos << ", " << dims << ", " << imagePos << ", " << blockSize << ");" << endl;
			break;
		}
		case DistributedImageStorageType::Sequence:
		{
			s << "writesequenceblock(\"" << uniqueName() << "\", \"" << currentWriteTarget() << "\", " << filePos << ", " << dims << ", " << imagePos << ", " << blockSize << ");" << endl;
			break;
		}
		default: throw ITLException("Invalid write target type.");
		}

		return s.str();
	}

	size_t DistributedImageBase::startConcurrentWrite(const std::vector<nn5::NN5Process>& processes)
	{
		if (currentWriteTargetType() == DistributedImageStorageType::NN5)
		{
			// TODO: NN5 compression defaults to LZ4
			return nn5::startConcurrentWrite(dimensions(), dataType(), currentWriteTarget(), nn5ChunkSize, nn5::NN5Compression::LZ4, processes);
		}
		return 0;
	}

	string DistributedImageBase::emitEndConcurrentWrite(const Vec3c& chunk) const
	{
		if (currentWriteTargetType() == DistributedImageStorageType::NN5)
		{
			return string("endconcurrentwrite(\"") + currentWriteTarget() + "\", " + toString(chunk) + ")";
		}
		else
		{
			return "";
		}
	}

	vector<Vec3c> DistributedImageBase::getChunksThatNeedEndConcurrentWrite() const
	{
		if (currentWriteTargetType() == DistributedImageStorageType::NN5)
		{
			return nn5::getChunksThatNeedEndConcurrentWrite(currentWriteTarget());
		}
		return vector<Vec3c>();
	}

    void DistributedImageBase::writeComplete()
    {
		if (currentWriteTargetType() == DistributedImageStorageType::NN5)
		{
			nn5::endConcurrentWrite(currentWriteTarget(), false);
		}

        // Temporary image corresponding to old read source is not needed anymore as it
        // is not up to date (unless read source and write target are the same).
        if(currentReadSource() != currentWriteTarget())
        {
            if(currentReadSource() == this->tempFilename1)
                fs::remove_all(this->tempFilename1);
            if(currentReadSource() == this->tempFilename2)
                fs::remove_all(this->tempFilename2);
        }
        
	    setReadSourceInternal(currentWriteTarget(), true);
    }
}
