#pragma once

#include <string>
#include <iostream>
#include <list>
#include <vector>

#include "math/aabox.h"

namespace itl2
{
	namespace io
	{
		struct DistributedImageProcess
		{
			AABoxc readBlock;
			AABoxc writeBlock;
		};

		/**
		Finds out if a chunk (given its bounding box) is 'safe' or not.
		Safe chunks can be written to without any synchronization or post-processing of the results.
		*/
		inline bool isChunkSafe(const AABoxc& box, const std::vector<DistributedImageProcess>& processes)
		{
			std::vector<size_t> readerIndices;
			std::vector<size_t> writerIndices;
			for (size_t n = 0; n < processes.size(); n++)
			{
				const auto& process = processes[n];
				if (process.readBlock.overlapsExclusive(box))
					readerIndices.push_back(n);
				if (process.writeBlock.overlapsExclusive(box))
					writerIndices.push_back(n);
			}

			if (writerIndices.empty())
			{
				// No writers, the chunk is never written to, so it is safe.
				return true;
			}
			if (writerIndices.size() > 1)
			{
				// Multiple writers, the chunk is unsafe as the writers can write simultaneously.
				return false;
			}
			// One writer.
			if (readerIndices.empty())
			{
				// No readers, one writer, the chunk is safe.
				return true;
			}
			if (readerIndices.size() > 1)
			{
				// Multiple readers, one writer, the chunk is not safe as it can be read from and written to simultaneously.
				return false;
			}
			// One reader, one writer.
			if (readerIndices[0] == writerIndices[0])
			{
				// Reader and writer are the same process.
				// The chunk is safe as the reader/writer process should control its possibly overlapping reads and writes internally.
				return true;
			}
			// Reader and writer are different processes.
			// The chunk is not safe as the reader and the writer might access the chunk simultaneously.
			return false;
		}
	}
}