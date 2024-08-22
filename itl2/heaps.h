#pragma once

#include <map>
#include <deque>
#include <queue>


namespace itl2
{
	/**
	Stores elements in a map where priority is the key.
	Groups multiple elements of the same priority to the same map entry.
	Good if there are small number of priorities but large number of elements per priority.
	*/
	template<typename item_t, typename priority_t> class BucketMap
	{
	private:
		std::map<priority_t, std::deque<item_t> > items;

	//public:
	//	const priority_t find(const item_t& value, const priority_t& notFound) const
	//	{
	//		for (const auto& pair : items)
	//		{
	//			if (std::find(pair.second.begin(), pair.second.end(), value) != pair.second.end())
	//				return pair.first;
	//		}

	//		return notFound;
	//	}

	//	void remove(const item_t& value)
	//	{
	//		for (auto& pair : items)
	//		{
	//			auto& it = std::find(pair.second.begin(), pair.second.end(), value);
	//			if (it != pair.second.end())
	//				pair.second.erase(it);
	//		}
	//	}

	public:

		/**
		Pushes the given element value to the priority queue.
		*/
		void push(const priority_t priority, const item_t& value)
		{
			items[priority].push_back(value);
		}

		/**
		Returns reference to the item on the top of the priority queue.
		*/
		const item_t& topItem() const
		{
			return items.rbegin()->second.front();
		}

		/**
		Returns priority of the item on the top of the priority queue.
		*/
		const priority_t topPriority() const
		{
			return items.rbegin()->first;
		}

		/**
		Returns all item that have the top priority.
		*/
		const std::deque<item_t>& topItems() const
		{
			return items.rbegin()->second;
		}

		/**
		Removes all items that have the top priority from the priority queue.
		*/
		void popTopItems()
		{
			items.erase(items.rbegin()->first);
		}

		/**
		Removes the top element from the priority queue.
		*/
		void pop()
		{
			auto& l = items.rbegin()->second;
			l.pop_front();
			if (l.empty())
				items.erase(items.rbegin()->first);
		}

		/**
		Checks if the container has no elements.
		*/
		bool empty() const
		{
			return items.empty();
		}

		size_t size() const
		{
			size_t s = 0;
			for (const auto& item : items)
				s += item.second.size();
			return s;
		}
	};

	/**
	Grouped heap
	When multiple items are pushed with the same priority, they are mapped to single element in the underlying queue
	*/
	template<typename item_t, typename priority_t> class GHeap
	{
	public:

		struct Entry
		{
			priority_t priority;
			size_t birthday;
			std::deque<item_t> items;

			bool operator < (const Entry& right) const
			{
				// Large priorities first, small birthdays first
				if (!NumberUtils<priority_t>::equals(priority, right.priority))
					return NumberUtils<priority_t>::lessThan(priority, right.priority);
				else
					return birthday > right.birthday;
			}
		};

	private:

		std::priority_queue<Entry> queue;

		/**
		Next birthday value to use.
		*/
		size_t nextBirthday;

		Entry lastPushEntry;

	public:

		/**
		Constructor.
		*/
		GHeap() :
			nextBirthday(1)
		{
			lastPushEntry.priority = 0;
			lastPushEntry.birthday = 0;
		}

		/**
		Pushes the given element value to the priority queue.
		*/
		void push(const priority_t priority, const item_t& value)
		{
			if (lastPushEntry.priority != priority)
			{
				// Push lastPushItems to the queue and make new 
				if (lastPushEntry.items.size() > 0)
				{
					// NOTE: This reveals that addition pattern is many times A, B, A, B, A, B... and so this implementation fails!
					//cout << "Pushing " << lastPushEntry.items.size() << " values with priority " << lastPushEntry.priority << endl;
					queue.push(lastPushEntry);
				}
				lastPushEntry.priority = priority;
				lastPushEntry.birthday = nextBirthday;
				lastPushEntry.items.clear();
				nextBirthday++;
			}

			lastPushEntry.items.push_back(value);
		}


		/**
		Returns reference to the top element in the priority queue. This element will be removed on a call to pop().
		Calling top on empty queue is undefined behaviour.
		*/
		const Entry& top() const
		{
			if (queue.empty() || lastPushEntry.priority >= queue.top().priority)
				return lastPushEntry;

			return queue.top();
		}

		/**
		Returns reference to the item on the top of the priority queue.
		*/
		const item_t& topItem() const
		{
			return top().items.front();
		}

		/**
		Returns priority of the item on the top of the priority queue.
		*/
		const priority_t topPriority() const
		{
			return top().priority;
		}

		/**
		Removes the top element from the priority queue.
		*/
		void pop()
		{
			if (queue.empty() || lastPushEntry.priority >= queue.top().priority)
			{
				lastPushEntry.items.pop_front();
				if (lastPushEntry.items.size() <= 0 && !queue.empty())
				{
					lastPushEntry = queue.top();
					queue.pop();
				}
			}
			else
			{
				const_cast<Entry&>(queue.top()).items.pop_front();
				if (queue.top().items.size() <= 0)
					queue.pop();
			}
		}

		/**
		Checks if the container has no elements.
		*/
		bool empty() const
		{
			return queue.empty();
		}


	};

	/**
	Hierarchical heap.
	Divides expected priority range into number of bins, and each bin is a separate (smaller) priority queue.
	*/
	template<typename item_t, typename priority_t> class HHeap
	{
	public:

		struct Entry
		{
			priority_t priority;
			size_t birthday;
			item_t item;

			bool operator < (const Entry& right) const
			{
				// Large priorities first, small birthdays first
				if (!NumberUtils<priority_t>::equals(priority, right.priority))
					return NumberUtils<priority_t>::lessThan(priority, right.priority);
				else
					return birthday > right.birthday;
			}
		};


	private:

		std::vector<std::priority_queue<Entry> > queues;

		/**
		Expected minimum and maximum priority values.
		*/
		priority_t minp, maxp;

		/**
		Next birthday value to use.
		*/
		size_t nextBirthday;

		/**
		Index of the highest-priority non-empty queue.
		*/
		size_t topQueueIndex;

		/**
		Calculates index of bucket for given priority.
		*/
		size_t bucket_index(priority_t priority) const
		{
			coord_t i = (coord_t)floor(((double)priority - (double)minp) / ((double)maxp - (double)minp) * queues.size());
			if (i < 0)
				i = 0;
			else if (intuitive::gt(i, queues.size() - 1))
				i = (coord_t)queues.size() - 1;
			return i;
		}

	public:
		/**
		Constructor.
		@param minp Minimum expected priority value.
		@param maxp Maximum expected priority value.
		@param bucketCount Count of buckets.
		*/
		HHeap(priority_t minp, priority_t maxp, size_t bucketCount = 256) :
			minp(minp),
			maxp(maxp),
			nextBirthday(0),
			topQueueIndex(0)
		{
			queues.reserve(bucketCount);
			for (size_t n = 0; n < bucketCount; n++)
				queues.push_back(priority_queue<Entry>());
		}

		/**
		Pushes the given element value to the priority queue.
		*/
		void push(const priority_t priority, const item_t& value)
		{
			size_t bi = bucket_index(priority);
			queues[bi].push({ priority, nextBirthday, value });
			nextBirthday++;
			if (bi > topQueueIndex)
				topQueueIndex = bi;
		}


		/**
		Returns reference to the top element in the priority queue. This element will be removed on a call to pop().
		Calling top on empty queue is undefined behaviour.
		*/
		const Entry& top() const
		{
			return queues[topQueueIndex].top();
		}

		/**
		Returns reference to the item on the top of the priority queue.
		*/
		const item_t& topItem() const
		{
			return queues[topQueueIndex].top().item;
		}

		/**
		Returns priority of the item on the top of the priority queue.
		*/
		priority_t topPriority() const
		{
			return queues[topQueueIndex].top().priority;
		}

		/**
		Returns birthday of the item on the top of the priority queue.
		*/
		size_t topBirthday() const
		{
			return queues[topQueueIndex].top().birthday;
		}

		/**
		Removes the top element from the priority queue.
		*/
		void pop()
		{
			queues[topQueueIndex].pop();
			while (topQueueIndex > 0 && queues[topQueueIndex].empty())
				topQueueIndex--;
		}

		/**
		Checks if the container has no elements.
		*/
		bool empty() const
		{
			return queues[topQueueIndex].empty();
		}

		/**
		Calculates measure of non-uniformity of item counts in buckets.
		Zero corresponds to totally uniform item counts.
		*/
		double nonUniformity() const
		{
			std::vector<double> counts;
			counts.reserve(queues.size());
			for (const auto& q : queues)
				counts.push_back((double)q.size());
			return stddev(counts);
		}

		/**
		Calculates total count of elements in the heap.
		Complexity is linear in bucket count.
		*/
		size_t size() const
		{
			size_t count = 0;
			for (const auto& q : queues)
				count += q.size();
			return count;
		}
	};

	namespace tests
	{
		inline void hheap();
	}
}