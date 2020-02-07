#pragma once

#include <vector>
//#include <atomic>

#include "ompatomic.h"

namespace itl2
{
	/**
	Disjoint-set forest that contains all sets in some range [0, max].
	*/
	class IndexForest
	{
	private:
		struct Element
		{
			//std::atomic<size_t> parent;
			OmpAtomic<size_t> parent;
			size_t rank;

			Element(size_t parent) :
				parent(parent), 
				rank(0)
			{
			}

			Element(const Element& other)
			{
				//parent.store(other.parent);
				parent.set(other.parent.get());
				rank = other.rank;
			}
		};

		std::vector<Element> items;

	public:

		/**
		Creates forest containing no sets.
		*/
		IndexForest()
		{

		}

		/**
		Creates forest that contains elements in range [0, maxValue[.
		*/
		IndexForest(size_t maxValue)
		{
			initialize(maxValue);
		}

		/**
		Creates forest that contains elements in range [0, maxValue[.
		Erases old forest, if any.
		*/
		void initialize(size_t maxValue)
		{
			items.clear();
			items.reserve(maxValue);
			for (size_t n = 0; n < maxValue; n++)
				items.push_back(Element(n));
			items.shrink_to_fit();
		}

		/**
		Finds set where x is contained.
		Multiple threads can call find_set concurrently, and one of the threads can
		call union_sets at the same time.
		*/
		size_t find_set(size_t x)
		{
			Element& element = items[x];
			//std::atomic<size_t>& parent = element.parent;
			OmpAtomic<size_t>& parent = element.parent;
			if (parent != x)
			{
				parent = find_set(parent);
			}
			return parent;
		}

		/**
		Merges sets containing keys x and y.
		Multiple threads can call find_set concurrently, and one of the threads can
		call union_sets at the same time.
		*/
		void union_sets(size_t x0, size_t y0)
		{
			size_t setX = find_set(x0);
			size_t setY = find_set(y0);
			if (setX != setY)
			{
				Element& elementX = items[setX];
				Element& elementY = items[setY];
				size_t& rankX = elementX.rank;
				size_t& rankY = elementY.rank;
				if (rankX > rankY)
				{
					elementY.parent = setX;
				}
				else
				{
					elementX.parent = setY;
					if (rankX == rankY)
						rankY++;
				}
			}
		}
	};

}