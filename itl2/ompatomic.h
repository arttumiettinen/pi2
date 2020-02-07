#pragma once

namespace itl2
{
	/**
	Provides atomic read and write access to a single value using OpenMP.
	*/
	template<typename T> class OmpAtomic
	{
	private:
		volatile T value;

	public:

		OmpAtomic() : value(T())
		{
		}

		OmpAtomic(const T& value) : value(value)
		{
		}

		OmpAtomic(const OmpAtomic&) = delete;
		OmpAtomic& operator=(OmpAtomic const&) = delete;

		T get() const
		{
			T local;

#if defined(_MSC_VER)
			local = value;
#else
			#pragma omp atomic read
			local = value;
#endif

			return local;
		}

		void set(const T& val)
		{
#if defined(_MSC_VER)
			value = val;
#else
			#pragma omp atomic write
			value = val;
#endif
		}

		operator T() const
		{
			return get();
		}

		void operator=(const T& val)
		{
			set(val);
		}
	};

}