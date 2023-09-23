#include <cstdint>
#include <iostream>

#ifndef AP_H
#define AP_H

template <size_t T> class ap_int {
public:
	ap_int<T> ()
	{
	}
	ap_int<T> (int64_t value)
	{
		this->value = value;
	}
	template<size_t T2>
	ap_int<T> (const ap_int<T2> & other)
	{
		// this->value = other.value;
		// return;

		if (T == 128)
			this->value = other.value;
		else if (other.value < 0)
			this->value = other.value | (~((__int128(1UL) << T) - 1));
		else
			this->value = other.value & ((__int128(1UL) << T) - 1);
	}
	~ap_int (){}

	template <size_t T2>
	ap_int<T> & operator= (const ap_int<T2> & other)
	{
		// this->value = other.value;
		// return *this;

		if (T == 128)
			this->value = other.value;
		else if (other.value < 0)
			this->value = other.value | (~((__int128(1UL) << T) - 1));
		else
			this->value = other.value & ((__int128(1UL) << T) - 1);
		return *this;
	}

	ap_int<T> & operator= (const int64_t& value)
	{
		// this->value = value;
		// return *this;

		if (T == 128)
			this->value = value;
		else if (value < 0)
			this->value = value | (~((__int128(1UL) << T) - 1));
		else
			this->value = value & ((__int128(1UL) << T) - 1);
		return *this;
	}

	template <size_t T2>
	ap_int<T> & operator+= (const ap_int<T2>& other)
	{
		// this->value = (this->value + other.value);
		// return *this;

		if (T == 128)
			this->value = (this->value + other.value);
		else if (other.value < 0)
			this->value = (this->value + other.value) | (~((__int128(1UL) << T) - 1));
		else
			this->value = (this->value + other.value) & ((__int128(1UL) << T) - 1);
		return *this;
	}
	ap_int<T> operator-()
	{
		ap_int<T> ret;
		ret.value = this->value * -1;
		return ret;
	}
	template <size_t T2>
	ap_int<T> operator+(const ap_int<T2> & other)
	{
		return (this->value + other.value);
	}
	template <size_t T2>
	ap_int<T> operator-(ap_int<T2> const & other)
	{
		return (this->value - other.value);
	}
	ap_int<T> operator<<(int const & other)
	{
		return (this->value << other);
	}
	ap_int<T> operator>>(int const & other)
	{
		return (this->value >> other);
	}
	template <size_t T2>
	ap_int<T> operator*(ap_int<T2> const & other)
	{
		return (this->value * other.value);
	}
	template <size_t T2>
	ap_int<T> operator/(ap_int<T2> const & other)
	{
		return (this->value * other.value);
	}
	// friend std::ostream& operator<<(std::ostream& out, const ap_int& v)
	// {
	// 	int64_t one = (v.value >> 64);
	// 	int64_t fractional = (v.value - (integer << TF)) / (1UL << TF);
	// 	out << integer << "." << fractional;
	// 	return out;
	//
	// 	return out;
	// }

	float to_fixed (int frac_bits)
	{
		return this->value / float (1UL<<frac_bits);
	}

	__int128 value;
private:
};

#endif // AP_H

