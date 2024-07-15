#pragma once
#include "RationalPolynom.hpp"

namespace twmath{
	template<class T, size_t num_size, size_t den_size>
	struct DiscreteTransferFunction{
		RationalPolynom<T, num_size, den_size> rp;
		T sample_time;
		
		constexpr DiscreteTransferFunction(const DiscreteTransferFunction&) = default;
		constexpr DiscreteTransferFunction(const RationalPolynom<T, num_size, den_size>& other, const T& sample_time) : rp(other), sample_time(sample_time){}
		constexpr DiscreteTransferFunction(const Polynom<T, num_size>& num, const Polynom<T, den_size>& den, const T& sample_time) : rp(num, den), sample_time(sample_time){}
		constexpr DiscreteTransferFunction(const T (&num)[num_size], const T (&den)[den_size], const T& sample_time) : rp(num, den), sample_time(sample_time){}
		constexpr DiscreteTransferFunction(const T& num, const T (&den)[den_size], const T& sample_time) : rp(num, den), sample_time(sample_time){}
		constexpr DiscreteTransferFunction(const T (&num)[num_size], const T& den, const T& sample_time) : rp(num, den), sample_time(sample_time){}
		
		constexpr DiscreteTransferFunction& operator=(const DiscreteTransferFunction&) = default;
		constexpr DiscreteTransferFunction& operator=(const RationalPolynom<T, num_size, den_size>& other){this->rp = other; return *this;}
		constexpr DiscreteTransferFunction& operator=(const Polynom<T, num_size>& num){this->rp = num; return *this;}
		constexpr DiscreteTransferFunction& operator=(const T (&num)[num_size]){this->rp = num; return *this;}
		constexpr DiscreteTransferFunction& operator=(const T& num){this->rp = num; return *this;}
	};

	/*Add errors for calculations with different sample times*/
	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator + (const DiscreteTransferFunction<Tl,lnum,lden>& l, const DiscreteTransferFunction<Tr, rnum, rden>& r){return DiscreteTransferFunction(l.rp + r.rp);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator + (const DiscreteTransferFunction<Tl,lnum,lden>& l, const Tr& r){return DiscreteTransferFunction(l.rp + r);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator + (const Tl& l, const DiscreteTransferFunction<Tr, rnum, rden>& r){return DiscreteTransferFunction(l + r.rp);}

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator - (const DiscreteTransferFunction<Tl,lnum,lden>& l, const DiscreteTransferFunction<Tr, rnum, rden>& r){return DiscreteTransferFunction(l.rp - r.rp);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator - (const DiscreteTransferFunction<Tl,lnum,lden>& l, const Tr& r){return DiscreteTransferFunction(l.rp - r);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator - (const Tl& l, const DiscreteTransferFunction<Tr, rnum, rden>& r){return DiscreteTransferFunction(l - r.rp);}

	template<class Tl, size_t lnum, size_t lden>
	constexpr auto operator - (const DiscreteTransferFunction<Tl,lnum,lden>& p){return DiscreteTransferFunction(-p.rp);}

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator * (const DiscreteTransferFunction<Tl,lnum,lden>& l, const DiscreteTransferFunction<Tr, rnum, rden>& r){return DiscreteTransferFunction(l.rp * r.rp);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator * (const DiscreteTransferFunction<Tl,lnum,lden>& l, const Tr& r){return DiscreteTransferFunction(l.rp * r);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator * (const Tl& l, const DiscreteTransferFunction<Tr, rnum, rden>& r){return DiscreteTransferFunction(l * r.rp);}

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator / (const DiscreteTransferFunction<Tl,lnum,lden>& l, const DiscreteTransferFunction<Tr, rnum, rden>& r){return DiscreteTransferFunction(l.rp / r.rp);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator / (const DiscreteTransferFunction<Tl,lnum,lden>& l, const Tr& r){return DiscreteTransferFunction(l.rp / r);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator / (const Tl& l, const DiscreteTransferFunction<Tr, rnum, rden>& r){return DiscreteTransferFunction(l / r.rp);}

	template<class Stream, class T, size_t num_size, size_t den_size>
	Stream& print_pretty(Stream& stream, const DiscreteTransferFunction<T, num_size, den_size>& dtf, const char* name=nullptr, const char* indent=""){
		return print_pretty(stream, dtf.rp, name, "z", indent);
	}
}