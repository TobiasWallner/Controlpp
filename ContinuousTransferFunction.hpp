#pragma once
#include "RationalPolynom.hpp"

namespace twmath{
	template<class T, size_t num_size, size_t den_size>
	struct ContinuousTransferFunction{
		RationalPolynom<T, num_size, den_size> rp;
		
		constexpr ContinuousTransferFunction(const ContinuousTransferFunction&) = default;
		constexpr ContinuousTransferFunction(const RationalPolynom<T, num_size, den_size>& other) : rp(other){}
		constexpr ContinuousTransferFunction(const Polynom<T, num_size>& num, const Polynom<T, den_size>& den) : rp(num, den){}
		constexpr ContinuousTransferFunction(const T (&num)[num_size], const T (&den)[den_size]) : rp(num, den){}
		constexpr ContinuousTransferFunction(const T& num, const T (&den)[den_size]) : rp(num, den){}
		constexpr ContinuousTransferFunction(const T (&num)[num_size], const T& den) : rp(num, den){}
		
		constexpr ContinuousTransferFunction& operator=(const ContinuousTransferFunction&) = default;
		constexpr ContinuousTransferFunction& operator=(const RationalPolynom<T, num_size, den_size>& other){this->rp = other; return *this;}
		constexpr ContinuousTransferFunction& operator=(const Polynom<T, num_size>& num){this->rp = num; return *this;}
		constexpr ContinuousTransferFunction& operator=(const T (&num)[num_size]){this->rp = num; return *this;}
		constexpr ContinuousTransferFunction& operator=(const T& num){this->rp = num; return *this;}
	};


	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator + (const ContinuousTransferFunction<Tl,lnum,lden>& l, const ContinuousTransferFunction<Tr, rnum, rden>& r){return ContinuousTransferFunction(l.rp + r.rp);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator + (const ContinuousTransferFunction<Tl,lnum,lden>& l, const Tr& r){return ContinuousTransferFunction(l.rp + r);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator + (const Tl& l, const ContinuousTransferFunction<Tr, rnum, rden>& r){return ContinuousTransferFunction(l + r.rp);}

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator - (const ContinuousTransferFunction<Tl,lnum,lden>& l, const ContinuousTransferFunction<Tr, rnum, rden>& r){return ContinuousTransferFunction(l.rp - r.rp);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator - (const ContinuousTransferFunction<Tl,lnum,lden>& l, const Tr& r){return ContinuousTransferFunction(l.rp - r);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator - (const Tl& l, const ContinuousTransferFunction<Tr, rnum, rden>& r){return ContinuousTransferFunction(l - r.rp);}

	template<class Tl, size_t lnum, size_t lden>
	constexpr auto operator - (const ContinuousTransferFunction<Tl,lnum,lden>& p){return ContinuousTransferFunction(-p.rp);}

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator * (const ContinuousTransferFunction<Tl,lnum,lden>& l, const ContinuousTransferFunction<Tr, rnum, rden>& r){return ContinuousTransferFunction(l.rp * r.rp);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator * (const ContinuousTransferFunction<Tl,lnum,lden>& l, const Tr& r){return ContinuousTransferFunction(l.rp * r);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator * (const Tl& l, const ContinuousTransferFunction<Tr, rnum, rden>& r){return ContinuousTransferFunction(l * r.rp);}

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator / (const ContinuousTransferFunction<Tl,lnum,lden>& l, const ContinuousTransferFunction<Tr, rnum, rden>& r){return ContinuousTransferFunction(l.rp / r.rp);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator / (const ContinuousTransferFunction<Tl,lnum,lden>& l, const Tr& r){return ContinuousTransferFunction(l.rp / r);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator / (const Tl& l, const ContinuousTransferFunction<Tr, rnum, rden>& r){return ContinuousTransferFunction(l / r.rp);}

	template<class Stream, class T, size_t num_size, size_t den_size>
	Stream& print_pretty(Stream& stream, const ContinuousTransferFunction<T, num_size, den_size>& ctf, const char* name=nullptr, const char* indent=""){
		return print_pretty(stream, ctf.rp, name, "s", indent);
	}
}