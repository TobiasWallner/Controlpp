#pragma once

#include "Polynom.hpp"

namespace twmath{
	template<class T, size_t num_size=1, size_t den_size=1>
	struct RationalPolynom{
		Polynom<T, num_size> num;
		Polynom<T, den_size> den;
		
		constexpr RationalPolynom(const RationalPolynom&) = default;
		constexpr RationalPolynom(const Polynom<T, num_size>& num, const Polynom<T, den_size>& den) : num(num), den(den){}
		constexpr RationalPolynom(const T (&num)[num_size], const T (&den)[den_size]) : num(num), den(den){}
		constexpr RationalPolynom(const T& num, const T (&den)[den_size]) : num(num), den(den){}
		constexpr RationalPolynom(const T (&num)[num_size], const T& den) : num(num), den(den){}
		
		constexpr RationalPolynom& operator=(const RationalPolynom&) = default;
		constexpr RationalPolynom& operator=(const Polynom<T, num_size>& num){this->num = num; this->den = 1; return *this;}
		constexpr RationalPolynom& operator=(const T (&num)[num_size]){this->num = num; this->den = 1; return *this;}
		constexpr RationalPolynom& operator=(const T& num){this->num = num; this->den = 1; return *this;}
		
	};

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator + (const RationalPolynom<Tl,lnum,lden>& l, const RationalPolynom<Tr, rnum, rden>& r){return RationalPolynom(l.num * r.den + l.den * r.num, l.den * r.den);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator + (const RationalPolynom<Tl,lnum,lden>& l, const Tr& r){return RationalPolynom(l.num + l.den * r, l.den);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator + (const Tl& l, const RationalPolynom<Tr, rnum, rden>& r){return RationalPolynom(r.den * l + r.num, r.den);}

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator - (const RationalPolynom<Tl,lnum,lden>& l, const RationalPolynom<Tr, rnum, rden>& r){return RationalPolynom(l.num * r.den - l.den * r.num, l.den * r.den);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator - (const RationalPolynom<Tl,lnum,lden>& l, const Tr& r){return RationalPolynom(l.num - l.den * r, l.den);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator - (const Tl& l, const RationalPolynom<Tr, rnum, rden>& r){return RationalPolynom(r.den * l - r.num, r.den);}

	template<class Tl, size_t lnum, size_t lden>
	constexpr auto operator - (const RationalPolynom<Tl,lnum,lden>& p){return RationalPolynom(-p.num, p.den);}

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator * (const RationalPolynom<Tl,lnum,lden>& lhs, const RationalPolynom<Tr, rnum, rden>& rhs){return RationalPolynom(lhs.num * rhs.num, lhs.den * rhs.den);}
	template<class Tl, size_t lnum, size_t lden, class Tr>
	constexpr auto operator * (const RationalPolynom<Tl,lnum,lden>& lhs, const Tr& rhs){return RationalPolynom(lhs.num * rhs, lhs.den);}
	template<class Tl, class Tr, size_t rnum, size_t rden>
	constexpr auto operator * (const Tl& lhs, const RationalPolynom<Tr, rnum, rden>& rhs){return RationalPolynom(lhs * rhs.num, rhs.den);}

	template<class Tl, size_t lnum, size_t lden, class Tr, size_t rnum, size_t rden>
	constexpr auto operator / (const RationalPolynom<Tl,lnum,lden>& lhs, const RationalPolynom<Tr, rnum, rden>& rhs){return RationalPolynom(lhs.num * rhs.den, lhs.den * rhs.num);}
	template<class Tl, size_t lnum, size_t ldenom, class Tr>
	constexpr auto operator / (const RationalPolynom<Tl,lnum,ldenom>& lhs, const Tr& rhs){return RationalPolynom(lhs.num, lhs.den * rhs);}
	template<class Tl, class Tr, size_t rnum, size_t rdenom>
	constexpr auto operator / (const Tl& lhs, const RationalPolynom<Tr,rnum,rdenom>& rhs){return RationalPolynom(lhs * rhs.den, rhs.num);}

	template<class Stream, class T, size_t num_size, size_t denom_size>
	Stream& operator << (Stream& stream, const RationalPolynom<T,num_size,denom_size>& rp){return stream << "(" << rp.num << ") / (" << rp.den << ")";}

	template<class Stream, class T, size_t num_size, size_t denom_size>
	Stream& print_pretty (Stream& stream, const RationalPolynom<T,num_size,denom_size>& rp, const char* name = nullptr, const char* varname = "x", const char* indentation = ""){
		std::stringstream str_num;
		std::stringstream str_den;
		
		print(str_num, rp.num, varname);
		print(str_den, rp.den, varname);
		
		const size_t name_len = (name == nullptr) ? 0 : std::strlen(name)+3;
		const char* assignment = (name == nullptr) ? "" : " = ";
		
		if (str_num.str().size() > str_den.str().size()){
			size_t frac_size = str_num.str().size();
			stream << indentation;
			for(size_t i = 0; i < name_len; i++) stream << ' ';
			stream << str_num.str() << '\n';
			
			stream << indentation << name << assignment;
			for(size_t i = 0; i < frac_size; ++i) stream << '-';
			stream << '\n';
			
			stream << indentation;
			size_t offset = (str_num.str().size() - str_den.str().size())/2;
			for(size_t i = 0; i < offset + name_len; ++i) stream << ' ';
			stream << str_den.str() << '\n';
		}else{
			size_t frac_size = str_den.str().size();
			stream << indentation;
			size_t offset = (str_den.str().size() - str_num.str().size())/2;
			for(size_t i = 0; i < offset + name_len; ++i) stream << ' ';
			stream << str_num.str() << '\n';
			
			stream << indentation << name << assignment;
			for(size_t i = 0; i < frac_size; ++i) stream << '-';
			stream << '\n';
			
			stream << indentation;
			for(size_t i = 0; i < name_len; i++) stream << ' ';
			stream << str_den.str() << '\n';
		}
		return stream;
	}
}