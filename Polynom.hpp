#pragma once

#include "Vector.hpp"

namespace twmath{
	template<class T, size_t N>
	struct Polynom {
		twmath::Vector<T, N> vec;
		
		constexpr Polynom() = default;
		constexpr Polynom(const Polynom&) = default;
		constexpr Polynom(const T& v) : Polynom(&v, &v+1){}
		
		template<class Iterator>
		constexpr Polynom(Iterator first, Iterator last){this->assign(first, last);}
		
		template<class Iterator>
		constexpr Polynom(Iterator first){this->assign(first);}
		
		constexpr Polynom(const T (&array)[N]) : vec(array){}
		constexpr Polynom& operator=(const T (&array)[N]){vec = array; return *this;}
		
		template<class Ta, size_t Na>
		constexpr Polynom(const Polynom<Ta, Na>& p) : vec(p.vec){}
		
		template<class Ta, size_t Na>
		constexpr Polynom& operator= (const Polynom<Ta, Na>& p){this->vec = p.vec; return *this;}
		
		constexpr Polynom(const Vector<T, N>& v) : vec(v){}
		
		template<class Ta, size_t Na>
		constexpr Polynom(const Vector<Ta, Na>& v) : vec(v){}
		
		template<class Ta, size_t Na>
		constexpr Polynom& operator= (const Vector<Ta, Na>& v){this->vec = v; return *this;}
		
		template<class Ta> 
		constexpr decltype(std::declval<T>() * std::declval<Ta>()) eval(Ta x) const {
			decltype(std::declval<T>() * std::declval<Ta>()) result = static_cast<decltype(std::declval<T>() * std::declval<Ta>())>(this->at(0));
			if(1 < N){
				result += x * this->at(1);
			}
			for(size_t i = 2; i < N; ++i){
				x *= x;
				result += x * this->at(i);
			}
			return result;
		}

		template<class Ta>  
		constexpr decltype(std::declval<T>() * std::declval<Ta>()) operator() (const Ta& x) const {return this->eval(x);}
		
		constexpr T& at(size_t i){return vec.at(i);}
		constexpr const T& at(size_t i)const{return vec.at(i);}
		constexpr T& operator[](size_t i){return vec.at(i);}
		constexpr const T& operator[](size_t i)const{return vec.at(i);}
		
		constexpr T* begin(){return vec.begin();}
		constexpr const T* begin()const{return vec.begin();}
		constexpr const T* cbegin()const{return vec.begin();}
		
		constexpr T* end(){return vec.end();}
		constexpr const T* end()const{return vec.end();}
		constexpr const T* cend()const{return vec.end();}
		
		constexpr size_t size() const {return vec.size();}

		constexpr T& front(){return *this->begin();}
		constexpr const T& front()const{return *this->begin();}
		
		constexpr T& back(){return *(this->end()-1);}
		constexpr const T& back()const{return *(this->end()-1);}

		constexpr size_t order() const {
			size_t i = vec.size();
			while (i > 0 && vec[i-1] == static_cast<T>(0)) --i;
			return i - 1;
		}
	
		template<class Iterator>
		constexpr Polynom& assign(Iterator first, const Iterator last){
			auto this_itr = this->begin();
			for (;  this_itr != this->end() && first != last; ++this_itr, ++first) *this_itr = *first;
			for (; this_itr != this->end(); ++this_itr) *this_itr = T(0);
			return *this;
		}
		
		template<class Iterator>
		constexpr Polynom& assign(Iterator first){vec.assign(first); return *this;}
		
	};

	template<class Tl, size_t Nl, class Tr, size_t Nr> 
	constexpr Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), (Nl > Nr) ? Nl : Nr> operator+ (const Polynom<Tl, Nl>& lhs, const Polynom<Tr, Nr>& rhs){
		Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), (Nl > Nr) ? Nl : Nr> result;
		size_t i = 0;
		for(; i < Nl && i < Nr; ++i) result[i] = lhs[i] + rhs[i];
		for(; i < Nl; ++i) result[i] = lhs[i];
		for(; i < Nr; ++i) result[i] = rhs[i];
		return result;
	}

	template<class Tl, class Tr, size_t N> 
	constexpr Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), N> operator+ (const Tl& lhs, const Polynom<Tr, N>& rhs){
		Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), N> result;
		result[0] = lhs + rhs[0];
		std::copy(rhs.begin()+1, rhs.end(), result.begin()+1);
		return result;
	}

	template<class Tl, class Tr, size_t N> 
	constexpr Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), N> operator+ (const Polynom<Tl, N>& lhs, const Tr& rhs){
		Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), N> result;
		result[0] = lhs[0] + rhs;
		std::copy(lhs.begin()+1, lhs.end(), result.begin()+1);
		return result;
	}

	template<class Tl, size_t Nl, class Tr, size_t Nr> 
	constexpr Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), (Nl > Nr) ? Nl : Nr> operator- (const Polynom<Tl, Nl>& lhs, const Polynom<Tr, Nr>& rhs){
		Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), (Nl > Nr) ? Nl : Nr> result;
		size_t i = 0;
		for(; i < Nl && i < Nr; ++i) result[i] = lhs[i] - rhs[i];
		for(; i < Nl; ++i) result[i] = lhs[i];
		for(; i < Nr; ++i) result[i] = -rhs[i];
		return result;
	}

	template<class Tl, class Tr, size_t N> 
	constexpr Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), N> operator- (const Tl& lhs, const Polynom<Tr, N>& rhs){
		Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), N> result;
		result[0] = lhs - rhs[0];
		for(size_t i = 1; i < result.size(); ++i) result[i] = -rhs;
		return result;
	}
	
	template<class Tr, size_t N> 
	constexpr Polynom<decltype(-std::declval<Tr>()), N> operator- (const Polynom<Tr, N>& rhs){
		return -rhs.vec;
	}

	template<class Tl, class Tr, size_t N> 
	constexpr Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), N> operator- (Polynom<Tl, N> lhs, const Tr& rhs){
		Polynom<decltype(std::declval<Tl>() + std::declval<Tr>()), N> result;
		result[0] = lhs[0] - rhs;
		std::copy(lhs.begin()+1, lhs.end(), result.begin()+1);
		return result;
	}

	template<class Tl, class Tr, size_t Nl, size_t Nr> 
	constexpr Polynom<decltype(std::declval<Tl>() * std::declval<Tr>()), Nl + Nr - 1> operator* (const Polynom<Tl, Nl>& lhs, const Polynom<Tr, Nr>& rhs){
		return twmath::conv(lhs.vec, rhs.vec);
	}

	template<class Tl, class Tr, size_t Nl> 
	constexpr Polynom<decltype(std::declval<Tl>() * std::declval<Tr>()), Nl> operator* (const Polynom<Tl, Nl>& lhs, const Tr& rhs){
		return Polynom<decltype(std::declval<Tl>() * std::declval<Tr>()), Nl>(lhs.vec * rhs);
	}

	template<class Tl, class Tr, size_t Nr> 
	constexpr Polynom<decltype(std::declval<Tl>() * std::declval<Tr>()), Nr> operator* (const Tl& lhs, const Polynom<Tr, Nr>& rhs){
		return Polynom<decltype(std::declval<Tl>() * std::declval<Tr>()), Nr>(lhs * rhs.vec);
	}

	template<class T, size_t N1, size_t N2>
	struct PolynomDivRest{
		Polynom<T, N1> div;
		Polynom<T, N2> mod;
	};

	template<class Tl, size_t Nl, class Tr, size_t Nr>
	constexpr PolynomDivRest<decltype(std::declval<Tl>() / std::declval<Tr>()), Nl, Nr - 1> div_mod (const Polynom<Tl, Nl>& l, const Polynom<Tr, Nr>& r){
		using ValueType = decltype(std::declval<Tl>() / std::declval<Tr>());
		PolynomDivRest<ValueType, Nl, Nr - 1> result;
		Polynom<ValueType, Nl> a(l);
		
		long ai = static_cast<long>(a.order());
		const auto r_order = static_cast<long>(r.order());
		long divi = ai - r_order;
		for (size_t i = divi+1; i < result.div.size(); ++i) 
			result.div[i] = static_cast<ValueType>(0);

		for(; divi >= 0; --ai, (void)--divi) {
			ValueType d = static_cast<ValueType>(a[ai] / r[r_order]);
			result.div[divi] = d;
			
			long aii = ai-1;
			long rii = static_cast<long>(r_order)-1;
			for(; rii >= 0; --aii, --rii)
				a[aii] -= r[rii] * d;

		}
		result.mod.assign(a.begin(), a.begin()+r.size()-1);
		return result;
	}
	
	template<class Tl, size_t Nl, class Tr, size_t Nr>
	constexpr Polynom<decltype(std::declval<Tl>() / std::declval<Tr>()), Nl> div (const Polynom<Tl, Nl>& l, const Polynom<Tr, Nr>& r){return div_mod<Tl, Nl, Tr, Nr>(l, r).div;}
	
	template<class Tl, size_t Nl, class Tr>
	constexpr Polynom<decltype(std::declval<Tl>() / std::declval<Tr>()), Nl> div (const Polynom<Tl, Nl>& l, const Tr& r){return l.vec / r;}
	
	template<class Tl, size_t Nl, class Tr, size_t Nr>
	constexpr Polynom<decltype(std::declval<Tl>() / std::declval<Tr>()), Nl> operator/ (const Polynom<Tl, Nl>& l, const Polynom<Tr, Nr>& r){return div_mod<Tl, Nl, Tr, Nr>(l, r).div;}
	
	template<class Tl, size_t Nl, class Tr>
	constexpr Polynom<decltype(std::declval<Tl>() / std::declval<Tr>()), Nl> operator/ (const Polynom<Tl, Nl>& l, const Tr& r){return l.vec / r;}
	
	template<class Tl, size_t Nl, class Tr, size_t Nr>
	constexpr Polynom<decltype(std::declval<Tl>() / std::declval<Tr>()), Nr - 1> div (const Polynom<Tl, Nl>& l, const Polynom<Tr, Nr>& r){return div_mod<Tl, Nl, Tr, Nr>(l, r).mod;}
	
	template<class Tl, size_t Nl, class Tr, size_t Nr>
	constexpr Polynom<decltype(std::declval<Tl>() / std::declval<Tr>()), Nr - 1> operator% (const Polynom<Tl, Nl>& l, const Polynom<Tr, Nr>& r){return div_mod<Tl, Nl, Tr, Nr>(l, r).mod;}

	template<class T, size_t N, class Ta> 
	constexpr decltype(std::declval<T>() * std::declval<Ta>()) eval (const Polynom<T, N>& f, const T& x){return f.eval(x);}

	template<class Stream, class T, size_t N>
	Stream& print (Stream& stream, const Polynom<T, N>& v, const char* symbol = "x"){
		auto itr = v.begin();
		size_t i = 0;
		bool is_first = true;
		bool is_second = true;
		for(; itr != v.end(); ++itr, ++i){
			if (*itr != T(0)) {
				if (is_first) {
					stream << *itr;
					is_first = false;
				}
				else if (is_second) {
					if (*itr == T(1)) {
						stream << " + " << symbol;
					}
					else {
						stream << " + " << *itr << symbol;
					}
					is_second = false;
				}
				else {
					if (*itr == T(1)) {
						stream << " + " << symbol << "^" << i;
					}
					else {
						stream << " + " << *itr << symbol << "^" << i;
					}
					
				}
			}
		}
		return stream;
	}

	template<class Stream, class T, size_t N>
	Stream& operator << (Stream& stream, const Polynom<T, N>& v){return print(stream, v, "x");}

}