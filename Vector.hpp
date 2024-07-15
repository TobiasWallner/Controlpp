#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <utility>
#include <type_traits>

#include "BaseTypeArithmetic.hpp"

/*
	Implicit conversion between templated classes??? -> template parameter deduction failed.
*/

/*
	To disable boundary checks:
		Define: DISABLE_SMATH_BOUNDARY_CHECKS
		
	To handle errors, define one of the following
		SMATH_THROW 				(default)
		SMATH_EXIT 					(calls exit(-1) to exit the program)
		SMATH_TRAP 					(traps the execution in a while loop)
		SMATH_CUSTOM_ERROR_HANDLER 	(calls the user defined function: 'void twmath_custom_error_handler(const char* error_message)')
		
	Set the error output stream by defining (only for SMATH_EXIT, SMATH_TRAP or SMATH_RETURN_ZERO):
		SMATH_CERR	(default is std::cerr)
*/

#define TWMATH_ENABLE_IF(condition) typename std::enable_if<(condition), int>::type = 0

#ifdef DISABLE_SMATH_BOUNDARY_CHECKS
	#define twmath_assert(condition, message)
#else
	
	#ifndef SMATH_CERR
		#define SMATH_CERR std::cout
	#endif

	#if defined(SMATH_CUSTOM_ERROR_HANDLER)
		void twmath_custom_error_handler(const char* error_message);
		#define twmath_assert(condition, message) if(!(condition)){std::stringstream s; s << message; twmath_custom_error_handler(s.str().c_str());}
	#elif defined(SMATH_TRAP)
		#define twmath_assert(condition, message) if(!(condition)){SMATH_CERR << message; while(true){};}
	#elif defined(SMATH_EXIT)
		#define twmath_assert(condition, message) if(!(condition)){SMATH_CERR << message; exit(-1);}
	#else
		#include <exception>
		#define twmath_assert(condition, message) if(!(condition)){std::stringstream s; s << message; throw std::runtime_error(s.str().c_str());}
	#endif

#endif

namespace twmath{
	
	template<class T, size_t N>
	struct Vector{
		public:
		T _values[N];

		
		Vector() = default;
		constexpr Vector(const T& v){for(T& elem : *this) elem = v;}
		constexpr Vector(const T (&array)[N]){
			this->assign(array);
		}
		inline Vector& operator=(const T (&array)[N]){
			this->assign(array);
			return *this;
		}
		
		constexpr Vector(const Vector&) = default;
		inline Vector& operator= (const Vector<T, N>& other){
			this->assign(other.begin());
			return *this;
		}
		
		template<class Ta>
		constexpr Vector(const Vector<Ta, N>& v){this->assign(v.begin());}
		
		template<class Ta>
		constexpr Vector& operator=(const Vector<Ta, N>& v){return this->assign(v.begin());}
		
		constexpr Vector(const T* first, const T* last){this->assign(first, last);}
		
		template<class Ta, size_t Na>
		constexpr Vector(const Vector<Ta, Na>& v){this->assign(v.begin(), v.end());}
		
		template<class Ta, size_t Na>
		constexpr Vector& operator=(const Vector<Ta, Na>& v){this->assign(v.begin(), v.end());}
		
		template<size_t sub_size>
		constexpr Vector<T, sub_size> subvector(ptrdiff_t offset) const {
			offset = (offset >= 0) ? offset : offset + this->size();
			twmath_assert(sub_size + offset <= this->size(), 
					"subvector would read from elements past the end of this:'\n'"
					<< "failed because of: " << "sub_size + offset > this->size(): "
					<< "(" << sub_size << " + " << offset << " > " << this->size() << ")");
			Vector<T, sub_size> result;
			size_t result_i = 0;
			size_t this_i = offset;
			for(; result_i < sub_size && this_i < N; ++result_i, (void)++this_i){
				result.at(result_i) = this->at(this_i);
			}
			return result;
		}
		
		template<size_t sub_size>
		constexpr Vector& assign(const Vector<T, sub_size>& subvector, ptrdiff_t offset){
			offset = (offset >= 0) ? offset : offset + this->size();
			twmath_assert(subvector.size() + offset <= this->size(), 
					"subvector would write to elements past the end of this:'\n'"
					<< "failed because of: " << "subvector.size() + offset > this->size(): "
					<< "(" << subvector.size() << " + " << offset << " > " << this->size() << ")");
			size_t result_i = 0;
			size_t this_i = offset;
			for(; result_i < subvector.size(); ++result_i, (void)++this_i){
				this->at(this_i) = subvector.at(result_i);
			}
			return *this;
		}
		
		template<class Itr>
		constexpr Vector& assign(Itr first, const Itr last){
			auto this_first = this->begin();
			const auto this_end = this->end();
			for(; first != last; ++this_first, (void)++first){
				twmath_assert(this_first != this_end, "Out of range assignment");
				*this_first = *first;
			}
			return *this;
		}
		
		template<class Itr>
		constexpr Vector& assign(Itr first, const Itr last, T rest){
			auto this_first = this->begin();
			const auto this_end = this->end();
			for(; this_first != this_end && first != last; ++this_first, (void)++first){
				*this_first = *first;
			}
			for(; this_first != this_end; ++this_first){
				*this_first = rest;
			}
			return *this;
		}
		
		constexpr Vector& assign(const T (&array)[N]){
			for(size_t i = 0; i < N; ++i) this->at(i) = array[i];
			return *this;
		}
		
		template<class Itr>
		constexpr Vector& assign(Itr first){
			auto this_first = this->begin();
			const auto this_end = this->end();
			for(; this_first != this_end; ++this_first, (void)++first){
				*this_first = *first;
			}
			return *this;
		}
		
		constexpr T& operator[] (ptrdiff_t i) {return this->at(i);}
		constexpr const T& operator[](ptrdiff_t i) const {return this->at(i);}
		
		constexpr T& at(ptrdiff_t i) {
			i = (i >= 0) ? i : i + this->size();
			twmath_assert(i < this->size(), "Out of range element access. Accessed element '" << i << "', but vector has size '" << this->size() << "'.");
			return this->_values[i];
		}
		
		constexpr const T& at(ptrdiff_t i) const {
			i = (i >= 0) ? i : i + this->size();
			twmath_assert(i < this->size(), "Out of range element access. Accessed element '" << i << "', but vector has size '" << this->size() << "'.");
			return this->_values[i];
		}
		
		constexpr T* begin(){return this->_values;}
		constexpr const T* begin()const{return this->_values;}
		constexpr const T* cbegin()const{return this->_values;}
		
		constexpr T* end(){return this->_values + this->size();}
		constexpr const T* end()const{return this->_values + this->size();}
		constexpr const T* cend()const{return this->_values + this->size();}
		
		constexpr size_t size()const{return N;}
		
		constexpr T& front(){return *this->begin();}
		constexpr const T& front()const{return *this->begin();}
		
		constexpr T& back(){return *(this->end()-1);}
		constexpr const T& back()const{return *(this->end()-1);}
		
		template<size_t M>
		constexpr Vector<T, M> subvector(ptrdiff_t offset){
			offset = (offset >= 0) ? offset : offset + this->size();
			return Vector<T, M>(this->begin() + offset, this->end());
		}
		
		constexpr Vector& fill(const T& value){
			for(T& elem : *this) elem = value;
			return *this;
		}
		
		constexpr Vector& fill(const T& value, ptrdiff_t first, ptrdiff_t last){
			first = (first >= 0) ? first : first + this->size();
			last = (last >= 0) ? last : last + this->size();
			for(size_t i = first; i < this->size() && i != last; ++i){
				this->at(i) = value;
			}
			return *this;
		}
		
		constexpr Vector& set_zero(){
			for(T& elem : *this) elem = T(0);
			return *this;
		}
		
		constexpr Vector& set_zero(ptrdiff_t first, ptrdiff_t last){
			first = (first >= 0) ? first : first + this->size();
			last = (last >= 0) ? last : last + this->size();
			for(T* itr = this->begin() + first; itr < this->end() && itr < this->begin() + last; ++itr){
				*itr = T(0);
			}
			return *this;
		}
		
		constexpr Vector& set_one(){
			for(T& elem : *this) elem = T(1);
			return *this;
		}
		
		constexpr Vector& set_one(size_t first, size_t last){
			for(T* itr = this->begin() + first; itr < this->end() && itr < this->begin() + last; ++itr){
				*itr = T(1);
			}
			return *this;
		}
		
		template <class U = T, typename std::enable_if<(std::is_same<T, U>::value && N == 1), int>::type = 0>
		constexpr U to_scalar() const {
			return U(*this->begin());
		}
		
		// conversion to scalar
		template <class U = T, typename std::enable_if<(std::is_same<T, U>::value && N == 1), int>::type = 0>
		constexpr operator U () const {
			return U(*this->begin());
		}
	};
	
	
	template<class Tl, class Tr, size_t N, class F>
	constexpr Vector<F, N> for_each (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs, F (*f)(const Tl&, const Tr&)){
		Vector<F, N> result;
		for(size_t i = 0; i < N; ++i) result.at(i) = f(lhs.at(i), rhs.at(i));
		return result;
	}
	template<class Tl, class Tr, size_t N, class F>
	constexpr Vector<F, N> for_each (const Tl& lhs, const Vector<Tr, N>& rhs, F (*f)(const Tl&, const Tr&)){
		Vector<F, N> result;
		for(size_t i = 0; i < N; ++i) result.at(i) = f(lhs, rhs.at(i));
		return result;
	}
	template<class Tl, class Tr, size_t N, class F>
	constexpr Vector<F, N> for_each (const Vector<Tl, N>& lhs, const Tr& rhs, F (*f)(const Tl&, const Tr&)){
		Vector<F, N> result;
		for(size_t i = 0; i < N; ++i) result.at(i) = f(lhs.at(i), rhs);
		return result;
	}
	template<class T, size_t N, class F>
	constexpr Vector<F, N> for_each (const Vector<T, N>& v, F (*f)(const T&)){
		Vector<F, N> result;
		for(size_t i = 0; i < N; ++i) result.at(i) = f(v.at(i));
		return result;
	}
	
	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() == std::declval<Tr>()), N> equal (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() == std::declval<Tr>()), N> equal (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() == std::declval<Tr>()), N> equal (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() == std::declval<Tr>()), N> operator== (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() == std::declval<Tr>()), N> operator== (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() == std::declval<Tr>()), N> operator== (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::equal<Tl, Tr>);}
	
	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() != std::declval<Tr>()), N> not_equal (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::not_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() != std::declval<Tr>()), N> not_equal (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::not_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() != std::declval<Tr>()), N> not_equal (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::not_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() != std::declval<Tr>()), N> operator!= (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::not_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() != std::declval<Tr>()), N> operator!= (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::not_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() != std::declval<Tr>()), N> operator!= (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::not_equal<Tl, Tr>);}

	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() < std::declval<Tr>()), N> less (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::less<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() < std::declval<Tr>()), N> less (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::less<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() < std::declval<Tr>()), N> less (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::less<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() < std::declval<Tr>()), N> operator< (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::less<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() < std::declval<Tr>()), N> operator< (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::less<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() < std::declval<Tr>()), N> operator< (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::less<Tl, Tr>);}
	
	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() <= std::declval<Tr>()), N> less_equal (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::less_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() <= std::declval<Tr>()), N> less_equal (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::less_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() <= std::declval<Tr>()), N> less_equal (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::less_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() <= std::declval<Tr>()), N> operator<= (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::less_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() <= std::declval<Tr>()), N> operator<= (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::less_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() <= std::declval<Tr>()), N> operator<= (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::less_equal<Tl, Tr>);}
	
	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() > std::declval<Tr>()), N> greater (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::greater<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() > std::declval<Tr>()), N> greater (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::greater<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() > std::declval<Tr>()), N> greater (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::greater<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() > std::declval<Tr>()), N> operator> (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::greater<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() > std::declval<Tr>()), N> operator> (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::greater<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() > std::declval<Tr>()), N> operator> (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::greater<Tl, Tr>);}
	
	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() >= std::declval<Tr>()), N> greater_equal (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::greater_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() >= std::declval<Tr>()), N> greater_equal (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::greater_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() >= std::declval<Tr>()), N> greater_equal (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::greater_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() >= std::declval<Tr>()), N> operator>= (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::greater_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() >= std::declval<Tr>()), N> operator>= (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::greater_equal<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() >= std::declval<Tr>()), N> operator>= (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::greater_equal<Tl, Tr>);}
	
	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() + std::declval<Tr>()), N> add (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() + std::declval<Tr>()), N> add (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() + std::declval<Tr>()), N> add (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() + std::declval<Tr>()), N> operator+ (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() + std::declval<Tr>()), N> operator+ (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() + std::declval<Tr>()), N> operator+ (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() - std::declval<Tr>()), N> sub (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() - std::declval<Tr>()), N> sub (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() - std::declval<Tr>()), N> sub (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() - std::declval<Tr>()), N> operator- (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() - std::declval<Tr>()), N> operator- (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() - std::declval<Tr>()), N> operator- (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() * std::declval<Tr>()), N> mul (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() * std::declval<Tr>()), N> mul (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() * std::declval<Tr>()), N> mul (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr auto operator* (const Vector<Tl, N>& lhs, const Tr& rhs){return mul(lhs, rhs);}
	template<class Tl, class Tr, size_t N>
	constexpr auto operator* (const Tl& lhs, const Vector<Tr, N>& rhs){return mul(lhs, rhs);}
	
	
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() / std::declval<Tr>()), N> div (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() / std::declval<Tr>()), N> div (const Vector<Tl, N>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr Vector<decltype(std::declval<Tl>() / std::declval<Tr>()), N> div (const Tl& lhs, const Vector<Tr, N>& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	template<class Tl, class Tr, size_t N>
	constexpr auto operator/ (const Vector<Tl, N>& lhs, const Tr& rhs){return div(lhs, rhs);}
	template<class Tl, class Tr, size_t N>
	constexpr auto operator/ (const Tl& lhs, const Vector<Tr, N>& rhs){return div(lhs, rhs);}
	
	template<class T, size_t N>constexpr Vector<decltype(-std::declval<T>()), N> negate (const Vector<T, N>& a){return for_each(a, twmath_base::negate<T>);}
	template<class T, size_t N>constexpr Vector<decltype(-std::declval<T>()), N> operator- (const Vector<T, N>& a){return for_each(a, twmath_base::negate<T>);}
	
	template<class T, size_t N>
	constexpr T reduce(const Vector<T, N>& v, T (*f)(const T&, const T&), T init = T(0)){
		for (const T& elem : v) init = f(init, elem);
		return init;
	}
	
	template<class T, size_t N>
	constexpr T sum(const Vector<T, N>& v){return reduce(v, twmath_base::add<T>, T(0));}
	
	template<class T, size_t N>
	constexpr T prod(const Vector<T, N>& v){return reduce(v, twmath_base::mul<T>, T(1));}
	
	template<class T> T min(const T& lhs, const T& rhs){return (lhs < rhs) ? lhs : rhs;}
	template<class T> T max(const T& lhs, const T& rhs){return (lhs > rhs) ? lhs : rhs;}
	
	template<class T, size_t N>
	constexpr T min(const Vector<T, N>& v){return reduce(v, min, v[0]);}
	
	template<class T, size_t N>
	constexpr T max(const Vector<T, N>& v){return reduce(v, max, v[0]);}
	
	template<class T, size_t N>
	constexpr T mean(const Vector<T, N>& v){return sum(v)/N;}
	
	template<class T, size_t N>
	constexpr T square(const Vector<T, N>& v){return mul(v, v);}
	
	template<class Tl, class Tr, size_t N>
	constexpr auto dot(const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return sum(mul(lhs, rhs));}
	template<class Tl, class Tr, size_t N>
	constexpr auto operator* (const Vector<Tl, N>& lhs, const Vector<Tr, N>& rhs){return dot(lhs, rhs);}
	
	template<class T, size_t N>
	constexpr T var(const Vector<T, N>& v){return sum(square(v - mean(v))) / N;}
	
	template<class T, size_t N>
	constexpr T stddev(const Vector<T, N>& v){return sqrt(var(v));}
	
	
	
	template<class T,size_t N>
	constexpr Vector<T, N> sqr_norm(const Vector<T, N>& v){return dot(v, v);}
	
	template<class T,size_t N>
	constexpr Vector<T, N> norm(const Vector<T, N>& v){return sqrt(sqr_norm(v));}
	
	template<class T,size_t N>
	constexpr Vector<T, N> normalize(const Vector<T, N>& v){return v / norm(v);}
	
	template<class T,size_t N>
	constexpr Vector<T, N> reverse(const Vector<T, N>& v){
		Vector<T, N> result;
		for (size_t i = 0, j = N-1; i < N; ++i, --j) 
			result.at(i) = v.at(j);
		return result;
	}
	
	template<class T,size_t N>
	constexpr Vector<T, N> rotate_r(const Vector<T, N>& v, size_t n){
		Vector<T, N> result;
		size_t i = 0;
		size_t j = n;
		for (; j < N; ++i, (void)++j) result.at(j) = v.at(i);
		j = 0;
		for (; i < N; ++i, (void)++j) result.at(j) = v.at(i);
		return result;
	}
	
	template<class T,size_t N>
	constexpr Vector<T, N> rotate_l(const Vector<T, N>& v, size_t n){
		Vector<T, N> result;
		size_t i = 0;
		size_t j = N - n;
		for (; j < N; ++i, (void)++j) result.at(j) = v.at(i);
		j = 0;
		for (; i < N; ++i, (void)++j) result.at(j) = v.at(i);
		return result;
	}
	
	template<class Tl, class Tr>
	constexpr Vector<decltype(std::declval<Tl>() * std::declval<Tr>()), 3> cross(const Vector<Tl, 3>& lhs, const Vector<Tr, 3>& rhs){
		Vector<decltype(std::declval<Tl>() * std::declval<Tr>()), 3> result({
			lhs.at(1) * rhs.at(2) - lhs.at(2) * rhs.at(1),
			lhs.at(2) * rhs.at(0) - lhs.at(0) * rhs.at(2),
			lhs.at(0) * rhs.at(1) - lhs.at(1) * rhs.at(0),
		});
		return result;
	}
	
	template<class T, size_t N>
	constexpr Vector<T, N> linspace(const T& start, const T& end){
		Vector<T, N> result;
		T step = (end - start) / (N + 1);
		for(size_t i = 0; i < N; ++i)
			result[i] = start + step * i;
		return result;
	}
	
	template<class T, size_t N>
	constexpr Vector<T, N> fill(const T& n){
		Vector<T, N> result;
		for(size_t i = 0; i < N; ++i) result.at(i) = n;
		return result;
	}
	
	template<class T, size_t N>
	constexpr Vector<T, N> zeros(){return fill<T, N>(0);}
	
	template<class T, size_t N>
	constexpr Vector<T, N> ones(){return fill<T, N>(1);}
	
	template<class T, size_t N>
	constexpr Vector<T, N> logspace(const T& start, const T& end, size_t num){
		using namespace std;
		return for_each(linspace(log2(start), log2(end), num), exp2);
	}
	
	template<class T, size_t N>
	constexpr auto softmax(const Vector<T, N>& v){
		using namespace std;
		return normalize(for_each(v, exp));
	}
	
	template<class Tl, class Tr, size_t Nl, size_t Nr>
	constexpr Vector<decltype(std::declval<Tl>() * std::declval<Tr>()), Nl + Nr - 1> conv (const Vector<Tl, Nl>& lhs, const Vector<Tr, Nr>& rhs){
		using ValueType = decltype(std::declval<Tl>() * std::declval<Tr>());
		Vector<ValueType, Nl + Nr - 1> result;
		size_t resi = 0;
		size_t lhs_first = 0;
		size_t rhs_first = 0;
		size_t lhs_last = 1;
		size_t rhs_last = 1;
		for(; resi < Nl + Nr - 1; ++resi){
			size_t lhsi = lhs_last-1;
			size_t rhsi = rhs_first;
			ValueType sum(0);
			for(; rhsi < rhs_last; ++rhsi, --lhsi) sum += lhs.at(lhsi) * rhs.at(rhsi);
			result.at(resi) = sum;
			
			rhs_first += (lhs_last == Nl);
			lhs_first += (rhs_last == Nr);
			lhs_last += (lhs_last != Nl);
			rhs_last += (rhs_last != Nr);
			
		}
		return result;
	}
	
	template<class Tl, class Tr, size_t Nl, size_t Nr>
	constexpr Vector<decltype(std::declval<Tl>() * std::declval<Tr>()), Nl + Nr - 1> corr (const Vector<Tl, Nl>& lhs, const Vector<Tr, Nr>& rhs){
		using ValueType = decltype(std::declval<Tl>() * std::declval<Tr>());
		Vector<ValueType, Nl + Nr - 1> result;
		size_t resi = 0;
		size_t lhs_first = 0;
		size_t rhs_first = 0;
		size_t lhs_last = 1;
		size_t rhs_last = 1;
		for(; resi < Nl + Nr - 1; ++resi){
			size_t lhsi = lhs_last-1;
			size_t rhsi = rhs_first;
			ValueType sum(0);
			for(; rhsi < rhs_last; ++rhsi, --lhsi){
				sum += lhs[Nl-1-lhsi] * rhs[rhsi];
			}
			result.at(resi) = sum;
			
			rhs_first += (lhs_last == Nl);
			lhs_first += (rhs_last == Nr);
			lhs_last += (lhs_last != Nl);
			rhs_last += (rhs_last != Nr);
			
		}
		return result;
	}
	
	template<class Stream, class T, size_t N>
	Stream& operator << (Stream& stream, const Vector<T, N>& v){
		auto itr = v.begin();
		if(itr != v.end()){
			stream << *itr;
			++itr;
		}
		for(; itr != v.end(); ++itr){
			stream << ", " << *itr;
		}
		return stream;
	}
}

