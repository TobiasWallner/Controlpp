#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <utility>
#include <type_traits>

#include "BaseTypeArithmetic.hpp"
#include "definitions.hpp"
#include "VectorToken.hpp"

namespace twmath{
	
	template<class T, size_t N>
	struct Vector : public StaticVectorToken{
		public:
		using ValueType = T;
		
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
		
		static constexpr size_t ssize(){return N;}
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
	
	template<class T, size_t N> struct value_type<Vector<T, N>>{using type = T;};

}

