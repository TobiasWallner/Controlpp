#pragma once

#include "BaseTypeArithmetic.hpp"

#include "definitions.hpp"


namespace twmath{

	template<class T, size_t N>
	struct VectorRef : public StaticVectorToken{
		T* const _ptr = nullptr;
		const size_t _increment = 1;
		
	public:
	
		using ValueType = T;
	
		VectorRef(T* ptr, size_t increment)
			: _ptr(ptr)
			, _increment(increment){}
			
		inline VectorRef& operator=(const T (&array)[N]){
			this->assign(array);
			return *this;
		}
		
		constexpr VectorRef(const VectorRef&) = default;
		inline VectorRef& operator= (const VectorRef<T, N>& other){
			this->assign(other.begin());
			return *this;
		}
		
		template<class Ta>
		constexpr VectorRef(const VectorRef<Ta, N>& v){this->assign(v.begin());}
		
		template<class Ta>
		constexpr VectorRef& operator=(const VectorRef<Ta, N>& v){return this->assign(v.begin());}
		
		constexpr VectorRef(const T* first, const T* last){this->assign(first, last);}
		
		template<class Ta, size_t Na>
		constexpr VectorRef(const VectorRef<Ta, Na>& v){this->assign(v.begin(), v.end());}
		
		template<class Ta, size_t Na>
		constexpr VectorRef& operator=(const VectorRef<Ta, Na>& v){this->assign(v.begin(), v.end());}
		
		template<size_t sub_size>
		constexpr VectorRef<T, sub_size> sub_vector(size_t offset, size_t stride) const {
			return VectorRef(&this->at(offset), this->_increment * stride);
		}
		
		template<bool b>
		constexpr VectorRef& assign(const VectorRef<T, N>& v){
			for(size_t i = 0; i < N; ++i){this->at(i) = v.at(i);}
			return *this;
		}
		
		template<class Itr>
		constexpr VectorRef& assign(Itr first, const Itr last){
			for(size_t i = 0; i < N && first != last; ++i, (void)++first){this->at(i) = *first;}
			return *this;
		}
		
		constexpr VectorRef& assign(const T (&array)[N]){
			for(size_t i = 0; i < N; ++i) this->at(i) = array[i];
			return *this;
		}
		
		template<class Itr>
		constexpr VectorRef& assign(Itr first){
			for(size_t i = 0; i < N; ++i, (void)++first){this->at(i) = *first;}
			return *this;
		}
		
		constexpr T& operator[] (ptrdiff_t i) {return this->at(i);}
		constexpr const T& operator[](ptrdiff_t i) const {return this->at(i);}
		
		constexpr T& at(ptrdiff_t i) {
			i = (i >= 0) ? i : i + this->size();
			twmath_assert(i < this->size(), "Out of range element access. Accessed element '" << i << "', but vector has size '" << this->size() << "'.");
			return *(this->_ptr + i * this->_increment);
		}
		
		constexpr const T& at(ptrdiff_t i) const {
			i = (i >= 0) ? i : i + this->size();
			twmath_assert(i < this->size(), "Out of range element access. Accessed element '" << i << "', but vector has size '" << this->size() << "'.");
			return *(this->_ptr + i * this->_increment);
		}
		
		static constexpr size_t ssize(){return N;}
		constexpr size_t size()const{return N;}
		
		constexpr VectorRef& fill(const T& value){
			for(T& elem : *this) elem = value;
			return *this;
		}
		
		constexpr VectorRef& fill(const T& value, ptrdiff_t first, ptrdiff_t last){
			first = (first >= 0) ? first : first + this->size();
			last = (last >= 0) ? last : last + this->size();
			for(size_t i = first; i < this->size() && i != last; ++i){
				this->at(i) = value;
			}
			return *this;
		}
		
		constexpr VectorRef& set_zero(){
			for(T& elem : *this) elem = T(0);
			return *this;
		}
		
		constexpr VectorRef& set_one(){
			for(T& elem : *this) elem = T(1);
			return *this;
		}
	};
	
	template<class T, size_t N> struct value_type<VectorRef<T, N>>{using type = T;};
	
}