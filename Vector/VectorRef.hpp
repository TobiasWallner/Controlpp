#pragma once

#include <type_traits>

#include "BaseTypeArithmetic.hpp"
#include "definitions.hpp"
#include "VectorTraits.hpp"
#include "ConstVectorRef.hpp"

namespace twmath{

	template<class T, size_t N>
	struct VectorRef : public StaticVectorToken{
		T* const _ptr = nullptr;
		const size_t _increment = 1;
		
	public:
	
		using ValueType = T;
	
		constexpr VectorRef(T* ptr, size_t increment=1)
			: _ptr(ptr)
			, _increment(increment){}
			
		constexpr VectorRef& operator=(const T (&array)[N]){
			this->assign(array);
			return *this;
		}
		
		constexpr VectorRef(const VectorRef&) = default;
		constexpr VectorRef& operator= (const VectorRef& v){return this->assign(v);}
		
		template<class V, TWMATH_ENABLE_IF(is_static_vector_v<V>), TWMATH_ENABLE_IF(V::ssize() == N)>
		inline VectorRef& operator= (const V& v){return this->assign(v);}
		
		template<class V, TWMATH_ENABLE_IF(is_dynamic_vector_v<V>)>
		inline VectorRef& operator= (const V& v){return this->assign(v);}
		
		template<size_t sub_size>
		constexpr VectorRef<T, sub_size> sub_vector(size_t offset=0, size_t stride=1) {
			twmath_assert(offset + sub_size * stride < N, "Error: Vector.sub_vector(size_t, size_t): subvector would access elemens not contained in origin");
			return VectorRef<T, sub_size>(&this->at(offset), this->_increment * stride);
		}
		
		template<size_t sub_size>
		constexpr ConstVectorRef<T, sub_size> sub_vector(size_t offset=0, size_t stride=1) const {
			twmath_assert(offset + sub_size * stride < N, "Error: Vector.sub_vector(size_t, size_t): subvector would access elemens not contained in origin");
			return ConstVectorRef<T, sub_size>(&this->at(offset), this->_increment * stride);
		}
		
		template<size_t sub_size>
		constexpr ConstVectorRef<T, sub_size> const_sub_vector(size_t offset=0, size_t stride=1) const {
			twmath_assert(offset + sub_size * stride < N, "Error: Vector.sub_vector(size_t, size_t): subvector would access elemens not contained in origin");
			return ConstVectorRef<T, sub_size>(&this->at(offset), this->_increment * stride);
		}
		
		template<class V, TWMATH_ENABLE_IF(is_static_vector_v<V>), TWMATH_ENABLE_IF(V::ssize() == N)>
		constexpr VectorRef& assign(const V& v){
			for(size_t i = 0; i < N; ++i){this->at(i) = v.at(i);}
			return *this;
		}
		
		template<class V, TWMATH_ENABLE_IF(is_dynamic_vector_v<V>)>
		constexpr VectorRef& assign(const V& v){
			twmath_assert(v.size() == this->size(), "Error: VectorRef: Cannot assign vector of size '" << v.size() << "' to vector of size '" << this->size() << "'.");
			for(size_t i = 0; i < N; ++i) this->at(i) = v.at(i);
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
		
		template<class Itr, TWMATH_ENABLE_IF(!is_vector_v<Itr>)>
		constexpr VectorRef& assign(Itr first){
			for(size_t i = 0; i < N; ++i, (void)++first){this->at(i) = *first;}
			return *this;
		}
		
		template<class Integer, TWMATH_ENABLE_IF(std::is_integral_v<Integer>)>
		constexpr T& operator[](Integer i) {return this->at(i);}
		
		template<class Integer, TWMATH_ENABLE_IF(std::is_integral_v<Integer>)>
		constexpr const T& operator[](Integer i) const {return this->at(i);}
		
		constexpr T& at(size_t i){
			twmath_assert(i < this->size(), "Out of range element access. Accessed element '" << i << "', but vector has size '" << this->size() << "'.");
			return this->_ptr[i * this->_increment];
		}
		
		template<class UInteger, TWMATH_ENABLE_IF(std::is_unsigned_v<UInteger>)>
		constexpr T& at(UInteger i){return this->at(static_cast<size_t>(i));}
		
		template<class SInteger, TWMATH_ENABLE_IF(std::is_signed_v<SInteger>)>
		constexpr T& at(SInteger i){
			size_t ui = (i >= 0) ? i : i + this->size();
			return this->at(ui);
		}
		
		constexpr const T& at(size_t i) const {
			twmath_assert(i < this->size(), "Out of range element access. Accessed element '" << i << "', but vector has size '" << this->size() << "'.");
			return this->_ptr[i * this->_increment];
		}
		
		template<class UInteger, TWMATH_ENABLE_IF(std::is_unsigned_v<UInteger>)>
		constexpr const T& at(UInteger i) const {return this->at(static_cast<size_t>(i));}
		
		template<class SInteger, TWMATH_ENABLE_IF(std::is_signed_v<SInteger>)>
		constexpr const T& at(SInteger i) const {
			size_t ui = (i >= 0) ? i : i + this->size();
			return this->at(ui);
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
			for(size_t i = 0; i < this->size(); ++i) this->at(i) = T(0);
			return *this;
		}
		
		constexpr VectorRef& set_one(){
			for(size_t i = 0; i < this->size(); ++i) this->at(i) = T(1);
			return *this;
		}
		
		// ---------------------- conversion to scalar ---------------------------
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && N == 1))>
		constexpr operator U& () {return this->at(0);}
		
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && N == 1))>
		constexpr operator const U& () const {return this->at(0);}
		
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && N == 1))>
		constexpr operator U () const {return U(this->at(0));}
		
	};
	
	template<class T, size_t N> struct value_type<VectorRef<T, N>>{using type = T;};
	
}