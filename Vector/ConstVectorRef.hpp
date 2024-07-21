#pragma once

#include <type_traits>

#include "BaseTypeArithmetic.hpp"
#include "definitions.hpp"
#include "VectorTraits.hpp"


namespace twmath{

	template<class T, size_t N>
	struct ConstVectorRef : public StaticVectorToken{
		const T* const _ptr = nullptr;
		const size_t _increment = 1;
		
	public:
	
		using ValueType = T;
	
		constexpr ConstVectorRef(T const * ptr, size_t increment=1)
			: _ptr(ptr)
			, _increment(increment){}
			
		constexpr ConstVectorRef& operator=(const T (&array)[N]){
			this->assign(array);
			return *this;
		}
		
		constexpr ConstVectorRef(const ConstVectorRef&) = default;
		
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
		
		template<class Integer, TWMATH_ENABLE_IF(std::is_integral_v<Integer>)>
		constexpr const T& operator[](Integer i) const {return this->at(i);}
		
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

		// ---------------------- conversion to scalar ---------------------------
		
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && N == 1))>
		constexpr operator const U& () const {return this->at(0);}
		
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && N == 1))>
		constexpr operator U () const {return U(this->at(0));}

	};
	
	template<class T, size_t N> struct value_type<ConstVectorRef<T, N>>{using type = T;};
	
}