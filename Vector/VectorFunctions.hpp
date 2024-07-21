#pragma once

#include <type_traits>

#include "MatrixTraits.hpp"

namespace twmath{
	
	// ----------------------- for each -----------------------
	
	template<class Vl, class Vr, class F, 
		TWMATH_ENABLE_IF(is_static_vector_v<Vl> && is_static_vector_v<Vr> && (Vl::ssize() == Vr::ssize()))>
	constexpr Vector<F, Vl::ssize()> for_each (const Vl& lhs, const Vr& rhs, F (*f)(const value_type_t<Vl>&, const value_type_t<Vr>&)){
		Vector<F, Vl::ssize()> result;
		for(size_t i = 0; i < Vl::ssize(); ++i) result.at(i) = f(lhs.at(i), rhs.at(i));
		return result;
	}
	template<class Tl, class Vr, class F,
		TWMATH_ENABLE_IF(!is_vector_v<Tl> && !is_matrix_v<Tl> && is_static_vector_v<Vr>)>
	constexpr Vector<F, Vr::ssize()> for_each (const Tl& lhs, const Vr& rhs, F (*f)(const Tl&, const value_type_t<Vr>&)){
		Vector<F, Vr::ssize()> result;
		for(size_t i = 0; i < Vr::ssize(); ++i) result.at(i) = f(lhs, rhs.at(i));
		return result;
	}
	template<class Vl, class Tr, class F,
		TWMATH_ENABLE_IF(is_static_vector_v<Vl> && !is_vector_v<Tr> && !is_matrix_v<Tr>)>
	constexpr Vector<F, Vl::ssize()> for_each (const Vl& lhs, const Tr& rhs, F (*f)(const value_type_t<Vl>&, const Tr&)){
		Vector<F, Vl::ssize()> result;
		for(size_t i = 0; i < Vl::ssize(); ++i) result.at(i) = f(lhs.at(i), rhs);
		return result;
	}
	template<class V, class F,
		TWMATH_ENABLE_IF(is_static_vector_v<V>)>
	constexpr Vector<F, V::ssize()> for_each (const V& v, F (*f)(const value_type_t<V>&)){
		Vector<F, V::ssize()> result;
		for(size_t i = 0; i < V::ssize(); ++i) result.at(i) = f(v.at(i));
		return result;
	}
	
	// ======================= comparisons =======================
	// ----------------------- equal / operator== -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto equal (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::equal<value_type_t<Vl>, value_type_t<Vr>>);
	}
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto operator== (const Vl& lhs, const Vr& rhs){
		return equal(lhs, rhs);
	}
	
	// ----------------------- not_equal / operator!= -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto not_equal (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::not_equal<value_type_t<Vl>, value_type_t<Vr>>);
	}
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto operator!= (const Vl& lhs, const Vr& rhs){
		return not_equal(lhs, rhs);
	}
	
	// ----------------------- less / operator< -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto less (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::less<value_type_t<Vl>, value_type_t<Vr>>);
	}
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto operator< (const Vl& lhs, const Vr& rhs){
		return less(lhs, rhs);
	}
	
	// ----------------------- less_equal / operator<= -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto less_equal (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::less_equal<value_type_t<Vl>, value_type_t<Vr>>);}
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto operator<= (const Vl& lhs, const Vr& rhs){
		return less_equal(lhs, rhs);
	}
	
	// ----------------------- greater / operator> -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto greater (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::greater<value_type_t<Vl>, value_type_t<Vr>>);
	}
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto operator> (const Vl& lhs, const Vr& rhs){
		return greater(lhs, rhs);
	}
	
	// ----------------------- greater_equal / operator>= -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto greater_equal (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::greater_equal<value_type_t<Vl>, value_type_t<Vr>>);
	}
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto operator>= (const Vl& lhs, const Vr& rhs){
		return greater_equal(lhs, rhs);
	}
	
	// ======================= binary arithmetic operators =======================
	// ----------------------- add / operator+ -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto add (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::add<value_type_t<Vl>, value_type_t<Vr>>);
	}
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto operator+ (const Vl& lhs, const Vr& rhs){
		return add(lhs, rhs);
	}
	
	// ----------------------- sub / operator- -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto sub (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::sub<value_type_t<Vl>, value_type_t<Vr>>);
	}
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto operator- (const Vl& lhs, const Vr& rhs){
		return sub(lhs, rhs);
	}
	
	// ----------------------- mul / operator* -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto mul (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::mul<value_type_t<Vl>, value_type_t<Vr>>);
	}
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> != is_vector_v<Vr>) && !is_matrix_v<Vl> && !is_matrix_v<Vr>)> 
	constexpr auto operator* (const Vl& lhs, const Vr& rhs){return mul(lhs, rhs);}
	
	// ----------------------- div / operator/ -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> || is_vector_v<Vr>) && !(is_matrix_v<Vl> || is_matrix_v<Vr>))> 
	constexpr auto div (const Vl& lhs, const Vr& rhs){
		return for_each(lhs, rhs, twmath_base::div<value_type_t<Vl>, value_type_t<Vr>>);
	}
	template<class Vl, class Vr, TWMATH_ENABLE_IF((is_vector_v<Vl> != is_vector_v<Vr>) && !is_matrix_v<Vl> && !is_matrix_v<Vr>)> 
	constexpr auto operator/ (const Vl& lhs, const Vr& rhs){return div(lhs, rhs);}
	
	// ----------------------- negate -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	constexpr auto negate (const V& a){return for_each(a, twmath_base::negate<value_type_t<V>>);}
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)>
	constexpr auto operator- (const V& a){return for_each(a, twmath_base::negate<value_type_t<V>>);}
	
	// ----------------------- reduce -----------------------
	
	template<class V, class T=value_type_t<V>, TWMATH_ENABLE_IF(is_vector_v<V>)>
	constexpr value_type_t<V> reduce(const V& v, T (*function)(const T&, const value_type_t<V>&), T init = T()){
		for(size_t i = 0; i < v.size(); ++i) init = function(init, v.at(i));
		return init;
	}
	
	// ----------------------- sum -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	constexpr value_type_t<V> sum(const V& v){return reduce<V, value_type_t<V>>(v, twmath_base::add<value_type_t<V>>, value_type_t<V>());}
	
	// ----------------------- prod -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	constexpr value_type_t<V> prod(const V& v){return reduce<V, value_type_t<V>>(v, twmath_base::mul<value_type_t<V>>, value_type_t<V>(1));}	
	
	// ----------------------- min -----------------------
	
	template<class T, TWMATH_ENABLE_IF(!is_vector_v<T>)> 
	constexpr T min(const T& lhs, const T& rhs){return (lhs < rhs) ? lhs : rhs;}	
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	constexpr value_type_t<V> min(const V& v){return reduce<V, value_type_t<V>>(v, min, v.at(0));}
	
	// ----------------------- max -----------------------
	
	template<class T, TWMATH_ENABLE_IF(!is_vector_v<T>)>
	constexpr T max(const T& lhs, const T& rhs){return (lhs > rhs) ? lhs : rhs;}
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	constexpr value_type_t<V> max(const V& v){return reduce<V, value_type_t<V>>(v, max, v[0]);}
	
	// ----------------------- mean -----------------------
	
	template<class V, class R=float, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	constexpr R mean(const V& v){return static_cast<R>(sum(v))/static_cast<R>(v.size());}
	
	// ----------------------- square -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	constexpr auto square(const V& v){return mul(v, v);}
	
	// ----------------------- dot / operator* -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF(is_static_vector_v<Vl>), TWMATH_ENABLE_IF(is_static_vector_v<Vr>), TWMATH_ENABLE_IF(Vl::ssize() == Vr::ssize())>  
	constexpr auto dot(const Vl& lhs, const Vr& rhs){return sum(mul(lhs, rhs));}
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF(is_vector_v<Vl>), TWMATH_ENABLE_IF(is_vector_v<Vr>)>
	constexpr auto operator* (const Vl& lhs, const Vr& rhs){return dot(lhs, rhs);}
	
	// ----------------------- var -----------------------
	
	template<class V, class R=float, TWMATH_ENABLE_IF(is_vector_v<V>)>
	constexpr R var(const V& v){return static_cast<R>(sum(square(v - mean(v)))) / static_cast<R>(v.size());}
	
	// ----------------------- stddev -----------------------
	
	template<class V, class R=float, TWMATH_ENABLE_IF(is_vector_v<V>)>
	constexpr R stddev(const V& v){return sqrt(var<V, R>(v));}
	
	// ----------------------- sqr_norm -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)>
	constexpr auto sqr_norm(const V& v){return dot(v, v);}
	
	// ----------------------- pow -----------------------
	
	template<class V, class T, TWMATH_ENABLE_IF(is_vector_v<V>)>
	constexpr auto pow(const V& v, const T& p){return for_each(v, p, twmath_base::pow);}
	
	// ----------------------- norm -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)>
	constexpr auto norm(const V& v){return sqrt(sqr_norm(v));}
	
	template<class V, class Tp, TWMATH_ENABLE_IF(is_vector_v<V>)>
	constexpr auto norm(const V& v, const Tp& p){return pow(sum(pow(v, p)), 1./p);}
	
	// ----------------------- normalize -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)>
	constexpr auto normalize(const V& v){return v / norm(v);}
	
	// ----------------------- reverse -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_static_vector_v<V>)>
	constexpr Vector<value_type_t<V>, V::ssize()> reverse(const V& v){
		Vector<value_type_t<V>, V::ssize()> result;
		for (size_t i = 0, j = V::ssize()-1; i < V::ssize(); ++i, --j) result.at(i) = v.at(j);
		return result;
	}
	
	// ----------------------- rotate_r -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_static_vector_v<V>)>
	constexpr Vector<value_type_t<V>, V::ssize()> rotate_r(const V& v, size_t n){
		Vector<value_type_t<V>, V::ssize()> result;
		size_t i = 0;
		size_t j = n;
		for (; j < V::ssize(); ++i, (void)++j) result.at(j) = v.at(i);
		j = 0;
		for (; i < V::ssize(); ++i, (void)++j) result.at(j) = v.at(i);
		return result;
	}
	
	// ----------------------- rotate_l -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_static_vector_v<V>)>
	constexpr Vector<value_type<V>, V::ssize()> rotate_l(const V& v, size_t n){
		Vector<value_type<V>, V::ssize()> result;
		size_t i = 0;
		size_t j = V::ssize() - n;
		for (; j < V::ssize(); ++i, (void)++j) result.at(j) = v.at(i);
		j = 0;
		for (; i < V::ssize(); ++i, (void)++j) result.at(j) = v.at(i);
		return result;
	}
	
	// ----------------------- cross -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF(is_static_vector_v<Vl>), TWMATH_ENABLE_IF(is_static_vector_v<Vr>), TWMATH_ENABLE_IF(Vl::ssize()==3), TWMATH_ENABLE_IF(Vr::ssize()==3)>
	constexpr Vector<decltype(std::declval<value_type_t<Vl>>() * std::declval<value_type_t<Vr>>()), 3> cross(const Vl& lhs, const Vr& rhs){
		Vector<decltype(std::declval<value_type_t<Vl>>() * std::declval<value_type_t<Vr>>()), 3> result({
			lhs.at(1) * rhs.at(2) - lhs.at(2) * rhs.at(1),
			lhs.at(2) * rhs.at(0) - lhs.at(0) * rhs.at(2),
			lhs.at(0) * rhs.at(1) - lhs.at(1) * rhs.at(0),
		});
		return result;
	}
	
	// ----------------------- linspace -----------------------
	
	template<class T, size_t N>
	constexpr Vector<T, N> linspace(const T& start, const T& end){
		Vector<T, N> result;
		T step = (end - start) / (N + 1);
		for(size_t i = 0; i < N; ++i)
			result[i] = start + step * i;
		return result;
	}
	
	// ----------------------- fill -----------------------
	
	template<class T, size_t N>
	constexpr Vector<T, N> fill(const T& n){
		Vector<T, N> result;
		for(size_t i = 0; i < N; ++i) result.at(i) = n;
		return result;
	}
	
	// ----------------------- generators -----------------------
	
	template<class T, size_t N>
	constexpr Vector<T, N> zeros(){return fill<T, N>(0);}
	
	template<class T, size_t N>
	constexpr Vector<T, N> ones(){return fill<T, N>(1);}
	
	template<class T, size_t N>
	constexpr Vector<T, N> logspace(const T& start, const T& end, size_t num){
		using namespace std;
		return for_each(linspace(log2(start), log2(end), num), exp2);
	}
	
	// ----------------------- standart functions ---------------------
	// ## basic operations
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto abs(const V& v){return for_each(v, twmath_base::abs);}
	
	// ## exponential functions
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto exp(const V& v){return for_each(v, twmath_base::exp);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto exp2(const V& v){return for_each(v, twmath_base::exp2);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto exp10(const V& v){return for_each(v, twmath_base::exp10);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto expm1(const V& v){return for_each(v, twmath_base::expm1);}
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto log(const V& v){return for_each(v, twmath_base::log);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto log2(const V& v){return for_each(v, twmath_base::log2);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto log10(const V& v){return for_each(v, twmath_base::log10);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto logp1(const V& v){return for_each(v, twmath_base::logp1);}
	
	// ## trigonometric functions
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto sin(const V& v){return for_each(v, twmath_base::sin);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto cos(const V& v){return for_each(v, twmath_base::cos);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto tan(const V& v){return for_each(v, twmath_base::tan);}
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto asin(const V& v){return for_each(v, twmath_base::asin);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto acos(const V& v){return for_each(v, twmath_base::acos);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto atan(const V& v){return for_each(v, twmath_base::atan);}
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto sinh(const V& v){return for_each(v, twmath_base::sinh);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto cosh(const V& v){return for_each(v, twmath_base::cosh);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto tanh(const V& v){return for_each(v, twmath_base::tanh);}
	
	// ## power
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto cube(const V& v){ return for_each(v, twmath_base::cube);}
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto sqrt(const V& v){return for_each(v, twmath_base::sqrt);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto cbrt(const V& v){return for_each(v, twmath_base::cbrt);}
	
	// ## rounding
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto ceil(const V& v){return for_each(v, twmath_base::ceil);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto floor(const V& v){return for_each(v, twmath_base::floor);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto trunk(const V& v){return for_each(v, twmath_base::trunk);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto round(const V& v){return for_each(v, twmath_base::round);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto rint(const V& v){return for_each(v, twmath_base::rint);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto llrint(const V& v){return for_each(v, twmath_base::llrint);}


	// ## classification
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline bool is_inf(const V& v){return for_each(v, twmath_base::is_inf);}
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline bool is_nan(const V& v){return for_each(v, twmath_base::is_nan);}
	
	// ----------------------- softmax -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_vector_v<V>)> 
	inline auto softmax(const V& v){ return normalize(exp(v));}
	
	// ----------------------- conv -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF(is_static_vector_v<Vl>), TWMATH_ENABLE_IF(is_static_vector_v<Vr>)>
	constexpr Vector<decltype(std::declval<value_type_t<Vl>>() * std::declval<value_type_t<Vr>>()), Vl::ssize() + Vr::ssize() - 1> conv (const Vl& lhs, const Vr& rhs){
		using ValueType = decltype(std::declval<value_type_t<Vl>>() * std::declval<value_type_t<Vr>>());
		Vector<ValueType, Vl::ssize() + Vr::ssize() - 1> result;
		size_t resi = 0;
		size_t lhs_first = 0;
		size_t rhs_first = 0;
		size_t lhs_last = 1;
		size_t rhs_last = 1;
		for(; resi < Vl::ssize() + Vr::ssize() - 1; ++resi){
			size_t lhsi = lhs_last-1;
			size_t rhsi = rhs_first;
			ValueType sum(0);
			for(; rhsi < rhs_last; ++rhsi, --lhsi) sum += lhs.at(lhsi) * rhs.at(rhsi);
			result.at(resi) = sum;
			
			rhs_first += (lhs_last == Vl::ssize());
			lhs_first += (rhs_last == Vr::ssize());
			lhs_last += (lhs_last != Vl::ssize());
			rhs_last += (rhs_last != Vr::ssize());
			
		}
		return result;
	}
	
	// ----------------------- corr -----------------------
	
	template<class Vl, class Vr, TWMATH_ENABLE_IF(is_static_vector_v<Vl>), TWMATH_ENABLE_IF(is_static_vector_v<Vr>)>
	constexpr Vector<decltype(std::declval<value_type_t<Vl>>() * std::declval<value_type_t<Vr>>()), Vl::ssize() + Vr::ssize() - 1> corr (const Vl& lhs, const Vr& rhs){
		using ValueType = decltype(std::declval<value_type_t<Vl>>() * std::declval<value_type_t<Vr>>());
		Vector<ValueType, Vl::ssize() + Vr::ssize() - 1> result;
		size_t resi = 0;
		size_t lhs_first = 0;
		size_t rhs_first = 0;
		size_t lhs_last = 1;
		size_t rhs_last = 1;
		for(; resi < Vl::ssize() + Vr::ssize() - 1; ++resi){
			size_t lhsi = lhs_last-1;
			size_t rhsi = rhs_first;
			ValueType sum(0);
			for(; rhsi < rhs_last; ++rhsi, --lhsi){
				sum += lhs[Vl::ssize()-1-lhsi] * rhs[rhsi];
			}
			result.at(resi) = sum;
			
			rhs_first += (lhs_last == Vl::ssize());
			lhs_first += (rhs_last == Vr::ssize());
			lhs_last += (lhs_last != Vl::ssize());
			rhs_last += (rhs_last != Vr::ssize());
			
		}
		return result;
	}
	
	// ----------------------- print / operator << -----------------------
	
	template<class Stream, class V, TWMATH_ENABLE_IF(is_vector_v<V>)>
	Stream& print (Stream& stream, const V& v, const char* seperator = ", "){
		if(v.size() != 0) stream << v.at(0);
		for(size_t i = 1; i < v.size(); ++i) stream << seperator << v.at(i);
		return stream;
	}
	
	template<class Stream, class V, TWMATH_ENABLE_IF(is_vector_v<V>)>
	Stream& operator << (Stream& stream, const V& v){return print(stream, v, ", ");}
	
}