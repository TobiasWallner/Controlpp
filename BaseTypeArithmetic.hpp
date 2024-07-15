#pragma once

#include <cstddef>
#include <cstdint>

#include <utility>
#include <cmath>


namespace twmath_base{	
	// # binary operators 
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() == std::declval<Tr>()) equal (const Tl& lhs, const Tr& rhs){return lhs == rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() != std::declval<Tr>()) not_equal (const Tl& lhs, const Tr& rhs){return lhs != rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() < std::declval<Tr>()) less (const Tl& lhs, const Tr& rhs){return lhs < rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() <= std::declval<Tr>()) less_equal (const Tl& lhs, const Tr& rhs){return lhs <= rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() > std::declval<Tr>()) greater (const Tl& lhs, const Tr& rhs){return lhs > rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() >= std::declval<Tr>()) greater_equal (const Tl& lhs, const Tr& rhs){return lhs >= rhs;}
	
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() + std::declval<Tr>()) add (const Tl& lhs, const Tr& rhs){return lhs + rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() - std::declval<Tr>()) sub (const Tl& lhs, const Tr& rhs){return lhs - rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() * std::declval<Tr>()) mul (const Tl& lhs, const Tr& rhs){return lhs * rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() / std::declval<Tr>()) div (const Tl& lhs, const Tr& rhs){return lhs / rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() % std::declval<Tr>()) rest (const Tl& lhs, const Tr& rhs){return lhs % rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() << std::declval<Tr>()) lshift (const Tl& lhs, const Tr& rhs){return lhs << rhs;}
	template<class Tl, class Tr> constexpr decltype(std::declval<Tl>() >> std::declval<Tr>()) rshift (const Tl& lhs, const Tr& rhs){return lhs >> rhs;}
	
	// # unary operators
	template<class T>constexpr decltype(-std::declval<T>()) negate (const T& a){return -a;}
	
	// ## basic operations
	template<class T>inline T abs(const T& n){using namespace std;return abs(n);}
	
	// ## exponential functions
	template<class T>inline T exp(const T& n){using namespace std;return exp(n);}
	template<class T>inline T exp2(const T& n){using namespace std;return exp2(n);}
	template<class T>inline T exp10(const T& n){using namespace std;return exp10(n);}
	template<class T>inline T expm1(const T& n){using namespace std;return expm1(n);}
	
	template<class T>inline T log(const T& n){using namespace std;return log(n);}
	template<class T>inline T log2(const T& n){using namespace std;return log2(n);}
	template<class T>inline T log10(const T& n){using namespace std;return log10(n);}
	template<class T>inline T logp1(const T& n){using namespace std;return logp1(n);}
	
	// ## trigonometric functions
	template<class T>inline T sin(const T& n){using namespace std;return sin(n);}
	template<class T>inline T cos(const T& n){using namespace std;return cos(n);}
	template<class T>inline T tan(const T& n){using namespace std;return tan(n);}
	
	template<class T>inline T asin(const T& n){using namespace std;return asin(n);}
	template<class T>inline T acos(const T& n){using namespace std;return acos(n);}
	template<class T>inline T atan(const T& n){using namespace std;return atan(n);}
	
	template<class T>inline T sinh(const T& n){using namespace std;return sinh(n);}
	template<class T>inline T cosh(const T& n){using namespace std;return cosh(n);}
	template<class T>inline T tanh(const T& n){using namespace std;return tanh(n);}
	
	// ## power
	template<class T>inline T square(const T& n){using namespace std;return n * n;}
	template<class T>inline T cube(const T& n){using namespace std;return n * n * n;}
	
	template<class T>inline T sqrt(const T& n){using namespace std;return sqrt(n);}
	template<class T>inline T cbrt(const T& n){using namespace std;return cbrt(n);}
	
	// ## rounding
	template<class T>inline T ceil(const T& n){using namespace std;return ceil(n);}
	template<class T>inline T floor(const T& n){using namespace std;return floor(n);}
	template<class T>inline T trunk(const T& n){using namespace std;return trunk(n);}
	template<class T>inline T round(const T& n){using namespace std;return round(n);}
	template<class T>inline int32_t rint(const T& n){using namespace std;return lrint(n);}
	template<class T>inline int64_t rint(const T& n){using namespace std;return llrint(n);}


	// ## classification
	template<class T>inline bool is_inf(const T& n){using namespace std;return isinf(n);}
	template<class T>inline bool is_nan(const T& n){using namespace std;return isnan(n);}
}