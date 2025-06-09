#pragma once


namespace control{
	
	// approximates the exponential of the scalar with the tailor expansion:
	// e^x = I + x + x^2 / 2 + ... + x^n / n!
	// the minimum number of n is 2. If n is set lower than 2, then 2 increments will be calculated regardless
	template<class T>
	constexpr T exp_taylor(const T& x, size_t n = 4){
		T x_pow = x * x;
		T factorial = T(2);
		T result = 1 + x + x_pow / factorial;
		for(size_t i = 3; i <= n; ++i){
			x_pow = x_pow * x;
			factorial *= i;
			result = result + x_pow / factorial;
		}
		return result;
	}
	
	// improves the accuracy of the taylor expansion using e^x = (e^(x/s))^s
	// accuracy is increased because the taylor expansion is more precise for smaller values of x.
	// e^x ~ taylor_n(x/(2^s))^(2^s)
	template<class T>
	constexpr T exp_taylor_scale(const T& x, size_t taylor_order = 4, size_t scaling = 4){
		T t = exp_taylor(x/(1 << scaling), taylor_order);
		for(size_t i = 0; i < scaling; ++i) t = t * t;
		return t;
	}

}