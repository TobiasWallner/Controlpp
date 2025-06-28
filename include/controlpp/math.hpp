#pragma once

#include <Eigen/Core>

namespace controlpp{
	
	/**
	 * \brief Exponential function with a taylor approximation
	 * 
	 * \f[
	 * 	\exp{x} = I + x + x^2 / 2 + ... + x^n / n!
	 * \f]
	 * 
	 * The minimum number of n is 2. If n is set lower than 2, then 2 increments will be calculated regardless
	 * 
	 * \tparam T The value type
	 * \param x The value taken to the exponent
	 * \param n The order of the taylor approximation (default: 3)
	 * 
	 * \result the exponentiated value
	 * 
	 * \see controlpp::exp_taylor_scale
	 * \see controlpp::mexp
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	constexpr Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> exp_taylor(const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& x, size_t n){
		using Matrix = Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>;
		const Matrix I = Matrix::Identity();
		Matrix x_pow = x * x;
		T factorial = T(2);
		Matrix result = I + x + x_pow / factorial;
		for(size_t i = 3; i <= n; ++i){
			x_pow = x_pow * x;
			factorial *= i;
			result = result + x_pow / factorial;
		}
		return result;
	}

	/**
	 * \brief Calculates the matrix exponent \f$ \exp{\mathbf{\M}} \f$
	 * 
	 * Uses a scaled taylor approximation for the exponential.
	 * 
	 * Uses the following relationship:
	 * 
	 * \f[
	 * 	exp{x} = exp{x/s*s} = \left( exp{x/s} \right) ^ {s}
	 * \f]
	 * 
	 *  to improve accuracy, by scaling the value first, 
	 * allowing for smaller taylor orders with increased accuracy.
	 * 
	 * \tparam T The value type of the matrix elements
	 * \tparam Rows The number of rows of the matrix
	 * \tparam Cols The number of columns of the matrix
	 * 
	 * \param M The matrix
	 * \param taylor_order The order of the taylor polynomial used to approximate the exponential function 
	 * \param scaling The scaling factor used to improve the accuracy of the exponential function
	 * 
	 * \return The resulting exponentiated matrix
	 * 
	 * \see controlpp::exp_taylor
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	constexpr Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> exp_taylor_scaled(const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& M, size_t taylor_order = 4, size_t scaling = 10){
		using Matrix = Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>;
		T s = static_cast<T>(1 << scaling);
		Matrix scaled_M = M/s;
		Matrix t = exp_taylor(scaled_M, taylor_order);
		for(size_t i = 0; i < scaling; ++i) t = t * t;
		return t;
	}

	/**
	 * \brief Calculates the matrix exponent \f$ \exp{\mathbf{\M}} \f$
	 * 
	 * Uses a scaled taylor approximation for the exponential.
	 * 
	 * \tparam T The value type of the matrix elements
	 * \tparam Rows The number of rows of the matrix
	 * \tparam Cols The number of columns of the matrix
	 * 
	 * \param M The matrix
	 * \param taylor_order The order of the taylor polynomial used to approximate the exponential function 
	 * \param scaling The scaling factor used to improve the accuracy of the exponential function
	 * 
	 * \return The resulting exponentiated matrix
	 * 
	 * \see controlpp::exp_taylor
	 * \see controlpp::exp_taylor_scaled
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	constexpr Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> mexp(const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& M, size_t taylor_order = 4){
		// use the absolute maxima as a (fast) normation factor to calculate the scaling
		size_t norm = M.cwiseAbs().maxCoeff();

		// a fast log2 to get the s^n scaling factor
		size_t n = std::bit_width(norm);

		// actual exponent calculation
		return exp_taylor_scaled(M, taylor_order, n);
	}

	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	constexpr Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> identity_like([[maybe_unused]]const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& unused){
		return Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>::Identity();
	}

	template<class T, int LSize, int RSize>
	Eigen::Matrix<T, LSize + RSize, LSize + RSize> join_to_diagonal(const Eigen::Vector<T, LSize>& l, const Eigen::Vector<T, RSize>& r){
		Eigen::Matrix<T, LSize + RSize, LSize + RSize> result;
		result.diagonal().head(LSize) = l;
		result.diagonal().tail(RSize) = r;
		return result;
	}

	template<class T, int LSize, int RSize>
	Eigen::Vector<T, LSize + RSize> join_to_vector(const Eigen::Vector<T, LSize>& l, const Eigen::Vector<T, RSize>& r){
		Eigen::Vector<T, LSize + RSize> result;
		result.head(LSize) = l;
		result.tail(RSize) = r;
		return result;
	}
}