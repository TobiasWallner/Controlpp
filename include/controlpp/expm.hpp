#pragma once

/**
 * @file expm.hpp   
 * @brief Includes functions for matrix exponentiation
 * 
 */


#include "math.hpp"

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
	 * \see controlpp::expm_taylor_scale
	 * \see controlpp::expm
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> expm_taylor(
			const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& x, 
			int n
	){
		Eigen::Matrix<T, Rows, Cols> result;
		const auto I = Eigen::Matrix<T, Rows, Cols>::Identity();
		result.setZero();
		for(int i = n; i > 0; --i){
			Eigen::Matrix<T, Rows, Cols> new_result = (result + I) * x / static_cast<T>(i);
			result = new_result;
		}
		result += I;
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
	 * \see controlpp::expm_taylor
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> expm_taylor_scaled(
			const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& M, 
			int taylor_order = 8, 
			int scaling = 10
	){
		using Matrix = Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>;
		T s = static_cast<T>(1 << scaling);
		Matrix scaled_M = M/s;
		Matrix t = expm_taylor(scaled_M, taylor_order);
		for(int i = 0; i < scaling; ++i){
			const Matrix temp = t * t;
			t = temp;
		};
		return t;
	}


	/**
	 * @brief Approximates \f$\exp{A}\f$ using a pade fraction
	 * 
	 * @throws `std::invalid_argument` if the order is greater than 31
	 * 
	 * @tparam T The value type of the matrix entries (usually `float` or `double`)
	 * @tparam Rows The number of rows
	 * @tparam Cols The number of columns
	 * @tparam Options Matrix options (See: [Store Orders](https://libeigen.gitlab.io/eigen/docs-nightly/group__TopicStorageOrders.html))
	 * @param A The input matrix
	 * @param Order The order of the numerator and denominator of the pade fraction
	 * @return An approximation of the exponential \f$\exp{A}\f$
	 */
	template<class T, int N, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, N, N> expm_pade(
			const Eigen::Matrix<T, N, N, Options, MaxRows, MaxCols>& A,
			int Order = 5
	){
		if(Order >= 32){
			throw std::invalid_argument("expm_pade: Order must be less than 32");
		}
		const Eigen::Matrix<T, N, N> I = Eigen::Matrix<T, N, N>::Identity();
		const Eigen::Matrix<T, N, N> A2 = A * A;

		// even odd split of polynomials: allows to reuse results
		Eigen::Matrix<T, N, N> even;
		Eigen::Matrix<T, N, N> odd;
		even.setZero();
		odd.setZero();

		// calculate pade parameters
		T buffer[32];
		pade_params(buffer, Order+1, Order, Order);

		// horner chains to evaluate the polynomials
		for(int k = Order; k >= 0; --k){
			const bool is_first_iteration = (k == Order);
			const T pk = static_cast<T>(buffer[k]);
			if((k & 1) == 0){
				// even
				if(!is_first_iteration) even *= A2; // skipp on the first iteration
				even += I * pk;
			}else{
				// odd
				if(!is_first_iteration) odd *= A2; // skipp on the first iteration
				odd += I * pk;
			}
		}
		odd *= A;

		const Eigen::Matrix<T, N, N> num = even + odd;
		const Eigen::Matrix<T, N, N> den = even - odd;

		const Eigen::Matrix<T, N, N> result = den.partialPivLu().solve(num);
		
		return result;
	}

	/**
	 * @brief Applies scaling and squaring to the pade approximation of the matrix exponential
	 * @tparam T The value type of the matrix entries (usually `float` or `double`)
	 * @tparam N The size of the square matrix
	 * @tparam Options The storage order of the matrix (see Eigen documentation)
	 * @tparam MaxRows The maximum number of rows for the matrix (default is `N`)
	 * @tparam MaxCols The maximum number of columns for the matrix (default is `N`)
	 * @param M Matrix to exponentiate
	 * @param order Order of the Pade approximation (default is 7)
	 * @param scaling Scaling factor for the matrix (default is 7)
	 * @return The exponentiated matrix \f$\exp{M}\f$ using the scaled Pade approximation
	 */
	template<class T, int N, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, N, N> expm_pade_scaled(
			const Eigen::Matrix<T, N, N, Options, MaxRows, MaxCols>& M, 
			int order,
			int scaling
	){
		Eigen::Matrix<T, N, N> scaled_M = M * std::ldexp(T(1), -scaling);
		Eigen::Matrix<T, N, N> t = expm_pade(scaled_M, order);
		for(int i = 0; i < scaling; ++i){
			const Eigen::Matrix<T, N, N> temp = t * t;
			t = temp;
		};
		return t;
	}

	/**
	 * @brief Pade for the scaled pade matrix exponential
	 */
	struct ExpmPadeParams{
		int order; ///< The order of the pade approximation
		int scaling; ///< The scaling and squaring factor for the matrix
	};

	/**
	 * @brief Determines the order and scaling factor for the scaled pade approximation of the matrix exponential
	 * 
	 * Higham, Nicholas J., and Desmond J. Higham. "A new scaling and squaring algorithm for the matrix exponential." 
	 * 
	 * @tparam T The value type of the matrix entries (usually `float` or `double`)
	 * @tparam MaxRows The maximum number of rows for the matrix (default is `N`)
	 * @tparam MaxCols The maximum number of columns for the matrix (default is `N`)
	 * @tparam Options The storage order of the matrix (see Eigen documentation)
	 * @tparam N The size of the square matrix
	 * @param M The matrix for which to determine the parameters
	 * @return The parameters for the scaled pade approximation, including the order and scaling factor
	 */
	template<class T, int N, int Options, int MaxRows, int MaxCols>
	ExpmPadeParams expm_pade_params(const Eigen::Matrix<T, N, N, Options, MaxRows, MaxCols>& M){
		const T norm = M.cwiseAbs().colwise().sum().maxCoeff();

		// Theta values from Higham
		constexpr T theta3  = T(1.495585217958292e-2);
		constexpr T theta5  = T(2.539398330063230e-1);
		constexpr T theta7  = T(9.504178996162932e-1);
		constexpr T theta9  = T(2.097847961257068);
		constexpr T theta13 = T(5.371920351148152);

		if(norm <= theta3){
			return {3, 0};
		}
		
		if(norm <= theta5){
			return {5, 0};
		}
		
		if(norm <= theta7){
			return {7, 0};
		}
		
		if(norm <= theta9){
			return {9, 0};
		}
		
		if(norm <= theta13){
			return {13, 0};
		}
		
		const double ratio = static_cast<double>(norm / theta13);
		int scaling = 0;
		if (ratio > 1.0) {
			scaling = static_cast<int>(std::ceil(std::log2(ratio)));
			if (scaling < 0) scaling = 0;
		}
		return {13, scaling};
	
	}

	/**
	 * @brief Matrix exponential using a scaled Pade approximation with automatic order and scaling factor determination
	 * @tparam T The value type of the matrix entries (usually `float` or `double`)
	 * @tparam N The size of the square matrix
	 * @tparam Options The storage order of the matrix (see Eigen documentation)
	 * @tparam MaxRows The maximum number of rows for the matrix (default is `N`)
	 * @tparam MaxCols The maximum number of columns for the matrix (default is `N`)
	 * @param M The matrix to exponentiate
	 * @return The exponentiated matrix \f$\exp{M}\f$ using the scaled Pade approximation with automatically determined order and scaling factor
	 */
	template<class T, int N, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, N, N> expm_pade_scaled(
			const Eigen::Matrix<T, N, N, Options, MaxRows, MaxCols>& M
	){
		const ExpmPadeParams params = expm_pade_params(M);
		return expm_pade_scaled(M, params.order, params.scaling);
	}

	/**
	 * \brief Calculates the matrix exponent \f$ \exp{\mathbf{\M}} \f$
	 * 
	 * Uses a scaled pade approximation for the exponential.
	 * 
	 * The order and scaling factor are automatically determined based on the norm of the matrix.
	 * 
	 * \tparam T The value type of the matrix elements
	 * \tparam Rows The number of rows of the matrix
	 * \tparam Cols The number of columns of the matrix
	 * 
	 * \param M The matrix
	 * 
	 * \return The resulting exponentiated matrix
	 * 
	 * \see controlpp::expm_pade
	 * \see controlpp::expm_pade_scaled
	 */
	template<class T, int N, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, N, N> expm(const Eigen::Matrix<T, N, N, Options, MaxRows, MaxCols>& M){
		// actual exponent calculation
		return expm_pade_scaled(M);
	}
}