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

	/**
	 * \brief returns the minor matrix excluding the provided column and row
	 * \param A The source matrix
	 * \param ex_row The row to be excluded
	 * \param ex_col The column to be excluded
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, Rows-1, Cols-1> minor(const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& A, size_t ex_row, size_t ex_col){
		Eigen::Matrix<T, Rows-1, Cols-1> result;
		size_t A_row = 0;
		size_t A_col = 0;
		for(size_t result_row = 0; result_row < (Rows-1); ++result_row, (void)++A_row){
			for(size_t result_col = 0; result_col < (Cols-1); ++result_col, (void)++A_col){
				A_row += (A_row == ex_row) ? 1 : 0;
				A_col += (A_col == ex_col) ? 1 : 0;
				result(result_row, result_col) = A(A_row, A_col);
			}
		}
		return result;
	}

	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, Rows, Cols> adj(const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& A){
		Eigen::Matrix<T, Rows, Cols> result;
		for(int row = 0; row < Rows; ++row){
			const bool row_sign = (row & 1 == 1); // even are positive, odd are negative
			for(int col = 0; col < Cols; ++col){
				const bool col_sign = (col & 1 == 1); // even are positive, odd are negative
				const bool value_sign = row_sign != col_sign; // sign of the value (both positive/negative: value is positive) (different signs: value is negative)
				const auto min = minor(A, row, col);
				const T det = min.determinant();
				if(value_sign){
					// negative sign
					result(col, row) = -det;
				}else{
					// positive sign
					result(col, row) = det;
				}
			}
		}
		return result;
	}

	/**
	 * \brief Creates a companion matrix from a vector
	 * 
	 * given the vector:
	 * 
	 * \f[
	 * \vec{v} = \left[a_0, a_1, a_2, \hdots, a_n \right]
	 * \f]
	 * 
	 * creates the companion matrix:
	 * 
	 * \f[
	 * \mathbf{C} = \begin{bmatrix}
	 * 		0 		& 1 		& 0 		& \cdots	& 0 		\\
	 * 		0 		& 0 		& 1 		& \cdots 	& 0 		\\
	 * 		\vdots 	& \vdots 	& \ddots 	& \ddots 	& \vdots 	\\
	 * 		0 		& 0 		& \cdots 	& 0 		& 1 		\\
	 * 		b_0 	& b_1 		& \cdots 	& b_{n-2} 	& b_{n-1}
	 * \end{bmatrix}
	 * \f]
	 * 
	 * where \f$b_i\f$ is:
	 * 
	 * \f[
	 * b_i = \frac{a_i}{a_n}
	 * \f]
	 * 
	 * \returns a companion matrix
	 */
	template<class T, int N>
	Eigen::Matrix<T, N-1, N-1> companion(const Eigen::Vector<T, N>& v){
		Eigen::Matrix<T, N-1, N-1> result;
		result.template block<N-2, N-2>(0, 1) = Eigen::Matrix<T, N-2, N-2>::Identity();
		result.col(0).setZero();
		result.row(N-2) = -v.head(N-1) / v(N-1);
		return result;
	}

	template<class T, int N>
	Eigen::Vector<T, N> real(const Eigen::Vector<std::complex<T>, N>& v){
		Eigen::Vector<T, N> result;
		for(int i = 0; i < N; ++i){
			result(i) = std::real(v(i));
		}
		return result;
	}

	template<class T, int N>
	Eigen::Vector<T, N> imag(const Eigen::Vector<std::complex<T>, N>& v){
		Eigen::Vector<T, N> result;
		for(int i = 0; i < N; ++i){
			result(i) = std::imag(v(i));
		}
		return result;
	}

}