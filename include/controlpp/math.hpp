#pragma once

#include <Eigen/Core>
#include <Eigen/Eigenvalues> 

namespace controlpp{
	
	/**
	 * \brief Calculates the product of all numbers in the closed open range [from, to)
	 * \returns An integer
	 */
	constexpr unsigned long product_over(unsigned long from, unsigned long to){
		unsigned long product = 1;
		for(; from < to; ++from){
			product *= from;
		}
		return product;
	}

	/**
	 * \brief Power function for integral exponents
	 */
	template<class T>
	constexpr T pow(const T& base, const int& exp){
		T result = static_cast<T>(1);
		if(exp >= 0){
			for(int i = 0; i < exp; ++i){
				result *= base;
			}
		}else{
			for(int i = 0; i < -exp; ++i){
				result /= base;
			}
		}
		return result;
	}

	/**
	 * \brief Unwinds modulo jumps
	 * 
	 * Checks wheather jumps in occur that are larger than the threshold,
	 * and if so the modulo is added or subtracted n times so that the step \f$y_{k+1} - y_{k}\f$ is within the bounds [-threshold, +threshold].
	 * 
	 * \param y The vector to be unwrapped
	 * \param modulo The modulo to unwrap. Will be added/subtracted n times to minimize variance.
	 * \returns An eigen vector, the size of the input with the unwrapped values
	 */
	template<class T, int N>
	Eigen::Vector<T, N> unwrap(const Eigen::Vector<T, N>& y, const T& modulo){
		Eigen::Vector<T, N> result;
		if constexpr (N == Eigen::Dynamic) result.resize(y.size());
		T offset = 0;
		result(0) = y(0);
		for(int i = 1; i < y.size(); ++i){
			while(((y(i) + offset) - result(i-1)) > (modulo/2)){
				offset -= modulo;
			}
			while(((y(i) + offset) - result(i-1)) < (-modulo/2)){
				offset += modulo;
			}
			result(i) = y(i) + offset;
		}
		return result;
	}

	/**
	 * \brief Unwinds phase jumps of \f$2 \pi\f$ in radiants
	 * \param phase A vector op phases (rad)
	 * \returns A vector with the corrected phases (rad)
	 * \see unwrap
	 * \see unwrap_deg
	 */
	template<class T, int N>
	Eigen::Vector<T, N> unwrap_rad(const Eigen::Vector<T, N>& phases){
		const T modulo = std::numbers::pi_v<T> * static_cast<T>(2);
		return unwrap(phases, modulo);
	}

	/**
	 * \brief Unwinds phase jumps of \f$2 \pi\f$ in degrees
	 * \param phase A vector op phases (deg)
	 * \returns A vector with the corrected phases (deg)
	 * \see unwrap
	 * \see unwrap_rad
	 */
	template<class T, int N>
	Eigen::Vector<T, N> unwrap_deg(const Eigen::Vector<T, N>& phases){
		const T modulo = 360;
		return unwrap(phases, modulo);
	}

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
	constexpr Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> exp_taylor(
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
	 * \see controlpp::exp_taylor
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	constexpr Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> exp_taylor_scaled(
			const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& M, 
			int taylor_order = 8, 
			int scaling = 10
	){
		using Matrix = Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>;
		T s = static_cast<T>(1 << scaling);
		Matrix scaled_M = M/s;
		Matrix t = exp_taylor(scaled_M, taylor_order);
		for(int i = 0; i < scaling; ++i) t = t * t;
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
	constexpr Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> mexp(const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& M, unsigned int taylor_order = 8){
		// use the absolute maxima as a (fast) normation factor to calculate the scaling
		int norm = M.cwiseAbs().maxCoeff();

		// a fast log2 to get the s^n scaling factor
		int n = std::bit_width(static_cast<unsigned int>(norm))+1;

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
			const bool row_sign = ((row & 1) == 1); // even are positive, odd are negative
			for(int col = 0; col < Cols; ++col){
				const bool col_sign = ((col & 1) == 1); // even are positive, odd are negative
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

	/**
	 * \brief Solves the continuous time Lyapunov equation
	 * 
	 * This function computes the stabilizing symmetric solution of the lyapunov equation:
     * 
	 * \f[
	 * A^\top X + X A + Q = 0
	 * \f]
	 * 
	 * where A and Q are parameter and X is the matrix being solved for.
	 * 
	 * ----
	 * 
	 * Soves the Lyapunov equation by:
	 * 
	 * 1. Building the Hamilton matrix
	 * 2. Compute its stable eigenvectors
	 * 3. Re-Partitions the eigenvectors
	 * 4. Recovers X from the partitions
	 * 
	 * \param A State matrix (\f$n \times n\f$).
	 * \param Q State weighting matrix (\f$n \times n\f$, symmetric positive semidefinite).
	 *
	 * \tparam T Scalar type (e.g., `double`, `float`).
	 * \tparam NStates Dimension of the matrices
	 * 
	 * \returns X the solution to the Lyapunov equation as an Eigen::Matrix with the dimensions `NStates x NStates`.
	 */
	template<class T, int NStates,
			int AOpt, int AMaxR, int AMaxC,
			int QOpt, int QMaxR, int QMaxC
  	>
	constexpr Eigen::Matrix<T, NStates, NStates> lyapunov_solver(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NStates, QOpt, QMaxR, QMaxC>& Q
	){
		// 1. Build the hamilton matrix
		Eigen::Matrix<T, 2*NStates, 2*NStates> H;
		H.topLeftCorner(NStates, NStates) = A;
		H.topRightCorner(NStates, NStates).setZero();
		H.bottomLeftCorner(NStates, NStates) = - Q;
		H.bottomRightCorner(NStates, NStates) = -A.transpose();

		// 2. Compute stable eigenvectors
		Eigen::ComplexEigenSolver<Eigen::Matrix<T, 2*NStates, 2*NStates>> ces;
		ces.compute(H);
		const auto& H_eigvals = ces.eigenvalues();
		const auto& H_eigvecs = ces.eigenvectors();

		// in a 2n hamilton matrix are exactly n stable ones
		Eigen::Matrix<std::complex<T>, 2*NStates, NStates> StableEigenVecs;
		int si = 0; // stable eigenvalue iterator
		int ei = 0; // eigen value iterator
		for(; (si < NStates) && (ei < 2*NStates); ++ei){
			if(H_eigvals(ei).real() < static_cast<T>(0)){// only stable ones
				StableEigenVecs.col(si) = H_eigvecs.col(ei);
				++si;
			}
		}

		// 3. Repartition
		const auto TopPartition = StableEigenVecs.template block<NStates, NStates>(0, 0);
		const auto BottomPartition = StableEigenVecs.template block<NStates, NStates>(NStates, 0);

		// 4. Recover Solution (BottomPartition * TopPartition^-1)
		const Eigen::Matrix<std::complex<T>, NStates, NStates> X = TopPartition.transpose().partialPivLu().solve(BottomPartition.transpose()).transpose();

		// The result is real, make sure it is real because it should be
		const Eigen::Matrix<T, NStates, NStates> realX = X.real();

        // force the result X to be symetric to combat small numerical errors
        const Eigen::Matrix<T, NStates, NStates> result = static_cast<T>(0.5) * (realX + realX.transpose());
		return result;
	}

	/**
	 * \brief Solves the continuous time riccati equation (CARE)
	 * 
	 * This function computes the stabilizing symmetric solution of the CARE:
     * 
	 * \f[
	 * A^\top X + X A - (X B + N) R^{-1} (B^\top X + N) + Q = 0
	 * \f]
	 * 
	 * where \f$A, B, C, D\f$ are system matrices: 
	 * 
	 * \f[
	 * \dot{x} = A x + B u
	 * y = C x + D u
	 * \f]
	 * 
	 * with the system states \f$x\f$, inputs \f$u\f$ and outputs \f$y\f$,
	 * 
	 * as well as \f$Q\f$, \f$R\f$ the state and input weighting matrices.
	 * 
	 * ----
	 * 
	 * Soves the Riccatiy equation by:
	 * 
	 * 1. Building the Hamilton matrix (`controlpp::create_hamilton()`)
	 * 2. Compute its stable eigenvectors
	 * 3. Re-Partitions the eigenvectors
	 * 4. Recovers X from the partitions
	 * 
	 * \param A State matrix (\f$n \times n\f$).
	 * \param B Input matrix (\f$n \times m\f$).
	 * \param R Input weighting matrix (\f$m \times m\f$, symmetric positive definite).
	 * \param Q State weighting matrix (\f$n \times n\f$, symmetric positive semidefinite).
	 * \param N N Cross-term weighting matrix (\f$n × m\f$).
	 *
	 * \tparam T Scalar type (e.g., `double`, `float`).
	 * \tparam NStates Number of states.
	 * \tparam NInputs Number of control inputs.
	 * \tparam NOutputs Number of outputs.
	 * 
	 * @see Zhou, Doyle, and Glover (1996), *Robust and Optimal Control*.
	 */
	template<class T, int NStates, int NInputs,
			int AOpt, int AMaxR, int AMaxC,
			int BOpt, int BMaxR, int BMaxC,
			int ROpt, int RMaxR, int RMaxC,
			int QOpt, int QMaxR, int QMaxC,
			int NOpt, int NMaxR, int NMaxC
  	>
	constexpr Eigen::Matrix<T, NStates, NStates> care_solver(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NInputs, BOpt, BMaxR, BMaxC>& B,
		const Eigen::Matrix<T, NInputs, NInputs, ROpt, RMaxR, RMaxC>& R,
		const Eigen::Matrix<T, NStates, NStates, QOpt, QMaxR, QMaxC>& Q,
		const Eigen::Matrix<T, NStates, NInputs, NOpt, NMaxR, NMaxC>& N
	){
		// 1. Build the hamilton matrix
        const Eigen::Matrix<T, NInputs, NStates> R_inv_DT_C = R.ldlt().solve((N.transpose()));
		const Eigen::Matrix<T, NStates, NStates> A_ = A - B * R_inv_DT_C;

		Eigen::Matrix<T, 2*NStates, 2*NStates> H;
		H.topLeftCorner(NStates, NStates) = A_;
		H.topRightCorner(NStates, NStates) =  - B * R.ldlt().solve(B.transpose());
		H.bottomLeftCorner(NStates, NStates) = - (Q - N * R_inv_DT_C);
		H.bottomRightCorner(NStates, NStates) = -A_.transpose();


		// 2. Compute stable eigenvectors
		Eigen::ComplexEigenSolver<Eigen::Matrix<T, 2*NStates, 2*NStates>> ces;
		ces.compute(H);
		const auto& H_eigvals = ces.eigenvalues();
		const auto& H_eigvecs = ces.eigenvectors();

		// in a 2n hamilton matrix are exactly n stable ones
		Eigen::Matrix<std::complex<T>, 2*NStates, NStates> StableEigenVecs;
		int si = 0; // stable eigenvalue iterator
		int ei = 0; // eigen value iterator
		for(; (si < NStates) && (ei < 2*NStates); ++ei){
			if(H_eigvals(ei).real() <= 0){// only stable ones
				StableEigenVecs.col(si) = H_eigvecs.col(ei);
				++si;
			}
		}

		// 3. Repartition
		const auto TopPartition = StableEigenVecs.template block<NStates, NStates>(0, 0);
		const auto BottomPartition = StableEigenVecs.template block<NStates, NStates>(NStates, 0);

		// 4. Recover Solution (BottomPartition * TopPartition^-1)
		const Eigen::Matrix<std::complex<T>, NStates, NStates> X = TopPartition.transpose().partialPivLu().solve(BottomPartition.transpose()).transpose();

		// The result is real, make sure it is real because it should be
		const Eigen::Matrix<T, NStates, NStates> realX = X.real();

        // force the result X to be symetric to combat small numerical errors
        const Eigen::Matrix<T, NStates, NStates> result = static_cast<T>(0.5) * (realX + realX.transpose());
		return result;
	}

}