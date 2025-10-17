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
	 * \see controlpp::mexp_taylor_scale
	 * \see controlpp::mexp
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> mexp_taylor(
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
	 * \see controlpp::mexp_taylor
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> mexp_taylor_scaled(
			const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& M, 
			int taylor_order = 8, 
			int scaling = 10
	){
		using Matrix = Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>;
		T s = static_cast<T>(1 << scaling);
		Matrix scaled_M = M/s;
		Matrix t = mexp_taylor(scaled_M, taylor_order);
		for(int i = 0; i < scaling; ++i){
			const Matrix temp = t * t;
			t = temp;
		};
		return t;
	}

	/**
	 * \brief Calculates the numerator parameters of the pade approximation
     * 
     * Calculates the parameter:
     * 
     * \f[
     * P_k = \frac{(m + n - k)! m!}{(m + n)! k! (m - k)!}
     * \f]
     * 
     * For a Pade approximation like:
     * 
     * \f[
     * A = \frac{\sum_{k = 0}^{m} P_k s^k}{\sum_{k = 0}^{n} Q_k s^k} 
     * \f]
     * 
     * \param m Is the order of the numerator
     * \param n Is the order of the denominator
     * \param k Is the order of the parameter
     * 
     * \tparam The result type of the division
     * 
     * \returns The pade parameter of the numerator at the k-th order
     * 
     * \see pade_den_param
     * \see https://ris.utwente.nl/ws/portalfiles/portal/134422804/Some_remarks_on_Pade-approximations.pdf
	 */
    template<class T = double>
	constexpr T pade_num_param(unsigned long m, unsigned long n, unsigned long k){
		const unsigned long long num1 = product_over(k+1, m+1);
		const unsigned long long den1 = product_over(m + n - k + 1, m + n + 1);
		const unsigned long long den2 = product_over(1, m-k + 1);
		const T result = static_cast<T>(num1) / static_cast<T>(den1 * den2);
		return result;
	}

    /**
	 * \brief Calculates the numerator parameters of the pade approximation
     * 
     * Calculates the parameter:
     * 
     * \f[
     * Q_k = \frac{(n + m - k)! n!}{(n + m)! k! (n - k)!}
     * \f]
     * 
     * For a Pade approximation like:
     * 
     * \f[
     * A = \frac{\sum_{k = 0}^{m} P_k s^k}{\sum_{k = 0}^{n} Q_k s^k} 
     * \f]
     * 
     * \param m Is the order of the numerator
     * \param n Is the order of the denominator
     * \param k Is the order of the parameter
     * 
     * \tparam T the result type of the function (used for the final division)
     * 
     * \returns The pade parameter of the numerator at the k-th order
     * 
     * \see pade_num_param
     * \see https://ris.utwente.nl/ws/portalfiles/portal/134422804/Some_remarks_on_Pade-approximations.pdf
	 */
    template<class T = double>
	constexpr T pade_den_param(unsigned long m, unsigned long n, unsigned long k){
		return pade_num_param(n, m, k);
	}

	/**
	 * @brief Approximates \f$\exp{A}\f$ using a pade fraction
	 * @tparam T The value type of the matrix entries (usually `float` or `double`)
	 * @tparam Rows The number of rows
	 * @tparam Cols The number of columns
	 * @tparam Options Matrix options (See: [Store Orders](https://libeigen.gitlab.io/eigen/docs-nightly/group__TopicStorageOrders.html))
	 * @param A The input matrix
	 * @param Order The order of the numerator and denominator of the pade fraction
	 * @return An approximation of the exponential \f$\exp{A}\f$
	 */
	template<class T, int N, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, N, N> mexp_pade(
			const Eigen::Matrix<T, N, N, Options, MaxRows, MaxCols>& A,
			int Order = 5
	){
		const Eigen::Matrix<T, N, N> I = Eigen::Matrix<T, N, N>::Identity();
		const Eigen::Matrix<T, N, N> A2 = A * A;

		// even odd split of polynomials: allows to reuse results
		Eigen::Matrix<T, N, N> even;
		Eigen::Matrix<T, N, N> odd;
		even.setZero();
		odd.setZero();

		// horner chains to evaluate the polynomials
		for(int k = Order; k >= 0; --k){
			const bool is_first_iteration = (k == Order);
			const T pk = pade_num_param<T>(Order, Order, k);
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
		
		std::cout << std::endl;

		const Eigen::Matrix<T, N, N> num = even + odd;
		const Eigen::Matrix<T, N, N> den = even - odd;

		const Eigen::Matrix<T, N, N> result = den.partialPivLu().solve(num);
		
		return result;
	}

	template<class T, int N, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, N, N> mexp_pade_scaled(
			const Eigen::Matrix<T, N, N, Options, MaxRows, MaxCols>& M, 
			int order = 5,
			int scaling = 5
	){
		const T s = static_cast<T>(1ULL << scaling);
		Eigen::Matrix<T, N, N> scaled_M = M/s;
		Eigen::Matrix<T, N, N> t = mexp_pade(scaled_M, order);
		for(int i = 0; i < scaling; ++i){
			const Eigen::Matrix<T, N, N> temp = t * t;
			t = temp;
		};
		return t;
	}

	/**
	 * \brief Calculates the matrix exponent \f$ \exp{\mathbf{\M}} \f$
	 * 
	 * Uses a scaled pade approximation for the exponential.
	 * 
	 * \tparam T The value type of the matrix elements
	 * \tparam Rows The number of rows of the matrix
	 * \tparam Cols The number of columns of the matrix
	 * 
	 * \param M The matrix
	 * \param order The order of the pade polynomial used to approximate the exponential function 
	 * \param scaling The scaling factor used to improve the accuracy of the exponential function
	 * 
	 * \return The resulting exponentiated matrix
	 * 
	 * \see controlpp::mexp_pade
	 * \see controlpp::mexp_pade_scaled
	 */
	template<class T, int N, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, N, N> mexp(const Eigen::Matrix<T, N, N, Options, MaxRows, MaxCols>& M, unsigned int order = 3){
		// actual exponent calculation
		const T maxColNorm = M.colwise().norm().maxCoeff();
		const T maxRowNorm = M.rowwise().norm().maxCoeff();
		const T maxNorm = std::max(maxColNorm, maxRowNorm);
		const T one = static_cast<T>(1);
		int scaling = static_cast<int>(std::log2(std::max(maxNorm, one))) + 1;
		return mexp_pade_scaled(M, order, scaling);
	}

	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols> identity_like([[maybe_unused]]const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& unused){
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
	 * @brief Solves the hamilton matrix
	 * 
	 * For example when solveing the CARE (continuous time ricatti equation).
	 * 
	 * The solution of the hamiltoin matrix is understood as the result of the operations:
	 * 
	 * 1. Compute its stable eigenvectors
	 * 2. Re-Partitions the eigenvectors
	 * 3. Recovers X (ricatti result) from the partitions
	 * 
	 * @tparam T The value type of the parameters/matrix elements (ususally `float` or `double`)
	 * @tparam N The size of the hamilton matrix
	 * @param H The hamilton matrix
	 * @returns The solution of the hamilton matrix
	 * 
	 * @see lyapunov_solver
	 * @see care_solver
	 */
	template<class T, int N>
	requires(N % 2 == 0)
	Eigen::Matrix<T, N/2, N/2> hamilton_solver(
		const Eigen::Matrix<T, N, N> H
	){
		// 2. Compute stable eigenvectors
		Eigen::ComplexEigenSolver<Eigen::Matrix<T, N, N>> ces;
		ces.compute(H);
		const auto& H_eigvals = ces.eigenvalues();
		const auto& H_eigvecs = ces.eigenvectors();

		// in a 2n hamilton matrix are exactly n stable ones
		Eigen::Matrix<std::complex<T>, N, N/2> StableEigenVecs;
		int si = 0; // stable eigenvalue iterator
		int ei = 0; // eigen value iterator
		for(; (si < N/2) && (ei < N); ++ei){
			if(H_eigvals(ei).real() <= 0){// only stable ones
				StableEigenVecs.col(si) = H_eigvecs.col(ei);
				++si;
			}
		}

		// 3. Repartition
		const auto TopPartition = StableEigenVecs.template block<N/2, N/2>(0, 0);
		const auto BottomPartition = StableEigenVecs.template block<N/2, N/2>(N/2, 0);

		// 4. Recover Solution (BottomPartition * TopPartition^-1)
		const Eigen::Matrix<std::complex<T>, N/2, N/2> X = TopPartition.transpose().partialPivLu().solve(BottomPartition.transpose()).transpose();

		// The result is real, make sure it is real because it should be
		const Eigen::Matrix<T, N/2, N/2> realX = X.real();

        // force the result X to be symetric to combat small numerical errors
        const Eigen::Matrix<T, N/2, N/2> result = static_cast<T>(0.5) * (realX + realX.transpose());
		return result;
	}

	/**
	 * @brief Solves the hamilton matrix
	 * 
	 * For example when solveing the CARE (continuous time ricatti equation).
	 * 
	 * The solution of the hamiltoin matrix is understood as the result of the operations:
	 * 
	 * 1. Compute its stable eigenvectors
	 * 2. Re-Partitions the eigenvectors
	 * 3. Recovers X (ricatti result) from the partitions
	 * 
	 * @tparam T The value type of the parameters/matrix elements (ususally `float` or `double`)
	 * @tparam N The size of the hamilton matrix
	 * @param H The hamilton matrix
	 * @returns The solution of the hamilton matrix
	 * 
	 * @see lyapunov_solver
	 * @see care_solver
	 */
	template<class T, int N>
	Eigen::Matrix<T, N/2, N/2> symplectic_solver(
		const Eigen::Matrix<T, N, N> S
	){
		// 2. Compute stable eigenvectors
		Eigen::ComplexEigenSolver<Eigen::Matrix<T, N, N>> ces;
		ces.compute(S);
		const auto& eigvals = ces.eigenvalues();
		const auto& eigvecs = ces.eigenvectors();

		// in a 2n hamilton matrix are exactly n stable ones
		Eigen::Matrix<std::complex<T>, N, N/2> StableEigenVecs;
		int si = 0; // stable eigenvalue iterator
		int ei = 0; // eigen value iterator
		for(; (si < N/2) && (ei < N); ++ei){
			if(std::abs(eigvals(ei)) < 1){// only stable ones
				StableEigenVecs.col(si) = eigvecs.col(ei);
				++si;
			}
		}

		// 3. Repartition
		const auto TopPartition = StableEigenVecs.template block<N/2, N/2>(0, 0);
		const auto BottomPartition = StableEigenVecs.template block<N/2, N/2>(N/2, 0);

		// 4. Recover Solution (BottomPartition * TopPartition^-1)
		const Eigen::Matrix<std::complex<T>, N/2, N/2> X = TopPartition.transpose().partialPivLu().solve(BottomPartition.transpose()).transpose();

		// The result is real, make sure it is real because it should be
		const Eigen::Matrix<T, N/2, N/2> realX = X.real();

        // force the result X to be symetric to combat small numerical errors
        const Eigen::Matrix<T, N/2, N/2> result = static_cast<T>(0.5) * (realX + realX.transpose());
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
	 * 
	 * \see hamilton_solver
	 */
	template<class T, int NStates,
			int AOpt, int AMaxR, int AMaxC,
			int QOpt, int QMaxR, int QMaxC
  	>
	Eigen::Matrix<T, NStates, NStates> lyapunov_solver(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NStates, QOpt, QMaxR, QMaxC>& Q
	){
		// 1. Build the hamilton matrix
		Eigen::Matrix<T, 2*NStates, 2*NStates> H;
		H.topLeftCorner(NStates, NStates) = A;
		H.topRightCorner(NStates, NStates).setZero();
		H.bottomLeftCorner(NStates, NStates) = - Q;
		H.bottomRightCorner(NStates, NStates) = -A.transpose();
		return hamilton_solver(H);
	}

	/**
     * \brief Calculates the energy (Gramian) scaling via Lyapunov equations
	 * 
	 * Calculates the scaling (similarity) matrix T so that the similar system:
	 * 
	 * \f[
	 * A' = T A T^{-1}, \quad B' = T B, \quad C' = C T^{-1}, \quad D' = D
	 * \f]
	 * 
	 * has equal state energy.
	 * 
	 * This is especially useful, if a system has very small and very large entries 
	 * of very small and large frequency components/dynamics.
	 * 
     * \tparam ValueType The value type of the matrices. (E.g.: `double`, `float`)
     * \tparam NStates The number of states
     * \tparam NInputs The number of inputs
     * \tparam NOutputs The number of outputs
	 * \param A The system state matrix (NStates x NStates)
	 * \param B The system input matrix (NStates x NInputs)
	 * \param C The system output matrix (NOutputs x NStates)
	 * \param D The system passthrough matrix (NOutputs x NInputs)
     * \return The tuple of matrices T and its inverse \f$T^{-1}\f$ [T, T_inverse]
     */
    template<class ValueType, int NStates, int NInputs, int NOutputs>
    std::tuple<Eigen::Matrix<ValueType, NStates, NStates>, Eigen::Matrix<ValueType, NStates, NStates>> energy_scaling_matrix(
		const Eigen::Matrix<ValueType, NStates, NStates> A,
		const Eigen::Matrix<ValueType, NStates, NInputs> B,
		const Eigen::Matrix<ValueType, NOutputs, NStates> C
	){
		// Calculate the two gramians Wc and Wo
		const Eigen::Matrix<ValueType, NStates, NStates> Bsqr = B * B.transpose();
		const Eigen::Matrix<ValueType, NStates, NStates> Wc = lyapunov_solver(A.transpose(), Bsqr);

		const Eigen::Matrix<ValueType, NStates, NStates> Csqr = C.transpose() * C;
		const Eigen::Matrix<ValueType, NStates, NStates> Wo = lyapunov_solver(A, Csqr);

		// force symetry and add some to the diagonal for numerical stability
		const ValueType eps = std::numeric_limits<ValueType>::epsilon();

		Eigen::Matrix<ValueType, NStates, NStates> Wc_sym = (Wc + Wc.transpose()) * static_cast<ValueType>(0.5);
		{
			ValueType jitter = std::max(ValueType(1), Wc_sym.diagonal().cwiseAbs().maxCoeff());
			Wc_sym.diagonal().array() += jitter * eps;
		}
		Eigen::Matrix<ValueType, NStates, NStates> Wo_sym = (Wo + Wo.transpose()) * static_cast<ValueType>(0.5);
		{
			ValueType jitter = std::max(ValueType(1), Wo_sym.diagonal().cwiseAbs().maxCoeff());
			Wo_sym.diagonal().array() += jitter * eps;
		}

		// calculate the cholesky factors using LLT
		const Eigen::LLT<Eigen::Matrix<ValueType, NStates, NStates>> Wc_colesky(Wc_sym);
		const Eigen::LLT<Eigen::Matrix<ValueType, NStates, NStates>> Wo_colesky(Wo_sym);

		// partition into triangle forms
		const Eigen::Matrix<ValueType, NStates, NStates> Rc = Wc_colesky.matrixU();
		const Eigen::Matrix<ValueType, NStates, NStates> Ro = Wo_colesky.matrixL();

		// calculate the SVD
		const Eigen::Matrix<ValueType, NStates, NStates> M = Ro * Rc;
		const Eigen::JacobiSVD<Eigen::Matrix<ValueType, NStates, NStates>> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
		const Eigen::Matrix<ValueType, NStates, NStates> U  = svd.matrixU();
		const Eigen::Matrix<ValueType, NStates, 1> s  = svd.singularValues();
		const Eigen::Matrix<ValueType, NStates, NStates> V = svd.matrixV();

		
		const Eigen::Matrix<ValueType, NStates, 1> s_sqrt = (s.array().max(eps)).sqrt().matrix();
		const Eigen::Matrix<ValueType, NStates, NStates> S_sqrt = s_sqrt.asDiagonal();		

		// Finally!!! we can calculate the resulting transformation matrices
		const Eigen::Matrix<ValueType, NStates, NStates> T = Rc.transpose().template triangularView<Eigen::Lower>().solve((S_sqrt * V.transpose()).transpose()).transpose();
		const Eigen::Matrix<ValueType, NStates, NStates> T_inverse = Ro.template triangularView<Eigen::Lower>().solve(U * S_sqrt);

		return {T, T_inverse};
    }

	

	/**
	 * \brief Solves the continuous time riccati equation (CARE)
	 * 
	 * This function computes the stabilizing symmetric solution of the CARE:
     * 
	 * \f[
	 * A^\top X + X A - (X B) R^{-1} (B^\top X) + Q = 0
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
	 *
	 * \tparam T Scalar type (e.g., `double`, `float`).
	 * \tparam NStates Number of states.
	 * \tparam NInputs Number of control inputs.
	 * 
	 * @see Zhou, Doyle, and Glover (1996), *Robust and Optimal Control*.
	 * @see hamilton_solver
	 */
	template<class T, int NStates, int NInputs,
			int AOpt, int AMaxR, int AMaxC,
			int BOpt, int BMaxR, int BMaxC,
			int ROpt, int RMaxR, int RMaxC,
			int QOpt, int QMaxR, int QMaxC
  	>
	Eigen::Matrix<T, NStates, NStates> care_solver(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NInputs, BOpt, BMaxR, BMaxC>& B,
		const Eigen::Matrix<T, NStates, NStates, QOpt, QMaxR, QMaxC>& Q,
		const Eigen::Matrix<T, NInputs, NInputs, ROpt, RMaxR, RMaxC>& R
	){
		Eigen::Matrix<T, 2*NStates, 2*NStates> H;
		H.topLeftCorner(NStates, NStates) = A;
		H.topRightCorner(NStates, NStates) =  -B * R.ldlt().solve(B.transpose());
		H.bottomLeftCorner(NStates, NStates) = -Q;
		H.bottomRightCorner(NStates, NStates) = -A;
		return hamilton_solver(H);
	}

	/**
	 * \brief Solves the discrete time riccati equation (DARE)
	 * 
	 * This function computes the stabilizing symmetric solution of the DARE:
     * 
	 * \f[
	 * A^\top X A - A^\top X B (R + B^\top X B)^{-1} B^\top X A + Q = 0
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
	 *
	 * \tparam T Scalar type (e.g., `double`, `float`).
	 * \tparam NStates Number of states.
	 * \tparam NInputs Number of control inputs.
	 * 
	 * @see hamilton_solver
	 */
	template<class T, int NStates, int NInputs,
			int AOpt, int AMaxR, int AMaxC,
			int BOpt, int BMaxR, int BMaxC,
			int ROpt, int RMaxR, int RMaxC,
			int QOpt, int QMaxR, int QMaxC
  	>
	Eigen::Matrix<T, NStates, NStates> dare_solver(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NInputs, BOpt, BMaxR, BMaxC>& B,
		const Eigen::Matrix<T, NStates, NStates, QOpt, QMaxR, QMaxC>& Q,
		const Eigen::Matrix<T, NInputs, NInputs, ROpt, RMaxR, RMaxC>& R
	){
		const Eigen::Matrix<T, NStates, NStates> brb = -B * R.llt().solve(B.transpose());
		const Eigen::Matrix<T, NStates, NStates> aq = A.transpose().colPivHouseholderQr().solve(Q);

		Eigen::Matrix<T, 2*NStates, 2*NStates> S;
		S.topLeftCorner(NStates, NStates) = A - brb * aq;
		S.topRightCorner(NStates, NStates) =  - A.colPivHouseholderQr().solve(brb.transpose()).transpose();
		S.bottomLeftCorner(NStates, NStates) = aq;
		S.bottomRightCorner(NStates, NStates) = -A;
		return symplectic_solver(S);
	}

	/**
	 * \brief Solves the continuous time riccati equation (CARE) with cross coupling
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
	 * 
	 * @see Zhou, Doyle, and Glover (1996), *Robust and Optimal Control*.
	 * @see hamilton_solver
	 */
	template<class T, int NStates, int NInputs,
			int AOpt, int AMaxR, int AMaxC,
			int BOpt, int BMaxR, int BMaxC,
			int ROpt, int RMaxR, int RMaxC,
			int QOpt, int QMaxR, int QMaxC,
			int NOpt, int NMaxR, int NMaxC
  	>
	Eigen::Matrix<T, NStates, NStates> care_solver(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NInputs, BOpt, BMaxR, BMaxC>& B,
		const Eigen::Matrix<T, NStates, NStates, QOpt, QMaxR, QMaxC>& Q,
		const Eigen::Matrix<T, NInputs, NInputs, ROpt, RMaxR, RMaxC>& R,
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
		return solve_hamilton(H);
	}

}