#pragma once

/**
 * @file dare.hpp
 * @brief This file contains the implementation of the Discrete-time Algebraic Riccati Equation (DARE) solver.
 * @author Tobias Wallner
 * @copyright 2026 Tobias Wallner
 */

#include <complex>
#include <iostream>

#include <Eigen/Dense>

namespace controlpp{

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
		StableEigenVecs.setZero();
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
	 * @returns The solution of the DARE as an Eigen::Matrix with the dimensions `NStates x NStates`.
	 * 
	 * @see symplectic_solver
	 */
	template<class T, int NStates, int NInputs,
			int AOpt, int AMaxR, int AMaxC,
			int BOpt, int BMaxR, int BMaxC,
			int ROpt, int RMaxR, int RMaxC,
			int QOpt, int QMaxR, int QMaxC
  	>
	Eigen::Matrix<T, NStates, NStates> dare_eigenvector_solver(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NInputs, BOpt, BMaxR, BMaxC>& B,
		const Eigen::Matrix<T, NStates, NStates, QOpt, QMaxR, QMaxC>& Q,
		const Eigen::Matrix<T, NInputs, NInputs, ROpt, RMaxR, RMaxC>& R
	){
		const Eigen::Matrix<T, NStates, NStates> brb = B * R.llt().solve(B.transpose());
		Eigen::ColPivHouseholderQR<Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>> At_qr(A.transpose());
		const Eigen::Matrix<T, NStates, NStates> aq = At_qr.solve(Q);
		Eigen::Matrix<T, 2*NStates, 2*NStates> S;
		S.topLeftCorner(NStates, NStates) = A + brb * aq;
		S.topRightCorner(NStates, NStates) = -A.colPivHouseholderQr().solve(brb.transpose()).transpose();
		S.bottomLeftCorner(NStates, NStates) = - aq;
		S.bottomRightCorner(NStates, NStates) = At_qr.solve(identity_like(A));
		return symplectic_solver(S);
	}

	/**
	 * @brief Builds the symplectic pencil matrix H for the discrete algebraic riccati equation (DARE)
	 * 
	 * To solve the DARE:
     * 
	 * \f[
	 * A^\top X A - A^\top X B (R + B^\top X B)^{-1} B^\top X A + Q = 0
	 * \f]
	 * 
	 * we can construct the symplectic pencil matrix H, which is used in the solution process.
	 * 
	 * the pencil is given by:
	 * 
	 * \f[
	 * H - \lambda J
	 * \f]
	 * 
	 * This function constructs the matrix H matrix:
	 * 
	 * \f[
	 * \begin{bmatrix}
	 *    A & 0 & B \\
	 *   -Q & I & 0 \\
	 *    0 & 0 & R
	 * \end{bmatrix}
	 * \f]
	 * 
	 * @see build_dare_symplectic_pencil_J
	 * 
	 * @tparam T 
	 * @tparam NStates 
	 * @tparam NInputs 
	 * @param A 
	 * @param B 
	 * @param Q 
	 * @param R 
	 * @return The symplectic pencil matrix H for the DARE as an Eigen::Matrix with the dimensions `(2*NStates + NInputs) x (2 NStates + NInputs)`.
	 */
	template<class T, int NStates, int NInputs,
			int AOpt, int AMaxR, int AMaxC,
			int BOpt, int BMaxR, int BMaxC,
			int ROpt, int RMaxR, int RMaxC,
			int QOpt, int QMaxR, int QMaxC
  	>
	Eigen::Matrix<T, 2*NStates + NInputs, 2*NStates + NInputs> build_dare_symplectic_pencil_H(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NInputs, BOpt, BMaxR, BMaxC>& B,
		const Eigen::Matrix<T, NStates, NStates, QOpt, QMaxR, QMaxC>& Q,
		const Eigen::Matrix<T, NInputs, NInputs, ROpt, RMaxR, RMaxC>& R
	){
		constexpr int H_rows = 2*NStates + NInputs;
		constexpr int H_cols = 2*NStates + NInputs;

		// build the matrix:
		// [  A   0   B ]
		// [ -Q   I   0 ]
		// [  0   0   R ]
		Eigen::Matrix<T, H_rows, H_cols> H;
		H.setZero();
		H.template block<NStates, NStates>(0, 0) = A;
		H.template block<NStates, NInputs>(0, 2 * NStates) = B;
		H.template block<NStates, NStates>(NStates, 0) = -Q;
		H.template block<NStates, NStates>(NStates, NStates).setIdentity();
		H.template block<NInputs, NInputs>(2 * NStates, 2 * NStates) = R;

		return H;
	}

    /**
	 * @brief Builds the symplectic pencil matrix J for the discrete algebraic riccati equation (DARE)
	 * 
	 * To solve the DARE:
     * 
	 * \f[
	 * A^\top X A - A^\top X B (R + B^\top X B)^{-1} B^\top X A + Q = 0
	 * \f]
	 * 
	 * we can construct the symplectic pencil matrix J, which is used in the solution process.
	 * 
	 * the pencil is given by:
	 * 
	 * \f[
	 * H - \lambda J
	 * \f]
	 * 
	 * This function constructs the matrix J matrix:
	 * 
	 * \f[
	 * \begin{bmatrix}
	 *   I & 0 & B \\
	 *   0 & A^T & 0 \\
	 *   0 & -B^T & 0
	 * \end{bmatrix}
	 * \f]
	 * 
	 * @see build_dare_symplectic_pencil_H
	 * 
	 * @tparam T 
	 * @tparam NStates 
	 * @tparam NInputs 
	 * @param A 
	 * @param B 
	 * @param Q 
	 * @param R 
	 * @return The symplectic pencil matrix J for the DARE as an Eigen::Matrix with the dimensions `(2 NStates + NInputs) x (2 NStates + NInputs)`.
	 */
	template<class T, int NStates, int NInputs,
			int AOpt, int AMaxR, int AMaxC,
			int BOpt, int BMaxR, int BMaxC
  	>
	Eigen::Matrix<T, 2*NStates + NInputs, 2*NStates + NInputs> build_dare_symplectic_pencil_J(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NInputs, BOpt, BMaxR, BMaxC>& B
	){	
		constexpr int J_rows = NStates + NStates + NInputs;
		constexpr int J_cols = NStates + NStates + NInputs;
		
		// build the matrix:
		// [ I   0    B]
		// [ 0   A^T  0]
		// [ 0  -B^T  0]
		Eigen::Matrix<T, J_rows, J_cols> J;
		J.setZero();
		J.template block<NStates, NStates>(0, 0).setIdentity();
		J.template block<NStates, NInputs>(0, 2 * NStates) = B;
		J.template block<NStates, NStates>(NStates, NStates) = A.transpose();
		J.template block<NInputs, NStates>(2 * NStates, NStates) = -B.transpose();

		return J;
	}

    /**
	 * @brief Determines the block size of a quasi-upper-triangular matrix S at a given block start index.
	 * 
	 * Asserts: That the block_start is within the valid range of the matrix S.
	 * 
	 * @tparam T The value type of the matrix elements (e.g., `float`, `double`).
	 * @tparam Rows The number of rows in the matrix S.
	 * @tparam Cols The number of columns in the matrix S.
	 * @tparam Options The options for the Eigen matrix (e.g., storage order).
	 * @tparam MaxRows The maximum number of rows in the matrix S.
	 * @tparam MaxCols The maximum number of columns in the matrix S.
	 * @param S The quasi-upper-triangular matrix to analyze.
	 * @param block_start The starting index of the block to analyze.
	 * @return The size of the block (1 or 2) at the specified block start index.
	 * @throws std::runtime_error If both S and T are zero, indicating an indeterminate pencil block.
	 */
	template<class T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	int pencil_block_size(const Eigen::Matrix<T, Rows, Cols, Options, MaxRows, MaxCols>& S, int block_start){
        if(block_start + 1 < S.rows()){
            // determine the blocksize (1x1 or 2x2)
            const T sub_diag = S(block_start + 1, block_start);
            
            // zero threshold for numerical stability
            const T local_diagonal = std::abs(S(block_start, block_start)) + std::abs(S(block_start + 1, block_start + 1));
            const T epsilon  = std::numeric_limits<T>::epsilon();
            const T zero_threshold = epsilon  * std::max(local_diagonal, T(1)) * T(1024);
			if(std::abs(sub_diag) > zero_threshold){
				return 2;
			}else{
				return 1;
			}
		}else{
			return 1;
		}
	}

	/**
	 * @brief Computes normalized characteristic-polynomial invariants for a 2x2 generalized eigenvalue problem.
	 *
	 * Given two 2x2 matrices S and T, the generalized eigenvalues are the roots
	 * of the characteristic equation
	 *
	 * @f[
	 *     \det(S - \lambda T) = c_0 + c_1\lambda + c_2\lambda^2.
	 * @f]
	 *
	 * The polynomial coefficients are
	 *
	 * @f[
	 *     c_0 = \det(S),
	 * @f]
	 *
	 * @f[
	 *     c_1 = -\left(
	 *         S_{00}T_{11} + S_{11}T_{00}
	 *         - S_{01}T_{10} - S_{10}T_{01}
	 *     \right),
	 * @f]
	 *
	 * @f[
	 *     c_2 = \det(T).
	 * @f]
	 *
	 * Assuming @f$c_2 \ne 0@f$, this function returns the normalized
	 * coefficients
	 *
	 * @f[
	 *     \begin{bmatrix}
	 *         c_0/c_2 \\
	 *         c_1/c_2
	 *     \end{bmatrix}.
	 * @f]
	 *
	 * These normalized coefficients are characteristic-polynomial invariants:
	 * they uniquely determine the monic polynomial associated with the
	 * generalized eigenvalues, including when the eigenvalues form a
	 * complex-conjugate pair.
	 *
	 * They are useful for verifying that two 2x2 matrix pairs have the same
	 * generalized eigenvalues without comparing diagonal entry ratios, which
	 * are not generally the individual eigenvalues of a real 2x2 Schur block.
	 *
	 * Under orthogonal equivalence transformations
	 *
	 * @f[
	 *     S' = U^\top S V,\qquad T' = U^\top T V,
	 * @f]
	 *
	 * where U and V are orthogonal, the determinant coefficients acquire
	 * the same common factor, @f$\det(U^\top)\det(V)@f$. Consequently, their
	 * normalized ratios remain unchanged whenever @f$\det(T)\ne 0@f$.
	 *
	 * @param S The 2x2 matrix in the generalized eigenvalue problem.
	 * @param T The 2x2 matrix paired with S.
	 *
	 * @return A vector containing @f$c_0/c_2@f$ and @f$c_1/c_2@f$.
	 *
	 * @note Requires @f$\det(T)\ne 0@f$. If T is singular or nearly singular,
	 *       the normalized coefficients may be undefined or numerically
	 *       unreliable.
	 */
	template<class Type>
	Eigen::Vector<Type, 2> generalized_charpoly_invariants(
    	const Eigen::Matrix<Type, 2, 2>& S,
    	const Eigen::Matrix<Type, 2, 2>& T
	){
		// det(S - lambda*T) = c0 + c1*lambda + c2*lambda^2
		const double c0 = S.determinant();

		const double c1 =
			-(S(0, 0) * T(1, 1)
			+ S(1, 1) * T(0, 0)
			- S(0, 1) * T(1, 0)
			- S(1, 0) * T(0, 1));

		const double c2 = T.determinant();

		return Eigen::Vector<Type, 2>(c0 / c2, c1 / c2);
	}

    /**
	 * @brief Determines if a given pencil block (Sb, Tb) is stable, meaning all eigenvalues of the pencil are inside the unit circle.
	 * 
	 * Asserts that:
	 *  - The matrices Sb and Tb have the same dimensions.
	 *  - The matrices Sb and Tb are either 1x1 or 2x2.
	 * 
	 * @tparam T The value type of the matrix elements (e.g., `float`, `double`).
	 * @tparam Options The options for the Eigen matrix (e.g., storage order).
	 * @tparam MaxRows The maximum number of rows in the matrices Sb and Tb.
	 * @tparam MaxCols The maximum number of columns in the matrices Sb and Tb.
	 * @tparam N The size of the matrices Sb and Tb (must be either 1 or 2).
	 * @param Sb The matrix representing the S part of the pencil block. (quasi-upper-triangular from the QZ decomposition)
	 * @param Tb The matrix representing the T part of the pencil block. (upper-triangular from the QZ decomposition)
	 * @return The boolean value indicating whether the pencil block is stable (true) or not (false).
	 * @throws std::runtime_error If both Sb and Tb have zero eigenvalues, indicating an indeterminate pencil block.
	 */
	template<class DerivedS, class DerivedT>
	bool is_stable_pencil_block(
		const Eigen::MatrixBase<DerivedS>& Sb,
		const Eigen::MatrixBase<DerivedT>& Tb
	){
		assert(Sb.rows() == Sb.cols());
		assert(Tb.rows() == Tb.cols());
		assert(Sb.rows() == Tb.rows());
		assert(Sb.rows() >= 1 && Sb.rows() <= 2);

		if(Sb.rows() == 1){
			// blocksize is 1

			// check if the eigenvalue is stable (inside the unit circle)
			const auto s = Sb(0, 0);
			const auto t = Tb(0, 0);

			// avoid division by zero
			const bool stable = (std::abs(s) < std::abs(t));
			return stable;
		}else{
			// blocksize is 2

			// solve det(S - lambda T) = 0 for the 2x2 block to find the eigenvalues
			// build the quadratic polynomial coefficients
			//
			// 
			// det | (S00 - lambda T00)   (S01 - lambda T01) |
			//     | (S10)                (S11 - lambda T11) |
			const auto a = Tb(0, 0) * Tb(1, 1);
			const auto b = Sb(1, 0) * Tb(0, 1) - Sb(0, 0) * Tb(1, 1) - Sb(1, 1) * Tb(0, 0);
			const auto c = Sb(0, 0) * Sb(1, 1) - Sb(0, 1) * Sb(1, 0);

            if(a == 0){
                // The characteristic polynomial is not quadratic.
                // At least one generalized eigenvalue is infinite
                // (assuming the pencil is regular), hence unstable.
                return false;
            }

            //const auto discriminant = b * b - 4 * a * c;
            const bool discriminant_positive = (b * b) >= (4 * a * c);
            if(discriminant_positive){
                // real eigenvalues

                // use the Jury criterion for real roots
                const bool cond_1 = (a * (a + b + c)) > 0;
                const bool cond_2 = (a * (a - b + c)) > 0;
                const bool cond_3 = (a * (a - c)) > 0;

                return cond_1 && cond_2 && cond_3;

            }else{
                // complex eigenvalues

                // check that the eigenvalues 
                // lambda1,2 = (-b +- sqrt(b^2 - 4ac)) / 2a 
                // are stable (inside the unit circle)
                // 
                // for stability, we need |lambda| < 1
                // which is equivalent to checking |lambda|^2 < 1 
                // which is equivalent to checking lambda * conj(lambda) < 1
                // and we know that lambda1 = conj(lambda2) for complex eigenvalues
                // so this is equivalent to checking |lambda_1| * |lambda_2| < 1 
                // therefore:
                // lambda1 * lambda2 = (-b + sqrt(b^2 - 4ac)/(2a) * (-b - sqrt(b^2 - 4ac)/(2a) = c/a
                // so we need to check if c/a < 1
                return std::abs(c) < std::abs(a);
            }
		}
	}

	template<class T, std::size_t N>
	struct QZSwap{
		Eigen::Matrix<T, N, N> V;
		Eigen::Matrix<T, N, N> U;
	};

	/**
	 * @brief Caclulates the V matrix for a qt 1x1 block swap
	 * 
	 * for upper triangular block (1) matrices
	 * 
	 * \f[
	 * S_b' =
	 * 	\begin{bmatrix}
	 * 		s_{22}' & *\\
	 * 		0 & s_{11}'
	 * 	\end{bmatrix},
	 * \qquad
	 * T_b' =
	 * 	\begin{bmatrix}
	 * 		t_{22}' & *\\
	 * 		0 & t_{11}'
	 * 	\end{bmatrix}.
	 * f\]
	 * 
	 * finds the V matrix such that:
	 * 
	 * f\[
	 * 	S_b' = U^\top S_b V,
	 * 	\qquad
	 * 	T_b' = U^\top T_b V,
	 * f\]
	 * 
	 * swaps the eigenvalus of the matrices while remaining the
	 * generalised Shur form for both of them
	 * 
	 * @param s Upper triangular generalised Schur S matrix from the QZ decomposition
	 * @param t Upper triangular generalised Schur T matrix from the QZ decomposition
	 * 
	 * @tparam T The value type. Usually float or double.
	 */
	template<class T>
	QZSwap<T, 2> qz_1x1_swap_matrix(const Eigen::Matrix<T, 2, 2>& s, const Eigen::Matrix<T, 2, 2>& t){
		// find the right generalised eigenvector
		// for the second eigenvalue: (sb - lambda tb) v = 0

		// the second row provides: (s11 - lambda1 t11) v1 = 0
		// remember: lambda1 = s11 / t11
		// --> (s11 - s11 / t11 * t11) = 0
		// --> always satisfied for every v1
		
		// The first row: (s00 - lambda1 t00) v1 + (s01 - lambda1 t01) v1 = 0
		const T v0 = - (s(0,1) * t(1,1) - s(1,1) * t(0,1));
		const T v1 = (s(0,0) * t(1,1) - s(1,1) * t(0,0));

		const Eigen::Vector<T, 2> v(v0, v1);

		// normalise it to unit length
		const Eigen::Vector<T, 2> v_u = v.normalized();

		// construct the full orthogonal V matrix
		const Eigen::Matrix<T, 2, 2> V({
			{v_u(0), -v_u(1)},
			{v_u(1), v_u(0)}
		});

		// calculate the right transformed matrices
		const Eigen::Vector<T, 2> a = s * V.col(0);
		// const Eigen::Matrix<T, 2, 2> B = t * B; // ignored because its parallel

		// now satisfy U^T A.col(0) = [*, 0] and similarly for B
		// u10 a0 + u11 a1 = 0
		// u10 b0 + u11 b1 = 0
		
		const T u11 = a(0);
		const T u10 = -a(1);
		const Eigen::Vector<T, 2> u_c1(u10, u11); // bottom row vector of U
		
		// norm it
		const Eigen::Vector<T, 2> u_n = u_c1.normalized();

		// build the U matrix
		const Eigen::Matrix<T, 2, 2> U({
			{u_n(1), u_n(0)},
			{-u_n(0), u_n(1)}
		});

		const QZSwap<T, 2> result{
			.V = V,
			.U = U
		};

		return result;
	}

	/**
	 * @brief Creates the matrices to swap eigenvalues of an adjacent 1x1 with a 2x2 block in Shur form
	 * 
	 * Where T and S contain the generalised 1x1 and 2x2 blocks in shur form
	 * f\[
	 *	S=
	 *	\begin{bmatrix}
	 *		a & b_1 & b_2\\
	 *		0 & d & e\\
	 *		0 & f & g
	 *	\end{bmatrix},
	 *	\qquad
	 *	T=
	 *	\begin{bmatrix}
	 *		\alpha & \beta_1 & \beta_2\\
	 *		0 & \delta & \epsilon\\
	 *		0 & 0 & \gamma
	 *	\end{bmatrix}.
	 *	f\]
	 * 
	 * @param S Generalised Schur matrix S
	 * @param T Generalised Schur matrix T
	 * 
	 * @throws `std::runtime_error` on: Singular Sylvester systems, non-finite sylvester solutions, degenerate systems, qz decomposition failure
	 */
	template<class T>
	QZSwap<T, 3> qz_1x2_swap_matrix(const Eigen::Matrix<T, 3, 3>& Sb, const Eigen::Matrix<T, 3, 3>& Tb){

		const T eps = std::numeric_limits<T>::epsilon();
		
		// find right invariant subspace
		// partition input matrices:
		//	S_b = 	a C
		//			0 A
		//
		//	T_b =   alpha F 
		//          0     B

		// Local diagonal blocks and upper-right coupling rows.
		const T a = Sb(0, 0);
		const T alpha = Tb(0, 0);
		const Eigen::Matrix<T, 2, 2> A = Sb.template block<2, 2>(1, 1);
		const Eigen::Matrix<T, 2, 2> B = Tb.template block<2, 2>(1, 1);
		const Eigen::Matrix<T, 1, 2> C = Sb.template block<1, 2>(0, 1);
		const Eigen::Matrix<T, 1, 2> F = Tb.template block<1, 2>(0, 1);

		// Solve:
		//   a*R - L*A = -C
		// alpha*R - L*B = -F
		//
		// Unknown ordering: [R0, R1, L0, L1]^T.

		Eigen::Matrix<T, 4, 4> M({
			{a,     T(0), -A(0, 0), -A(1, 0)},
			{T(0),  a,    -A(0, 1), -A(1, 1)},
			{alpha, T(0), -B(0, 0), -B(1, 0)},
			{T(0),  alpha,-B(0, 1), -B(1, 1)}
		});

		Eigen::Vector<T, 4> r(
			-C(0, 0),
			-C(0, 1),
			-F(0, 0),
			-F(0, 1)
		);

		// Scale the system to reduce overflow/underflow risk.

		const T scale = M.cwiseAbs().maxCoeff();

		if (scale == T(0)) {
			throw std::runtime_error("Controlpp: QZ 1x2 swap: degenerate Sylvester system");
		}

		M /= scale;
    	r /= scale;

		Eigen::FullPivLU<Eigen::Matrix<T, 4, 4>> lu(M);

		if (lu.rank() < 4) {
			throw std::runtime_error("Controlpp: QZ 1x2 swap: singular Sylvester system");
		}

		const Eigen::Vector<T, 4> z = lu.solve(r);

		if (!z.allFinite()) {
			throw std::runtime_error("QZ 1x2 swap: non-finite Sylvester solution");
		}

		const T residual = (M * z - r).norm();
    	const T residual_scale = M.norm() * z.norm() + r.norm();

		if (residual > eps * residual_scale) {
			throw std::runtime_error("Controlpp: QZ 1x2 swap: Sylvester solve residual too large");
		}

		const Eigen::Matrix<T, 1, 2> R = z.template head<2>().transpose();
    	const Eigen::Matrix<T, 1, 2> L = z.template tail<2>().transpose();

		// Build the right and left deflating-subspace bases.
		const Eigen::Matrix<T, 3, 2> X2({
			{R(0), R(1)},
			{T(1), T(0)},
			{T(0), T(1)}
		});
		const Eigen::Matrix<T, 3, 2> Y2({
			{L(0), L(1)},
			{T(1), T(0)},
			{T(0), T(1)}
		});

		Eigen::HouseholderQR<Eigen::Matrix<T, 3, 2>> qrX(X2);
    	Eigen::HouseholderQR<Eigen::Matrix<T, 3, 2>> qrY(Y2);

		const Eigen::Matrix<T, 3, 3> V = qrX.householderQ();
    	const Eigen::Matrix<T, 3, 3> U = qrY.householderQ();

		// Apply the orthogonal transformations.
		const Eigen::Matrix<T, 3, 3> S1 = U.transpose() * Sb * V;
		const Eigen::Matrix<T, 3, 3> T1 = U.transpose() * Tb * V;

		// Restore generalized Schur form in the leading 2x2 pair.
		const Eigen::Matrix<T, 2, 2> Slead = S1.template block<2, 2>(0, 0);
		const Eigen::Matrix<T, 2, 2> Tlead = T1.template block<2, 2>(0, 0);

		Eigen::RealQZ<Eigen::Matrix<T, 2, 2>> localQZ;
    	localQZ.compute(Slead, Tlead);

		if (localQZ.info() != Eigen::Success) {
        	throw std::runtime_error("Controlpp: QZ 1x2 swap: local 2x2 QZ decomposition failed");
		}

		Eigen::Matrix<T, 3, 3> Ulead = Eigen::Matrix<T, 3, 3>::Identity();
    	Eigen::Matrix<T, 3, 3> Vlead = Eigen::Matrix<T, 3, 3>::Identity();
		
		Ulead.template block<2, 2>(0, 0) = localQZ.matrixQ();
    	Vlead.template block<2, 2>(0, 0) = localQZ.matrixZ().transpose();

		const Eigen::Matrix<T, 3, 3> U_restored = U * Ulead;
    	const Eigen::Matrix<T, 3, 3> V_restored = V * Vlead;

		const QZSwap<T, 3> result{
			.V = V_restored,
			.U = U_restored
		};

		return result;
	}

	/**
	* @brief Creates the matrices to swap an adjacent 2x2 block with a 1x1
	* block in generalized Schur form.
	*
	* The input matrices are partitioned as:
	*
	* \[
	* S_b=\begin{bmatrix}A&C\\0&a\end{bmatrix},
	* \qquad
	* T_b=\begin{bmatrix}B&F\\0&\alpha\end{bmatrix},
	* \]
	* 
	* The resulting transformations move the 1x1 block to the leading
	* position and the 2x2 block to the trailing position.
	*
	* @param Sb Generalized Schur matrix S.
	* @param Tb Generalized Schur matrix T.
	*
	* @throws std::runtime_error on a singular Sylvester system, non-finite
	* solution, excessive solve residual, or local QZ decomposition failure.
	*/
	template<class T>
	QZSwap<T, 3> qz_2x1_swap_matrix(
		const Eigen::Matrix<T, 3, 3>& Sb, 
		const Eigen::Matrix<T, 3, 3>& Tb)
	{
		const T eps = std::numeric_limits<T>::epsilon();

		// Partition the input matrices:
		//
		// S = [ A  C ]       T = [ B  F ]
		//     [ 0  a ]           [ 0  alpha ]

		const Eigen::Matrix<T, 2, 2> A = Sb.template block<2, 2>(0, 0);
		const Eigen::Matrix<T, 2, 2> B = Tb.template block<2, 2>(0, 0);
		const Eigen::Matrix<T, 2, 1> C = Sb.template block<2, 1>(0, 2);
		const Eigen::Matrix<T, 2, 1> F = Tb.template block<2, 1>(0, 2);
		const T a = Sb(2, 2);
		const T alpha = Tb(2, 2);

		// Solve the coupled Sylvester equations:
		//
		// A*R - L*a     = -C
		// B*R - L*alpha = -F
		//
		// Unknown ordering: [R0, R1, L0, L1]^T.

		Eigen::Matrix<T, 4, 4> M({
			{A(0, 0), A(0, 1), -a,     T(0)},
			{A(1, 0), A(1, 1), T(0),  -a},
			{B(0, 0), B(0, 1), -alpha, T(0)},
			{B(1, 0), B(1, 1), T(0),  -alpha}
		});

		Eigen::Matrix<T, 4, 1> rhs(-C(0), -C(1), -F(0), -F(1));

		if (!M.allFinite() || !rhs.allFinite()) {
			throw std::runtime_error("Controlpp: QZ 2x1 swap: non-finite Sylvester system");
		}

		// Scale the system to reduce overflow/underflow risk.
		const T scale = std::max(M.cwiseAbs().maxCoeff(), rhs.cwiseAbs().maxCoeff());

		if (scale == T(0)) {
			throw std::runtime_error("Controlpp: QZ 2x1 swap: degenerate Sylvester system");
		}

		M /= scale;
		rhs /= scale;

		Eigen::FullPivLU<Eigen::Matrix<T, 4, 4>> lu(M);

		if (lu.rank() < 4) {
			throw std::runtime_error("Controlpp: QZ 2x1 swap: singular Sylvester system");
		}

		const Eigen::Matrix<T, 4, 1> z = lu.solve(rhs);

		if (!z.allFinite()) {
			throw std::runtime_error("Controlpp: QZ 2x1 swap: non-finite Sylvester solution");
		}

		const T residual = (M * z - rhs).norm();
		const T residual_scale = M.norm() * z.norm() + rhs.norm();

		if (residual > eps * residual_scale) {
			throw std::runtime_error("Controlpp: QZ 2x1 swap: Sylvester solve residual too large");
		}

		const Eigen::Matrix<T, 2, 1> R = z.template head<2>();
		const Eigen::Matrix<T, 2, 1> L = z.template tail<2>();

		// Build the right and left vectors associated with the 1x1 block.
		const Eigen::Matrix<T, 3, 1> X1({R(0),R(1),T(1)});
		const Eigen::Matrix<T, 3, 1> Y1({L(0),L(1),T(1)});

		// Construct orthogonal transformations whose first columns span
		// the right and left vectors, respectively.
		Eigen::HouseholderQR<Eigen::Matrix<T, 3, 1>> qrX(X1);
		Eigen::HouseholderQR<Eigen::Matrix<T, 3, 1>> qrY(Y1);

		const Eigen::Matrix<T, 3, 3> V = qrX.householderQ();
		const Eigen::Matrix<T, 3, 3> U = qrY.householderQ();

		// Apply the orthogonal transformations.
		const Eigen::Matrix<T, 3, 3> S1 = U.transpose() * Sb * V;
		const Eigen::Matrix<T, 3, 3> T1 = U.transpose() * Tb * V;

		// Restore generalized Schur form in the trailing 2x2 pair.
		const Eigen::Matrix<T, 2, 2> Strail = S1.template block<2, 2>(1, 1);
		const Eigen::Matrix<T, 2, 2> Ttrail = T1.template block<2, 2>(1, 1);

		Eigen::RealQZ<Eigen::Matrix<T, 2, 2>> localQZ;
		localQZ.compute(Strail, Ttrail);

		if (localQZ.info() != Eigen::Success) {
			throw std::runtime_error("Controlpp: QZ 2x1 swap: local 2x2 QZ decomposition failed");
		}

		// Embed the local QZ transformations in the trailing 2x2 block.
		Eigen::Matrix<T, 3, 3> Utrail = Eigen::Matrix<T, 3, 3>::Identity();
		Eigen::Matrix<T, 3, 3> Vtrail = Eigen::Matrix<T, 3, 3>::Identity();
		Utrail.template block<2, 2>(1, 1) = localQZ.matrixQ();
		Vtrail.template block<2, 2>(1, 1) = localQZ.matrixZ().transpose();

		const Eigen::Matrix<T, 3, 3> U_restored = U * Utrail;
		const Eigen::Matrix<T, 3, 3> V_restored = V * Vtrail;

		const QZSwap<T, 3> result{.V = V_restored, .U = U_restored};

		return result;
	}


	/**
	* @brief Computes orthogonal transformations that swap two adjacent 2x2
	*        diagonal blocks of a generalized real Schur matrix pair.
	*
	* Given a 4x4 matrix pair in generalized real Schur form,
	*
	* @f[
	* S_b = \begin{bmatrix} A & C \\ 0 & D \end{bmatrix},
	* \qquad
	* T_b = \begin{bmatrix} B & F \\ 0 & E \end{bmatrix},
	* @f]
	*
	* where @f$A,B,D,E,C,F \in \mathbb{R}^{2\times2}@f$, computes orthogonal
	* matrices @f$U,V \in \mathbb{R}^{4\times4}@f$ such that
	*
	* @f[
	* S_{\mathrm{new}} = U^\top S_b V,\qquad
	* T_{\mathrm{new}} = U^\top T_b V,
	* @f]
	*
	* and the generalized eigenvalues associated with the diagonal blocks
	* @f$(D,E)@f$ are moved to the leading diagonal block, while those
	* associated with @f$(A,B)@f$ are moved to the trailing diagonal block.
	*
	* The transformation is constructed by solving a coupled generalized
	* Sylvester equation to determine the invariant subspaces to exchange,
	* constructing orthogonal transformations using QR factorization, and
	* applying a local generalized real Schur decomposition to restore the
	* required block structure.
	*
	* The returned transformations preserve the generalized eigenvalues of
	* the matrix pair. The swap may be numerically unreliable or impossible
	* when the spectra of the two diagonal blocks overlap or are too close,
	* since the associated Sylvester system may be singular or ill-conditioned.
	*
	* @tparam T Scalar type used for the matrix calculations. Usually `double` or `float`.
	*
	* @param Sb The 4x4 matrix in generalized real Schur form, with two
	*           adjacent 2x2 diagonal blocks.
	* @param Tb The 4x4 matrix paired with @p Sb, in generalized real Schur
	*           form.
	*
	* @return A QZSwap containing the orthogonal left transformation @c U
	*         and right transformation @c V, satisfying
	*         @f$S_{\mathrm{new}}=U^\top S_bV@f$ and
	*         @f$T_{\mathrm{new}}=U^\top T_bV@f$.
	*
	* @throws std::runtime_error If the Sylvester system cannot be solved
	*         reliably, produces non-finite values, or the local generalized
	*         Schur decomposition fails.
	*
	* @note This routine operates on a local 4x4 matrix pair. For a larger
	*       generalized Schur decomposition, apply the returned transformations
	*       to the corresponding rows and columns of the full matrices and
	*       update the accumulated orthogonal transformation matrices
	*       consistently.
	*/
	template <class T>
	QZSwap<T, 4> qz_2x2_swap_matrix(
		const Eigen::Matrix<T, 4, 4>& Sb,
		const Eigen::Matrix<T, 4, 4>& Tb)
	{
		#define print_var(var) {std::cout << "--------\n" << #var << ":\n" << var << std::endl;}

		using Mat2 = Eigen::Matrix<T, 2, 2>;
		using Mat4 = Eigen::Matrix<T, 4, 4>;
		using Mat8 = Eigen::Matrix<T, 8, 8>;
		using Vec8 = Eigen::Matrix<T, 8, 1>;

		const T eps = std::numeric_limits<T>::epsilon();
		print_var(eps);

		// Partition the two matrix pencils:
		//
		// Sb = [ A  C ]       Tb = [ B  F ]
		//      [ 0  D ]            [ 0  E ]
		//
		const Mat2 A = Sb.template block<2, 2>(0, 0); print_var(A);
		const Mat2 C = Sb.template block<2, 2>(0, 2); print_var(C);
		const Mat2 D = Sb.template block<2, 2>(2, 2); print_var(D);

		const Mat2 B = Tb.template block<2, 2>(0, 0); print_var(B);
		const Mat2 F = Tb.template block<2, 2>(0, 2); print_var(F);
		const Mat2 E = Tb.template block<2, 2>(2, 2); print_var(E);

		// Solve:
		//     A*R - L*D = -C
		//     B*R - L*E = -F
		//
		// Unknown vector:
		//     z = [R00, R10, R01, R11,
		//          L00, L10, L01, L11]^T
		//
		// This is column-major vectorization.

		const Mat8 M({
			{ A(0, 0), 0,       A(0, 1), 0,       -D(0, 0), -D(1, 0), 0,        0        },
			{ A(1, 0), 0,       A(1, 1), 0,       0,        0,        -D(0, 0), -D(1, 0) },
			{ 0,       A(0, 0), 0,       A(0, 1), -D(0, 1), -D(1, 1), 0,        0        },
			{ 0,       A(1, 0), 0,       A(1, 1), 0,        0,        -D(0, 1), -D(1, 1) },

			{ B(0, 0), 0,       B(0, 1), 0,       -E(0, 0), -E(1, 0), 0,        0        },
			{ B(1, 0), 0,       B(1, 1), 0,       0,        0,        -E(0, 0), -E(1, 0) },
			{ 0,       B(0, 0), 0,       B(0, 1), -E(0, 1), -E(1, 1), 0,        0        },
			{ 0,       B(1, 0), 0,       B(1, 1), 0,        0,        -E(0, 1), -E(1, 1) }
		});

		const Vec8 rhs(
			-C(0, 0), -C(1, 0), -C(0, 1), -C(1, 1),
			-F(0, 0), -F(1, 0), -F(0, 1), -F(1, 1)
		);
		print_var(M);
		print_var(rhs);

		if (!M.allFinite() || !rhs.allFinite()) {
			throw std::runtime_error("qz_2x2_swap: non-finite Sylvester system");
		}

		// Scale the system before solving.
		const T scale = std::max(M.cwiseAbs().maxCoeff(), rhs.cwiseAbs().maxCoeff());
		print_var(scale);
		if (!(scale > T(0)) || !std::isfinite(scale)) {
			throw std::runtime_error("qz_2x2_swap: invalid Sylvester system scale");
		}

		const Mat8 Ms = M / scale; print_var(Ms);
		const Vec8 bs = rhs / scale; print_var(bs);

		Eigen::FullPivLU<Mat8> lu(Ms);

		if (lu.rank() < 8) {
			throw std::runtime_error("qz_2x2_swap: singular Sylvester system");
		}

		const Vec8 z = lu.solve(bs); print_var(z);

		if (!z.allFinite()) {
			throw std::runtime_error("qz_2x2_swap: non-finite Sylvester solution");
		}

		const T residual = (Ms * z - bs).norm(); print_var(residual);
		const T residual_scale = std::max(T(1), Ms.norm() * z.norm() + bs.norm()); print_var(residual_scale);

		if (residual > T(100) * eps * residual_scale) {
			throw std::runtime_error("qz_2x2_swap: inaccurate Sylvester solution");
		}

		const Mat2 R({
			{z(0), z(1)},
			{z(2), z(3)}
		});

		const Mat2 L({
			{z(4), z(5)},
			{z(6), z(7)}
		});
		print_var(R);
		print_var(L);

		// Construct bases for the subspaces to move to the front.
		const Eigen::Matrix<T, 4, 2> X({
			{R(0, 0), R(0, 1)},
			{R(1, 0), R(1, 1)},
			{1,       0      },
			{0,       1      }
		});

		const Eigen::Matrix<T, 4, 2> Y({
			{L(0, 0), L(0, 1)},
			{L(1, 0), L(1, 1)},
			{1,       0      },
			{0,       1      }
		});
		print_var(X);
		print_var(Y);

		// Orthogonal completions.
		Eigen::HouseholderQR<Eigen::Matrix<T, 4, 2>> qrX(X);
		Eigen::HouseholderQR<Eigen::Matrix<T, 4, 2>> qrY(Y);

		const Mat4 V = qrX.householderQ();
		const Mat4 U = qrY.householderQ();
		print_var(V);
		print_var(U);

		// Apply the initial orthogonal transformations.
		const Mat4 S1 = U.transpose() * Sb * V;
		const Mat4 T1 = U.transpose() * Tb * V;
		print_var(S1);
		print_var(T1);

		// Re-triangularize the leading 2x2 pair.
		// This preserves real generalized Schur form.
		const Mat2 Slead = S1.template block<2, 2>(0, 0);
		const Mat2 Tlead = T1.template block<2, 2>(0, 0);
		print_var(Slead);
		print_var(Tlead);

		Eigen::RealQZ<Mat2> localQZ;
		localQZ.compute(Slead, Tlead);

		if (localQZ.info() != Eigen::Success) {
			throw std::runtime_error("qz_2x2_swap: local RealQZ failed");
		}

		// Eigen's convention is S = Q * S_qz * Z.
		// Therefore the right transform for U^T*S*V is Z^T.
		const Mat4 Ulead({
			{localQZ.matrixQ()(0, 0), localQZ.matrixQ()(0, 1), 0, 0},
			{localQZ.matrixQ()(1, 0), localQZ.matrixQ()(1, 1), 0, 0},
			{0, 0, 1, 0},
			{0, 0, 0, 1}
		});

		const Mat4 Vlead({
			{localQZ.matrixZ()(0, 0), localQZ.matrixZ()(1, 0), 0, 0},
			{localQZ.matrixZ()(0, 1), localQZ.matrixZ()(1, 1), 0, 0},
			{0, 0, 1, 0},
			{0, 0, 0, 1}
		});

		print_var(Ulead);
		print_var(Vlead);


		// Apply the leading-block QZ transformations.
		const Mat4 Urestored = U * Ulead;
		const Mat4 Vrestored = V * Vlead;

		const Mat4 S2 = Urestored.transpose() * Sb * Vrestored;
		const Mat4 T2 = Urestored.transpose() * Tb * Vrestored;

		// Re-triangularize the trailing 2x2 pair.
		const Mat2 Strail = S2.template block<2, 2>(2, 2);
		const Mat2 Ttrail = T2.template block<2, 2>(2, 2);

		Eigen::RealQZ<Mat2> trailingQZ;
		trailingQZ.compute(Strail, Ttrail);

		if (trailingQZ.info() != Eigen::Success) {
			throw std::runtime_error("qz_2x2_swap: trailing RealQZ failed");
		}

		// Embed the trailing QZ transformations.
		const Mat4 Utrail({
			{1, 0, 0, 0},
			{0, 1, 0, 0},
			{0, 0, trailingQZ.matrixQ()(0, 0),
				trailingQZ.matrixQ()(0, 1)},
			{0, 0, trailingQZ.matrixQ()(1, 0),
				trailingQZ.matrixQ()(1, 1)}
		});

		const Mat4 Vtrail({
			{1, 0, 0, 0},
			{0, 1, 0, 0},
			{0, 0, trailingQZ.matrixZ()(0, 0),
				trailingQZ.matrixZ()(1, 0)},
			{0, 0, trailingQZ.matrixZ()(0, 1),
				trailingQZ.matrixZ()(1, 1)}
		});

		// Compose transformations in application order.
		const Mat4 Ufinal = Urestored * Utrail;
		const Mat4 Vfinal = Vrestored * Vtrail;

		return {.V = Vfinal, .U = Ufinal};
	}


    /**
	 * @brief Solves the discrete time riccati equation (DARE)
	 * 
	 * This function computes the stabilizing symmetric solution of the DARE:
     * 
	 * \f[
	 * A^\top X A - A^\top X B (R + B^\top X B)^{-1} B^\top X A + Q = 0
	 * \f]
	 * 
	 * 
	 * ------
	 * 
	 * This algorithm uses the ordered QZ decomposition with generalized Schur form to solve the DARE.
	 * 
	 * \tparam T Scalar type (e.g., `double`, `float`).
	 * \tparam NStates Number of states.
	 * \tparam NInputs Number of control inputs.
	 * 
	 * \param A State matrix (\f$n \times n\f$).
	 * \param B Input matrix (\f$n \times m\f$).
	 * \param R Input weighting matrix (\f$m \times m\f$, symmetric positive definite).
	 * \param Q State weighting matrix (\f$n \times n\f$, symmetric positive semidefinite).
	 * 
	 * @returns The solution of the DARE as an Eigen::Matrix with the dimensions `NStates x NStates`.
	 * @throws std::runtime_error if QZ failed to converge
	 */
	template<class T, int NStates, int NInputs,
			int AOpt, int AMaxR, int AMaxC,
			int BOpt, int BMaxR, int BMaxC,
			int ROpt, int RMaxR, int RMaxC,
			int QOpt, int QMaxR, int QMaxC
  	>
	Eigen::Matrix<T, NStates, NStates> dare_qz_solver(
		const Eigen::Matrix<T, NStates, NStates, AOpt, AMaxR, AMaxC>& A,
		const Eigen::Matrix<T, NStates, NInputs, BOpt, BMaxR, BMaxC>& B,
		const Eigen::Matrix<T, NStates, NStates, QOpt, QMaxR, QMaxC>& Q,
		const Eigen::Matrix<T, NInputs, NInputs, ROpt, RMaxR, RMaxC>& R
	){
		// 1. Build the symplectic pencil matrices H and J
		constexpr int PencilSize = 2 * NStates + NInputs;
		using Pencil = Eigen::Matrix<T, PencilSize, PencilSize>;

		const Pencil H = build_dare_symplectic_pencil_H(A, B, Q, R);

		const Pencil J = build_dare_symplectic_pencil_J(A, B);

		// 2. Solve the generalized eigenvalue problem H v = lambda J v
		Eigen::RealQZ<Pencil> qz;
		qz.compute(H, J);

		if(qz.info() != Eigen::Success){
			// QZ failed to converge
			throw std::runtime_error("Controlpp: dare_qz_solver(): QZ failed to converge.");
		}

		Pencil S = qz.matrixS(); // is quasi-upper-triangular
		Pencil Tmat = qz.matrixT(); // is upper-triangular
		Pencil Qz = qz.matrixQ(); // contains the left generalized Schur vectors
		Pencil Z = qz.matrixZ(); // contains the right generalized Schur vectors

		// 3. Reorder
		int block_start = 0;
        int stable_region = 0;
		while(block_start < PencilSize){
			
			const int block_size = pencil_block_size(S, block_start);

			// map to eigen::map for easier access
			auto Sb = S.block(block_start, block_start, block_size, block_size);
			auto Tb = Tmat.block(block_start, block_start, block_size, block_size);

			const bool is_stable = is_stable_pencil_block(Sb, Tb);

			// if the block is not stable, we need to swap it with a stable block
			// so that all stable blocks are in the top left corner of the matrix
			// @todo


			block_start += block_size;
		}

	}

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
		return dare_eigenvector_solver(A, B, Q, R);
	}

}