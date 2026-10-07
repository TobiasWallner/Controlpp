#pragma once

/**
 * @file dare.hpp
 * @brief This file contains the implementation of the Discrete-time Algebraic Riccati Equation (DARE) solver.
 * @author Tobias Wallner
 * @copyright 2026 Tobias Wallner
 */

#include <complex>

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