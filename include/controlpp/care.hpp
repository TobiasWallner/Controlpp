#pragma once

/**
 * @file care.hpp
 * @brief This file contains the implementation of the Continuous-time Algebraic Riccati Equation (CARE) solver.
 * @author Tobias Wallner
 * @copyright 2026 Tobias Wallner
 */

// std
#include <complex>

// linear algebra library
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
		return hamilton_solver(H);
	}
}