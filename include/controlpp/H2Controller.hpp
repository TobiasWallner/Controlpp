#pragma once

#include <controlpp/math.hpp>

namespace controlpp
{
    
    /**
	 * \brief Solves the continuous Riccati equation for the H2 control problem
	 * 
	 * Solves the following LQR Riccati equation for \f$X\f$:
	 * 
	 * \f[
	 * A^\top X + X A - (X B + C^\top D) R^{-1} (B^\top X + D^\top C) + C^\top C = 0
	 * \f]
	 * 
	 * for the system:
	 * 
	 * \f[
	 * \dot{x} = A x + B u\\
	 * \f]
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
	 * \param A 
	 * \param B2 
	 * \param C1
	 * \param D12 
	 * 
	 * \tparam T The data type of the matrix elements
	 * \tparam N The number of states
	 * \tparam M The number of imputs
	 * 
	 * \see controlpp::solve_continuous_h2_estimator_riccati()
	 */
	template<
		class T, int N, int M,
		int AOptions, int AMaxRows, int AMaxCols,
		int BOptions, int BMaxRows, int BMaxCols,
		int COptions, int CMaxRows, int CMaxCols,
		int DOptions, int DMaxRows, int DMaxCols
	>
	constexpr Eigen::Matrix<T, N, N> solve_continuous_h2_control_riccati(
		const Eigen::Matrix<T, N, N, AOptions, AMaxRows, AMaxCols>& A,
		const Eigen::Matrix<T, N, M, BOptions, BMaxRows, BMaxCols>& B2,
		const Eigen::Matrix<T, M, M, COptions, CMaxRows, CMaxCols>& C1,
		const Eigen::Matrix<T, N, N, DOptions, DMaxRows, DMaxCols>& D12
	){
		// 1. Build the hamilton matrix
		const auto Q = C1.transpose() * C1;
        const auto R = D12.transpose() * D12;

		const auto R_inv_DT_C = R.llt().solve(D12.transpose() * C1).eval();
		const auto A_ = (A - B2 * R_inv_DT_C).eval();

		Eigen::Matrix<T, 2*N, 2*N> H;
		H.topLeftCorner(N, N) = A_;
		H.topRightCorner(N, N) =  - B2 * R.llt().solve(B2.transpose());
		H.bottomLeftCorner(N, N) = - (Q - C1.transpose() * D12 * R_inv_DT_C);
		H.bottomRightCorner(N, N) = -A_.transpose();

		// 2. Compute stable eigenvectors
		Eigen::ComplexEigenSolver<Eigen::Matrix<T, 2*N, 2*N>> ces;
		ces.compute(H);
		const auto& H_eigvals = ces.eigenvalues();
		const auto& H_eigvecs = ces.eigenvectors();

		// in a 2n hamilton matrix are exactly n stable ones
		Eigen::Matrix<std::complex<T>, 2*N, N> StableEigenVecs;
		int si = 0; // stable eigenvalue iterator
		int ei = 0; // eigen value iterator
		for(; (si < N) && (ei < 2*N); ++ei){
			if(H_eigvals(ei).real() < 0){// only stable ones
				StableEigenVecs.col(si) = H_eigvecs.col(ei);
				++si;
			}
		}

		// 3. Repartition
		const auto TopPartition = StableEigenVecs.template block<N, N>(0, 0);
		const auto BottomPartition = StableEigenVecs.template block<N, N>(N, 0);

		// 4. Recover Solution (BottomPartition * TopPartition^-1)
		const Eigen::Matrix<std::complex<T>, N, N> X = TopPartition.transpose().partialPivLu().solve(BottomPartition.transpose()).transpose();

		// The result is real, make sure it is real because it should be
		const Eigen::Matrix<T, N, N> realX = X.real();

		return realX;
	}

	template<
		class T, int N, int M,
		int AOptions, int AMaxRows, int AMaxCols,
		int BOptions, int BMaxRows, int BMaxCols,
		int COptions, int CMaxRows, int CMaxCols,
		int DOptions, int DMaxRows, int DMaxCols
	>
	constexpr Eigen::Matrix<T, N, N> solve_continuous_h2_estimator_riccati(
		const Eigen::Matrix<T, N, N, AOptions, AMaxRows, AMaxCols>& A,
		const Eigen::Matrix<T, N, M, BOptions, BMaxRows, BMaxCols>& B1,
		const Eigen::Matrix<T, M, M, COptions, CMaxRows, CMaxCols>& C2,
		const Eigen::Matrix<T, N, N, DOptions, DMaxRows, DMaxCols>& D21
	){
		// 1. Build the hamilton matrix
		const auto W = B1 * B1.transpose();
        const auto S = D21 * D21.transpose();

		const auto S_inv_C2 = S.llt().solve(C2).eval();
		const auto A_ = (A - B1 * (D21.transpose() * S_inv_C2)).eval();

		Eigen::Matrix<T, 2*N, 2*N> H;
		H.topLeftCorner(N, N) = A_;
		H.topRightCorner(N, N) =  -(W - B1 * (D21.transpose() * S.llt().solve(D21 * B1.transpose())));
		H.bottomLeftCorner(N, N) = - C2.transpose() * S_inv_C2;
		H.bottomRightCorner(N, N) = -A_.transpose();

		// 2. Compute stable eigenvectors
		Eigen::ComplexEigenSolver<Eigen::Matrix<T, 2*N, 2*N>> ces;
		ces.compute(H);
		const auto& H_eigvals = ces.eigenvalues();
		const auto& H_eigvecs = ces.eigenvectors();

		// in a 2n hamilton matrix are exactly n stable ones
		Eigen::Matrix<std::complex<T>, 2*N, N> StableEigenVecs;
		int si = 0; // stable eigenvalue iterator
		int ei = 0; // eigen value iterator
		for(; (si < N) && (ei < 2*N); ++ei){
			if(H_eigvals(ei).real() < 0){// only stable ones
				StableEigenVecs.col(si) = H_eigvecs.col(ei);
				++si;
			}
		}

		// 3. Repartition
		const auto TopPartition = StableEigenVecs.template block<N, N>(0, 0);
		const auto BottomPartition = StableEigenVecs.template block<N, N>(N, 0);

		// 4. Recover Solution (BottomPartition * TopPartition^-1)
		const Eigen::Matrix<std::complex<T>, N, N> X = TopPartition.transpose().partialPivLu().solve(BottomPartition.transpose()).transpose();

		// The result is real, make sure it is real because it should be
		const Eigen::Matrix<T, N, N> realX = X.real();

		return realX;
	}

    /**
     * \brief An optimal controller for linear time-invariant (LTI) system. It minimizes the H2 norm of the the disturbance input gain to the output.
     * 
     * An H2 controller is part of the class of optimal controllers. 
     * Where the typical system model is a linear time-invaritan (LTI) system, that looks like:
     * 
     * \f[
     * \dot{x} = A x + B_1 w + B_2 u \\
     * z = C_1 x + D_{1} u \\
     * y = C_2 x + D_{2} w
     * \f]
     * 
     * Where:
     *  + inputs:
     *    - \f$w\f$ is the disturbance
     *    - \f$u\f$ are the inputs of the system
     * 
     *  + Outputs/States
     *    - \f$x\f$ describes the states of the system
     *    - \f$z\f$ is the performance output (this is the output of which we want to minimize the variance)
     *    - \f$y\f$ is the output that we can actually measure (what the controller sees)
     * 
     *  + Parameters
     *    - \f$A\f$ describes the dynamic of the system/plant
     *    - \f$B_1\f$ describes how the disturbance \f$w\f$ affects the systems states
     *    - \f$B_2\f$ describes how the inputs \f$u\f$ affect the system states
     *    - \f$C_1\f$ describes how the system states result in the output of the system that should be controlled
     *    - \f$C_2\f$ describes how the states generate the system output \f$y\f$ that we can measure
     *    - \f$D_1\f$ describes how the inputs directly affect the performance output
     *    - \f$D_2\f$ describes how the disturbance affects the measurement output \f$y\f$
     * 
     * ----
     * 
     * Note how \f$D_{11}\f$, and \f$D_{22}\f$ are assumed to be zero.
     * 
     * ----
     * 
     * The H2 controller solves a system such that the closed-loop transfer function from the disturbance \f$w\f$ 
     * to the performace output \f$z\f$ has a minimized variance.
     * 
     * The dynamic controller of such a system has the form of:
     * 
     * \f[
     * \dot{x_K} = A_K + x_K + B_K y \\
     * u = C_K x_K + D_K y
     * \f]
     * 
     * ----
     * 
     * The optimal solution can be found solveing the Riccati equations
     * for the state-feedback riccati equation:
     * 
     * \f[
     * A^\top X + X A - X B_2 R^{-1} B_2^\top X + Q = 0
     * \f]
     * 
     * with 
     *  - \f$Q = C_1^\top C_1\f$ the performance weight, aka. state cost matrix. It penalizes state deviations.
     *  - \f$R = D_1^\top D_1\f$ is the input weight, aka. control cost matrix. It penalizes control effort.
     * 
     * and the estimator riccati equation:
     * 
     * \f[
     * A^\top Y + Y A - Y C_2 S^{-1} C_2^\top Y + W = 0
     * \f]
     * 
     * where: 
     *  - \f$W = B_1 B_1^\top\f$ disturbance weight
     *  - \f$S = D_2 D_2^\top\f$ measurement noise weight. 
     * 
     * ----
     * 
     * X and Y from the riccati equations are then used to calculate the
     * 
     * optimal state-feedback gain:
     * 
     * \f[
     * F = -R^{-1} \left( B_2^\top X + D_1^\top C_1 \right)
     * \f]
     * 
     * and the optimal observer gian:
     * 
     * \f[
     * L = - \left( Y C_2^\top + B_1 D_2^\top \right) S^{-1}
     * \f]
     * 
     * ----
     * 
     * The H2 controller then is:
     * 
     * \f[
     * A_K = A + B_2 F + L C_2 \\
     * B_K = -L
     * C_K = F
     * D_K = 0
     * \f]
     * 
     * ----
     * 
     * The default H2 controller minimizes the total energy over all frequencies.
     * To emphasise the control of a certain frequency region one can define a weighting transfer function \f$W_z(s)\f$
     * that affects the performance output \f$z\f$:
     * 
     * \f[
     * z_w = W_z(s) z
     * \f]
     * 
     * I \f$W_z\f$ is large at a certain frequency the optimizer will try to reduce the gain in that band more 
     * aggressively than in the band that W is smaller.
     * 
     * \param A
     * \param B1
     * \param B2
     * \param C1
     * \param C2
     * \param D1
     * \param D2
     * 
     * 
     */
    template<class T
        int NStates, int NDisturbances, int NInputs, int NPerfOutputs, int NMeasOutputs,
        int AOptions,int B1Options, int B2Options, int C1Options, int D1Options, int C2Options, int D2Options,
    >
    ContinuousStateSpace<T, NStates, NMeasOutputs, NInputs> continous_H2_controller(
                const Eigen::Matrix<T, NStates, NStates, AOptions>& A,
                const Eigen::Matrix<T, NStates, NDisturbances, B1Options>& B1,
                const Eigen::Matrix<T, NStates, NInputs, B2Options>& B2,
                const Eigen::Matrix<T, NPerfOutputs, NStates, C1Options>& C1,
                const Eigen::Matrix<T, NPerfOutputs, NInputs, D1Options>& D12,
                const Eigen::Matrix<T, NMeasOutputs, NStates, C2Options>& C2,
                const Eigen::Matrix<T, NMeasOutputs, NDisturbances, D2Options>& D21)
        {   
            // Weighting Matrices
            const auto R = D12.transpose() * D12;
            const auto S = D21 * D21.transpose();

            // State feedback
            // A^\top X + X A - X B_2 R^{-1} B_2^\top X + Q = 0
            const auto X = controlpp::solve_continuous_h2_control_riccati(A, B2, C1, D12);

            // Estimator Feedback
            // A Y + Y A^\top - Y C_2^\top S^{-1} C_2 Y + W = 0
            const auto Y = controlpp::solve_continuous_h2_estimator_riccati(A, B1, C2, D21);

            // State gain
            // F = -R^{-1} \left( B_2^\top X + D_1^\top C_1 \right)
            const auto F = -R.llt().solve(B2.transpose() * X + D12.transpose() * C1);

            // Estimator gain
            // L = - \left( Y C_2^\top + B_1 D_2^\top \right) S^{-1}
            const auto M = Y * C2.transpose() + B1 * D21.transpose();
            const auto L = - S.llt().solve(M.transpose()).transpose();

            // Construct the H2 Controller
            const Eigen::Matrix<T, NStates, NStates> A_K = A + B2 * F + L * C2;
            const Eigen::Matrix<T, NStates, NMeasOutputs> B_K = -L;
            const Eigen::Matrix<T, NInputs, NStates> C_K = F;

            ContinuousStateSpace<T> result(A_K, B_K, C_K);
            return result;
        }


        /**
         * \brief Constructs a continuous H2 controller from a continuous state space plant model
         * 
         * \param 
         * 
         * \tparam T The value type of the plant and controller. Usually `double` or `float`.
         * \tparam NPlantOutputs The number of states of the plant. Also the number of states of the controller if no extra weighting functions are applied.
         * \tparam NPlantInputs The number of inputs of the plant.
         * 
         * \see controlpp::continous_H2_controller()
         */
        template<class T, int N, int M, int K>
        ContinuousStateSpace<T, NStates, NPlantOutputs, NPlantInputs> H2_controller(const ContinuousStateSpace<T, NPlantInputs, NPlantOutputs>& Gss){
            return continous_H2_controller(Gss.A(), )
        }

} // namespace controlpp
