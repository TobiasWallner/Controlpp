#pragma once

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>

// controlpp
#include <controlpp/math.hpp>
#include <controlpp/ContinuousStateSpace.hpp>

namespace controlpp
{
    template<class T, int NStates, int NInputs, int NPerfOutputs, int NMeasOutputs, int NDisturbances>
    constexpr ContinuousStateSpace<T, NStates, NMeasOutputs, NInputs> continous_H2_controller(
            const Eigen::Matrix<T, NStates, NStates>& A,
            const Eigen::Matrix<T, NStates, NDisturbances>& Bw,
            const Eigen::Matrix<T, NStates, NInputs>& Bu,
            const Eigen::Matrix<T, NPerfOutputs, NStates>& Cz,
            const Eigen::Matrix<T, NPerfOutputs, NInputs>& Duz,
            const Eigen::Matrix<T, NMeasOutputs, NStates>& Cy,
            const Eigen::Matrix<T, NMeasOutputs, NDisturbances>& Dwy,
            const Eigen::Matrix<T, NInputs, NInputs>& R,
            const Eigen::Matrix<T, NMeasOutputs, NMeasOutputs>& S
    ){   
        
        // Weighting Matrices
        const Eigen::Matrix<T, NStates, NStates> Q = Cz.transpose() * Cz;
        const Eigen::Matrix<T, NStates, NStates> W = Bw * Bw.transpose();


        // State feedback
        // A^\top X + X A - X B_u R^{-1} B_u^\top X + Q = 0
        const Eigen::Matrix<T, NStates, NInputs> N1 = Cz.transpose() * Duz;
        const Eigen::Matrix<T, NStates, NStates> X = controlpp::care_solver(A, Bu, R, Q, N1);

        // Estimator Feedback
        // A Y + Y A^\top - Y C_y^\top S^{-1} C_y Y + W = 0
        const Eigen::Matrix<T, NStates, NMeasOutputs> N2 = Bw * Dwy.transpose();
        const Eigen::Matrix<T, NStates, NStates> Y = controlpp::care_solver(A.transpose().eval(), Cy.transpose().eval(), S, W, N2);

        // State gain
        // F = -R^{-1} \left( B_u^\top X + D_{1u}^\top C_z \right)
        const Eigen::Matrix<T, NInputs, NStates> F = -R.llt().solve(Bu.transpose() * X + Duz.transpose() * Cz);

        // Estimator gain
        // L = - \left( Y C_y^\top + B_w D_{2w}^\top \right) S^{-1}
        const Eigen::Matrix<T, NStates, NMeasOutputs> M = Y * Cy.transpose() + Bw * Dwy.transpose();
        const Eigen::Matrix<T, NStates, NMeasOutputs> L = -S.llt().solve(M.transpose()).transpose();

        // Construct the H2 Controller
        const Eigen::Matrix<T, NStates, NStates> A_K = A + Bu * F + L * Cy;
        const Eigen::Matrix<T, NStates, NMeasOutputs> B_K = -L;
        const Eigen::Matrix<T, NInputs, NStates> C_K = F;

        ContinuousStateSpace<T, NStates, NMeasOutputs, NInputs> H2(A_K, B_K, C_K);
        return H2;
    }

    /**
     * \brief An optimal controller for linear time-invariant (LTI) system. It minimizes the H2 norm of the the disturbance input gain to the output.
     * 
     * An H2 controller is part of the class of optimal controllers. 
     * Where the typical system model is a linear time-invaritan (LTI) system, that looks like:
     * 
     * \f[
     * \dot{x} = A x + B_w w + B_u u \\
     * z = C_z x + D_{1u} u \\
     * y = C_y x + D_{2w} w
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
     *    - \f$B_w\f$ describes how the disturbance \f$w\f$ affects the systems states
     *    - \f$B_u\f$ describes how the inputs \f$u\f$ affect the system states
     *    - \f$C_z\f$ describes how the system states result in the output of the system that should be controlled
     *    - \f$C_y\f$ describes how the states generate the system output \f$y\f$ that we can measure
     *    - \f$D_{1u}\f$ describes how the inputs directly affect the performance output
     *    - \f$D_{2w}\f$ describes how the disturbance affects the measurement output \f$y\f$
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
     * A^\top X + X A - X B_u R^{-1} B_u^\top X + Q = 0
     * \f]
     * 
     * with 
     *  - \f$Q = C_z^\top C_z\f$ the performance weight, aka. state cost matrix. It penalizes state deviations.
     *  - \f$R = D_{1u}^\top D_{1u}\f$ is the input weight, aka. control cost matrix. It penalizes control effort.
     * 
     * and the estimator riccati equation:
     * 
     * \f[
     * A^\top Y + Y A - Y C_y S^{-1} C_y^\top Y + W = 0
     * \f]
     * 
     * where: 
     *  - \f$W = B_w B_w^\top\f$ disturbance weight
     *  - \f$S = D_{2w} D_{2w}^\top\f$ measurement noise weight. 
     * 
     * ----
     * 
     * X and Y from the riccati equations are then used to calculate the
     * 
     * optimal state-feedback gain:
     * 
     * \f[
     * F = -R^{-1} \left( B_u^\top X + D_{1u}^\top C_z \right)
     * \f]
     * 
     * and the optimal observer gian:
     * 
     * \f[
     * L = - \left( Y C_y^\top + B_w D_{2w}^\top \right) S^{-1}
     * \f]
     * 
     * ----
     * 
     * The H2 controller then is:
     * 
     * \f[
     * A_K = A + B_u F + L C_y \\
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
     * \param A System/plant dynamic
     * \param Bw disturbance input
     * \param Bu control input
     * \param Cz state to performance output
     * \param Cy state to measurement output
     * \param Duz throughput from control input to performance output 
     * \param Dwy throughput from disturbance to measurement output
     * \param r  Input weightin. penalizes control effort u. 
     *  - Large r: controller is gentle (low gain)
     *  - Small r: controller is aggressive (high gain)
     * \param s Measurement noise / estimator weighting. penalizes measurement trust:
     *  - Large s: controller assumes the "measurements are noisy". It does not trust the measurements and relies on the model.
     *  - Small s: controller assumes the "measurements are clean". It trusts the measurements more and the estimator reacts stronger to sensor signals.
     */
    template<class T, int NStates, int NInputs, int NPerfOutputs, int NMeasOutputs, int NDisturbances>
    requires(NInputs>1 && NMeasOutputs>1)
    constexpr ContinuousStateSpace<T, NStates, NMeasOutputs, NInputs> continous_H2_controller(
            const Eigen::Matrix<T, NStates, NStates>& A,
            const Eigen::Matrix<T, NStates, NDisturbances>& Bw,
            const Eigen::Matrix<T, NStates, NInputs>& Bu,
            const Eigen::Matrix<T, NPerfOutputs, NStates>& Cz,
            const Eigen::Matrix<T, NPerfOutputs, NInputs>& Duz,
            const Eigen::Matrix<T, NMeasOutputs, NStates>& Cy,
            const Eigen::Matrix<T, NMeasOutputs, NDisturbances>& Dwy,
            const Eigen::Vector<T, NInputs> r,
            const Eigen::Vector<T, NMeasOutputs> s
    ){   
        
        // Weighting Matrices
        const Eigen::Matrix<T, NInputs, NInputs> R = (Duz.transpose() * Duz + Eigen::DiagonalMatrix<T, NInputs>(r).toDenseMatrix()).eval();
        const Eigen::Matrix<T, NMeasOutputs, NMeasOutputs> S = (Dwy * Dwy.transpose() + Eigen::DiagonalMatrix<T, NMeasOutputs>(s).toDenseMatrix()).eval();

        return continous_H2_controller(A, Bw, Bu, Cz, Duz, Cy, Dwy, R, S);
    }

    /**
     * 
     * 
     */
    template<class T, std::convertible_to<T> U1, std::convertible_to<T> U2, int NStates, int NInputs, int NPerfOutputs, int NMeasOutputs, int NDisturbances>
    constexpr ContinuousStateSpace<T, NStates, NMeasOutputs, NInputs> continous_H2_controller(
            const Eigen::Matrix<T, NStates, NStates>& A,
            const Eigen::Matrix<T, NStates, NDisturbances>& Bw,
            const Eigen::Matrix<T, NStates, NInputs>& Bu,
            const Eigen::Matrix<T, NPerfOutputs, NStates>& Cz,
            const Eigen::Matrix<T, NPerfOutputs, NInputs>& Duz,
            const Eigen::Matrix<T, NMeasOutputs, NStates>& Cy,
            const Eigen::Matrix<T, NMeasOutputs, NDisturbances>& Dwy,
            const U1& r_scalar = 1,
            const U2& s_scalar = 1e-3
    ){

        // Weighting Vectors
        const Eigen::Vector<T, NInputs> r = Eigen::Vector<T, NInputs>::Constant(static_cast<T>(r_scalar));
        const Eigen::Vector<T, NMeasOutputs> s = Eigen::Vector<T, NMeasOutputs>::Constant(static_cast<T>(s_scalar));

        return continous_H2_controller(A, Bw, Bu, Cz, Duz, Cy, Dwy, r, s);
    }

    /**
     * \brief A continuous generalised plant model
     */
    template<class T, int NStates, int NInputs=1, int NPerfOutputs=1, int NMeasOutputs=1, int NDisturbances=1>
    struct ContinuousGeneralisedPlant{
        Eigen::Matrix<T, NStates, NStates> A;
        Eigen::Matrix<T, NStates, NDisturbances> Bw;
        Eigen::Matrix<T, NStates, NInputs> Bu;
        Eigen::Matrix<T, NPerfOutputs, NStates> Cz;
        Eigen::Matrix<T, NPerfOutputs, NInputs> Duz;
        Eigen::Matrix<T, NMeasOutputs, NStates> Cy;
        Eigen::Matrix<T, NMeasOutputs, NDisturbances> Dwy;
    };

    template<class T, int NStates, int NInputs=1, int NPerfOutputs=1, int NMeasOutputs=1, int NDisturbances=1>
    std::ostream& operator<<(std::ostream& stream, const ContinuousGeneralisedPlant<T, NStates, NInputs, NPerfOutputs, NMeasOutputs, NDisturbances>& G){
        stream << "A:\n" << G.A << "\n";
        stream << "Bw:\n" << G.Bw << "\n";
        stream << "Bu:\n" << G.Bu << "\n";
        stream << "Cz:\n" << G.Cz << "\n";
        stream << "Duz:\n" << G.Duz << "\n";
        stream << "Cy:\n" << G.Cy << "\n";
        stream << "Dwy:\n" << G.Dwy << "\n";
        return stream;
    };

    /**
     * \brief Constructs a continuous H2 controller from a continuous state space plant model
     * 
     * \param Gss Generalised plant in state space form
     * 
     * \tparam T The value type of the plant and controller. Usually `double` or `float`.
     * \tparam NPlantOutputs The number of states of the plant. Also the number of states of the controller if no extra weighting functions are applied.
     * \tparam NPlantInputs The number of inputs of the plant.
     * 
     * \see controlpp::continous_H2_controller()
     * \see controlpp::ContinuousGeneralisedPlant
     */
    template<class T, int NStates, int NInputs, int NPerfOutputs, int NMeasOutputs, int NDisturbances>   
    constexpr ContinuousStateSpace<T, NStates, NMeasOutputs, NInputs> H2_controller(const ContinuousGeneralisedPlant<T, NStates, NInputs, NPerfOutputs, NMeasOutputs, NDisturbances> & Gss){
        return continous_H2_controller(
            Gss.A, Gss.Bw, Gss.Bu, 
            Gss.Cz, Gss.Duz, 
            Gss.Cy, Gss.Dwy
        );
    }

} // namespace controlpp
