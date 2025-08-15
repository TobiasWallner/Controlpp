#pragma once

namespace controlpp
{
    

    /**
     * \brief An optimal controller for linear time-invariant (LTI) system. It minimizes the H2 norm of the the disturbance input gain to the output.
     * 
     * An H2 controller is part of the class of optimal controllers. 
     * Where the typical system model is a linear time-invaritan (LTI) system, that looks like:
     * 
     * \f[
     * \dot{x} = A x + B_1 w + B_2 u \\
     * z = C_1 x + D_{1} u \\
     * y = C_2 x + D_2 w
     * \f]
     * 
     * Where:
     *  - \f$A\f$ describes the dynamic of the system/plant
     *  - \f$x\f$ describes the states of the system
     *  - \f$B_1\f$ describes how the disturbance \f$w\f$ affects the systems states
     *  - \f$w\f$ is the disturbance
     *  - \f$B_2\f$ describes how the inputs \f$u\f$ affect the system states
     *  - \f$u\f$ are the inputs of the system
     *  - \f$z\f$ is the performance output (this is the output of which we want to minimize the variance)
     *  - \f$C_1\f$ describes how the system states result in the output of the system that should be controlled
     *  - \f$D_1\f$ describes how the inputs directly affect the performance output
     *  - \f$C_2\f$ describes how the states generate the system output \f$y\f$ that we can measure
     *  - \f$D2\f$ describes how the disturbance affects the measurement output \f$y\f$
     *  - \f$y\f$ is the output that we can actually measure (what the controller sees)
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
     */
    template<class T>
    class H2Controller{
        private:

        public:

        /**
         * \brief Constructs a H2 controller
         * 
         * \param A System/Plant dynamics matrix, with size \f$n \times n\f$ where n is the number of states.
         * It describes how the system states evolves based on the systems state.
         * \param B2 Control input matrix with size \f$n \times m\f$, where m is the size of the input/measurement vector.
         * It describes how the the control input alters the systems state.
         * \param C2 Measurement output matrix with size \f$p \times n\f$.
         * It describes how the systems states map to the outputs that can be measured and used by the estimator.
         * \param R Control cost matrix, with size \f$m \times m\f$.
         * It penalizes the magnitude of control effort.
         * \param Q State cost matrix, with size \f$n \times n\f$.
         * It penalizes the deviation of the state.
         * \param W Process noise covariance matrix, with size \f$n \times n\f$. Has to be symetric positive semidefinite!
         * It models the process disturbances affecting the state.
         */
        H2Controller(
            const Eigen::Matrix<T>& A, 
            const Eigen::Matrix<T>&B2, 
            const Eigen::Matrix<T>&C2, 
            const Eigen::Matrix<T>&R, 
            const Eigen::Matrix<T>&Q, 
            const Eigen::Matrix<T>&V, 
            const Eigen::Matrix<T>&W){

        }
    }

} // namespace controlpp
