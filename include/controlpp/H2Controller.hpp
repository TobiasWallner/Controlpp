#pragma once

namespace controlpp
{
    

    /**
     * \brief   
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
