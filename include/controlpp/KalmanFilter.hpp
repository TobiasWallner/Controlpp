#pragma once

// std
#include <concepts>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>

namespace controlpp
{
    
    /**
     * \brief Kalman filter
     * 
     * The Kalman filter is an efficient, recursive algorithm for estimating the internal state of a dynamic system from a series of noisy, 
     * incomplete measurements. It combines predictions from a mathematical model with real-time sensor data to produce an estimate that is 
     * statistically optimal under certain assumptions.
     * 
     * ----
     * 
     * 1. Predict state and error covariance
     * 
     * \f[
     * \hat{x}_{k} = \mathbf{A} x_{k-1}\\
     * \hat{\mathbf{P}}_{k} = \mathbf{A} \mathbf{P}_{k-1} \mathbf{A}^{\top} + \mathbf{Q}\\
     * \f]
     * 
     * 2. Calculate Kalman Gain
     * 
     * \f[
     * \mathbf{K}_{k} = \mathbf{P}_{k-1} \mathbf{H}^{\top} \left( \mathbf{H} \mathbf{P}_{k-1} \mathbf{H}^{\top} + \mathbf{R} \right)^{-1}
     * \f]
     * 
     * 3. Update state estimate and error covariance
     * 
     * \f[
     * x_{k} = \hat{x}_{k} + \mathbf{K}_{k} \left( z_k - \mathbf{H} \hat{x}_{k} \right)\\
     * \mathbf{P}_{k} = (I - \mathbf{K}_{k} \mathbf{H}) \hat{\mathbf{P}}_{k}
     * \f]
     * 
     * with:
     * 
     * - \f$\mathbf{A}\f$ The system matrix, aka. state transition matrix  
     *   You get it from the systems identification/modeling/transfer function/state space
     * 
     * - \f$\mathbf{B}\f$ The control input matrix  
     *   You get if from the systems identification/modeling/transfer function/state space
     * 
     * - \f$\mathbf{H}\f$ The observation matrix  
     *   The observation matrix maps the state vector x to the measuremnts z:  
     *   \f[
     *     z_k = \mathbf{H} x_k
     *   \f]
     * 
     * - \f$\mathbf{R}\f$ The measurement noise covariance matrix  
     *   Captures the uncertainty in the measurement z.
     *   
     *   Diagonal elements represent the variance (= squared standar deviation) of the measurement noise.  
     *   - Large diagonal elements \f$\Rightarrow\f$ higher trust in model  
     *   - Small diagonal elements \f$\Rightarrow\f$ higher trust in measurements  
     *  
     *   Off diagonal elements tell the filter how measurements are correlated, for example when
     *   two sensors depend on the same noise or are coupled.
     *   - Higher off-diagonal elements \f$\Rightarrow\f$ higher coupling
     *   - Zero off-diagonal elements \f$\Rightarrow\f$ no coupling
     * 
     * - \f$\mathbf{Q}\f$ The process noise covariance matrix  
     *   Models imperfections in the system model and assumes additive gaussian noise on the state prediction step.
     *   It is used to:
     *   - model uncertainties in the prediction
     *   - account for model simplifications 
     *   - unmodeled dynamics
     *   
     *   diagonal entries reflect how uncertain you are about the state variable:
     *   - higher diagonal elements \f$\Rightarrow\f$ more uncertainty. Model can't reliably predict this state, let measurements correct it.
     *   - lower diagonal elements \f$\Rightarrow\f$ less uncertainty. The filter will trust the models prediction more.
     *  
     * - \f$\mathbf{P}\f$ The error covariance matrix
     *   Represents the uncertainty of the Kalman's filters estimate x.  
     *   Diagonal elements represent the variance of the state variable x:
     *   - Large diagonal elements \f$\Rightarrow\f$ the filter is uncertain about the corresponding state
     *   - Small diagonal elements \f$\Rightarrow\f$ the filter is confident about the value of the corresponding state
     *  
     *   Off-Diagonal elements estimate correlations between states:
     *   - Non-Zero \f$\Rightarrow\f$ those states are correlated
     *   - Zero \f$\Rightarrow\f$ those states are uncorrelated
     * 
     *   For Initialisation `P_0` should reflect how much you trust your initial guess of `x_0`:
     *   - diagonal elements smaller 1: you are quite confident about the state
     *   - diagonal elements larger 1: you are unsure about the initial guess and sensor readings should correct it quite quickly
     *   - much larger than 1 (e.g.: 1000): The initial guess is completely unreliable the filter should deduce it all from measurements.
     * 
     * - \f$\hat{\mathbf{P}}\f$ The predicted error covariance matrix
     * - \f$x\f$ The estimated state vector
     * - \f$\hat{x}\f$ The predicted estimated state vector
     * 
     */
    template<class T, int NStates, int NMeasurements=1, int NInputs = 0>
    class KalmanFilter{
        private:

        Eigen::Matrix<T, NStates, NStates> _A;
        Eigen::Matrix<T, NMeasurements, NStates> _H;
        Eigen::Matrix<T, NStates, NInputs> _B;
        Eigen::Matrix<T, NStates, NStates> _P;
        Eigen::Matrix<T, NMeasurements, NMeasurements> _R;
        Eigen::Matrix<T, NStates, NStates> _Q;

        Eigen::Vector<T, NStates> _x;

        public:

        /**
         * \brief Initialises a Kalman filter
         * 
         * \param A The system matrix, aka. state transition matrix
         * \param B Control input matrix
         * \param H Observation matrix
         * \param R Measurement noise covariance matrix
         * \param Q Process noise covariance matrix
         * \param P_0 Starting error covariance matrix
         * \param x_0 Starting state vector
         */
        template<
            int AOptions_ = 0, int AMaxRows_ = NStates, int AMaxCols_ = NStates,
            int HOptions_ = 0, int HMaxRows_ = NMeasurements, int HMaxCols_ = NStates,
            int BOptions_ = 0, int BMaxRows_ = NStates, int BMaxCols_ = NInputs,
            int ROptions_ = 0, int RMaxRows_ = NMeasurements, int RMaxCols_ = NMeasurements,
            int POptions_ = 0, int PMaxRows_ = NStates, int PMaxCols_ = NStates,
            int QOptions_ = 0, int QMaxRows_ = NStates, int QMaxCols_ = NStates
        >
        KalmanFilter(
            const Eigen::Matrix<T, NStates, NStates, AOptions_, AMaxRows_, AMaxCols_>& A, 
            const Eigen::Matrix<T, NMeasurements, NStates, HOptions_, HMaxRows_, HMaxCols_>& H,
            const Eigen::Matrix<T, NStates, NInputs, BOptions_, BMaxRows_, BMaxCols_>& B = Eigen::Matrix<T, NStates, NInputs>::Zero(),
            const Eigen::Matrix<T, NMeasurements, NMeasurements, ROptions_, RMaxRows_, RMaxCols_>& R = Eigen::Matrix<T, NMeasurements, NMeasurements>::Identity(),
            const Eigen::Matrix<T, NStates, NStates, POptions_, PMaxRows_, PMaxCols_>& P_0 = Eigen::Matrix<T, NStates, NStates>::Identity()*static_cast<T>(1e+3),
            const Eigen::Matrix<T, NStates, NStates, QOptions_, QMaxRows_, QMaxCols_>& Q = Eigen::Matrix<T, NStates, NStates>::Identity()*static_cast<T>(1e-2),
            const Eigen::Vector<T, NStates>& x_0 = Eigen::Vector<T, NStates>::Zero()
            )
            : _A(A)
            , _H(H)
            , _B(B)
            , _P(P_0)
            , _R(R)
            , _Q(Q)
            , _x(x_0){}

        /// @brief Set the system matrix
        /// @param new_A the new system matrix
        /// @return a reference to self for method chaining
        template<int AOptions_ = 0, int AMaxRows_ = NStates, int AMaxCols_ = NStates>
        KalmanFilter& setA(const Eigen::Matrix<T, NStates, NStates, AOptions_, AMaxRows_, AMaxCols_>& new_A){this->_A = new_A; return *this;}

        /// @brief Set the system matrix
        /// @param new_A the new system matrix
        /// @return a reference to self for method chaining
        template<int AOptions_ = 0, int AMaxRows_ = NStates, int AMaxCols_ = NStates>
        KalmanFilter& setSystemMatrix(const Eigen::Matrix<T, NStates, NStates, AOptions_, AMaxRows_, AMaxCols_>& new_A){this->_A = new_A; return *this;}

        /// @brief Sets the observation matrix
        /// @param new_H the new observation matrix
        /// @return a reference to self for method chaining
        template<int HOptions_ = 0, int HMaxRows_ = NMeasurements, int HMaxCols_ = NStates>
        KalmanFilter& setH(const Eigen::Matrix<T, NMeasurements, NStates, HOptions_, HMaxRows_, HMaxCols_>& new_H){this->_H = new_H; return *this;}
        
        /// @brief Sets the observation matrix
        /// @param new_H the new observation matrix
        /// @return a reference to self for method chaining
        template<int HOptions_ = 0, int HMaxRows_ = NMeasurements, int HMaxCols_ = NStates>
        KalmanFilter& setObservationMatrix(const Eigen::Matrix<T, NMeasurements, NStates, HOptions_, HMaxRows_, HMaxCols_>& new_H){this->_H = new_H; return *this;}

        /// @brief Sets the control input matrix
        /// @param new_B the new control input matrix
        /// @return a reference to self for method chaining
        template<int BOptions_ = 0, int BMaxRows_ = NStates, int BMaxCols_ = NInputs>
        KalmanFilter& setB(const Eigen::Matrix<T, NStates, NInputs, BOptions_, BMaxRows_, BMaxCols_>& new_B){this->_B = new_B; return *this;}
        
        /// @brief Sets the control input matrix
        /// @param new_B the new control input matrix
        /// @return a reference to self for method chaining
        template<int BOptions_ = 0, int BMaxRows_ = NStates, int BMaxCols_ = NInputs>
        KalmanFilter& setControlInputMatrix(const Eigen::Matrix<T, NStates, NInputs, BOptions_, BMaxRows_, BMaxCols_>& new_B){this->_B = new_B; return *this;}

        /// @brief Sets the measurement covariance matrix
        /// @param new_R The new measurement covariance matrix
        /// @return a reference to self for method chaining
        template<int ROptions_ = 0, int RMaxRows_ = NMeasurements, int RMaxCols_ = NMeasurements>
        KalmanFilter& setR(const Eigen::Matrix<T, NMeasurements, NMeasurements, ROptions_, RMaxRows_, RMaxCols_>& new_R){this->_R = new_R; return *this;}
        
        /// @brief Sets the measurement covariance matrix
        /// @param new_R The new measurement covariance matrix
        /// @return a reference to self for method chaining
        template<int ROptions_ = 0, int RMaxRows_ = NMeasurements, int RMaxCols_ = NMeasurements>
        KalmanFilter& setMeasurementCovMatrix(const Eigen::Matrix<T, NMeasurements, NMeasurements, ROptions_, RMaxRows_, RMaxCols_>& new_R){this->_R = new_R; return *this;}

        /// @brief Sets the error covariance matrix
        /// @param new_P The new error covariance matrix
        /// @return a reference to self for method chaining
        template<int POptions_ = 0, int PMaxRows_ = NStates, int PMaxCols_ = NStates>
        KalmanFilter& setP(const Eigen::Matrix<T, NStates, NStates, POptions_, PMaxRows_, PMaxCols_>& new_P){this->_P = new_P; return *this;}
        
        /// @brief Sets the error covariance matrix
        /// @param new_P The new error covariance matrix
        /// @return a reference to self for method chaining
        template<int POptions_ = 0, int PMaxRows_ = NStates, int PMaxCols_ = NStates>
        KalmanFilter& setErrorCovMatrix(const Eigen::Matrix<T, NStates, NStates, POptions_, PMaxRows_, PMaxCols_>& new_P){this->_P = new_P; return *this;}
        
        /// @brief Sets the process noise covariance matrix 
        /// @param new_Q The new process noise covariance matrix
        /// @return a reference to self for method chaining
        template<int QOptions_ = 0, int QMaxRows_ = NStates, int QMaxCols_ = NStates>
        KalmanFilter& setQ(const Eigen::Matrix<T, NStates, NStates, QOptions_, QMaxRows_, QMaxCols_>& new_Q){this->_Q = new_Q; return *this;}

        /// @brief Sets the process noise covariance matrix 
        /// @param new_Q The new process noise covariance matrix
        /// @return a reference to self for method chaining
        template<int QOptions_ = 0, int QMaxRows_ = NStates, int QMaxCols_ = NStates>
        KalmanFilter& setProcessNoiseCovMatrix(const Eigen::Matrix<T, NStates, NStates, QOptions_, QMaxRows_, QMaxCols_>& new_Q){this->_Q = new_Q; return *this;}

        /// @brief Sets the state vector
        /// @param new_x The new state vector
        /// @return a reference to self for method chaining 
        template<int POptions_ = 0, int PMaxRows_ = NStates, int PMaxCols_ = NStates>
        KalmanFilter& setx(const Eigen::Vector<T, NStates>& new_x){this->_x = new_x; return *this;}

        /// @brief Sets the state vector
        /// @param new_x The new state vector
        /// @return a reference to self for method chaining 
        template<int POptions_ = 0, int PMaxRows_ = NStates, int PMaxCols_ = NStates>
        KalmanFilter& setStateVector(const Eigen::Vector<T, NStates>& new_x){this->_x = new_x; return *this;}

        /**
         * \brief predicts the next system state
         * 
         * Current best prediction of the next measurement based on the new control input
         * 
         * \param u Control input vector
         * \returns An estimate of the next system state
         */
        Eigen::Vector<T, NStates> makePrediction(const Eigen::Vector<T, NInputs>& u) const {
            const Eigen::Vector<T, NStates> x_p = _A * _x + _B * u;
            return x_p;
        }

        void add(const Eigen::Vector<T, NMeasurements>& z, const Eigen::Vector<T, NInputs>& u = Eigen::Vector<T, NInputs>::Zero()){
            // 1. Prediction
            const Eigen::Vector<T, NStates> x_p = makePrediction(u);
            const Eigen::Matrix<T, NStates, NStates> P_p = _A * _P * _A.transpose() + _Q;

            // 2. Update
            // Kalman factor: K = P_p * _H.transpose() * (H * P_p * H.transpose() + R).inverse();
            // U = P_p * _H.transpose()
            // S = H * P_p * H.transpose() + R
            // K = U * S.inverse();
            // K = (S.transpose().inverse() * U.transpose()).transpose()
            //
            // > note S is already symetric, so no need to transpose
            //
            // K = (S.inverse() * U.transpose()).transpose()
            const auto U = P_p * _H.transpose();
            const auto S = _H * P_p * _H.transpose() + _R;
            const Eigen::Matrix<T, NStates, NMeasurements> K = (S.ldlt().solve(U.transpose())).transpose();

            _x = x_p + K * (z - _H * x_p);

            const auto I = Eigen::Matrix<T, NStates, NStates>::Identity();
            _P = (I - K * _H) * P_p;
        }

        /// @overload
        template<std::same_as<T> U>
        requires (NMeasurements == 1 && NInputs != 1)
        void add(const U& z, const Eigen::Vector<U, NInputs>& u = Eigen::Vector<U, NInputs>::Zero()){
            const Eigen::Vector<double, 1> v_z(z);
            this->add(v_z, u);
        }

        /// @overload
        template<std::same_as<T> U>
        requires (NMeasurements != 1 && NInputs == 1)
        void add(const Eigen::Vector<U, NMeasurements>& z, const U& u = static_cast<T>(0)){
            const Eigen::Vector<double, 1> v_u(u);
            this->add(z, v_u);
        }

        /// @overload
        template<std::same_as<T> U>
        requires (NMeasurements == 1 && NInputs == 1)
        void add(const U& z, const U& u = static_cast<T>(0)){
            const Eigen::Vector<double, 1> v_z(z);
            const Eigen::Vector<double, 1> v_u(u);
            this->add(v_z, v_u);
        }

        /**
         * \brief returns the current best estimate of the system states
         * \returns A vector of the systmes states with the dimeansions: states
         */
        const Eigen::Vector<T, NStates>& estimate() const {
            return this->_x;
        }

        /**
         * \brief returns the error covariance matrix
         * \returns a matrix with the dimensions of: states x states
         */
        const Eigen::Matrix<T, NStates, NStates>& errorCovariance() const {
            return this->_P;
        }

        /// @brief Resets the state vector and error covariance matrix
        /// @param x_0 The new state vector after the reset
        /// @param P_0 The new error covariance matrix after the reset
        void reset(
            const Eigen::Vector<T, NStates>& x_0 = Eigen::Vector<T, NStates>::Zero(),
            const Eigen::Matrix<T, NStates, NStates>& P_0 = Eigen::Matrix<T, NStates, NStates>::Identity()*static_cast<T>(1e+3))
        {
            this->_P = P_0;
            this->_x = x_0;   
        }

    };

} // namespace controlpp
