#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

#include <controlpp/DiscreteTransferFunction.hpp>

namespace controlpp
{
    
    /**
     * \brief Solves the overdefined system \f$y = X p\f$ for p 
     * 
     * There have to be more measurements than parameters.
     * Meaning `XRows >= XCols` has to be true.
     * 
     * \param X the systems matrix that describes how the parameters p can be transformed into the measured output y
     * \param y the actual measured system output
     * 
     * \returns the approximated optimal solution for the parameter vector p
     */
    template<class T, int XRows, int XCols, int XOpt, int XMaxRows, int XMaxCols>
    requires((XRows >= XCols) || (XRows == Eigen::Dynamic) || (XCols == Eigen::Dynamic))
    Eigen::Vector<T, XCols> least_squares(const Eigen::Matrix<T, XRows, XCols, XOpt, XMaxRows, XMaxCols>& X, const Eigen::Vector<T, XRows>& y) noexcept {
        Eigen::Vector<T, XCols> result = X.colPivHouseholderQr().solve(y).eval();
        return result;
    }

    /**
     * \brief Calculates the recursive least square for online parameter estimation
     * 
     * Uses the recursive least squares algorithm:
     * 
     * \f[
     * k_j = \frac{P_{j-1} s_j}{q + s_{j}^{T} P_{j-1} s_j}
     * \f]
     * 
     * \f[
     * P_j = \left( P_{j-1} - k_{j} s_{j}^{T} P_{j-1}\right) \frac{1}{q}
     * \f]
     * 
     * \f[
     * p_j = p_{j-1} + k_{j} \left( y_j - s_{j}^{T} p_{j-1}\right)
     * \f]
     * 
     * ---
     * 
     * Note that for a number of outputs greater than one (`NOutputs>1`)
     * the model uses the 'Multi-Output System with shared parameter vector'.
     * 
     * This allows to use multiple sensors observing the same
     * system states and parameters to increase the estimation result.
     * 
     */
    template<class T, size_t NParams, size_t NMeasurements = 1>
    class ReccursiveLeastSquares{
        private:

        Eigen::Matrix<T, NParams, NParams> _cov; ///< previous covariance
        Eigen::Vector<T, NParams> _param; ///< previous parameter estimate
        double _memory;

        public:

        /**
         * \brief Creates a recursive least square object with start parameters/covariance and a memory factor
         * 
         * \param param_hint The start value of the parameter vector.
         * If there is no prior knowledge of the values, 0 is often a good choice.
         * 
         * \param cov_hint The start value of the covariance matrix. 
         * The covariance matrix is a measure of the uncertainty of the parameter vector.
         * As a starting point use the square of the standard deviation of the noise if known.
         * If there is no prior knowledge of the uncertainties setting it to a diagonal matrix with elements much greater than 1 is often a good choice
         * 
         * \param memory The value memory that determines how much the past determines the new estimate.
         * It has to be within the open-closed range: \f$(0, 1]\f$. 
         * Remembers more of the past with higher `memory` and forgets more with lower `memory`
         * - `memory` = 1: no forgetting, converges to standard least squares
         * - `memory` < 1: forgetting older values with an exponential decay
         * Often used values are between 0.8 and 0.98
         */
        inline ReccursiveLeastSquares(
            const Eigen::Vector<T, NParams>& param_hint = Eigen::Vector<T, NParams>().setZero(), 
            const Eigen::Matrix<T, NParams, NParams>& cov_hint = (Eigen::Matrix<T, NParams, NParams>::Identity() * T(1000)),
            double memory = 0.95)
            : _cov(cov_hint)
            , _param(param_hint)
            , _memory(memory)
            {
                if(memory <= T(0) || memory > T(1)){
                    throw std::invalid_argument("Error: ReccursiveLeastSquares::ReccursiveLeastSquares(): memory has to be in the open-closed range of: (0, 1]");
                }
            }

        /**
         * \brief Adds a new input output pair that updates the estimate
         * \param y The new system measurements/outputs
         * \param s The known system inputs
         * \returns The new parameter state estimate
         */
        void add(const T& y, const Eigen::Matrix<T, NMeasurements, NParams>& s) noexcept {
            const Eigen::Matrix<T, NMeasurements, NMeasurements> I = Eigen::Matrix<T, NMeasurements, NMeasurements>::Identity();

            // Gain
            const auto A = (this->_cov * s.transpose());
            const auto B = (this->_memory * I + s * this->_cov * s.transpose());
            
            /*
                Note that B is positive definite:
                    - memory > 0 ... is a positive scalar
                    - I ... is positive definite
                    - s * Cov * s^t ... is a squared therm which is positive definite

                thus the 'ldlt' solver can be used which is more numerically stable an also faster
            */

            // calculate: K = A * B^-1
            const Eigen::Vector<T, NParams> K = B.transpose().ldlt().solve(A.transpose()).transpose();

            // Update
            this->_param = this->_param + K * (y - s * this->_param);

            // Covariance
            this->_cov = (this->_cov - K * s * this->_cov) / this->_memory;
        }

        /**
         * \brief returns the current best estimate
         * \returns the parameter vector
         */
        inline const Eigen::Vector<T, NParams>& estimate() const noexcept {return this->_param;}

        /**
         * \brief returns the current covariance
         * 
         * acts as a measure of the uncertainty of the estimate (higher values signal higher uncertainty)
         * 
         * \returns the current covariance matrix
         */
        inline const Eigen::Matrix<T, NParams, NParams>& covariance() const noexcept {return this->_cov;}
    };

    /**
     * \brief Calculates the recursive least square for online parameter estimation
     * 
     * Solves the following system for \f$\vec{p}\f$ online one interation after another
     * 
     * \f[
     * y_k = \vec{s}_k^T \vec{p}_k
     * \f]
     * 
     * With the 
     * - measurement \f$y\f$, 
     * - the data vector \f$\vec{s}\f$ and 
     * - the parameter vector \f$\vec{p}\f$
     * 
     * ----
     * 
     * Uses the recursive least squares algorithm:
     * 
     * \f[
     * k_j = \frac{P_{j-1} s_j}{q + s_{j}^{T} P_{j-1} s_j}
     * \f]
     * 
     * \f[
     * P_j = \left( P_{j-1} - k_{j} s_{j}^{T} P_{j-1}\right) \frac{1}{q}
     * \f]
     * 
     * \f[
     * p_j = p_{j-1} + k_{j} \left( y_j - s_{j}^{T} p_{j-1}\right)
     * \f]
     */
    template<class T, size_t NParams>
    class ReccursiveLeastSquares<T, NParams, 1>{
        private:

        Eigen::Matrix<T, NParams, NParams> _cov; ///< previous covariance
        Eigen::Vector<T, NParams> _param; ///< previous parameter estimate
        double _memory;

        public:

        /**
         * \brief Creates a recursive least square object with start parameters/covariance and a memory factor
         * 
         * \param param_hint The start value of the parameter vector.
         * If there is no prior knowledge of the values, 0 is often a good choice.
         * 
         * \param cov_hint The start value of the covariance matrix. 
         * The covariance matrix is a measure of the uncertainty of the parameter vector.
         * As a starting point use the square of the standard deviation of the noise if known.
         * If there is no prior knowledge of the uncertainties setting it to a diagonal matrix with elements much greater than 1 is often a good choice
         * 
         * \param memory The value memory that determines how much the past determines the new estimate.
         * It has to be within the open-closed range: \f$(0, 1]\f$. 
         * Remembers more of the past with higher `memory` and forgets more with lower `memory`
         * - `memory` = 1: no forgetting, converges to standard least squares
         * - `memory` < 1: forgetting older values with an exponential decay
         * Often used values are between `0.9 `and `0.995`.
         */
        inline ReccursiveLeastSquares(
            const Eigen::Vector<T, NParams>& param_hint = Eigen::Vector<T, NParams>().setZero(), 
            const Eigen::Matrix<T, NParams, NParams>& cov_hint = (Eigen::Matrix<T, NParams, NParams>::Identity() * T(1000)),
            double memory = 0.99)
            : _cov(cov_hint)
            , _param(param_hint)
            , _memory(memory)
            {
                if(memory <= T(0) || memory > T(1)){
                    throw std::invalid_argument("Error: ReccursiveLeastSquares::ReccursiveLeastSquares(): memory has to be in the open-closed range of: (0, 1]");
                }
            }

        /**
         * \brief Adds a new input output pair that updates the estimate
         * \param y The new system measurements/outputs
         * \param s The known system inputs
         * \returns The new parameter state estimate
         */
        void add(const T& y, const Eigen::Vector<T, NParams>& s) noexcept {
            // Gain
            const Eigen::Vector<T, NParams> K = (this->_cov * s) / (this->_memory + s.transpose() * this->_cov * s);

            // Update
            this->_param = this->_param + K * (y - s.transpose() * this->_param);

            // Covariance
            this->_cov = (this->_cov - K * s.transpose() * this->_cov) / this->_memory;
        }

        /**
         * \brief returns the current best estimate
         * \returns the parameter vector
         */
        inline const Eigen::Vector<T, NParams>& estimate() const noexcept {return this->_param;}

        /**
         * \brief returns the current covariance
         * 
         * acts as a measure of the uncertainty of the estimate (higher values signal higher uncertainty)
         * 
         * \returns the current covariance matrix
         */
        inline const Eigen::Matrix<T, NParams, NParams>& covariance() const noexcept {return this->_cov;}
    };

    /**
     * \brief Estimates a discrete transfer function from online data points
     * 
     * Estimates time discrete transfer functions represented with positive powers of \f$z\f$ in the form of:
     * 
     * \f[
     * ARX = \frac{B(z)}{A(z)} u_k = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2} \cdots b_m z^{-m} }{a_0 + a_1 z^{-1} + a_2 z^{-2} \cdots 1 z^{-n}} u_k
     * \f]
     *
     * > Note: That \f$a_n\f$ has been set to 1
     * 
     * Estimates/Identifies model parameters online using the recursive least squares algorithm
     * 
     *  ----
     * 
     * The least-squares identification of the parameter for the identification task of the ARX model 
     * is unbiased and consistent if the stochastic disturbance satisfies the Yule-Walker equations 
     * of an autoregressive signal process with zero-mean white noise corresponding to the transfer 
     * function being identified.
     */
    template<class T, size_t NumOrder, size_t DenOrder, size_t Measurements=1>
    class DtfEstimator{
        private:

        ReccursiveLeastSquares<T, NumOrder+1 + DenOrder> rls;
        
        Eigen::Vector<T, NumOrder> uk = Eigen::Vector<T, NumOrder>::Zero();
        Eigen::Vector<T, DenOrder> neg_yk = Eigen::Vector<T, DenOrder>::Zero();
        T _memory;

        public:

        /**
         * \brief Initialises a discrete transfer function estimator with optional hints and memory/decay parameters
         * 
         * Note that \f$a_0\f$ will be set to 1 with an uncertainty of 0. 
         * 
         * This also means, that the `DenominatorUncertainty` starts at \f$a_1\f$ wheras the `NumeratorUncertainty` starts at \f$b_0\f$
         */
        DtfEstimator(
            const DiscreteTransferFunction<T, NumOrder, DenOrder>& hint = DiscreteTransferFunction<T, NumOrder, DenOrder>(Eigen::Vector<T, NumOrder+1>().setZero(), Eigen::Vector<T, DenOrder+1>().setZero()),
            const Eigen::Vector<T, NumOrder+1> NumeratorUncertainty = Eigen::Vector<T, NumOrder+1>().setOnes()*T(1000),
            const Eigen::Vector<T, DenOrder> DenominatorUncertainty = Eigen::Vector<T, DenOrder>().setOnes()*T(1000),
            const T& memory = 0.99)
            : rls(
                controlpp::join_to_vector<T, NumOrder+1, DenOrder>(hint.num().vector(), hint.den().vector().tail(DenOrder).eval()), 
                controlpp::join_to_diagonal<T, NumOrder+1, DenOrder>(NumeratorUncertainty, DenominatorUncertainty), 
                memory)
        {
            static_assert(NumOrder <= DenOrder, "The Discrete Transfer Function has to be propper. Meaning `NumOrder <= DenOrder` has to be true.");
            if(memory <= T(0) || memory > T(1)){
                throw std::invalid_argument("Error: DtfEstimator::DtfEstimator(): memory has to be in the open-closed range of: (0, 1]");
            }
        }

        /**
         * \brief Adds an interation step to the estimate
         * 
         * Advances the estimation by another input output pair
         * 
         * \param y the systems outpout value
         * \param u the systems input value
         * 
         * \returns the 
         */
        void add(const T& y, const T& u){
            // add to the recursive least squares solver
            Eigen::Vector<T, 1 + NumOrder + DenOrder> s;
            s(0) = u;
            s.segment(1, NumOrder) = this->uk;
            s.tail(DenOrder) = this->neg_yk;
            this->rls.add(y, s);
            
            // update uk
            std::reverse_copy(this->uk.data(), this->uk.data()+this->uk.size(), this->uk.data()+1);
            this->uk(0) = u;
            
            // update yk
            std::reverse_copy(this->neg_yk.data(), this->neg_yk.data()+this->neg_yk.size(), this->neg_yk.data()+1);
            this->neg_yk(0) = -y;
        }

        /**
         * \brief returns the current estimate
         * \returns a DiscreteTransferFunction that represents the current best estimate
         */
        DiscreteTransferFunction<T, NumOrder, DenOrder> estimate(){
            DiscreteTransferFunction<T, NumOrder, DenOrder> result;
            Eigen::Vector<T, NumOrder+1 + DenOrder> est = rls.estimate();
            result.num().vector() = est.head(NumOrder+1);
            result.den().vector()(0) = T(1);
            result.den().vector().tail(DenOrder) = est.tail(DenOrder);
            return result;
        }

    };
    
} // namespace controlpp
