#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

namespace controlpp
{
    
    /**
     * \brief Solves the overdefined system \f$y = X p\f$ for p 
     * 
     * 
     * \param X the systems matrix that describes how the parameters p can be transformed into the measured output y
     * \param y the actual measured system output
     * 
     * \returns the approximated optimal solution for the parameter vector p
     */
    template<class T, int XRows, int XCols, int XOpt, int XMaxRows, int XMaxCols>
    Eigen::Vector<T, XCols> least_squares(const Eigen::Matrix<T, XRows, XCols, XOpt, XMaxRows, XMaxCols>& X, const Eigen::Vector<T, XRows>& y){
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
        ReccursiveLeastSquares(
            const Eigen::Vector<T, NParams>& param_hint = Eigen::Vector<T, NParams>().setZero(), 
            const Eigen::Matrix<T, NParams, NParams>& cov_hint = (Eigen::Matrix<T, NParams, NParams>::Identity() * T(1000)),
            double memory = 0.95)
            : _cov(cov_hint)
            , _param(param_hint)
            , _memory(memory)
            {}

        /**
         * \brief Adds a new input output pair that updates the estimate
         * \param y The new system measurements/outputs
         * \param s The known system inputs
         * \returns The new parameter state estimate
         */
        void add(const T& y, const Eigen::Matrix<T, NMeasurements, NParams>& s){
            const Eigen::Matrix<T, NMeasurements, NMeasurements> I = Eigen::Matrix<T, NMeasurements, NMeasurements>::Identity();

            // Gain
            const auto A = (this->_cov * s.transpose());
            const auto B = (this->_memory * I + s * this->_cov * s.transpose());
            
            /*
                Note that B is positive definite:
                    - memory > 0 ... is a positive scalar
                    - I ... is positive definite
                    - s * Cov * s^t
            */

            // calculate: K = A * B^-1
            const Eigen::Vector<T, NParams> K = B.transpose().partialPivLu().solve(A.transpose()).transpose();

            // Update
            this->_param = this->_param + K * (y - s * this->_param);

            // Covariance
            this->_cov = (this->_cov - K * s * this->_cov) / this->_memory;
        }

        /**
         * \brief returns the current best estimate
         * \returns the parameter vector
         */
        const Eigen::Vector<T, NParams>& estimate() const {return this->_param;}

        /**
         * \brief returns the current covariance
         * 
         * acts as a measure of the uncertainty of the estimate (higher values signal higher uncertainty)
         * 
         * \returns the current covariance matrix
         */
        const Eigen::Matrix<T, NParams, NParams>& covariance() const {return this->_cov;}
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
         * Often used values are between 0.8 and 0.98
         */
        ReccursiveLeastSquares(
            const Eigen::Vector<T, NParams>& param_hint = Eigen::Vector<T, NParams>().setZero(), 
            const Eigen::Matrix<T, NParams, NParams>& cov_hint = (Eigen::Matrix<T, NParams, NParams>::Identity() * T(1000)),
            double memory = 0.95)
            : _cov(cov_hint)
            , _param(param_hint)
            , _memory(memory)
            {}

        /**
         * \brief Adds a new input output pair that updates the estimate
         * \param y The new system measurements/outputs
         * \param s The known system inputs
         * \returns The new parameter state estimate
         */
        void add(const T& y, const Eigen::Vector<T, NParams>& s){
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
        const Eigen::Vector<T, NParams>& estimate() const {return this->_param;}

        /**
         * \brief returns the current covariance
         * 
         * acts as a measure of the uncertainty of the estimate (higher values signal higher uncertainty)
         * 
         * \returns the current covariance matrix
         */
        const Eigen::Matrix<T, NParams, NParams>& covariance() const {return this->_cov;}
    };

    /**
     * \brief Estimates a discrete transfer function from online data points
     * 
     * \f[
     * ARX = \frac{B(z)}{A(z)} u_k = \frac{b_0 + b_1 z^{1} + b_2 z^{2} \cdots b_m z^{m} }{1 + a_1 z^{1} + a_2 z^{2} \cdots 1 z^{n}} u_k
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
    template<class T, size_t NumSize, size_t DenSize, class Measurements=1>
    class DTFEstimator{
        private:

        ReccursiveLeastSquares<T, NumSize + DenSize - 1> rls;

        Eigen::Vector<T, DenSize> uk = Eigen::Vector<T, DenSize>::Zero();
        Eigen::Vector<T, DenSize - 1> neg_yk = Eigen::Vector<T, DenSize - 1>::Zero();
        T _memory;

        public:

        /**
         * \brief Initialises a discrete transfer function estimator with optional hints and memory/decay parameters
         * 
         * Note that \f$a_0\f$ will be set to 1 with an uncertainty of 0. 
         * 
         * This also means, that the `DenominatorUncertainty` starts at \f$a_1\f$ wheras the `NumeratorUncertainty` starts at \f$b_0\f$
         */
        DTFEstimator(
            const DiscreteTransferFunction<T, NumSize, DenSize>& hint = DiscreteTransferFunction<T, NumSize, DenSize>(Eigen::Vector<T, NumSize>().setZero(), Eigen::Vector<T, DenSize>().setZero()),
            const Eigen::Vector<T, NumSize> NumeratorUncertainty = Eigen::Vector<T, NumSize>().setOnes()*T(1000),
            const Eigen::Vector<T, DenSize - 1> DenominatorUncertainty = Eigen::Vector<T, DenSize-1>().setOnes()*T(1000),
            const T& memory = 0.98)
            : rls(
                controlpp::join_to_vector<T, NumSize, DenSize-1>(hint.num().vector(), hint.den().vector().tail(DenSize - 1).eval()), 
                controlpp::join_to_diagonal<T, NumSize, DenSize-1>(NumeratorUncertainty, DenominatorUncertainty), 
                memory){}

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
        const void add(const T& y, const T& u){
            this->uk.head(this->uk.size()-1) = this->uk.tail(this->uk.size()-1);
            this->uk(this->uk.size()-1) = u;

            const auto s = join_to_vector<T, NumSize, DenSize-1>(this->uk.head(NumSize), (this->neg_yk).eval());
            this->rls.add(y, s);

            this->neg_yk.head(this->neg_yk.size()-1) = this->neg_yk.tail(this->neg_yk.size()-1);
            this->neg_yk(this->neg_yk.size()-1) = -y;
        }

        /**
         * \brief returns the current estimate
         * \returns a DiscreteTransferFunction that represents the current best estimate
         */
        DiscreteTransferFunction<T, NumSize, DenSize> estimate(){
            DiscreteTransferFunction<T, NumSize, DenSize> result;
            Eigen::Vector<T, NumSize + DenSize - 1> est = rls.estimate();
            result.num().vector() = est.head(NumSize);
            result.den().vector().head(DenSize-1) = est.tail(DenSize-1);
            result.den().vector()(DenSize-1) = T(1);
            return result;
        }

    };
    
} // namespace controlpp
