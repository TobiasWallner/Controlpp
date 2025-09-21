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
    Eigen::Vector<T, XCols> least_squares(const Eigen::Matrix<T, XRows, XCols, XOpt, XMaxRows, XMaxCols>& X, const Eigen::Vector<T, XRows>& y)  {
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
     * For a number of outputs greater than one (`NOutputs>1`)
     * the model uses the 'Multi-Output System with shared parameter vector'.
     * This allows to use multiple sensors observing the same
     * system states and parameters to increase the estimation result.
     * 
     * There is a default regularisation therm (default: 1e-9) added to the diagonal elements of the covariance updata
     * 
     * There is gain clamping (default: [-10, +10]) applied to the parameter update therem K.
     */
    template<class T, size_t NParams, size_t NMeasurements = 1>
    class ReccursiveLeastSquares{
        private:

        Eigen::Matrix<T, NParams, NParams> _cov;    ///< previous covariance
        Eigen::Vector<T, NParams> _param;           ///< previous parameter estimate
        Eigen::Vector<T, NParams> _K;
        T _memory = 0.98;
        T _cov_regularisation = 1e-9;
        T _gain_clamp = 10;
        

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
         * 
         */
        inline ReccursiveLeastSquares(
            const Eigen::Vector<T, NParams>& param_hint = Eigen::Vector<T, NParams>().setOne(), 
            T memory = 0.99)
            : _param(param_hint)
            , _memory(memory)
            {
                if(memory <= T(0) || memory > T(1)){
                    throw std::invalid_argument("Error: ReccursiveLeastSquares::ReccursiveLeastSquares(): memory has to be in the open-closed range of: (0, 1]");
                }
                this->_cov = (Eigen::Matrix<T, NParams, NParams>::Identity() * T(1000));
                this->_K.setZero();
            }

        /**
         * \brief Adds a new input output pair that updates the estimate
         * \param y The new system measurements/outputs
         * \param s The known system inputs
         * \returns The new parameter state estimate
         */
        void add(const T& y, const Eigen::Matrix<T, NMeasurements, NParams>& s)  {
            this->_cov.diagonal().array() += this->_cov_regularisation;

            // Gain
            const auto A = (this->_cov * s.transpose()).eval();
            const Eigen::Matrix<T, NMeasurements, NMeasurements> I_m = Eigen::Matrix<T, NMeasurements, NMeasurements>::Identity();
            auto B1 = (s * this->_cov * s.transpose()).eval();
            B1.diagonal().array() += this->_memory;
            const auto B = ((B1 + B1.transpose())/2).eval(); // force symetry

            // calculate: K = A * B^-1
            this->_K = B.transpose().llt().solve(A.transpose()).transpose().eval();

            // Limit the update gain
            for(int i = 0; i < this->_K.size(); ++i){
                this->_K.at(i) = std::clamp(this->_K.at(i), -this->_gain_clamp, this->_gain_clamp);
            }

            // Update
            this->_param += this->_K * (y - s * this->_param);

            // Covariance
            this->_cov -= this->_K * s * this->_cov;
            this->_cov /= this->_memory;
        }

        /**
         * \brief returns the current best estimate
         * \returns the parameter vector
         */
        inline const Eigen::Vector<T, NParams>& estimate() const  {return this->_param;}

        /**
         * \brief returns the current covariance
         * 
         * acts as a measure of the uncertainty of the estimate (higher values signal higher uncertainty)
         * 
         * \returns the current covariance matrix
         */
        [[nodiscard]] inline const Eigen::Matrix<T, NParams, NParams>& cov() const  {return this->_cov;}

        inline void set_cov(const Eigen::Matrix<T, NParams, NParams>& cov) {
            this->_cov = cov;
        }

        /**
         * @brief Sets the memory factor
         * 
         * The memory factor [0, 1], where:
         * - larger values: Old parameters and inputs are remembered for longer (slower changes, more robust to noise).
         * - smaller values: Old parameters are forgotten more quickly (faster changes)
         * 
         * @param memory The new memory factor
         */
        inline void set_memory(const T& memory){
            this->_memory = memory;
        }

        /**
         * @brief Returns the memory factor
         * 
         * The memory factor [0, 1], where:
         * - larger values: Old parameters and inputs are remembered for longer (slower changes, more robust to noise).
         * - smaller values: Old parameters are forgotten more quickly (faster changes)
         * 
         * @returns The current memory factor
         */
        [[nodiscard]] inline const T& memory() const {
            return this->_memory;
        }

        /**
         * \brief Limits the update gain K. 
         * 
         * Limits the update gain K from -gain_clamp to +gain_clamp. 
         * Prevents too fast updates and ill conditioned updates. For example from a lack of excitation variety
         * 
         * default is 10
         * 
         * \param gain_clamp The new gain clamp
         */
        inline void set_gain_clamp(const T& gain_clamp) {
            this->_gain_clamp = gain_clamp;
        }

        /**
         * @brief Returns the currect gain clamp factor
         * @return The currect gain clamp factor
         */
        [[nodiscard]] inline const T& gain_clamp(){
            return this->_gain_clamp;
        }



        /**
         * @brief Returns the gain K used in the parameter update
         * 
         * The gain K can be seen as a measure of uncertainty
         * 
         * @return The gain vector;
         */
        inline const Eigen::Vector<T, NParams>& gain() const  {return this->_K;}
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
        Eigen::Vector<T, NParams> _K;
        T _memory = 0.98;
        T _cov_regularisation = 1e-9;
        T _gain_clamp = 10;
        

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
         * 
         * \param cov_regularisation A value that will be added to the diagonal of the covariance matrix before each update
         * to prevent the covariance to be become too small, ill formed and unregular. This is mainly to increase numerical stability. 
         * 
         */
        inline ReccursiveLeastSquares(
            const Eigen::Vector<T, NParams>& param_hint = Eigen::Vector<T, NParams>().setOnes(), 
            const Eigen::Matrix<T, NParams, NParams>& cov_hint = (Eigen::Matrix<T, NParams, NParams>::Identity() * T(1000)),
            T memory = 0.99,
            T cov_regularisation = 1e-9
        )
            : _cov(cov_hint)
            , _param(param_hint)
            , _memory(memory)
            , _cov_regularisation(cov_regularisation)
        {
            if(memory <= T(0) || memory > T(1)){
                throw std::invalid_argument("Error: ReccursiveLeastSquares::ReccursiveLeastSquares(): memory has to be in the open-closed range of: (0, 1]");
            }
            this->_K.setZero();
        }

        inline void set_cov(const Eigen::Matrix<T, NParams, NParams>& cov){
            this->_cov = cov;
        }

        inline void set_param(const Eigen::Vector<T, NParams>& param){
            this->_param = param;
        }

        inline void set_memory(const T& memory){
            this->_memory = memory;
        }

        inline void set_gain_clamp(const T& gain_clamp){
            this->_gain_clamp = gain_clamp;
        }

        [[nodiscard]] inline const T& gain_clamp() const {
            return this->_gain_clamp;
        }

        /**
         * \brief Adds a new input output pair that updates the estimate
         * \param y The new system measurements/outputs
         * \param s The known system inputs
         * \returns The new parameter state estimate
         */
        void add(const T& y, const Eigen::Vector<T, NParams>& s)  {
            this->_cov.diagonal().array() += this->_cov_regularisation;

            // Gain
            const auto A = (this->_cov * s).eval();
            const T B = this->_memory + s.transpose() * this->_cov * s;
            this->_K = A / B;

            for(int i = 0; i < this->_K.size(); ++i){
                this->_K(i) = std::clamp(this->_K(i), -this->_gain_clamp, this->_gain_clamp);
            }

            // Update
            this->_param += this->_K * (y - s.transpose() * this->_param);

            this->_cov -= this->_K * s.transpose() * this->_cov;
            this->_cov /= this->_memory;
        }

        /**
         * \brief returns the current best estimate
         * \returns the parameter vector
         */
        inline const Eigen::Vector<T, NParams>& estimate() const  {return this->_param;}

        /**
         * \brief returns the current covariance
         * 
         * acts as a measure of the uncertainty of the estimate (higher values signal higher uncertainty)
         * 
         * \returns the current covariance matrix
         */
        inline const Eigen::Matrix<T, NParams, NParams>& cov() const  {return this->_cov;}

        /**
         * @brief Returns the gain used for the updata
         * 
         * The gain can be seen as a measurement of uncertainty
         * 
         * @return The gain used in the parameter update
         */
        inline const Eigen::Vector<T, NParams>& gain() const  {return this->_K;}
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
        
        /*
            TODO

            Possible cause for numerical instability:

            For early iterations, when uk and neg_yk are mostly zero, the system matrix s may be poorly conditioned (near-zero values), leading to numerical instability in the RLS update, especially in the gain calculation.

        */
        Eigen::Vector<T, NumOrder> uk = Eigen::Vector<T, NumOrder>::Zero();
        Eigen::Vector<T, DenOrder> neg_yk = Eigen::Vector<T, DenOrder>::Zero();
        T _memory = 0.98;

        public:

        /**
         * \brief Initialises a discrete transfer function estimator with optional hints and memory/decay parameters
         * 
         * Note that the hint needs to be a propper transfer function with `a_0 != 0`
         * 
         * This also means, that the `DenominatorUncertainty` starts at \f$a_1\f$ wheras the `NumeratorUncertainty` starts at \f$b_0\f$
         * 
         */
        DtfEstimator(
            DiscreteTransferFunction<T, NumOrder, DenOrder> hint = DiscreteTransferFunction<T, 0, 0>({static_cast<T>(1)}, {static_cast<T>(1)}),
            const Eigen::Vector<T, NumOrder+1> NumeratorUncertainty = Eigen::Vector<T, NumOrder+1>().setOnes()*T(1000),
            const Eigen::Vector<T, DenOrder> DenominatorUncertainty = Eigen::Vector<T, DenOrder>().setOnes()*T(1000),
            const T& memory = 0.99)
        {
            T a_0 = hint.den().at(0);
            if(a_0 != static_cast<T>(0)){
                hint.num() /= a_0;
                hint.den() /= a_0;
            }

            auto param = controlpp::join_to_vector<T, NumOrder+1, DenOrder>(hint.num().vector(), hint.den().vector().tail(DenOrder).eval());
            this->rls.set_param(param);

            auto cov = controlpp::join_to_diagonal<T, NumOrder+1, DenOrder>(NumeratorUncertainty, DenominatorUncertainty);
            this->rls.set_cov(cov);

            this->rls.set_memory(memory);

            static_assert(NumOrder <= DenOrder, "The Discrete Transfer Function has to be propper. Meaning `NumOrder <= DenOrder` has to be true.");
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
            std::copy_backward(this->uk.data(), this->uk.data()+this->uk.size(), this->uk.data()+1);
            this->uk(0) = u;
            
            // update yk
            std::copy_backward(this->neg_yk.data(), this->neg_yk.data()+this->neg_yk.size(), this->neg_yk.data()+1);
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

        const Eigen::Matrix<T, NumOrder+1 + DenOrder, NumOrder+1 + DenOrder> cov(){
            return this->rls.cov();
        }

        const Eigen::Vector<T, NumOrder+1 + DenOrder>& gain(){
            return this->rls.gain();
        }

        void set_gain_clamp(const T& g){
            this->rls.set_gain_clamp(g);
        }

        [[nodiscard]] const T& gain_clamp() const {
            return this->rls.gain_clamp();
        }

    };
    
} // namespace controlpp
