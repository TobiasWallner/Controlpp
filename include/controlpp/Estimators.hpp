#pragma once

#include <cassert>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <tl/expected.hpp>

#include <controlpp/DiscreteTransferFunction.hpp>
#include <controlpp/algorithm.hpp>

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

    enum class dft_estimate_error{
        data_ranges_different_lenth,
        data_range_too_small
    };

    inline std::string_view to_string(dft_estimate_error err){
        switch(err){
            case dft_estimate_error::data_ranges_different_lenth : return "Data ranges have different length, but should have equal size";
            case dft_estimate_error::data_range_too_small : return "Data ranges are too small for the system order";
            default : return "dft_estimate_error::error";
        }
    }

    inline std::ostream& operator<< (std::ostream& stream, dft_estimate_error err){
        return stream << controlpp::to_string(err);
    }

    /**
     * @brief Estimates a discrete time transfer function of a specific order for the input (u) and output (y) data pairs
     * 
     * Finds the optimal parameters that minimize the cost function:
     * \f[
     * \left( y - U p \right)^2
     * \f]
     * 
     * where 
     * - y: is the output of the system
     * - p: are the parameters of the transfer function (this is what we optimise)
     * - U: A matrix where each row contains the current and past inputs and outputs.
     * 
     * @tparam T The data type
     * @tparam NumOrder The numerator order for the resulting transfer function
     * @tparam DenOrder The denominator order for the resulting transfer function
     * @tparam N The size of the data
     * @param y The output data of the system
     * @param u The input data of the system
     * @param regularization adds an extra (\f$\lambda\f$) therm \f$ \lambda p^2 \f$ to the cost function that penalises large parameters (p).
     * @param hint adds an extra penalty to the cost function if parameters deviate from the hint \f$ \lambda (p - p_\text{hint})^2 \f$. Note that regularisation needs to be non-zero for the hint to take effect.
     * @returns A transfer function that best fits: `y = Tf * u`
     */
    template<class T, int NumOrder, int DenOrder>
    requires((NumOrder != Eigen::Dynamic) && (DenOrder != Eigen::Dynamic))
    tl::expected<DiscreteTransferFunction<T, NumOrder, DenOrder>, dft_estimate_error> 
    dft_estimate(
        const Eigen::Vector<T, Eigen::Dynamic>& u,
        const Eigen::Vector<T, Eigen::Dynamic>& y,
        const T& regularization = T(0),
        const DiscreteTransferFunction<T, NumOrder, DenOrder>& hint = DiscreteTransferFunction<T, NumOrder, DenOrder>({T(0)}, {T(1)})
    ){
        if(y.size() != u.size()){
            return tl::unexpected(dft_estimate_error::data_ranges_different_lenth);
        }

        constexpr int P = NumOrder + 1 + DenOrder;
        Eigen::Vector<T, P> s;
        s.setZero();
        
        // input (u) partition
        auto s_u = s.head(NumOrder + 1);

        // output (y partition)
        auto s_y = s.tail(DenOrder);

        const std::size_t start = std::max(s_u.size()-1, s_y.size());

        if(static_cast<std::size_t>(y.size()) <= start || static_cast<std::size_t>(u.size()) <= start){
            return tl::unexpected(dft_estimate_error::data_range_too_small);
        }

        // fill s_u and s_y with start values
        for(std::size_t k = 0; k < static_cast<std::size_t>(s_u.size()); ++k){
            s_u(k) = u(start - k);
        }
        for(std::size_t k = 0; k < static_cast<std::size_t>(s_y.size()); ++k){
            s_y(k) = -y(start - k - 1);
        }

        // sums
        Eigen::Matrix<T, P, P> S;
        S.setZero();
        
        Eigen::Vector<T, P> r;
        r.setZero();

        // first iteration out of the loop
        S.noalias() += s * s.transpose();
        r.noalias() += y(start) * s;

        // make an memory efficient sum to solve the least squares equation
        for(std::size_t k = start+1; (k < static_cast<std::size_t>(y.size())) && (k < static_cast<std::size_t>(u.size())); ++k){
            // push next input, output
            controlpp::shift_up(s_u.data(), s_u.data() + s_u.size(), u(k));
            controlpp::shift_up(s_y.data(), s_y.data() + s_y.size(), -y(k-1));

            S.noalias() += s * s.transpose();
            r.noalias() += y(k) * s;
        }

        // S should be symetric --> enforce symetry
        S = (S + S.transpose()) * T(0.5);

        // add regularisation to S
        S.diagonal().array() += regularization;

        // add hint penalty to r
        const T hint_a0 = hint.den(0);

        const Eigen::Vector<T, NumOrder + 1> hint_scaled_num = (hint.num().vector() / hint_a0);
        const Eigen::Vector<T, DenOrder + 1> hint_scaled_den = (hint.den().vector() / hint_a0);
        const Eigen::Vector<T, DenOrder> hint_scaled_den_tail = hint_scaled_den.template tail<DenOrder>();

        const Eigen::Vector<T, P> hint_v = controlpp::join_to_vector(hint_scaled_num, hint_scaled_den_tail);
        r += regularization * hint_v;

        // solve r = S * param
        // since S is symetric positive definite use llt
        std::cout << "S:\n" << S << std::endl;
        std::cout << "r:\n" << r.transpose() << std::endl;
        const Eigen::Vector<T, P> params = S.ldlt().solve(r);

        // partition result into numerator and denominator
        const Eigen::Vector<T, NumOrder + 1> num = params.head(NumOrder + 1);
        Eigen::Vector<T, DenOrder + 1> den;
        den(0) = T(1);
        den.tail(DenOrder) = params.tail(DenOrder);

        // construct resulting transfer function
        DiscreteTransferFunction<T, NumOrder, DenOrder> result(num, den);
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
        T _memory = 0.98;                           ///< memory factor (in literature: the forgetting factor, but i think higher value corresponds to higher memory makes more sense)
        T _cov_regularisation = 1e-9;               ///< regularisation for numerical stability
        T _gain_clamp = 10;                         ///< clamping the maximal gain for correction therms
        

        public:

        /**
         * \brief Creates a recursive least square object with start parameters/covariance and a memory factor
         * 
         * Implicitly sets the covariance matrix to a diagonal matrix where each diagonal element has a
         * value of 1000. Large diagonal elements mean that the algorithm is very uncertain about the parameters.
         * 
         * Assertions/Assumptions:
         * -----------------------
         * - 0 < memory <= 1
         * - cov_regularisation >= 0
         * 
         * \param param_hint The start value of the parameter vector.
         * If there is no prior knowledge of the values, 0 is often a good choice.
         * 
         * \param cov_hint The start value of the covariance matrix. 
         * The covariance matrix is a measure of the uncertainty of the parameter vector.
         * As a starting point use the square of the standard deviation of the noise if known.
         * If there is no prior knowledge of the uncertainties setting it to a diagonal matrix with elements much greater than 1 is often a good choice
         * 
         * \param memory The value memory that determines how much the past determines the new estimate. (In literature often called the: forgetting factor)
         * It has to be within the open-closed range: \f$(0, 1]\f$. 
         * Remembers more of the past with higher `memory` and forgets more with lower `memory`
         * - `memory` = 1: no forgetting, converges to standard least squares
         * - `memory` < 1: forgetting older values with an exponential decay  
         * Often used values are between 0.8 and 0.98
         * 
         * \param cov_regularisation Adds a small positive number to the diagonal of the covariance matrix P at every iteration for numerical stability:
         * - it prevents the covariance matrix from collapsing to zero
         * - reduces the risk of the matrix becoming ill conditioned
         * - improves robustness of the gain computation when excitation is weak  
         * Typical values are between `1e-12` and `1e-6` depending on numerical precision and signal scaling
         * 
         * \param gain_clamp Limits/clamps the elements of the recursive gain vector `K = P s / (λ + sᵀ P s)` to the range `[-gain_clamp, +gain_clamp]`. 
         * Prevents sudden large scale parameter updates caused by:
         * - poor excitation
         * - nearly singulare covariance matrices
         * - Very small gain denominators  
         * Important notes: clamping the gain breaks strict optimality and introduces a bias.
         * Typical values are between 5 and 50.
         * 
         */
        inline ReccursiveLeastSquares(
                const Eigen::Vector<T, NParams>& param_hint = Eigen::Vector<T, NParams>().setOne(), 
                T memory = 0.99,
                T cov_regularisation = 1e-9,
                T gain_clamp = 10
            )
            : _param(param_hint)
            , _memory(memory)
            , _cov_regularisation(cov_regularisation)
            , _gain_clamp(gain_clamp)
            {
                assert(0 < memory);
                assert(memory <= 1);
                assert(cov_regularisation > 0);
                
                this->_cov = (Eigen::Matrix<T, NParams, NParams>::Identity() * T(1000));
                this->_K.setZero();
            }

        /**
         * \brief Adds a new input output pair that updates the estimate
         * \param y The new system measurements/outputs
         * \param s The known system inputs
         * \returns The new parameter state estimate
         */
        void input(const T& y, const Eigen::Matrix<T, NMeasurements, NParams>& s)  {
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
         * @brief Set the value for the covariant regularisation
         * 
         * Often used values: `1e-6` or `1e-9`
         * 
         * @param reg The regularisation coefficient that will be added to the diagonal of the covariance matrix
         */
        inline void set_cov_regularisation(const T& reg){
            this->_cov_regularisation = reg;
        }

        [[nodiscard]] inline T cov_regularisation() const {
            return this->_cov_regularisation;
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
        void input(const T& y, const Eigen::Vector<T, NParams>& s)  {
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
            DiscreteTransferFunction<T, NumOrder, DenOrder> hint,
            const Eigen::Vector<T, NumOrder+1> NumeratorUncertainty,
            const Eigen::Vector<T, DenOrder> DenominatorUncertainty,
            const T& memory = 0.995)
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

        DtfEstimator(
            DiscreteTransferFunction<T, NumOrder, DenOrder> hint,
            const T& uncertainty = 1000,
            const T& memory = 0.995)
            : DtfEstimator(
                hint, 
                Eigen::Vector<T, NumOrder+1>().setOnes()*T(uncertainty), 
                Eigen::Vector<T, DenOrder>().setOnes()*T(uncertainty),
                memory){}

        DtfEstimator() 
            : DtfEstimator(
                DiscreteTransferFunction<T, 0, 0>({T(1)}, {T(1)}), // hint
                1000, // uncertainty
                0.995 // memory
            ){}


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
        void input(const T& y, const T& u){
            // add to the recursive least squares solver
            Eigen::Vector<T, 1 + NumOrder + DenOrder> s;
            s(0) = u;
            s.segment(1, NumOrder) = this->uk;
            s.tail(DenOrder) = this->neg_yk;
            this->rls.input(y, s);
            
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

        inline void set_cov_regularisation(const T& reg){
            this->rls.cov_regularisation(reg);
        }

        [[nodiscard]] inline T cov_regularisation() const {
            return this->rls.cov_regularisation();
        }

    };

    
      
    
} // namespace controlpp
