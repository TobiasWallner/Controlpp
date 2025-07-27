/**
 * \file TimeVariantControl.hpp
 * 
 * \brief Contains controllers for systems with runtime varying sample-times like high time jitters
 * 
 * This file contains controllers that have been analytically/symbolically calculated to allow for varying runtime sample-times.
 *
 */

namespace controlpp
{
    /**
     * \brief A P (Proportional) controller
     * 
     * has the following transfer function
     * 
     * \f[
     * y_k = k_p u_k
     * \f]
     * 
     * where:
     * 
     * - \f$y_k\f$ is the output of the controller of the data sample k
     * - \f$u_k\f$ is the k-th input sample
     * - \f$k_p\f$ is the controller gain
     * 
     * \tparam T the data type the controller uses, like `float` or `double`
     */
    template<class T>
    class P{
        public:
        using value_type = T;
        private:

        T kp_;

        public:

        /// @brief Default copy constructor
        constexpr P(const P&) = default;

        /**
         * \brief Construct a P (Proportional) controller with a constant gain
         * \param kp the proportional gain of the controller
         */
        constexpr P(const T& kp) : kp_(kp){}

        /**
         * \brief sets the proportional gain
         * \param kp the new proportional gain for the next sample value
         */
        constexpr void gain(const T& kp){this->kp_ = kp;}

        /**
         * \brief returns the current gain
         */
        [[nodiscard]] constexpr const T& gain() const {return this->kp_;}

        /**
         * \brief Adds a new sample
         * \param u the new sample input
         * \returns the new output based on the input
         */
        constexpr T input(const T& u) const {return u * this->kp_;}

        /**
         * \brief Adds a new sample
         * \param u the new sample input
         * \returns the new output based on the input
         */
        constexpr T operator() (const T& u) const {return this->input(u);}
    };

    /**
     * \brief An I (Integrator) controller
     * 
     * (With variable sample-time)
     * 
     * An I controller with the following transfer function:
     * 
     * \f[
     * I(s) = \frac{k_i}{s}
     * \f]
     * 
     * discretised with the tustin transformation:
     * 
     * \f[
     * s = \frac{2}{T_s} \frac{1 - z^{-1}}{1 + z^{-1}}
     * \f]
     * 
     * to the time series:
     * 
     * \f[
     * y_k = \frac{k_i T_s}{2} \left( u_k + u_{k-1} \right) + y_{k-1}
     * \f]
     * 
     * where:
     * - \f$k_i\f$ is the integral gain
     * - \f$T_s\f$ is the sample time interval between uk−1uk−1​ and ukuk​ (in seconds)
     * 
     * ---
     * 
     * allows varying sample-times \f$T_s\f$
     * 
     * \tparam T the data type the controller uses, like `float` or `double`
     * \see controlpp::IAntiWindup
     */
    template<class T>
    class I{
        public:
        using value_type = T;
        private:
        T ki_;
        T u_k1_ = static_cast<T>(0);
        T y_k1_ = static_cast<T>(0);

        public:

        /**
         * \brief Constructs an I (Integrator) controller
         * \param ki The integrator constant.
         */
        constexpr I(const T& ki) : ki_(ki){}

        /// @brief Default copy constructor 
        constexpr I(const I&) = default;

        /// @brief Sets the integrator gain
        /// @param ki The integrator gain for the next sample
        constexpr void gain(const T& ki) {return this->ki_ = ki;}

        /// @brief returns the integrator gain 
        [[nodiscard]] constexpr const T& gain() const {return this->ki_;}

        /// @brief Sets the previous control input
        /// @param u_k1 The previous control input
        constexpr void u_k1(const T& u_k1) {return this->u_k1_ = u_k1;}

        /// @brief Returns the previous control input
        [[nodiscard]] constexpr const T& u_k1() const {return this->u_k1_;}

        /// @brief Sets the previous control output
        /// @details Useful to set the initial integration state
        /// @param y_k1 The previous control output
        constexpr void y_k1(const T& y_k1) {return this->y_k1_ = y_k1;}

        /// @brief returns the previous control output
        [[nodiscard]] constexpr const T& y_k1() const {return this->y_k1_;}

        /**
         * \brief resets (clears) the internal states
         * 
         * resets the internal states to zero
         */
        constexpr void reset(){
            this->u_k1_ = static_cast<T>(0);
            this->y_k1_ = static_cast<T>(0);
        }

        /**
         * \brief Adds a new value input and sample time to calculate the next output
         * \details Also advances the internal states
         * \returns The control output
         */
        constexpr T input(const T& u, const T& Ts){
            const T y = this->ki_ * Ts / 2 * (u + this->u_k1_) + this->y_k1_;
            this->y_k1_ = y;
            this->u_k1_ = u;
            return y;
        }

        /**
         * \brief Adds a new value input and sample time to calculate the next output
         * \details Also advances the internal states
         * \returns The control output
         * \see controlpp::I::input(const T& u, const T& Ts)
         */
        constexpr T operator() (const T& u, const T& Ts){return this->input(u, Ts);}
    };

    /**
     * \brief An I (integrating) controller with anti windup protection (rate and value limits) 
     * 
     * (With variable sample-time)
     * 
     * Anti-windup prevents the controller’s integrator from accumulating excessive error when 
     * the actuator cannot keep up with the commanded output. Without anti-windup protection, 
     * the integrator continues to grow—even though the actuator is saturated and unable to reflect 
     * these changes—resulting in increasingly unrealistic control signals. Once the system reaches 
     * its target, the controller must then "unwind" or de-integrate the surplus error, which leads 
     * to a delayed response and potential overshoot. Anti-windup mitigates this by limiting or 
     * adjusting the integrator during saturation, improving stability and responsiveness.
     * 
     * Its an I controller with a transfer function like:
     * 
     * \f[
     * I(s) = \frac{k_i}{s}
     * \f]
     * 
     * discretised with the tustin transformation:
     * 
     * \f[
     * s = \frac{2}{T_s} \frac{1 - z^{-1}}{1 + z^{-1}}
     * \f]
     * 
     * but with a conditional time-series:
     * 
     * \f[
     * y_k = 
     *  \begin{cases}
     *      \if y_k > max \text{ then } max \\
     *      \if y_k < min \text{ then } max \\
     *      \if \frac{y_k - y_{k-1}}{Ts} > v_p \text{ then } v_p * Ts \\
     *      \if \frac{y_k - y_{k-1}}{Ts} < v_n \text{ then } v_n * Ts \\
     *      \text{else } \frac{k_i T_s}{2} \left( u_k + u_{k-1} \right) + y_{k-1}
     *  \end{cases}
     * \f]
     * 
     * where:
     * 
     * - \f$k_i\f$ is the integrator gain
     * - \f$T_s\f$ is the sample-time in seconds
     * - \f$u_k\f$ and \f$u_{k-1}\f$ are the current and previous control inputs
     * - \f$y_k\f$ and \f$y_{k-1}\f$ are the current and previous control outputs
     * - \f$max\f$ is the maximal output clamp at which positive integration stops
     * - \f$min\f$ is the minimal outpout clamp at which negative integration stops
     * - \f$v_p\f$ is the positive rate output speed clamp. (Limits integration speed)
     * - \f$v_n\f$ is the negative rate output speed clamp. (Limits integration speed)
     * 
     * \tparam T the data type the controller uses, like `float` or `double`
     * 
     * \see controlpp::I
     */
    template<class T>
    class IAntiWindup{
        public:
        using value_type = T;
        private:
        T ki_;
        T pv_;
        T nv_;
        T max_;
        T min_;
        
        T u_k1_ = static_cast<T>(0);
        T y_k1_ = static_cast<T>(0);

        public:

        /**
         * \brief Constructs an I (Integrator) controller with anti-windup controller
         * \param ki The integrator constant. Assumes: \f$k_i > 0\f$
         * \param pv The maximal positive velocity of the output. Assumes: \f$v_n < 0 < v_p\f$
         * \param nv The maximal negaive velocity of the output. Assumes: \f$v_n < 0 < v_p\f$
         * \param max The maximal output value. Assumes: \f$min < max\f$
         * \param min The minimal output value. Assumes: \f$min < max\f$
         */
        constexpr IAntiWindup(const T& ki, const T& pv, const T& nv, const T& max, const T& min) 
            : ki_(ki)
            , pv_(pv)
            , nv_(nv)
            , max_(max)
            , min_(min){}

        /// @brief Default copy constructor 
        constexpr IAntiWindup(const IAntiWindup&) = default;

        /// @brief Sets the integrator gain
        /// @param ki The integrator gain for the next sample;
        constexpr void gain(const T& ki) {return this->ki_ = ki;}

        /// @brief returns the integrator gain 
        [[nodiscard]] constexpr const T& gain() const {return this->ki_;}

        /// @brief Sets the previous control input
        /// @param u_k1 The previous control input
        constexpr void u_k1(const T& u_k1) {return this->u_k1_ = u_k1;}

        /// @brief Returns the previous control input
        [[nodiscard]] constexpr const T& u_k1() const {return this->u_k1_;}

        /// @brief Sets the previous control output
        /// @details Useful to set the initial integration state
        /// @param y_k1 The previous control output
        constexpr void y_k1(const T& y_k1) {return this->y_k1_ = y_k1;}

        /// @brief returns the previous control output
        [[nodiscard]] constexpr const T& y_k1() const {return this->y_k1_;}
        
        /// @brief Sets the positive rate limit of the I controller
        /// @param pv The new positive rate limit
        constexpr void prate_limit(const T& pv) {return this->pv_ = pv;}

        /// @brief returns the positive rate limit
        [[nodiscard]] constexpr const T& prate_limit() const {return this->pv_;}

        /// @brief Sets the negative rate limit of the I controller
        /// @param nv The new negative rate limit
        constexpr void nrate_limit(const T& nv) {return this->nv_ = nv;}

        /// @brief returns the negative rate limit
        [[nodiscard]] constexpr const T& nrate_limit() const {return this->nv_;}

        /// @brief Sets the maximal value limit
        /// @param max The new maximal value limit
        constexpr void max_limit(const T& max) {return this->max_ = max;}

        /// @brief returns the maximal value limit
        [[nodiscard]] constexpr const T& max_limit() const {return this->max_;}

        /// @brief Sets the minimal value limit
        /// @param min The new minimal value limit
        constexpr void min_limit(const T& min) {return this->min_ = min;}

        /// @brief returns the minimal value limit
        [[nodiscard]] constexpr const T& min_limit() const {return this->min_;}

        /**
         * \brief resets (clears) the internal states
         * 
         * resets the internal states to zero
         * 
         * \param u_k1 The previous input `u`
         * \param y_k1 The previous output `y`
         */
        constexpr void reset(){
            this->u_k1_ = static_cast<T>(0);
            this->y_k1_ = static_cast<T>(0);
        }

        /**
         * \brief Adds a new value input and sample time to calculate the next output
         * \details Also advances the internal states
         * \returns The control output
         */
        constexpr T input(const T& u, const T& Ts){
            const T v = std::clamp(this->ki_ / 2 * (u + this->u_k1_), this->nv_, this->pv_); 
            const T dy = v * Ts;
            const T y = std::clamp(dy + this->y_k1_, this->min_, this->max_);
            
            this->y_k1_ = y;
            this->u_k1_ = u;
            return y;
        }

        /**
         * \brief Adds a new value input and sample time to calculate the next output
         * \details Also advances the internal states
         * \returns The control output
         * \see controlpp::I::input(const T& u, const T& Ts)
         */
        constexpr T operator() (const T& u, const T& Ts){return this->input(u, Ts);}
    };

    /**
     * \brief A D (differntial) controller
     * 
     * (With variable sample-time)
     * 
     * An ideal D (differential) control element has the transfer function:
     * 
     * \f[
     * D(s) = k_d * s.
     * \f]
     * 
     * However, such a transfer function is:
     * 
     * 1. not realisable
     * 2. will amplify lots of noise
     * 
     * So a low pass filter is added:
     * 
     * \f[
     * D(s) = \frac{k_d s}{1 + \frac{s}{\omega}}
     * \f]
     * 
     * where:
     * 
     * * \f$k_d\f$ is the differential gain
     * * \f$omega\f$ is the low-pass bandwidth in radians per second
     * 
     * Discretised with the tustin transformation:
     * 
     * \f[
     * s = \frac{2}{T_s} \frac{1 - z^{-1}}{1 + z^{-1}}
     * \f]
     * 
     * with the sample time \f$T_s\f$ and the delay operator \f$z^{-1}\f$ to:
     * 
     * \f[
     * y_k = \frac{1}{\frac{\omega T_s}{2} + 1} \left[ \omega k_d \left( u_k - u_{k-1} \right) - \left( \frac{\omega T_s}{2} - 1 \right) y_{k-1} \right]
     * \f]
     * 
     * 
     * \tparam T The data type of the D controller (inputs/outputs/states). Usually `float` or `double`.
     */
    template<class T>
    class D{
        private:

        T k_d_; ///< differential gain
        T omega_; ///< low pass bandwidth

        T u_k1_; ///< previous control input
        T y_k1_; ///< previous control output

        public:

        /**
         * \brief Construct a D control element
         * \param k_d The differential gain
         * \param omega The bandwidth of the low-pass filter in radians per second
         */
        constexpr D(const T& k_d, const T& omega)
            : k_d_(k_d)
            , omega_(omega){}

        /// @brief default copy constructor 
        constexpr D(const D&) = default;

        /// @brief Sets the differential gain
        /// @param k_d The new differential gain
        constexpr void gain(const T& k_d){this->k_d_ = k_d;}

        /// @brief Returns the differential gain
        /// @return The differentialgain
        [[nodiscard]] constexpr const T& gain(){return this->k_d_;}

        /// @brief Sets the low-pass bandwidth
        /// @param omega The new low-pass bandwidth
        constexpr void omega(const T& omega){this->omega_ = omega;}

        /// @brief Returns the differential gain
        /// @return The differentialgain
        [[nodiscard]] constexpr const T& omega(){return this->omega_;}

        /// @brief Sets the previous control input
        /// @param u_k1 The previous control input
        constexpr void u_k1(const T& u_k1) {return this->u_k1_ = u_k1;}

        /// @brief Returns the previous control input
        [[nodiscard]] constexpr const T& u_k1() const {return this->u_k1_;}

        /// @brief Sets the previous control output
        /// @details Useful to set the initial integration state
        /// @param y_k1 The previous control output
        constexpr void y_k1(const T& y_k1) {return this->y_k1_ = y_k1;}

        /// @brief returns the previous control output
        [[nodiscard]] constexpr const T& y_k1() const {return this->y_k1_;}

        /**
         * \brief Input a new sample to the controller
         * \param u The control input
         * \param Ts The sample-time
         * \see D::operator()(const T& u, const T& Ts)
         */
        constexpr T input(const T& u, const T& Ts){
            // calculate new output
            const T w_Ts_half = this->omega_ * Ts / 2;
            const T a = this->omega_ * this->k_d_;
            const T b = (w_Ts_half-1);
            const T c = (w_Ts_half + 1)
            const T y = (a * (u - this->u_k1_) - b * this->y_k1_) / c;

            // update internal states
            this->y_k1_ = y;
            this->u_k1_ = u;

            return y;
        }

        /**
         * \brief Input a new sample to the controller
         * \param u The control input
         * \param Ts The sample-time
         * \see D::input(const T& u, const T& Ts)
         */
        constexpr T operator() (const T& u, const T& Ts){return this->input(u, Ts);}

        /**
         * \brief resets (clears) the internal states
         * 
         * resets the internal states (\f$u_{k-1}\f$ and \f$y_{k-1}\f$) to zero
         */
        constexpr void reset(){
            this->u_k1_ = static_cast<T>(0);
            this->y_k1_ = static_cast<T>(0);
        }
    };

    /**
     * \brief A low-pass element of order 1
     * 
     * (With varying smaple-time)
     * 
     * With the continuous transfer function:
     * 
     * \f[
     * LP1 = \frac{K}{1 + \frac{s}{\omega}}
     * \f]
     * 
     * Discretised with the tustin transformation:
     * 
     * \f[
     * s = \frac{2}{T_s} \frac{1 - z^{-1}}{1 + z^{-1}}
     * \f]
     * 
     * to the sample series:
     * 
     * \f[
     * y_k = \frac{K \omega T_s \left( u_k - u_{k-1} \right) - \left( 2 - \omega T_s \right) y_{k-1}}{2 + \omega T_s}
     * \f]
     * 
     * with:
     * 
     * - \f$K\f$ The gain of the low-pass
     * - \f$T_s\f$ The sample-time in seconds
     * - \f$\omega\f$ The -3dB bandwidth of the low-pass
     * - \f$u_k\f$ and \f$u_{k-1}\f$ the current and previous control input
     * - \f$y_k\f$ and \f$y_{k-1}\f$ the current and previous control output
     * 
     */
    template<class T>
    class LowPassO1{
        private:
        T K_;
        T omega_;

        T u_k1_;
        T y_k1_;
        public:
        
        /**
         * \brief Constructs a low-pass filter
         * \param K The gain of the filter
         * \param omega The -3dB bandwidth of the lowpass filter
         */
        constexpr LowPassO1(const T& K, const T& omega)
            : K_(K)
            , omega_(omega){}

        /// @brief Default copy constructor 
        constexpr LowPassO1(const LowPassO1&) = default;

        /// @brief Sets the gain of the filter
        /// @param K The new gain
        constexpr void gain(const T& K){this->K_ = K;}

        /// @brief Returns the gain (\f$K\f$) of the filter
        /// @return The gain (\f$K\f$) of the filter
        [[nodiscard]] constexpr const T& gain(){return this->K_;}

        /// @brief Sets the low-pass bandwidth
        /// @param omega The new low-pass bandwidth
        constexpr void omega(const T& omega){this->omega_ = omega;}

        /// @brief Returns the differential gain
        /// @return The differentialgain
        [[nodiscard]] constexpr const T& omega(){return this->omega_;}

        /// @brief Sets the previous control input
        /// @param u_k1 The previous control input
        constexpr void u_k1(const T& u_k1) {return this->u_k1_ = u_k1;}

        /// @brief Returns the previous control input
        [[nodiscard]] constexpr const T& u_k1() const {return this->u_k1_;}

        /// @brief Sets the previous control output
        /// @details Useful to set the initial integration state
        /// @param y_k1 The previous control output
        constexpr void y_k1(const T& y_k1) {return this->y_k1_ = y_k1;}

        /// @brief returns the previous control output
        [[nodiscard]] constexpr const T& y_k1() const {return this->y_k1_;}

        /**
         * \brief Input a new sample to the controller
         * \param u The control input
         * \param Ts The sample-time
         * \see operator()(const T& u, const T& Ts)
         */
        constexpr T input(const T& u, const T& Ts){
            // calculate control output:
            const T a = this->K_ * this->omega_ * Ts;
            const T b = 2 - this->omega_ * Ts;
            const T c = 2 + this->omega_ * Ts;
            const T y = (a * (u - this->u_k1_) - b * this->y_k1_) / c;

            // update states
            this->u_k1_ = u;
            this->y_k1_ = y;

            return y;
        }

        /**
         * \brief Input a new sample to the controller
         * \param u The control input
         * \param Ts The sample-time
         * \see input(const T& u, const T& Ts)
         */
        constexpr T operator() (const T& u, const T& Ts){return this->input(u, Ts);}

    };

    template<class T>
    class LowPassO2{
        private:

        T D_;
        T omega_

        public:

        constexpr T input(const T& u, const T& Ts){
            // calculate output
            const T omega2 = this->omega_ * this->omega_;
            const T Ts2 = Ts * Ts;

            const T d = 4 / (omega2 * Ts2) + 1;
            const T e = 4 * this->D_;

            const T a = d + e;
            const T b = 2 * d;
            const T c = d - e;

            const T y = (K (u + 2 * this->u_k1_ + this->u_k2_) - b * this->y_k1_ - c * this->y_k2_)/a;

            // update states
            this->y_k2_ = this->y_k1_;
            this->y_k1_ = y;
            this->u_k2_ = this->u_k1_;
            this->u_k1_ = u;

            return y;
        }
    };

    class PI{
        //TODO
    };

    class PIAntiWindup{
        //TODO
    };

    class PID{
        //TODO
    };

    class PIDAntiWindup{
        //TODO
    };

    class LeadLag{
        //TODO
    };

} // namespace controlpp
