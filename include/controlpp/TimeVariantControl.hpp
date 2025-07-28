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
        constexpr void kp(const T& kp){this->kp_ = kp;}

        /**
         * \brief returns the current gain
         */
        [[nodiscard]] constexpr const T& kp() const {return this->kp_;}

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

        /// @brief Namespace that contains sample-time variant controllers  
        namespace timevar{
        /**
         * \brief An I (Integrator) controller, with varying smaple-time
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
            constexpr void ki(const T& ki) {return this->ki_ = ki;}

            /// @brief returns the integrator gain 
            [[nodiscard]] constexpr const T& ki() const {return this->ki_;}

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
         * \brief An I (integrating) controller with anti windup protection (rate and value limits), with varying smaple-time
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
            T vn_;
            T vp_;
            T min_;
            T max_;
            
            T u_k1_ = static_cast<T>(0);
            T y_k1_ = static_cast<T>(0);

            public:

            /**
             * \brief Constructs an I (Integrator) controller with anti-windup controller
             * \param ki The integrator constant. Assumes: \f$k_i > 0\f$
             * \param vn The maximal negaive velocity of the output. Assumes: \f$v_n < 0 < v_p\f$
             * \param vp The maximal positive velocity of the output. Assumes: \f$v_n < 0 < v_p\f$
             * \param max The maximal output value. Assumes: \f$min < max\f$
             * \param min The minimal output value. Assumes: \f$min < max\f$
             */
            constexpr IAntiWindup(const T& ki, const T& vn, const T& vp, const T& min, const T& max) 
                : ki_(ki)
                , vp_(vp)
                , vn_(vn)
                , min_(min)
                , max_(max)
                {}

            /// @brief Default copy constructor 
            constexpr IAntiWindup(const IAntiWindup&) = default;

            /// @brief Sets the integrator gain
            /// @param ki The integrator gain for the next sample;
            constexpr void ki(const T& ki) {return this->ki_ = ki;}

            /// @brief returns the integrator gain 
            [[nodiscard]] constexpr const T& ki() const {return this->ki_;}

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
            constexpr void prate_limit(const T& pv) {return this->vp_ = pv;}

            /// @brief returns the positive rate limit
            [[nodiscard]] constexpr const T& prate_limit() const {return this->vp_;}

            /// @brief Sets the negative rate limit of the I controller
            /// @param nv The new negative rate limit
            constexpr void nrate_limit(const T& nv) {return this->vn_ = nv;}

            /// @brief returns the negative rate limit
            [[nodiscard]] constexpr const T& nrate_limit() const {return this->vn_;}

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
                const T v = std::clamp(this->ki_ / 2 * (u + this->u_k1_), this->vn_, this->vp_); 
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
         * \brief A D (differntial) controller, with varying smaple-time
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

            T u_k1_ = static_cast<T>(0); ///< previous control input
            T y_k1_ = static_cast<T>(0); ///< previous control output

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
            constexpr void kd(const T& k_d){this->k_d_ = k_d;}

            /// @brief Returns the differential gain
            /// @return The differentialgain
            [[nodiscard]] constexpr const T& kd(){return this->k_d_;}

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
                const T a = this->omega_ * this->k_d_;
                const T b = this->omega_ * Ts / 2 - 1;
                const T c = this->omega_ * Ts / 2 + 1;
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
         * \brief A PT1 (proportional time first order) filter (=low-pass element of order 1), with varying smaple-time
         * 
         * 
         * With the continuous transfer function:
         * 
         * \f[
         * \text{PT}_{1} = \frac{K}{1 + \frac{s}{\omega}}
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
        class PT1{
            private:
            T K_;
            T omega_;

            T u_k1_ = static_cast<T>(0);
            T y_k1_ = static_cast<T>(0);
            public:
            
            /**
             * \brief Constructs a low-pass filter
             * \param K The gain of the filter
             * \param omega The -3dB bandwidth of the lowpass filter
             */
            constexpr PT1(const T& K, const T& omega)
                : K_(K)
                , omega_(omega){}

            /// @brief Default copy constructor 
            constexpr PT1(const PT1&) = default;

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
                const T b = this->omega_ * Ts - 2;
                const T c = this->omega_ * Ts + 2;
                const T y = (a * (u + this->u_k1_) - b * this->y_k1_) / c;

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
        using LowPassO1 = PT1<T>;

        /**
         * \brief A PT2 (proportional time second order) filter (=low-pass element of order 1), with varying smaple-time
         * 
         * With the transfer function:
         * 
         * \f[
         * \text{PT}_{2} = \frac{K}{\frac{s^2}{\omega^2} + \frac{2 D s}{\omega} + 1}
         * \f]
         * 
         * and the discretisation using the tustin transformation:
         * 
         * \f[
         * s = \frac{2}{T_s} \frac{1 - z^{-1}}{1 + z^{-1}}
         * \f]
         * 
         * resulting in the controlers time series:
         * 
         * \f[
         * y_k = \frac{K \left( u_k + 2 u_{k-1} + u_{k-2} \right) - b y_{k-1} - c y_{k-2}}{a} \\
         * a = \frac{4}{\omega^2 + T_s^2} + 4 D + 1 \\
         * b = 2 - \frac{8}{\omega^2 T_s^2} \\
         * c = \frac{4}{\omega^2 + T_s^2} - 4 D + 1
         * \f]
         * 
         * with:
         * 
         * - \f$K\f$ the gain of the filter
         * - \f$T_s\f$ the sample-time in seconds 
         * - \f$\omega\f$ the characteristic frequency of the filter (bandwidth)
         * - \f$D\f$ the dampening factor
         *      - Aperiodic borderline case: D = 1 
         *      - Overdamped(creep-in) : D > 1 
         *      - Underdamped(exponential decreasing oscillations): 0 < D < 1 
         *      - Oscillating: D = 0
         * - \f$u_{k}\f$, \f$u_{k-1}\f$ and \f$u_{k-2}\f$ are current and previous control inputs
         * - \f$y_{k}\f$, \f$y_{k-1}\f$ and \f$y_{k-2}\f$ are current and previous control outputs
         * 
         * \tparam T Value type of the filter like `float` or `double`
         */
        template<class T>
        class PT2{
            private:
            T K_;
            T D_;
            T omega_;

            T u_k1_ = static_cast<T>(0);
            T u_k2_ = static_cast<T>(0);
            T y_k1_ = static_cast<T>(0);
            T y_k2_ = static_cast<T>(0);

            public:
            PT2(const PT2&) = default;

            /**
             * \brief constructs a PT2 element
             * \param K The gain of the filter
             * \param D Then dampening factor of the filter
             * \param omega The characteristic frequency of the filter (bandwidth)
             */
            PT2(const T& K, const T& D, const T& omega)
                : K_(K)
                , D_(D)
                , omega_(omega)
                {}

            /**
             * \brief Sets the gain of the filter
             * \param K The new gain for the filter
             */ 
            constexpr void gain(const T& K){this->K_ = K;}

            /**
             * \brief Returns the gain factor
             * \returns The gain of the filter
             */
            [[nodiscard]] constexpr const T& gain(){return this->K_;}

            /**
             * \brief Sets the dampening factor
             * 
             *  - Aperiodic borderline case: D = 1 
             *  - Overdamped(creep-in) : D > 1 
             *  - Underdamped(exponential decreasing oscillations): 0 < D < 1 
             *  - Oscillating: D = 0
             * 
             * \param D The new dampening factor
             */ 
            constexpr void dampening(const T& D){this->D_ = D;}

            /**
             * \brief Returns the dampening factor
             * 
             *  - Aperiodic borderline case: D = 1 
             *  - Overdamped(creep-in) : D > 1 
             *  - Underdamped(exponential decreasing oscillations): 0 < D < 1 
             *  - Oscillating: D = 0
             * 
             * \returns The dampening factor of the filter
             */
            [[nodiscard]] constexpr const T& dampening(){return this->D_;}

            /// @brief Sets the low-pass bandwidth
            /// @param omega The new low-pass bandwidth
            constexpr void omega(const T& omega){this->omega_ = omega;}

            /// @brief Returns the differential gain
            /// @return The differentialgain
            [[nodiscard]] constexpr const T& omega(){return this->omega_;}

            /// @brief Sets the previous \f$[k-1]\f$ control input
            /// @param u_k1 The previous \f$[k-1]\f$ control input
            constexpr void u_k1(const T& u_k1) {return this->u_k1_ = u_k1;}

            /// @brief Returns the previous \f$[k-1]\f$ control input
            [[nodiscard]] constexpr const T& u_k1() const {return this->u_k1_;}

            /// @brief Sets the previous \f$[k-2]\f$ control input
            /// @param u_k2 The previous \f$[k-2]\f$ control input
            constexpr void u_k2(const T& u_k2) {return this->u_k2_ = u_k2;}

            /// @brief Returns the previous \f$[k-2]\f$ control input
            [[nodiscard]] constexpr const T& u_k2() const {return this->u_k2_;}

            /// @brief Sets the previous \f$[k-1]\f$ control output
            /// @param y_k1 The previous \f$[k-1]\f$ control output
            constexpr void y_k1(const T& y_k1) {return this->y_k1_ = y_k1;}

            /// @brief returns the previous \f$[k-1]\f$ control output
            [[nodiscard]] constexpr const T& y_k1() const {return this->y_k1_;}

            /// @brief Sets the previous \f$[k-2]\f$ control output
            /// @param y_k2 The previous \f$[k-2]\f$ control output
            constexpr void y_k2(const T& y_k2) {return this->y_k2_ = y_k2;}

            /// @brief returns the previous \f$[k-2]\f$ control output
            [[nodiscard]] constexpr const T& y_k2() const {return this->y_k2_;}

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see operator()(const T& u, const T& Ts)
             */
            constexpr T input(const T& u, const T& Ts){
                // calculate output
                const T omega2 = this->omega_ * this->omega_;
                const T Ts2 = Ts * Ts;

                const T d = 4 / (omega2 * Ts2) + 1;
                const T e = 4 * this->D_;

                const T a = d + e;
                const T b = 2 * d;
                const T c = d - e;

                const T y = (this->K_ (u + 2 * this->u_k1_ + this->u_k2_) - b * this->y_k1_ - c * this->y_k2_) / a;

                // update states
                this->y_k2_ = this->y_k1_;
                this->y_k1_ = y;
                this->u_k2_ = this->u_k1_;
                this->u_k1_ = u;

                return y;
            }

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see input(const T& u, const T& Ts)
             */
            constexpr T operator()(const T& u, const T& Ts){return this->intput(u, Ts);}

            /**
             * \brief resets (clears) the internal states
             * resets the internal states (\f$u_{k-1}\f$, \f$u_{k-2}\f$, \f$y_{k-1}\f$ and \f$y_{k-2}\f$) to zero
             */
            constexpr void reset(){
                this->u_k1_ = static_cast<T>(0);
                this->u_k2_ = static_cast<T>(0);
                this->y_k1_ = static_cast<T>(0);
                this->y_k2_ = static_cast<T>(0);
            }
        };
        template<class T>
        using LowPassO2 = PT2<T>;

        /**
         * \brief A PI (Proportional Integral) controller with varying sample-times
         * 
         * With the transfer function:
         * 
         * \f[
         * \text{PI} = k_p + \frac{k_i}{s}
         * \f]
         * 
         * and the discretisation using the tustin transformation:
         * 
         * \f[
         * s = \frac{2}{T_s} \frac{1 - z^{-1}}{1 + z^{-1}}
         * \f]
         * 
         * resulting in the controlers time series:
         * 
         * \f[
         * y_k = \left( \frac{k_i T}{2} + k_p \right) u_k + \left( \frac{k_i T}{2} - k_p \right) u_{k-1} + y_{k-1}
         * \f]
         * 
         * with:
         * 
         * - \f$T_s\f$ the sample-time in seconds 
         * - \f$u_{k}\f$ and \f$u_{k-1}\f$ are current and previous control inputs
         * - \f$y_{k}\f$ and \f$y_{k-1}\f$ are current and previous control outputs
         * 
         * ----
         * 
         * The purely proportional component can amplify high-frequency noise.
         * While this does not improve control performance, it may lead to increased power consumption due to rapid steering actions.
         * To mitigate this, the proportional gain (P) can be replaced with a first-order low-pass filter (PT1), 
         * effectively transforming the PI controller into a PT1I controller.
         * 
         * \see controlpp::vartime::PT1I
         * 
         * \tparam T Value type of the filter like `float` or `double`
         */
        template<class T>
        class PI{
            private:

            T kp_;
            T ki_;

            T u_k1_ = static_cast<T>(0);
            T y_k1_ = static_cast<T>(0);

            public:

            /**
             * \param ki The integral gain
             * \param kp The proportional gain
             */
            PI(const T& kp, const T& ki)
                : kp_(kp)
                , ki_(ki)
                {}

            PI(const PI&) = default;

            /**
             * \brief sets the proportional gain
             * \param kp the new proportional gain for the next sample value
             */
            constexpr void kp(const T& kp){this->kp_ = kp;}

            /**
             * \brief returns the current gain
             */
            [[nodiscard]] constexpr const T& kp() const {return this->kp_;}

            /**
             * \brief sets the integral gain
             * \param ki the new proportional gain for the next sample value
             */
            constexpr void ki(const T& ki){this->ki_ = ki;}

            /**
             * \brief returns the current gain
             */
            [[nodiscard]] constexpr const T& ki() const {return this->ki_;}

            /// @brief Sets the previous \f$[k-1]\f$ control input
            /// @param u_k1 The previous \f$[k-1]\f$ control input
            constexpr void u_k1(const T& u_k1) {return this->u_k1_ = u_k1;}

            /// @brief Returns the previous \f$[k-1]\f$ control input
            [[nodiscard]] constexpr const T& u_k1() const {return this->u_k1_;}

            /// @brief Sets the previous \f$[k-2]\f$ control input
            /// @param u_k2 The previous \f$[k-2]\f$ control input
            constexpr void u_k2(const T& u_k2) {return this->u_k2_ = u_k2;}

            /// @brief Returns the previous \f$[k-2]\f$ control input
            [[nodiscard]] constexpr const T& u_k2() const {return this->u_k2_;}

            /// @brief Sets the previous \f$[k-1]\f$ control output
            /// @param y_k1 The previous \f$[k-1]\f$ control output
            constexpr void y_k1(const T& y_k1) {return this->y_k1_ = y_k1;}

            /// @brief returns the previous \f$[k-1]\f$ control output
            [[nodiscard]] constexpr const T& y_k1() const {return this->y_k1_;}

            /// @brief Sets the previous \f$[k-2]\f$ control output
            /// @param y_k2 The previous \f$[k-2]\f$ control output
            constexpr void y_k2(const T& y_k2) {return this->y_k2_ = y_k2;}

            /// @brief returns the previous \f$[k-2]\f$ control output
            [[nodiscard]] constexpr const T& y_k2() const {return this->y_k2_;}

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see operator()(const T& u, const T& Ts)
             */
            constexpr T input(const T& u, const T& Ts){
                // calculate new output
                const T p = ki_ * Ts / 2;

                const T a = p + kp_;
                const T b = p - kp_;

                const T y = a * u + b * u_k1_ + y_k1_;

                // update states
                this->y_k1_ = y;
                this->u_k1_ = u;

                return y;
            }

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see input(const T& u, const T& Ts)
             */
            constexpr T operator()(const T& u, const T& Ts){return this->intput(u, Ts);}

            /**
             * \brief resets (clears) the internal states
             * resets the internal states (\f$u_{k-1}\f$ and \f$y_{k-1}\f$) to zero
             */
            constexpr void reset(){
                this->u_k1_ = static_cast<T>(0);
                this->y_k1_ = static_cast<T>(0);
            }
        };

        /**
         * \brief A PI (Proportional Integral) controller with varying sample-times
         * 
         * Anti-windup prevents the controller’s integrator from accumulating excessive error when 
         * the actuator cannot keep up with the commanded output. Without anti-windup protection, 
         * the integrator continues to grow—even though the actuator is saturated and unable to reflect 
         * these changes—resulting in increasingly unrealistic control signals. Once the system reaches 
         * its target, the controller must then "unwind" or de-integrate the surplus error, which leads 
         * to a delayed response and potential overshoot. Anti-windup mitigates this by limiting or 
         * adjusting the integrator during saturation, improving stability and responsiveness.
         * 
         * Has the transfer function:
         * 
         * \f[
         * \text{PI} = k_p + \frac{k_i}{s}
         * \f]
         * 
         * and the discretisation using the tustin transformation:
         * 
         * \f[
         * s = \frac{2}{T_s} \frac{1 - z^{-1}}{1 + z^{-1}}
         * \f]
         * 
         * resulting in the controlers time series:
         * 
         * \f[
         * y_k = 
         *  \begin{cases}
         *      \if y_k > max \text{ then } max \\
         *      \if y_k < min \text{ then } max \\
         *      \if \frac{y_k - y_{k-1}}{Ts} > v_p \text{ then } v_p * Ts \\
         *      \if \frac{y_k - y_{k-1}}{Ts} < v_n \text{ then } v_n * Ts \\
         *      \text{else } \left( \frac{k_i T}{2} + k_p \right) u_k + \left( \frac{k_i T}{2} - k_p \right) u_{k-1} + y_{k-1}
         *  \end{cases}
         * \f]
         * 
         * 
         * 
         * with:
         * 
         * - \f$T_s\f$ the sample-time in seconds 
         * - \f$u_{k}\f$ and \f$u_{k-1}\f$ are current and previous control inputs
         * - \f$y_{k}\f$ and \f$y_{k-1}\f$ are current and previous control outputs
         * 
         * \see controlpp::vartime::PI
         * 
         * \tparam T Value type of the filter like `float` or `double`
         */
        template<class T>
        class PIAntiWindup{
            private:
            T kp_;
            T ki_;

            T vn_;
            T vp_;

            T min_;
            T max_;

            T u_k1_ = static_cast<T>(0);
            T y_k1_ = static_cast<T>(0);

            public:

            /// @brief Constructs a PI (proportional + integral) controller with anti-windup
            /// @param kp The proportional gain
            /// @param ki The integral gain.
            /// @param vn The minimal negative velocity of the control output. Assumes: \f$v_n < 0 < v_p\f$
            /// @param vp The maximal positive velocity of the control output. Assumes: \f$v_n < 0 < v_p\f$
            /// @param min The minimal value of the control output. Assumes: \f$min < max\f$
            /// @param max The maximal value of the control output. Assumes: \f$min < max\f$
            PIAntiWindup(const T& kp, const T& ki, const T& vn, const T& vp, const T& min, const T& max)
                : kp_(kp)
                , ki_(ki)
                , vn_(vn)
                , vp_(vp)
                , min_(min)
                , max_(max)
                {}

            PIAntiWindup(const PIAntiWindup&) = default;

            /**
             * \brief sets the proportional gain
             * \param kp the new proportional gain for the next sample value
             */
            constexpr void kp(const T& kp){this->kp_ = kp;}

            /**
             * \brief returns the current gain
             */
            [[nodiscard]] constexpr const T& kp() const {return this->kp_;}

            /**
             * \brief sets the integral gain
             * \param ki the new proportional gain for the next sample value
             */
            constexpr void ki(const T& ki){this->ki_ = ki;}

            /**
             * \brief returns the current gain
             */
            [[nodiscard]] constexpr const T& ki() const {return this->ki_;}

            /// @brief Sets the previous \f$[k-1]\f$ control input
            /// @param u_k1 The previous \f$[k-1]\f$ control input
            constexpr void u_k1(const T& u_k1) {return this->u_k1_ = u_k1;}

            /// @brief Returns the previous \f$[k-1]\f$ control input
            [[nodiscard]] constexpr const T& u_k1() const {return this->u_k1_;}

            /// @brief Sets the previous \f$[k-2]\f$ control input
            /// @param u_k2 The previous \f$[k-2]\f$ control input
            constexpr void u_k2(const T& u_k2) {return this->u_k2_ = u_k2;}

            /// @brief Returns the previous \f$[k-2]\f$ control input
            [[nodiscard]] constexpr const T& u_k2() const {return this->u_k2_;}

            /// @brief Sets the previous \f$[k-1]\f$ control output
            /// @param y_k1 The previous \f$[k-1]\f$ control output
            constexpr void y_k1(const T& y_k1) {return this->y_k1_ = y_k1;}

            /// @brief returns the previous \f$[k-1]\f$ control output
            [[nodiscard]] constexpr const T& y_k1() const {return this->y_k1_;}

            /// @brief Sets the previous \f$[k-2]\f$ control output
            /// @param y_k2 The previous \f$[k-2]\f$ control output
            constexpr void y_k2(const T& y_k2) {return this->y_k2_ = y_k2;}

            /// @brief returns the previous \f$[k-2]\f$ control output
            [[nodiscard]] constexpr const T& y_k2() const {return this->y_k2_;}

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see operator()(const T& u, const T& Ts)
             */
            constexpr T input(const T& u, const T& Ts){
                // calculate new output
                const T dy = (this->ki_ * Ts / 2 + this->kp_) * u + (this->ki_ * Ts / 2 - this->kp_) * this->u_k1_;
                const T v = dy / Ts;

                const T v_rate_clamped = std::clamp(v, this->vn_, this->vp_);
                const T dy_rate_clamped = v_rate_clamped * Ts;

                const T y = dy_rate_clamped + y_k1_;
                const T y_clamped = std::clamp(y, this->min_, this->max_);

                // update states
                this->y_k1_ = y_clamped;
                this->u_k1_ = u;

                return y_clamped;
            }

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see input(const T& u, const T& Ts)
             */
            constexpr T operator()(const T& u, const T& Ts){return this->intput(u, Ts);}

            /**
             * \brief resets (clears) the internal states
             * resets the internal states (\f$u_{k-1}\f$, \f$y_{k-1}\f$) to zero
             */
            constexpr void reset(){
                this->u_k1_ = static_cast<T>(0);
                this->y_k1_ = static_cast<T>(0);
            }
        };


        /**
         * \brief A Thamed PI controller (PI controller with a low pass filter) with variable sample-times
         * 
         * The purely proportional component of a PI controller can amplify high-frequency noise.
         * While this does not improve control performance, it may lead to increased power consumption due to rapid steering actions.
         * To mitigate this, the proportional gain (P) can be replaced with a first-order low-pass filter (PT1), 
         * effectively transforming the PI controller into a PT1I controller.
         * 
         * With the transfer function:
         * 
         * \f[
         * \text{PI} = \frac{s k_p + k_i}{s + \frac{s^2}{\omega}}
         * \f]
         * 
         * and the discretisation using the tustin transformation:
         * 
         * \f[
         * s = \frac{2}{T_s} \frac{1 - z^{-1}}{1 + z^{-1}}
         * \f]
         * 
         * resulting in the controlers time series:
         * 
         * \f[
         * y_k = (K * (b_0 * u_k + b_1 * u_{k-1} + b_2 * u_{k-2}) - a_1 * y_{k-1} - a_2 * y_{k-2}) / a_0
         * \f]
         * 
         * where:
         * 
         * \f[
         * K = \omega * T_s
         * a0 = 4 + 2 * \omega * T_s
         * a1 = -8
         * a2 = 4 - 2 * \omega * T_s
         * b0 = 2 * k_p + k_i * T_s
         * b1 = 2 * k_i * T_s
         * b2 = k_i * T_s - 2 * k_p
         * \f]
         * 
         * with:
         * 
         * - \f$T_s\f$ the sample-time in seconds 
         * - \f$u_{k}\f$ and \f$u_{k-1}\f$ are current and previous control inputs
         * - \f$y_{k}\f$ and \f$y_{k-1}\f$ are current and previous control outputs
         * - \f$k_p\f$ the proportional gain
         * - \f$k_i\f$ the integral gain
         * - \f$\omega\f$ the characteristic frequency of the low pass filter aimed to thame the P control at high frequencies
         * 
         * 
         * \see controlpp::vartime::PT1I
         * 
         * \tparam T Value type of the filter like `float` or `double`
         */
        template<class T>
        class ThamedPI{
            private:

            T kp_;
            T ki_;
            T omega_;

            T u_k1_ = static_cast<T>(0);
            T u_k2_ = static_cast<T>(0);
            T y_k1_ = static_cast<T>(0);
            T y_k2_ = static_cast<T>(0);

            public:

            constexpr ThamedPI(const ThamedPI&) = default;

            /**
             * \brief Constructs a PT1I (proportional time delayed + integral) controller
             * \param kp The proportional gain of the controller
             * \param ki The integral gain of the controller
             * \param omega The lowpass filter bandwidth
             */
            constexpr ThamedPI(const T& kp, const T& ki, const T& omega)
                : kp_(kp)
                , ki_(ki)
                , omega_(omega){}

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see operator()(const T& u, const T& Ts)
             */
            constexpr T input(const T& u, const T& Ts){
                // calculate filter parameters
                const T K = omega_ * Ts;

                const T a0 = 4 + 2 * omega_ * Ts;
                const T a1 = -8;
                const T a2 = 4 - 2 * omega_ * Ts;

                const T b0 = 2 * kp_ + ki_ * Ts;
                const T b1 = 2 * ki_ * Ts;
                const T b2 = ki_ * Ts - 2 * kp_;

                // calculate filter output
                const T y = (K * (b0 * u + b1 * u_k1_ + b2 * u_k2_) - a1 * y_k1_ - a2 * y_k2_) / a0;

                // update states
                u_k2_ = u_k1_;
                u_k1_ = u;
                y_k2_ = y_k1_;
                y_k1_ = y;

                return y;
            }

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see input(const T& u, const T& Ts)
             */
            constexpr T operator()(const T& u, const T& Ts){return this->intput(u, Ts);}

            /**
             * \brief resets (clears) the internal states
             * resets the internal states (\f$u_{k-1}\f$, \f$u_{k-2}\f$, \f$y_{k-1}\f$ and \f$y_{k-2}\f$) to zero
             */
            constexpr void reset(){
                this->u_k1_ = static_cast<T>(0);
                this->u_k2_ = static_cast<T>(0);
                this->y_k1_ = static_cast<T>(0);
                this->y_k2_ = static_cast<T>(0);
            }
        };

        class PT1IAntiWindup{
            //TODO
        };

        /**
         * \brief A PID (proportional + integral + differential) controller with variable sample-times
         * 
         * with the transfer function:
         * 
         * \f[
         * PID(s) = \left( k_p + \frac{k_i}{s} + k_d s \right)\frac{1}{1 + \frac{s}{\omega}}
         * \f]
         * 
         * The therm:
         * 
         * \f[
         * \frac{1}{1 + \frac{s}{\omega}}
         * \f]
         * 
         * has been added to:
         * 
         * 1. make the D controller realisable
         * 2. thame the D controller for high frequencies
         * 
         * 
         */
        template<class T>
        class PID{
            private:

            T kp_;
            T ki_;
            T kd_;
            T omega_;

            T y_k1_ = static_cast<T>(0);
            T y_k2_ = static_cast<T>(0);
            T u_k1_ = static_cast<T>(0);
            T u_k2_ = static_cast<T>(0);

            public:

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see operator()(const T& u, const T& Ts)
             */
            constexpr T input(const T& u, const T& Ts){
                // parameters
                const T a0 = 4 + 2 * omega_ * Ts;
                const T a1 = -8;
                const T a2 = 4 - 2 * omega_ * Ts;

                const T b0 = 4 * kd_ * omega_ + 2 * kp_ + ki_ * Ts;
                const T b1 = 2 * ki_ * Ts - 8 * kd_ * omega_;
                const T b2 = 4 * kd_ * omega_ - 2 * kp_ + ki_ * Ts;

                // control output
                const T y = (b0 * u + b1 * u_k1_ + b2 * u_k2_ - a1 * y_k1_ - a2 * y_k2_) / a0;

                // update states
                y_k2_ = y_k1_;
                y_k1_ = y;

                u_k2_ = u_k1_;
                u_k1_ = u;

                return y;
            }

            /**
             * \brief resets (clears) the internal states
             * resets the internal states (\f$u_{k-1}\f$, \f$u_{k-2}\f$, \f$y_{k-1}\f$ and \f$y_{k-2}\f$) to zero
             */
            constexpr T operator() (const T& u, const T& Ts){return this->input(u, Ts);}

            /**
             * \brief resets (clears) the internal states
             * resets the internal states (\f$u_{k-1}\f$, \f$u_{k-2}\f$, \f$y_{k-1}\f$ and \f$y_{k-2}\f$) to zero
             */
            constexpr void reset(){
                this->u_k1_ = static_cast<T>(0);
                this->u_k2_ = static_cast<T>(0);
                this->y_k1_ = static_cast<T>(0);
                this->y_k2_ = static_cast<T>(0);
            }

        };

        class PIDAntiWindup{
            //TODO
        };

        /**
         * \brief A Lead-Lag control element with varying sample-times
         * 
         * Has a continuous transfer function of:
         * 
         * \f[
         * LeadLag(s) = \frac{1 + \frac{s}{\omega_1}}{1 + \frac{s}{\omega_2}}
         * \f]
         * 
         * transformed with the tustin transformation to the following time series:
         * 
         * \f[
         * y_k = \frac{ \left( T_s \omega_2 + \frac{2 \omega_2}{\omega_1} \right) u_k + \left( T_s \omega_2 - \frac{2 \omega_2}{\omega_1} \right) u_{k-1} - \left( T_s \omega_2 - 2 \right)}{ T_s \omega_2 + 2 }
         * \f]
         * 
         * with:
         * 
         * - \f$T_s\f$ The sample-time in seconds
         * - \f$\omega_1\f$ The characteristic frequency of the phase lead
         * - \f$omega_2\f$ The characteristic frequency of the phase lag
         * - \f$u_{k}\f$ and \f$u_{k-1}\f$ are current and previous control inputs
         * - \f$y_{k}\f$ and \f$y_{k-1}\f$ are current and previous control outputs
         * 
         */
        template<class T>
        class LeadLag{
            private:

            T omega1_;
            T omega2_;

            T u_k1_ = static_cast<T>(0);
            T y_k1_ = static_cast<T>(0);

            public:

            LeadLag(const T& omega1, const T& omega2)
                : omega1_(omega1)
                , omega2_(omega2)
                {}

            LeadLag(const LeadLag&) = default;

            /// @brief Sets the low-pass bandwidth
            /// @param omega1 The new low-pass bandwidth
            constexpr void omega1(const T& omega1){this->omega1_ = omega1;}

            /// @brief Returns the differential gain
            /// @return The differentialgain
            [[nodiscard]] constexpr const T& omega1(){return this->omega1_;}

            /// @brief Sets the low-pass bandwidth
            /// @param omega2 The new low-pass bandwidth
            constexpr void omega2(const T& omega2){this->omega2_ = omega2;}

            /// @brief Returns the differential gain
            /// @return The differentialgain
            [[nodiscard]] constexpr const T& omega2(){return this->omega2_;}

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

            constexpr T input(const T& u, const T& Ts){
                // calculate filter output
                const T a = this->omega2_ * (Ts + 2 / this->omega1_);
                const T b = this->omega2_ * (Ts - 2 / this->omega1_);
                const T c = Ts * this->omega2_ - 2;
                const T d = Ts * this->omega2_ + 2;
                const T y = (a * u + b * this->u_k1_ - c * this->y_k1_) / d;

                // update state
                this->y_k1_ = y;
                this->u_k1_ = u;
                
                return y;
            }

            constexpr void reset(){
                this->u_k1_ = static_cast<T>(0);
                this->y_k1_ = static_cast<T>(0);
            }
        };

    } // vartime
} // namespace controlpp
