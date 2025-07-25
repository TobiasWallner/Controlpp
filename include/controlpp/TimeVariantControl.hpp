/**
 * \file TimeVariantControl.hpp
 * 
 * \brief Contains controlers for systems with runtime variing sample-times like high time jitters
 * 
 * This file contains controllers that have been analytically/symbolically calculated to allow for variing runtime sample-times.
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
         * \brief Construct a P (Proportional) controler with a constant gain
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
        constexpr const T& kp() const {return this->kp_;}

        /**
         * \brief Adds a new sample
         * \param u the new sample input
         * \returns the new output based on the input
         */
        constexpr T add(const T& u) const {return u * this->kp_;}

        /**
         * \brief Adds a new sample
         * \param u the new sample input
         * \returns the new output based on the input
         */
        constexpr T operator() (const T& u) const {return this->add(u);}
    };

    /**
     * \brief An I (Integrator) controller with variing sample times
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
     * - \f$T_s\f$ is the sample time
     * 
     * ---
     * 
     * allows variing sample-times \f$T_s\f$
     * 
     * \tparam T the data type the controller uses, like `float` or `double`
     * \see controlpp::AntiWindupI
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
         * \param ki The integrator constant
         */
        constexpr I(const T& ki) : ki_(ki){}

        /// @brief Default copy constructor 
        constexpr I(const I&) = default;

        /// @brief Sets the integrator gain
        /// @param ki The integrator gain for the next sample;
        constexpr void ki(const T& ki) {return this->ki_ = ki;}

        /// @brief returns the integrator gain 
        constexpr const T& ki() const {return this->ki_;}

        /// @brief Sets the previous control input
        /// @param u_k1 The previous control input
        constexpr void u_k1(const T& u_k1) {return this->u_k1_ = u_k1;}

        /// @brief Returns the previous control input
        constexpr const T& u_k1() const {return this->u_k1_;}

        /// @brief Sets the previous control output
        /// @details Useful to set the initial integration state
        /// @param y_k1 The previous control output
        constexpr void y_k1(const T& y_k1) {return this->y_k1_ = y_k1;}

        /// @brief returns the previous control output
        constexpr const T& y_k1() const {return this->y_k1_;}

        /**
         * \brief resets (clears) the internal states
         * 
         * resets the internal states to zero
         * 
         * \param u_k1 The previous input `u`
         * \param y_k1 The previous output `y`
         */
        constexpr void reset(const T& u_k1, const T& y_k1){
            this->u_k1_ = static_cast<T>(0);
            this->y_k1_ = static_cast<T>(0);
        }

        /**
         * \brief Adds a new value input and sample time to calculate the next output
         * \details Also advances the internal states
         * \returns The control output
         */
        constexpr T add(const T& u, const T& Ts){
            const T y = this->ki_ * Ts / 2 * (u + this->u_k1_) + this->y_k1_;
            this->y_k1_ = y;
            this->u_k1_ = u;
            return y;
        }

        /**
         * \brief Adds a new value input and sample time to calculate the next output
         * \details Also advances the internal states
         * \returns The control output
         * \see controlpp::I::add(const T& u, const T& Ts)
         */
        constexpr T operator() (const T& u, const T& Ts){return this->add(u, Ts);}
    };

    /**
     * \brief An I (integrating) controller with anti windup protection (rate and value limits) 
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
     * \see controlpp::I
     */
    class AntiWindupI{
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
         * \param ki The integrator constant
         * \param pv The maximal positive velocity of the output
         * \param nv The maximal negaive velocity of the output
         * \param max The maximal output value
         * \param min The minimal output value
         */
        constexpr AntiWindupI(const T& ki, const T& pv, const T& nv, const T& max, const T& min) 
            : ki_(ki)
            , pv_(pv)
            , nv_(nv)
            , max_(max)
            , min_(min){}

        /// @brief Default copy constructor 
        constexpr AntiWindupI(const AntiWindupI&) = default;

        /// @brief Sets the integrator gain
        /// @param ki The integrator gain for the next sample;
        constexpr void ki(const T& ki) {return this->ki_ = ki;}

        /// @brief returns the integrator gain 
        constexpr const T& ki() const {return this->ki_;}

        /// @brief Sets the previous control input
        /// @param u_k1 The previous control input
        constexpr void u_k1(const T& u_k1) {return this->u_k1_ = u_k1;}

        /// @brief Returns the previous control input
        constexpr const T& u_k1() const {return this->u_k1_;}

        /// @brief Sets the previous control output
        /// @details Useful to set the initial integration state
        /// @param y_k1 The previous control output
        constexpr void y_k1(const T& y_k1) {return this->y_k1_ = y_k1;}

        /// @brief returns the previous control output
        constexpr const T& y_k1() const {return this->y_k1_;}
        
        /// @brief Sets the positive rate limit of the I controller
        /// @param pv The new positive rate limit
        constexpr void prate_limit(const T& pv) {return this->pv_ = pv;}

        /// @brief returns the positive rate limit
        constexpr const T& prate_limit() const {return this->pv_;}

        /// @brief Sets the negative rate limit of the I controller
        /// @param nv The new negative rate limit
        constexpr void nrate_limit(const T& nv) {return this->nv_ = nv;}

        /// @brief returns the negative rate limit
        constexpr const T& nrate_limit() const {return this->nv_;}

        /// @brief Sets the maximal value limit
        /// @param max The new maximal value limit
        constexpr void max_limit(const T& max) {return this->max_ = max;}

        /// @brief returns the maximal value limit
        constexpr const T& max_limit() const {return this->max_;}

        /// @brief Sets the minimal value limit
        /// @param min The new minimal value limit
        constexpr void min_limit(const T& min) {return this->min_ = min;}

        /// @brief returns the minimal value limit
        constexpr const T& min_limit() const {return this->min_;}

        /**
         * \brief resets (clears) the internal states
         * 
         * resets the internal states to zero
         * 
         * \param u_k1 The previous input `u`
         * \param y_k1 The previous output `y`
         */
        constexpr void reset(const T& u_k1, const T& y_k1){
            this->u_k1_ = static_cast<T>(0);
            this->y_k1_ = static_cast<T>(0);
        }

        /**
         * \brief Adds a new value input and sample time to calculate the next output
         * \details Also advances the internal states
         * \returns The control output
         */
        constexpr T add(const T& u, const T& Ts){
            T y = this->ki_ * Ts / 2 * (u + this->u_k1_) + this->y_k1_;
            
            if(y > this->max_){
                y = this->max_;
            }else if(y < this->min_){
                y = this->min_;
            }else{
                const T rate = (y - y_k1_)/Ts;
                if(rate > this->pv_){
                    y = this->pv_ * Ts;
                }else if(rate < this->nv_){
                    y = this->nv_ * Ts;
                }
            }
            
            this->y_k1_ = y;
            this->u_k1_ = u;
            return y;
        }

        /**
         * \brief Adds a new value input and sample time to calculate the next output
         * \details Also advances the internal states
         * \returns The control output
         * \see controlpp::I::add(const T& u, const T& Ts)
         */
        constexpr T operator() (const T& u, const T& Ts){return this->add(u, Ts);}
    };

    class D{

    };

} // namespace controlpp
