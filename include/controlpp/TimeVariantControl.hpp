#pragma once
/**
 * \file TimeVariantControl.hpp
 * 
 * \brief Contains controllers for systems with runtime varying sample-times like high time jitters
 * 
 * This file contains controllers that have been analytically/symbolically calculated to allow for varying runtime sample-times.
 *
 */

#include <algorithm>
#include <limits>

#include <controlpp/ContinuousTransferFunction.hpp>
#include <controlpp/ContinuousStateSpace.hpp>
#include <controlpp/DiscreteTransferFunction.hpp>
#include <controlpp/DiscreteStateSpace.hpp>
#include <controlpp/transformations.hpp>

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
    class PControl{
        public:
        using value_type = T;
        private:

        T kp_;

        public:

        /// @brief Default copy constructor
        constexpr PControl(const PControl&) = default;

        /**
         * \brief Construct a P (Proportional) controller with a constant gain
         * \param kp the proportional gain of the controller
         */
        constexpr PControl(const T& kp) : kp_(kp){}

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
        constexpr T input(const T& u) {return u * this->kp_;}

        /**
         * \brief Adds a new sample
         * \param u the new sample input
         * \param Ts sample time (unused - there for compatibility)
         * \returns the new output based on the input
         */
        constexpr T input(const T& u, [[maybe_unused]]const T& Ts) {return this->input(u);}

        /**
         * \brief Adds a new sample
         * \param u the new sample input
         * \returns the new output based on the input
         */
        constexpr T operator() (const T& u) {return this->input(u);}

        /**
         * \brief Adds a new sample
         * \param u the new sample input
         * \param Ts sample time (unused - there for compatibility)
         * \returns the new output based on the input
         */
        constexpr T operator() (const T& u, [[maybe_unused]]const T& Ts) {return this->input(u);}

        /**
         * \brief resets the internal states of the filter
         * 
         * For P control elements this is a no-op.
         */
        constexpr void reset() {/*no operation*/}
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
        class IControl{
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
            constexpr IControl(const T& ki) : ki_(ki){}

            /// @brief Default copy constructor 
            constexpr IControl(const IControl&) = default;

            /// @brief Sets the integrator gain
            /// @param ki The integrator gain for the next sample
            constexpr void set_param(const T& ki) {return this->ki_ = ki;}

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
            T min_;
            T max_;
            T vn_;
            T vp_;
            
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
            constexpr IAntiWindup(const T& ki, const T& min, const T& max, const T& vn = std::numeric_limits<T>::lowest(), const T& vp = std::numeric_limits<T>::max()) 
                : ki_(ki)
                , min_(min)
                , max_(max)
                , vn_(vn)
                , vp_(vp)
                {}

            /// @brief Default copy constructor 
            constexpr IAntiWindup(const IAntiWindup&) = default;

            /// @brief Sets the integrator gain
            /// @param ki The integrator gain for the next sample;
            constexpr void set_param(const T& ki) {return this->ki_ = ki;}
            
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
        class DControl{
            private:
            ContinuousStateSpace<T, 1> css_;
            Eigen::Vector<T, 1> states_ = Eigen::Vector<T, 1>::Zero();

            public:

            /**
             * \brief Construct a D control element
             * \param k_d The differential gain
             * \param omega The bandwidth of the low-pass filter in radians per second
             */
            constexpr DControl(const T& k_d, const T& omega){
                this->set_params(k_d, omega);
            }

            /// @brief default copy constructor 
            constexpr DControl(const DControl&) = default;

            constexpr void set_params(const T& k_d, const T& omega){
                const ContinuousTransferFunction<T, 1, 1> tf(
                    Eigen::Vector<T, 2>(static_cast<T>(0), k_d), 
                    Eigen::Vector<T, 2>(static_cast<T>(1), static_cast<T>(1) / omega));
                this->css_ = to_state_space(tf);
            }

            constexpr T input(const T& u, const T& Ts){
                const DiscreteStateSpace dss = discretise_tustin(this->css_, Ts);
                const auto [x, y] = dss.eval(this->states_, u);
                this->states_ = x;
                return y;
            }

            constexpr T operator() (const T& u, const T& Ts) {return this->input(u, Ts);}

            constexpr void reset(){
                this->states_.setZero();
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
         * y_k = \frac{a  (u_k + u_{k-1}) - b  y_{k-1}}{c};
         * a = K  \omega_ * T_s
         * b = \omega * T_s - 2
         * c = \omega * T_s + 2
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
        class PT1Control{
            private:
            T _K = static_cast<T>(1);       ///< The gain of the filter
            T _omega = static_cast<T>(1);   ///< The characteristic frequency of the filter
            T _x = static_cast<T>(0);       ///< The internal state of the filter
        
            public:
            
            /**
             * \brief Constructs a low-pass filter
             * \param K The gain of the filter
             * \param omega The -3dB bandwidth of the lowpass filter
             */
            constexpr PT1Control(const T& K, const T& omega)
                : _K(K)
                , _omega(omega){}

            /// @brief Default copy constructor 
            constexpr PT1Control(const PT1Control&) = default;

            constexpr void set_gain(const T& g){this->_K = g;}
            [[nodiscard]] constexpr const T& gain() const {return this->_K;}

            constexpr void set_omega(const T& w){this->_omega = w;}
            [[nodiscard]] constexpr const T& omega() const {return this->_omega;}

            constexpr T input(const T& u, const T& Ts){
                // helper values
                const T p = (Ts * this->_omega + 2);
                const T n = (Ts * this->_omega - 2);

                // parameters
                // a_0 = 1
                const T a = n / p;
                const T b = this->_K * Ts * this->_omega / p;

                // calculate new state and output
                const T new_x = (-a) * this->_x + u;
                const T new_y = b * ((1 - a) * this->_x + u);

                // update state and output
                this->_x = new_x;
                return new_y;
            }

            constexpr T operator() (const T& u, const T& Ts) {return this->input(u, Ts);}

            constexpr void reset(){
                this->_x = static_cast<T>(0);
            }
        };

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
         * 
         * TODO: calculate the tustin discretised states space manuall and build the state space matrices directly to improve performance-
         */
        template<class T>
        class PT2Control{
            private:
            ContinuousStateSpace<T, 2> css_;
            Eigen::Vector<T, 2> states_ = Eigen::Vector<T, 2>::Zero();

            public:
            PT2Control(const PT2Control&) = default;

            /**
             * \brief constructs a PT2 element
             * \param K The gain of the filter
             * \param D Then dampening factor of the filter
             * \param omega The characteristic frequency of the filter (bandwidth)
             */
            PT2Control(const T& K, const T& D, const T& omega){
                this->set_params(K, D, omega);
            }

            constexpr void set_params(const T& K, const T& D, const T& omega){
                const ContinuousTransferFunction<T, 0, 2> tf(
                    Eigen::Vector<T, 1>(K), 
                    Eigen::Vector<T, 3>(static_cast<T>(1), static_cast<T>(2) * D / omega, static_cast<T>(1)/(omega * omega)));
                this->css_ = to_state_space(tf);
            }

            constexpr T input(const T& u, const T& Ts){
                const DiscreteStateSpace dss = discretise_tustin(this->css_, Ts);
                const auto [x, y] = dss.eval(this->states_, u);
                this->states_ = x;
                return y;
            }

            constexpr T operator() (const T& u, const T& Ts) {return this->input(u, Ts);}

            constexpr void reset(){
                this->states_.setZero();
            }
        };
        template<class T>
        using LowPassO2 = PT2Control<T>;

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
        class PIControl{
            private:
            ContinuousStateSpace<T, 1> css_;
            Eigen::Vector<T, 1> states_ = Eigen::Vector<T, 1>::Zero();

            public:

            /**
             * \param ki The integral gain
             * \param kp The proportional gain
             */
            PIControl(const T& kp, const T& ki){
                this->set_params(kp, ki);
            }

            PIControl(const PIControl&) = default;

            constexpr void set_params(const T& kp, const T& ki){
                const ContinuousTransferFunction<T, 1, 1> tf(
                    Eigen::Vector<T, 2>(ki, kp), 
                    Eigen::Vector<T, 2>(static_cast<T>(0), static_cast<T>(1)));
                this->css_ = to_state_space(tf);
            }

            constexpr T input(const T& u, const T& Ts){
                const DiscreteStateSpace dss = discretise_tustin(this->css_, Ts);
                const auto [x, y] = dss.eval(this->states_, u);
                this->states_ = x;
                return y;
            }

            constexpr T operator() (const T& u, const T& Ts) {return this->input(u, Ts);}

            constexpr void reset(){
                this->states_.setZero();
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
            PControl<T> P_;
            IAntiWindup<T> I_;

            T y_k1_ = static_cast<T>(0);

            public:

            /// @brief Constructs a PI (proportional + integral) controller with anti-windup
            /// @param kp The proportional gain
            /// @param ki The integral gain.
            /// @param vn The minimal negative velocity of the control output. Assumes: \f$v_n < 0 < v_p\f$
            /// @param vp The maximal positive velocity of the control output. Assumes: \f$v_n < 0 < v_p\f$
            /// @param min The minimal value of the control output. Assumes: \f$min < max\f$
            /// @param max The maximal value of the control output. Assumes: \f$min < max\f$
            PIAntiWindup(const T& kp, const T& ki, const T& min, const T& max, const T& vn = std::numeric_limits<T>::lowest(), const T& vp = std::numeric_limits<T>::max())
                : P_(kp)
                , I_(ki, min, max, vn, vp)
                {}

            PIAntiWindup(const PIAntiWindup&) = default;

            /**
             * \brief sets the proportional gain
             * \param kp the new proportional gain for the next sample value
             */
            constexpr void kp(const T& kp){this->P_.kp(kp);}

            /**
             * \brief returns the current gain
             */
            [[nodiscard]] constexpr const T& kp() const {return this->P_.kp();}

            /**
             * \brief sets the integral gain
             * \param ki the new proportional gain for the next sample value
             */
            constexpr void ki(const T& ki){this->I_.ki(ki);}

            /**
             * \brief returns the current gain
             */
            [[nodiscard]] constexpr const T& ki() const {return this->I_.ki_();}

            /// @brief Sets the positive rate limit of the I controller
            /// @param pv The new positive rate limit
            constexpr void prate_limit(const T& pv) {this->I_.prate_limit(pv);}

            /// @brief returns the positive rate limit
            [[nodiscard]] constexpr const T& prate_limit() const {return this->I_.prate_limit();}

            /// @brief Sets the negative rate limit of the I controller
            /// @param nv The new negative rate limit
            constexpr void nrate_limit(const T& nv) {this->I_.nrate_limit(nv);}

            /// @brief returns the negative rate limit
            [[nodiscard]] constexpr const T& nrate_limit() const {return this->I_.nrate_limit();}

            /// @brief Sets the maximal value limit
            /// @param max The new maximal value limit
            constexpr void max_limit(const T& max) {this->I_.max_limit(max);}

            /// @brief returns the maximal value limit
            [[nodiscard]] constexpr const T& max_limit() const {return this->I_.max_limit();}

            /// @brief Sets the minimal value limit
            /// @param min The new minimal value limit
            constexpr void min_limit(const T& min) {this->I_.min_limit(min);}

            /// @brief returns the minimal value limit
            [[nodiscard]] constexpr const T& min_limit() const {return this->I_.min_limit();}

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see operator()(const T& u, const T& Ts)
             */
            constexpr T input(const T& u, const T& Ts){
                // calculate new output
                T y_raw = I_.input(u, Ts) + P_.input(u);
                T rate_raw = (y_raw - y_k1_)/Ts;
                T rate = std::clamp(rate_raw, nrate_limit(), prate_limit());
                T y_rate_clamped = rate * Ts + y_k1_;
                T y = std::clamp(y_rate_clamped, min_limit(), max_limit());

                // update states
                this->y_k1_ = y;

                return y;
            }

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see input(const T& u, const T& Ts)
             */
            constexpr T operator()(const T& u, const T& Ts){return this->input(u, Ts);}

            /**
             * \brief resets (clears) the internal states
             * resets the internal states (\f$u_{k-1}\f$, \f$y_{k-1}\f$) to zero
             */
            constexpr void reset(){
                this->I_.reset();
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
        class PT1IControl{
            private:

            PT1Control<T> PT1_;
            IControl<T> I_;

            public:

            constexpr PT1IControl(const PT1IControl&) = default;

            /**
             * \brief Constructs a PT1I (proportional time delayed + integral) controller
             * \param kp The proportional gain of the controller
             * \param ki The integral gain of the controller
             * \param omega The lowpass filter bandwidth
             */
            constexpr PT1IControl(const T& kp, const T& ki, const T& omega)
                : PT1_(kp, omega)
                , I_(ki)
                {}

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see operator()(const T& u, const T& Ts)
             */
            constexpr T input(const T& u, const T& Ts){
                // calculate filter parameters
                const T result = PT1_(u, Ts) + I_(u, Ts);
                return result;
            }

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see input(const T& u, const T& Ts)
             */
            constexpr T operator()(const T& u, const T& Ts){return this->input(u, Ts);}

            /**
             * \brief resets (clears) the internal states
             * resets the internal states (\f$u_{k-1}\f$, \f$u_{k-2}\f$, \f$y_{k-1}\f$ and \f$y_{k-2}\f$) to zero
             */
            constexpr void reset(){
                PT1_.reset();
                I_.reset();
            }
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
        class PIDControl{
            private:
            PControl<T> P_;
            IControl<T> I_;
            DControl<T> D_;

            public:

            /**
             * \brief constructs a PID controller form gain parameters
             * \param kp The proportional gain
             * \param ki The integral gain
             * \param kd The differential gain
             * \param omega The thaming frequency
             */
            PIDControl(const T& kp, const T& ki, const T& kd, const T& omega)
                : P_(kp)
                , I_(ki)
                , D_(kd, omega){}

            /**
             * \brief Constructa a PID controller using the alpha tuning method
             * 
             * TODO: insert image with explanation of how the alpha tuning works
             * 
             * \param alpha Design parameter
             *  - Large alpha: 
             *      - increases phase margin but decreases gain at low frequencies
             *      - more robust but less performance
             *  - Small alpha:
             *      - decreses phase margin but increses gain at low frequencies
             *      - less robust but more performant
             * 
             * \param mag_G The magnitude of the plant that should be controlled at the desired cutoff frequency
             * 
             * \param omega_c The desired cutoff frequency
             */
            static PIDControl from_alpha_tune(const T& alpha, const T& omega_c, const T& mag_G){
                const T kp = 1 / (mag_G * alpha);
                const T ki = omega_c / (mag_G * alpha * alpha * alpha);
                const T kd = 1 / (omega_c * mag_G);
                const T omega = omega_c * alpha;
                PIDControl result(kp, ki, ki, omega);
                return result;
            }

            /**
             * \brief Input a new sample to the controller
             * \param u The control input
             * \param Ts The sample-time
             * \see operator()(const T& u, const T& Ts)
             */
            constexpr T input(const T& u, const T& Ts){
                const T result = P_(u, Ts) + I_(u, Ts) + D_(u, Ts);
                return result;
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
                this->I_.reset();
                this->D_.reset();
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
        class LeadLagControl{
            private:
            ContinuousStateSpace<T, 1> css_;
            Eigen::Vector<T, 1> states_ = Eigen::Vector<T, 1>::Zero();

            public:

            constexpr LeadLagControl(const T& omega1, const T& omega2){
                this->set_params(omega1, omega2);
            }

            constexpr void set_params(const T& omega1, const T& omega2){
                const ContinuousTransferFunction<T, 1, 1> tf(
                    Eigen::Vector<T, 2>(static_cast<T>(1), static_cast<T>(1) / omega1), 
                    Eigen::Vector<T, 2>(static_cast<T>(1), static_cast<T>(1) / omega2));
                
                this->css_ = to_state_space(tf);
            }

            constexpr T input(const T& u, const T& Ts){
                const DiscreteStateSpace dss = discretise_tustin(this->css_, Ts);
                const auto [x, y] = dss.eval(this->states_, u);
                this->states_ = x;
                return y;
            }

            constexpr T operator() (const T& u, const T& Ts) {return this->input(u, Ts);}

            constexpr void reset(){
                this->states_.setZero();
            }
        };

        /**
         * \brief A notch element with variable sample-time
         * 
         * Has the following continuous transfer function:
         * 
         * \f[
         * N(s) = 
         *   \frac{\omega^2 + 2 G_\text{min} W s \omega + s^2}
         *   {\omega^2 + 2 W s \omega + s^2}
         * \f]
         * 
         * where:
         *  - \f$G_\text{min}\f$ is the value at the notch
         *  - \f$W\f$ is the width of the notch
         * 
         */
        template<class T>
        class NotchControl{
            private:
            ContinuousStateSpace<T, 2> css_;
            Eigen::Vector<T, 2> states_ = Eigen::Vector<T, 2>::Zero();

            public:

            /**
             * \brief Creates a notch filter
             * \param w The width of the noth. 
             *      - \f$w > 0\f$: Larger value \f$\Rightarrow\f$ wider notch
             *      - \f$w = 0\f$: Needle
             * \param g_min The value at the notch peak
             * \param omega the frequency (rad/s) of the notch peak
             */
            constexpr NotchControl(const T& w, const T& g_min, const T& omega){
                this->set_params(w, g_min, omega);
            }

            constexpr void set_params(const T& w, const T& g_min, const T& omega){
                const ContinuousTransferFunction<T, 2, 2> tf(
                    Eigen::Vector<T, 3>((omega * omega), static_cast<T>(2) * g_min * w * omega, static_cast<T>(1)), 
                    Eigen::Vector<T, 3>((omega * omega), static_cast<T>(2) * w * omega, static_cast<T>(1)));
                
                this->css_ = to_state_space(tf);
            }

            constexpr T input(const T& u, const T& Ts){
                const DiscreteStateSpace dss = discretise_tustin(this->css_, Ts);
                const auto [x, y] = dss.eval(this->states_, u);
                this->states_ = x;
                return y;
            }

            constexpr T operator() (const T& u, const T& Ts) {return this->input(u, Ts);}

            constexpr void reset(){
                this->states_.setZero();
            }
        };

    } // vartime
} // namespace controlpp
