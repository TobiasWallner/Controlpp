#pragma once

/**
 * @file
 * @author Tobias Wallner
 * @brief Contains prefilters and trajectory planners for control systems
 */

#include <controlpp/math.hpp>

/**
 * @brief First order pre filter
 * 
 * Uses constant speed to approach a target point
 * 
 * @tparam T The used value type for the calculations, like `float` or `double`
 */
template<class T = float>
class PreFilterO1{
    private:
        T Ts_; ///< The sample time of the prefilter
        T pos_; ///< The current position of the prefilter
        T max_speed_; ///< The maximum speed of the prefilter
        bool reached_ = false; ///< Flag indicating if the prefilter has reached the target position

    public:

        /**
         * @brief Constructs a first order prefilter with a given sample time, maximum speed and initial position
         * @param Ts The sample time of the prefilter
         * @param max_speed The maximum speed of the prefilter
         * @param initial_pos The initial position of the prefilter (default: 0)
         * @throws `std::invalid_argument` if Ts <= 0 or max_speed <= 0
         */
        constexpr PreFilterO1(const T& Ts, const T& max_speed, const T& initial_pos = static_cast<T>(0))
            : Ts_(Ts)
            , pos_(initial_pos)
            , max_speed_(max_speed)
        {
            if(Ts_ <= static_cast<T>(0)){
                throw std::invalid_argument("Sample time must be greater than zero");
            }
            if(max_speed_ <= static_cast<T>(0)){
                throw std::invalid_argument("Max speed must be greater than zero");
            }
        }

        /**
         * @brief Inputs a new target position and returns the next position of the prefilter
         * @param target The target position to approach
         * @return The next position of the prefilter
         */
        constexpr T input(const T& target){
            const T delta = target - this->pos_;
            const T required_speed = delta / this->Ts_;
            T speed = 0;
            if(required_speed < -this->max_speed_){
                this->reached_ = false;
                speed = -this->max_speed_;
            }else if(required_speed > this->max_speed_){
                this->reached_ = false;
                speed = this->max_speed_;
            }else{
                this->reached_ = true;
                speed = required_speed;
            }
            const T new_pos = this->pos_ + speed * this->Ts_;
            this->pos_ = new_pos;
            return new_pos;
        }

        /**
         * @brief Inputs a new target position and returns the next position of the prefilter
         * @param target The target position to approach
         * @return The next position of the prefilter
         */
        constexpr T operator()(const T& target){return this->input(target);}

        /**
         * @brief Resets the prefilter to a new position
         * @param new_pos The new position to reset the prefilter to
         */
        constexpr void reset(const T& new_pos = static_cast<T>(0)){ this->pos_ = new_pos; }

        /**
         * @brief Returns the current position of the prefilter
         * @return The current position of the prefilter
         */
        constexpr T pos() const { return this->pos_; }

        /**
         * @brief Checks if the prefilter has reached the target position
         * @return Returns true if the prefilter has reached the target position, false otherwise
         */
        constexpr bool reached() const { return this->reached_; }
};

template<class T = float>
class PreFilterO2{
    private:
        T Ts_ = 0;
        T pos_ = 0;
        T vel_ = 0;
        T acc_ = 0;

        T max_acc_ = 0;
        T max_vel_ = 0;
        bool reached_ = false;

    public:

        constexpr PreFilterO2(const T& Ts, const T& max_acc, const T& max_vel, const T& initial_pos = static_cast<T>(0), const T& initial_vel = static_cast<T>(0))
            : Ts_(Ts)
            , pos_(initial_pos)
            , vel_(initial_vel)
            , max_acc_(max_acc)
            , max_vel_(max_vel)
        {
            if(Ts_ <= static_cast<T>(0)){
                throw std::invalid_argument("Sample time must be greater than zero");
            }
            if(max_acc_ <= static_cast<T>(0)){
                throw std::invalid_argument("Max acceleration must be greater than zero");
            }
            if(max_vel_ <= static_cast<T>(0)){
                throw std::invalid_argument("Max velocity must be greater than zero");
            }
        }

        /**
         * @brief Calculates the next position of the prefilter based on the target position and velocity and advances the internal states
         * 
         * The prefilter will calculate the next position based on the target position and velocity, the current position and velocity, and the maximum acceleration and velocity. 
         * The prefilter will try to reach the target position and velocity as fast as possible without exceeding the maximum acceleration and velocity. 
         * The prefilter will also try to finish exactly in finite time at the targeted position and velocity
         * 
         * Clamps the target_velocity to the maximum velocity and the acceleration to the maximum acceleration.
         * 
         * @param target_pos The target position to reach
         * @param target_vel The target velocity to reach (default: 0)
         * @return The next position of the prefilter
         */
        constexpr T input(const T& target_pos, const T& target_vel = static_cast<T>(0)){
            // target parabola 1:
            const T a_t1 = this->max_acc_ / T(2);
            const T d_t1 = target_pos - (target_vel * target_vel) / (T(4) * a_t1);

            // target parabola 2:
            const T a_t2 = -this->max_acc_ / T(2);
            const T d_t2 = target_pos - (target_vel * target_vel) / (T(4) * a_t2);

            // source parabola 1: 
            const T a_s1 = this->max_acc_ / T(2);
            const T d_s1 = this->pos_ - (this->vel_ * this->vel_) / (T(4) * a_s1);

            // source parabola 2: 
            const T a_s2 = this->max_acc_ / T(2);
            const T d_s2 = this->pos_ - (this->vel_ * this->vel_) / (T(4) * a_s2);

            // 1 1 --> no need to check --> there is no solution to the kiss
            
            // 1 2
            const std::optional<T> kiss_12 = controlpp::parabola_kiss(a_s1, d_s1, a_t2, d_t2);

            // 2 1
            const std::optional<T> kiss_21 = controlpp::parabola_kiss(a_s2, d_s2, a_t1, d_t1);

            // 2 2 --> no need to check --> there is no solution to the kiss

            // find the correct kiss
            bool k12 = true;
            if((kiss_12.has_value() == true) && (kiss_21.has_value() == false)){
                k12 = true;
            }else if((kiss_12.has_value() == false) && (kiss_21.has_value() == true)){
                k12 = false;
            }else /*if((kiss_12.has_value() == true) && (kiss_21.has_value() == true))*/{
                if(this->pos_ <= kiss_12.value() && kiss_12.value() <= target_pos){
                    k12 = true;
                }else /*if(this->pos_ <= kiss_21.value() && kiss_21.value() <= target_pos)*/{
                    k12 = false;
                }
                //else{
                    // error: no solution
                //}
            }
            //else{
                // error: no solution
            //}

            // get the right combination of polynomials
            // selected source parabol
            const T a_s = k12 ? a_s1 : a_s2;
            const T d_s = k12 ? d_s1 : d_s2;

            // selected target parabola
            const T a_t = k12 ? a_t1 : a_t2;
            const T d_t = k12 ? d_t1 : d_t2;

            // needed distance between the parabola for a kiss
            const T rho = k12 ? kiss_12.value() : kiss_21.value(); //< distance between the parabolas
            const T kiss = -(rho * a_t) / (a_s - a_t); //< point of the kiss measured from the center of the source parabola

            // time point of the source point
            const T t_s = this->vel_ / (2 * a_s);

            const T t = t_s + this->Ts_;
            T new_pos = 0;
            if(t < kiss){
                // next sample is on the source parabola
                new_pos = a_s * t * t + d_s;
            }else{
                // next sample is on the target parabola
                new_pos = a_t * (t - rho) * (t - rho) + d_t;
            }

            // calculate required speed and limit velocity
            const T required_speed = (new_pos - this->pos_) / this->Ts_;
            const T clipped_speed = std::clamp(required_speed, -this->max_vel_, this->max_vel_);
            new_pos = clipped_speed * this->Ts_;

            // update state
            this->acc_ = (clipped_speed - this->vel_) / this->Ts_;
            this->vel_ = clipped_speed;
            this->pos_ = new_pos;

            return new_pos;
        }

        constexpr T operator()(const T& target){return this->input(target);}

        /**
         * @brief Resets the prefilter to a new position and velocity
         * @param new_pos The new position to reset the prefilter to. Defaults to 0.
         * @param new_vel The new velocity to reset the prefilter to. Defaults to 0.
         */
        constexpr void reset(const T& new_pos = static_cast<T>(0), const T& new_vel = static_cast<T>(0)){ 
            this->pos_ = new_pos; 
            this->vel_ = new_vel; 
            this->acc_ = static_cast<T>(0);
        }

        /**
         * @brief Returns the current position of the prefilter
         * @return Returns the current position of the prefilter
         */
        constexpr T pos() const { return this->pos_; }

        /**
         * @brief Returns the current velocity of the prefilter
         * @return Returns the current velocity of the prefilter
         */
        constexpr T vel() const { return this->vel_; }

        /**
         * @brief Returns the current acceleration of the prefilter
         * @return Returns the current acceleration of the prefilter
         */
        constexpr T acc() const { return this->acc_; }

        /**
         * @brief Checks if the prefilter has reached the target position and velocity
         * @return Returns true if the prefilter has reached the target position and velocity, false otherwise
         */
        constexpr bool reached() const { return this->reached_; }
}