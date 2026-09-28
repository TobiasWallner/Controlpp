#pragma once

/**
 * @file
 * @author Tobias Wallner
 * @brief Contains prefilters and trajectory planners for control systems
 */

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
            const t_vel = std::clamp(target_vel, -this->max_vel_, this->max_vel_);
            // calculate the stopping position based on the current velocity and the maximum acceleration
            const T direction = (this->vel_ >= static_cast<T>(0)) ? static_cast<T>(1) : static_cast<T>(-1);
            const T v_2 = t_vel * t_vel;
            const T v0_2 = this->vel_ * this->vel_;
            const T stopping_distance = (v_2 - v0_2) / (static_cast<T>(2) * this->max_acc_);
            const T stopping_pos = this->pos_ + stopping_distance * direction
            const T next_stopping_pos = stopping_pos + this->vel_ * this->Ts_;
            const T distance_to_target = target_pos - this->pos_;

            // maximal distance that can be convered in one step from the target position with the target velocity and the maximal acceleration
            const T eps_step = std::abs((v_2 - v0_2) / (static_cast<T>(2) * this->max_acc_)) + std::abs(this->vel_ * this->Ts_) + std::abs(t_vel * this->Ts_);

            // maximal velocity that can be reached in one step from the target velocity
            const T eps_vel = std::abs(this->vel_ - t_vel) + this->max_acc_ * this->Ts_;

            // check if we are close enough to the target position and velocity
            // check before to avoid overshoot but mostly division by zero in the next step
            if(std::abs(distance_to_target) <= eps_step && std::abs(this->vel_ - t_vel) <= eps_vel){
                this->pos_ = target_pos;
                this->vel_ = t_vel;
                this->reached_ = true;
                return this->pos_;
            }
            this->reached_ = false;

            // calculate the new acceleration based on the current position, velocity, and target position and velocity
            T new_accel = T(0);
            if((this->pos_ < target_pos && next_stopping_pos >= target_pos) || (this->pos_ > target_pos && next_stopping_pos <= target_pos)){
                if(distance_to_target == static_cast<T>(0)){
                    new_accel = static_cast<T>(0);
                }else{
                    // decelerate with exact acceleration
                    new_accel = (v_2 - v0_2) / (static_cast<T>(2) * distance_to_target);
                }
            }else if(stopping_pos < target_pos){
                new_accel = this->max_acc_;
            }else{
                new_accel = -this->max_acc_;
            }

            // clamp the velocity to the maximum velocity
            const T new_vel = std::clamp(this->vel_ + new_accel * this->Ts_, -this->max_vel_, this->max_vel_);
            new_accel = (new_vel - this->vel_) / this->Ts_; // recalculate the acceleration based on the clamped velocity

            // calculate the new position 
            const T new_pos = this->pos_ + new_vel * this->Ts_;
            
            // update the internal states
            this->acc_ = new_accel;
            this->pos_ = new_pos;
            this->vel_ = new_vel;

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