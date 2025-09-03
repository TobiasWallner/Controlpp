#pragma once
/**
 * \file Generator.hpp
 * 
 * This file contains generators that create sampled outputs like:
 * - sine and cosine waves (TODO)
 * - rectangle waves (TODO)
 * - triangle waves (TODO)
 * - ramp waves (TODO)
 */

#include <Eigen/Dense>

#include <controlpp/math.hpp>

namespace controlpp
{
    /**
     * \brief A lite weight incremental sine generator
     * 
     * A lite weight sine generator that calculates incremental results.
     * This allows the use of a fast 2x2 Matrix Vector multiplication,
     * which is much faster than the standard `sin` and `cos` computations.
     * 
     * Thus this class is well suited for generating sine waves that are 
     * written to an output like an DAC (Digital to Analog Converter).
     * 
     * \tparam T The value type that the SineGenerator uses. Typically `float`, `double` or a custom fixpoint implementation.
     */
    template<class T>
    class SineGenerator{
        private:
            T frequency_;
            T sample_time_;
            Eigen::Vector<T, 2> x_; ///< internal states
            Eigen::Matrix<T, 2, 2> A_; ///< dynamic matrix

        public:

            /**
             * \brief Constructs a sine generator 
             * \details The initial position (internal state) is at `sin() == 0` and `cos() == 1`
             * \param frequency The frequency of the sine wave
             * \param sample_time The time that passes in between outputs (calls to `next()`)
             * \see controlpp::SineGenerator::next()
             */
            SineGenerator(const T& frequency, const T& sample_time)
                : frequency_(frequency)
                , sample_time_(sample_time)
                , x_(static_cast<T>(0), static_cast<T>(1))
            {
                // continuous dynamic matrix
                Eigen::Matrix<T, 2, 2> A_continuous({
                    {0, frequency},
                    {-frequency, 0}
                });

                // discretised dynamic matrix
                Eigen::Matrix<T, 2, 2> ATs = A_continuous * sample_time;
                this->A_ = controlpp::mexp(ATs);
            }

            /**
             * \brief Returns the current sine output
             * \details Does not advance the state
             * \returns The sine result
             */
            T sin() const {return this->x_(0);}

            /**
             * \brief Returns the current cosine output
             * \details Does not advance the state
             * \returns The cosine result
             */
            T cos() const {return this->x_(1);}

            /**
             * \brief calculates the next output
             * \details Advances the state
             * \returns An eigen vector containing the sine (at index 0) and cosine (at index 1)
             */
            Eigen::Vector<T, 2> next(){
                Eigen::Vector<T, 2> new_x = this->A_ * this->x_;
                this->x_ = new_x;
                return new_x;
            }

            /**
             * \brief Resets the Sine generator to the initial position
             * \details The initial position (internal state) is at `sin() == 0` and `cos() == 1`
             */
            void reset(){
                this->x_(0) = static_cast<T>(0);
                this->x_(1) = static_cast<T>(1);
            }

    };
} // namespace controlpp

