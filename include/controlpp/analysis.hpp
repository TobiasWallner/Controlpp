#pragma once
/**
 * \file analysis.hpp
 * 
 * \brief Functions to analyse filters/controllers/transfer functions/state spaces
 * 
 */

 // std
#include <complex>
#include <numbers>

// eigen
#include <Eigen/Dense>

// controlpp
#include <controlpp/DiscreteStateSpace.hpp>
#include <controlpp/DiscreteTransferFunction.hpp>
#include <controlpp/DiscreteFilter.hpp>
#include <controlpp/TimeSeries.hpp>

namespace controlpp
{
    /**
     * \brief calculates the step response of a system
     * 
     * \param dss A discrete state space model of the system
     * \param Ts The sampling frequency
     * \param simulation_time The time to be simulated
     * 
     * \returns A time series of dynamically sized `Eigen::Vector`s
     */
    template<class T, int NStates>
    TimeSeries<T> step(const DiscreteStateSpace<T, NStates, 1, 1>& dss, double Ts, double simulation_time){
        DssFilter filter(dss);
        
        const int expected_size = simulation_time/Ts;
        
        TimeSeries<T> timeseries(expected_size);
        
        T time = static_cast<T>(0);
        int i = 0;
        for(; i < expected_size; time += Ts, (void)++i){
            const T value = filter(static_cast<T>(1));
            timeseries.times(i) = time;
            timeseries.values(i) = value;
        }

        return timeseries;
    }

    /**
     * \brief Calculates the slowest (lowest) and fastest (highest) frequencies of a continuous transfer function
     * 
     * Example:
     * 
     * \code{.cpp}
     * const auto [slowest, fastest] = slowest_fastest_frequencies(tf);
     * \endcode
     * 
     * \param tf A continuous time transfer function
     * \param alternative The value to be returned if the transfer function has no dynamics and frequencies to be evaluated.
     * \returns A tuple `[slowest_frequency, fastest_frequency]` containing the fastest and slowest frequencies or the alternative.
     */
    template<class T, int NumOrder, int DenOrder>
    std::tuple<T, T> slowest_fastest_frequencies(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, T alternative = static_cast<T>(1)){
        // calculate the poles and zeros
        const Eigen::Vector<T, NumOrder> z = zeros(tf).array().abs();
        const Eigen::Vector<T, DenOrder> p = poles(tf).array().abs();

        // find the minimum and maximum of the zeros that is a non-zero number
        T max_z = std::numeric_limits<T>::lowest();
        T min_z = std::numeric_limits<T>::max();
        bool valid_z = false;
        for(int i = 0; i < z.size(); ++i){
            if(z(i) != static_cast<T>(0)){
                min_z = (z(i) < min_z) ? z(i) : min_z;
                max_z = (z(i) > max_z) ? z(i) : max_z;
                valid_z = true;
            }
        }

        // find the minimum and maximum of the poles that is a non-zero number
        T max_p = std::numeric_limits<T>::lowest();
        T min_p = std::numeric_limits<T>::max();
        bool valid_p = false;
        for(int i = 0; i < p.size(); ++i){
            if(p(i) != static_cast<T>(0)){
                min_p = (p(i) < min_p) ? p(i) : min_p;
                max_p = (p(i) > max_p) ? p(i) : max_p;
                valid_p = true;
            }
        }

        // get the fastest and slowest frequency if valid
        T fastest_frequency = alternative;
        T slowest_frequency = alternative;
        if(valid_p && valid_z){
            slowest_frequency = (min_p < min_z) ? min_p : min_z;
            fastest_frequency = (max_p > max_z) ? max_p : max_z;
        }else if(valid_p && !valid_z){
                slowest_frequency = min_p;
                fastest_frequency = max_p;
        }else if(!valid_p && valid_z){
                slowest_frequency = min_z;
                fastest_frequency = max_z;
        }

        // return result as tuple
        return std::tuple<T, T>(slowest_frequency, fastest_frequency);
    }

    template<class T, int NumOrder, int DenOrder>
    TimeSeries<T> step(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, const T& sample_time, const T& simulation_time){
        const auto css = to_state_space(tf);
        const auto dss = discretise_zoh(css, sample_time);
        return step(dss, sample_time, simulation_time);
    }

    template<class T, int NumOrder, int DenOrder>
    TimeSeries<T> step(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf){
        const auto [slowest_freq, fastest_freq] = slowest_fastest_frequencies(tf);
        const T sample_time = static_cast<T>(0.05) / fastest_freq;
        const T simulation_time = static_cast<T>(20) / slowest_freq;
        return step(tf, sample_time, simulation_time);
    }

    

} // namespace controlpp
