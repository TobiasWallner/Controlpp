#pragma once
/**
 * \file analysis.hpp
 * 
 * \brief Functions to analyse filters/controllers/transfer functions/state spaces
 * 
 */

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
     * \brief suggests the simulation time and the sample time for a continuous transfer function to get a nice simulation for e.g.: a step response analysis
     * \returns a tuple returning [sample_time, simulation_time]
     */
    template<class T, int NumOrder, int DenOrder>
    std::tuple<T, T> suggest_simulation_times (const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf){
        const Eigen::Vector<T, NumOrder> z = zeros(tf).array().abs();
        const Eigen::Vector<T, DenOrder> p = poles(tf).array().abs();

        // get the smallest non-zero value
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

        // use the larger simulation time
        const T simulation_time_factor = static_cast<T>(20);
        const T smaple_time_factor = static_cast<T>(0.05);
        T simulation_time = static_cast<T>(1);
        T sample_time = static_cast<T>(0.1);
        // calculate the expected simulation time
        if(valid_p && valid_z){
            // if there are only integrative elements (poles at zero), default to 1s
            if(min_p < min_z){
                simulation_time = static_cast<T>(simulation_time_factor) / min_p;
            }else{
                simulation_time = static_cast<T>(simulation_time_factor) / min_z;
            }

            if(max_p > max_z){
                sample_time = static_cast<T>(smaple_time_factor) / max_p;
            }else{
                sample_time = static_cast<T>(smaple_time_factor) / max_z;
            }
        }else if(valid_p && !valid_z){
                simulation_time = static_cast<T>(simulation_time_factor) / min_p;
                sample_time = static_cast<T>(smaple_time_factor) / max_p;
        }else if(!valid_p && valid_z){
                simulation_time = static_cast<T>(simulation_time_factor) / min_z;
                sample_time = static_cast<T>(smaple_time_factor) / max_z;
        }

        return std::tuple<T, T>(sample_time, simulation_time);
    }

    template<class T, int NumOrder, int DenOrder>
    TimeSeries<T> step(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, const T& sample_time, const T& simulation_time){
        const auto css = to_state_space(tf);
        const auto dss = s_to_z(css, sample_time);
        return step(dss, sample_time, simulation_time);
    }

    template<class T, int NumOrder, int DenOrder>
    TimeSeries<T> step(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf){
        const auto [sample_time, simulation_time] = suggest_simulation_times(tf);
        return step(tf, sample_time, simulation_time);
    }
} // namespace controlpp
