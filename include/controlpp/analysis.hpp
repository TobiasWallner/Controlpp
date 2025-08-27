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

namespace controlpp
{

    /**
     * \brief calculates the step response of a system
     * 
     * \param dss A discrete state space model of the system
     * \param Ts The sampling frequency
     * \param simulation_time The time to be simulated
     * 
     * \returns A tuple [time, value] that holds the result of the step analysis
     */
    template<class T, int NStates>
    std::tuple<std::vector<T>, std::vector<T>> step(const DiscreteStateSpace<T, NStates, 1, 1>& dss, double Ts, double simulation_time){
        DssFilter filter(dss);
        
        std::vector<T> times;
        std::vector<T> values;
        
        int expected_size = simulation_time/Ts + 1;
        times.reserve(expected_size);
        values.reserve(expected_size);
        
        for(T time = static_cast<T>(0); time < simulation_time; time += Ts){
            const T value = filter(static_cast<T>(1));
            times.push_back(time);
            values.push_back(value);
        }

        return std::tuple<std::vector<T>, std::vector<T>>(std::move(times), std::move(values));
    }
} // namespace controlpp
