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

    /**
     * \brief contains a bodeplot of frequencies (Hz), magnitudes (dB) and phases (deg)
     */
    template<class T, int N = Eigen::Dynamic>
    struct Bode{
        Eigen::Vector<T, N> frequencies;
        Eigen::Vector<T, N> magnitudes;
        Eigen::Vector<T, N> phases;
    };

    /**
     * \brief Prints a bode plot to an output stream as a `.csv` file
     * \param stream The stream to be printed to
     * \param bode The bode container with the frequencies, magnitudes and phases
     * \returns A reference to the stream object for operation chaining.
     * \see Bode
     */
    template<class T, int N = Eigen::Dynamic>
    std::ostream& operator<< (std::ostream& stream, const Bode<T, N>& bode){
        stream << "Frequencies (Hz), Magnitudes (dB), Phases (deg)" << std::endl;
        const int n = std::min({bode.frequencies.size(), bode.magnitudes.size(), bode.phases.size()});
        for(int i = 0; i < n; ++i){
            stream << bode.frequencies(i) << ", " << bode.magnitudes(i) << ", " << bode.phases(i);
            if(i < n-1) stream << "\n";
        }
        return stream;
    }

    /**
     * @brief Calculates the bode response for a pre defined frequency vector
     * @tparam T The value type of the transfer function
     * @tparam NumOrder The numerator order
     * @tparam DenOrder The denominator order
     * @tparam NSize The number of elements in the frequency vector (may also be `Eigen::Dynamic`)
     * @param tf The continuous transfer function to analyse
     * @param freqs_Hz The frequency vector at which to evaluate the transfer function
     * @returns The bode result as a `Bode` struct
     * @see ContinuousTransferFunction
     * @see Bode
     */
    template<class T, int NumOrder, int DenOrder, int NSize>
    Bode<T, NSize> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const Eigen::Vector<T, NSize>& freqs_Hz
    ){
        const Eigen::Vector<std::complex<T>, NSize> complex_magnitudes = tf.eval_frequencies_Hz(freqs_Hz);
        
        Bode<T, NSize> result;
        result.frequencies = freqs_Hz;
        result.magnitudes = complex_magnitudes.array().abs().log10() * static_cast<T>(20);
        result.phases = complex_magnitudes.array().arg() * static_cast<T>(180 / std::numbers::pi_v<T>);
        result.phases = unwrap_deg(result.phases);
        
        return result;
    }

    /**
     * \brief Calculates the bode response of a transfer function
     * \param slowest_freq_Hz The slowest/lowest frequency in Hz from which to calculate frequency responses
     * \param fastest_freq_Hz The fastest/highest frequency in Hz to which to calculate the frequency response
     * \param samples_per_decade The number of samples per decade of frequencies to be calculated
     * \returns A Bode containing the frequencies (Hz), magnitues (dB) and phases (deg)
     * \see Bode
     * \see ContinuousTransferFunction::eval_frequencies
     * \see bode(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, const int samples_per_decade)
     */
    template<class T, int NumOrder, int DenOrder, std::convertible_to<T> T1, std::convertible_to<T> T2>
    Bode<T, Eigen::Dynamic> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const T1& slowest_freq_Hz, 
            const T2& fastest_freq_Hz, 
            const int samples_per_decade=100
    ){
        const T decades = std::log10(fastest_freq_Hz) - std::log10(slowest_freq_Hz);
        const T samples = samples_per_decade * decades;

        const Eigen::Vector<T, Eigen::Dynamic> freqs_Hz = Eigen::Vector<T, Eigen::Dynamic>::LinSpaced(samples, std::log(slowest_freq_Hz), std::log(fastest_freq_Hz)).array().exp();
        return bode<T, NumOrder, DenOrder, Eigen::Dynamic>(tf, freqs_Hz);
    }

    /**
     * \brief Calculates the bode response of a transfer function. 
     * 
     * Infers the frequency range of the bode plot from the transfer function.
     * 
     * \param tf The transfer function
     * \param samples_per_decade The number of samples per decade of frequencies to be calculated
     * \returns A Bode containing the frequencies (Hz), magnitues (dB) and phases (deg)
     * \see Bode
     * \see bode(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, const T& slowest_freq_Hz, const T& fastest_freq_Hz, const int samples_per_decade)
     * \see ContinuousTransferFunction::eval_frequencies
     */
    template<class T, int NumOrder, int DenOrder>
    Bode<T, Eigen::Dynamic> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const int samples_per_decade=100
    ){
        const auto [slowest_freq_rad, fastest_freq_rad] = slowest_fastest_frequencies(tf);
        const T frequency_from_Hz = slowest_freq_rad / (static_cast<T>(10 * 2) * std::numbers::pi_v<T>);
        const T frequency_to_Hz = fastest_freq_rad * static_cast<T>(10 / 2) / std::numbers::pi_v<T>;
        return bode(tf, frequency_from_Hz, frequency_to_Hz, samples_per_decade);
    }

} // namespace controlpp
