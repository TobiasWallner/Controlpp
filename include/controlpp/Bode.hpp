#pragma once

//std
#include <numbers>
#include <cassert>
#include <concepts>

// eigen
#include <Eigen/Dense>

// controlpp
#include "ContinuousTransferFunction.hpp"

namespace controlpp{

    /**
     * \brief Frequency response data
     * 
     * This class holds frequency response data and can be used to either plot 
     * a transfer function or to do data driven design.
     * 
     * It holds frequencies and complplex magnitudes and allows arithmetic calculations
     * like transfer functions do.
     * 
     * \tparam T The data type of the class. Typically `float`, `double` or a custom fixpoint type
     * \tparam N The size of the bode data entries if known beforehand or `Eigen::Dynamic`
     */
    template<class T = double, int N = Eigen::Dynamic>
    class Bode{
        private:
            Eigen::Vector<T, N> freqs_; // frequencies in Hz
            Eigen::Vector<std::complex<T>, N> values_; // complex magnitudes
    
        public:
            Bode() = default;

            /**
             * @brief Constructs a bode from frequencies and complex magnitudes
             * @param freqs_Hz Frequency vector in Hz
             * @param mags Complex Magnitudes Vector
             */
            Bode(Eigen::Vector<T, N> freqs_Hz, Eigen::Vector<std::complex<T>, N> mags)
                : freqs_(freqs_Hz)
                , values_(mags){
                    assert(this->freqs_.size() == this->values_.size());
                }

            /**
             * \brief Returns a const-reference to the frequency vector in Hz
             * \returns const-reference to an eigen vector
             */
            const Eigen::Vector<T, N>& frequencies() const {
                return this->freqs_;
            }

            /**
             * @brief sets the frequencies assuming the input vector is in Hz
             * @see set_frequencies_Hz
             */
            void set_frequencies(const Eigen::Vector<T, N>& frequs){
                this->freqs_ = frequs;
            }

            /**
             * \brief Returns a const-reference to the complex magnitued vector
             * \returns const-reference to an complex magnitude vector
             */
            const Eigen::Vector<std::complex<T>, N>& values() const {
                return this->values_;
            }

            /**
             * @brief Creates a vector containing the absolute magnitudes
             * 
             * Calculates the absolute magnitudes from the complex magnitudes 
             */
            Eigen::Vector<T, N> magnitudes() const {
                return this->values_.array().abs();
            }

            /**
             * @brief Creates a vector of magnitudes in dB
             * 
             * Calculates the magnitudes from the complex magnitudes 
             */
            Eigen::Vector<T, N> magnitudes_dB() const {
                return static_cast<T>(20) * this->values_.array().abs().log10();
            }

            /**
             * @brief Creates a vector of phases in rad
             * 
             * Calculates the phases from the complex magnitudes
             */
            Eigen::Vector<T, N> phases() const {
                Eigen::Vector<T, N> result = this->values_.array().arg();
                return unwrap(result);
            }

            /**
             * @brief Creates a vector of phases in degree
             * 
             * Calculates the phases from the complex magnitudes
             */
            Eigen::Vector<T, N> phases_deg() const {
                Eigen::Vector<T, N> result = this->values_.array().arg() * static_cast<T>(180) / std::numbers::pi_v<T>;
                return unwrap_deg(result);
            }

            /**
             * @brief Sets the magnitude vector of complex magnitudes
             * @param mags The magnitude vector of complex magnitudes
             * @see set_magnitudes_and_phases
             * @see set_magnitudes_dB_and_phases_deg
             */
            void set_magnitudes(Eigen::Vector<std::complex<T>, N>& mags) {
                this->values_ = mags;
            }

            /**
             * @brief Sets the magnitude data from a vector of magnitudes and phases
             * @param mags A vector of absolute magnitudes (not dB)
             * @param phases A vector of absolute phases in rad (not degree)
             * @see set_magnitudes_dB_and_phases_deg
             */
            void set_magnitudes_and_phases(Eigen::Vector<T, N>& mags, Eigen::Vector<T, N>& phases) {
                auto real = mags.array() * phases.array().cos();
                auto imag = mags.array() * phases.array().sin();
                this->values_.real() = real;
                this->values_.imag() = imag;
            }

            /**
             * @brief Sets the magnitude data from a vector of magnitudes and phases
             * @param mags_dB A vector of absolute magnitudes in dB
             * @param phases_deg A vector of absolute phases in degree
             */
            void set_magnitudes_dB_and_phases_deg(Eigen::Vector<T, N>& mags_dB, Eigen::Vector<T, N>& phases_deg) {
                auto mags = (mags_dB.array() / T(20)).exp10();
                auto phases = phases_deg.array() * std::numbers::pi_v<T> / static_cast<T>(180);
                auto real = mags * phases.cos();
                auto imag = mags * phases.sin();
                this->values_.real() = real;
                this->values_.imag() = imag;
            }

            size_t size() const {return this->freqs_.size();}

            
    };

    // operator +
    //-----------------

    template<class T, int N>
    Bode<T, N> operator+ (const Bode<T, N>& l, const Bode<T, N>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies() && "Both operands need to have the same frequency vector");
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() + r.magnitudes().array();
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, std::convertible_to<T> T2>
    Bode<T, N> operator+ (const Bode<T, N>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() + r;
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, std::convertible_to<T> T2>
    Bode<T, N> operator+ (const T2& l, const Bode<T, N>& r){
        const Eigen::Vector<std::complex<T>, N> result = l + r.magnitudes().array();
        return Bode<T, N>(r.frequencies(), result);
    }

    template<class T, int N>
    Bode<T, N> operator+ (const Bode<T, N>& b){
        const Eigen::Vector<std::complex<T>, N> result = +b.magnitudes().array();
        return Bode<T, N>(b.frequencies(), result);
    }

    template<class T, int N, int NumOrder, int DenOrder>
    Bode<T, N> operator+ (const Bode<T, N>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, N> r_mags = r.eval_frequencies_Hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() + r_mags.array();
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, int NumOrder, int DenOrder>
    Bode<T, N> operator+ (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const Bode<T, N>& r){
        const Eigen::Vector<std::complex<T>, N> l_mags = l.eval_frequencies_Hz(r.frequencies());
        const Eigen::Vector<std::complex<T>, N> result = l_mags.array() + l.magnitudes().array();
        return Bode<T, N>(r.frequencies(), result);
    }

    // operator -
    //-----------------

    template<class T, int N>
    Bode<T, N> operator- (const Bode<T, N>& l, const Bode<T, N>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies() && "Both operands need to have the same frequency vector");
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() - r.magnitudes().array();
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, std::convertible_to<T> T2>
    Bode<T, N> operator- (const Bode<T, N>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() - r;
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, std::convertible_to<T> T2>
    Bode<T, N> operator- (const T2& l, const Bode<T, N>& r){
        const Eigen::Vector<std::complex<T>, N> result = l - r.magnitudes().array();
        return Bode<T, N>(r.frequencies(), result);
    }

    template<class T, int N>
    Bode<T, N> operator- (const Bode<T, N>& b){
        const Eigen::Vector<std::complex<T>, N> result = -b.magnitudes().array();
        return Bode<T, N>(b.frequencies(), result);
    }

    template<class T, int N, int NumOrder, int DenOrder>
    Bode<T, N> operator- (const Bode<T, N>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, N> r_mags = r.eval_frequencies_Hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() - r_mags.array();
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, int NumOrder, int DenOrder>
    Bode<T, N> operator- (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const Bode<T, N>& r){
        const Eigen::Vector<std::complex<T>, N> l_mags = l.eval_frequencies_Hz(r.frequencies());
        const Eigen::Vector<std::complex<T>, N> result = l_mags.array() - l.magnitudes();
        return Bode<T, N>(r.frequencies(), result);
    }

    // operator *
    //-----------------

    template<class T, int N>
    Bode<T, N> operator* (const Bode<T, N>& l, const Bode<T, N>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies() && "Both operands need to have the same frequency vector");
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes() * r.magnitudes();
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, std::convertible_to<T> T2>
    Bode<T, N> operator* (const Bode<T, N>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() * r;
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, std::convertible_to<T> T2>
    Bode<T, N> operator* (const T2& l, const Bode<T, N>& r){
        const Eigen::Vector<std::complex<T>, N> result = l * r.magnitudes().array();
        return Bode<T, N>(r.frequencies(), result);
    }

    template<class T, int N, int NumOrder, int DenOrder>
    Bode<T, N> operator* (const Bode<T, N>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, N> r_mags = r.eval_frequencies_Hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() * r_mags.array();
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, int NumOrder, int DenOrder>
    Bode<T, N> operator* (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const Bode<T, N>& r){
        const Eigen::Vector<std::complex<T>, N> l_mags = l.eval_frequencies_Hz(r.frequencies());
        const Eigen::Vector<std::complex<T>, N> result = l_mags.array() * l.magnitudes().array();
        return Bode<T, N>(r.frequencies(), result);
    }

    // operator /
    //-----------------

    template<class T, int N>
    Bode<T, N> operator/ (const Bode<T, N>& l, const Bode<T, N>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies() && "Both operands need to have the same frequency vector");
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() / r.magnitudes().array();
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, std::convertible_to<T> T2>
    Bode<T, N> operator/ (const Bode<T, N>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() / r;
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, std::convertible_to<T> T2>
    Bode<T, N> operator/ (const T2& l, const Bode<T, N>& r){
        const Eigen::Vector<std::complex<T>, N> result = l / r.magnitudes().array();
        return Bode<T, N>(r.frequencies(), result);
    }

    template<class T, int N, int NumOrder, int DenOrder>
    Bode<T, N> operator/ (const Bode<T, N>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, N> r_mags = r.eval_frequencies_Hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, N> result = l.magnitudes().array() / r_mags.array();
        return Bode<T, N>(l.frequencies(), result);
    }

    template<class T, int N, int NumOrder, int DenOrder>
    Bode<T, N> operator/ (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const Bode<T, N>& r){
        const Eigen::Vector<std::complex<T>, N> result = l.eval_frequencies_Hz(r.frequencies()) / l.magnitudes();
        return Bode<T, N>(r.frequencies(), result);
    }

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
        const auto frequs = bode.frequencies();
        const auto mags = bode.magnitudes_dB();
        const auto phases = bode.phases_deg();
        const int n = std::min({frequs.size(), mags.size(), phases.size()});
        for(int i = 0; i < n; ++i){
            stream << frequs(i) << ", " << mags(i) << ", " << phases(i);
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
     * @param freqs_Hz The frequency vector at which to evaluate the transfer function in Hz
     * @returns The bode result as a `Bode` struct
     * @see ContinuousTransferFunction
     * @see Bode
     */
    template<class T, int NumOrder, int DenOrder, int N>
    Bode<T, N> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const Eigen::Vector<T, N>& freqs_Hz
    ){
        const Eigen::Vector<std::complex<T>, N> complex_magnitudes = tf.eval_frequencies_Hz(freqs_Hz);
        const Bode<T, N> result(freqs_Hz, complex_magnitudes);
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

}