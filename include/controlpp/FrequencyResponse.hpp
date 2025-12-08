#pragma once

//std
#include <numbers>
#include <cassert>
#include <concepts>
#include <complex>

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
     * 
     * 
     * \tparam T The data type of the class. Typically `float`, `double` or a custom fixpoint type
     */
    template<class T = double>
    class FrequencyResponse{
        private:
            Eigen::Vector<T, Eigen::Dynamic> freqs_; // frequencies in Hz
            Eigen::Vector<std::complex<T>, Eigen::Dynamic> values_; // complex magnitudes
    
        public:
            FrequencyResponse() = default;

            /**
             * @brief Constructs a bode from frequencies and complex magnitudes
             * @param freqs_Hz Ordered (ascending) frequency vector in Hz
             * @param mags Complex Magnitudes Vector that corresponds to the frequencies
             */
            FrequencyResponse(const Eigen::Vector<T, Eigen::Dynamic>& freqs_Hz, const Eigen::Vector<std::complex<T>, Eigen::Dynamic>& values)
                : freqs_(freqs_Hz)
                , values_(values)
            {
                // assert same size
                assert(this->freqs_.size() == this->values_.size());

                // assert only positive frequencies
                assert([&](){
                    for(size_t i = 0; i < freqs_.size(); ++i){
                        if(freqs_(i) < 0){
                            return false;
                        }
                    }
                    return true;
                });

                // assert frequencies are ascending
                assert([&](){
                    for(size_t i = 1; i < freqs_.size(); ++i){
                        if(freqs_(i) < freqs_(i-1)){
                            return false;
                        }
                    }
                    return true;
                });
            }

            FrequencyResponse(Eigen::Vector<T, Eigen::Dynamic>&& freqs_Hz, Eigen::Vector<std::complex<T>, Eigen::Dynamic>&& values)
                : freqs_(std::move(freqs_Hz))
                , values_(std::move(values))
            {
                // assert same size
                assert(this->freqs_.size() == this->values_.size());

                // assert only positive frequencies
                assert([&](){
                    for(size_t i = 0; i < freqs_.size(); ++i){
                        if(freqs_(i) < 0){
                            return false;
                        }
                    }
                    return true;
                });

                // assert frequencies are ascending
                assert([&](){
                    for(size_t i = 1; i < freqs_.size(); ++i){
                        if(freqs_(i) < freqs_(i-1)){
                            return false;
                        }
                    }
                    return true;
                });
            }

            FrequencyResponse(const Eigen::Vector<T, Eigen::Dynamic>& freqs_Hz, Eigen::Vector<std::complex<T>, Eigen::Dynamic>&& values)
                : freqs_(freqs_Hz)
                , values_(std::move(values))
            {
                // assert same size
                assert(this->freqs_.size() == this->values_.size());

                // assert only positive frequencies
                assert([&](){
                    for(size_t i = 0; i < freqs_.size(); ++i){
                        if(freqs_(i) < 0){
                            return false;
                        }
                    }
                    return true;
                });

                // assert frequencies are ascending
                assert([&](){
                    for(size_t i = 1; i < freqs_.size(); ++i){
                        if(freqs_(i) < freqs_(i-1)){
                            return false;
                        }
                    }
                    return true;
                });
            }

            FrequencyResponse(Eigen::Vector<T, Eigen::Dynamic>&& freqs_Hz, const Eigen::Vector<std::complex<T>, Eigen::Dynamic>& values)
                : freqs_(std::move(freqs_Hz))
                , values_(values)
            {
                // assert same size
                assert(this->freqs_.size() == this->values_.size());

                // assert only positive frequencies
                assert([&](){
                    for(size_t i = 0; i < freqs_.size(); ++i){
                        if(freqs_(i) < 0){
                            return false;
                        }
                    }
                    return true;
                });

                // assert frequencies are ascending
                assert([&](){
                    for(size_t i = 1; i < freqs_.size(); ++i){
                        if(freqs_(i) < freqs_(i-1)){
                            return false;
                        }
                    }
                    return true;
                });
            }

            /**
             * \brief Returns a const-reference to the frequency vector in Hz
             * \returns const-reference to an eigen vector
             */
            const Eigen::Vector<T, Eigen::Dynamic>& frequencies() const {
                return this->freqs_;
            }

            /**
             * @brief sets the frequencies assuming the input vector is in Hz
             * @see set_frequencies_Hz
             */
            void set_frequencies(const Eigen::Vector<T, Eigen::Dynamic>& frequs){
                this->freqs_ = frequs;
            }

            /**
             * \brief Returns a const-reference to the complex magnitued vector
             * \returns const-reference to an complex magnitude vector
             */
            const Eigen::Vector<std::complex<T>, Eigen::Dynamic>& values() const {
                return this->values_;
            }

            Eigen::Vector<T, Eigen::Dynamic> real() const {
                return this->values_.real();
            }

            Eigen::Vector<T, Eigen::Dynamic> imag() const {
                return this->values_.imag();
            }

            /**
             * @brief Creates a vector containing the absolute magnitudes
             * 
             * Calculates the absolute magnitudes from the complex magnitudes 
             */
            Eigen::Vector<T, Eigen::Dynamic> magnitudes() const {
                return this->values_.array().abs();
            }

            /**
             * @brief Creates a vector of magnitudes in dB
             * 
             * Calculates the magnitudes from the complex magnitudes 
             */
            Eigen::Vector<T, Eigen::Dynamic> magnitudes_dB() const {
                return static_cast<T>(20) * this->values_.array().abs().log10();
            }

            /**
             * @brief Creates a vector of phases in rad
             * 
             * Calculates the phases from the complex magnitudes
             */
            Eigen::Vector<T, Eigen::Dynamic> phases() const {
                Eigen::Vector<T, Eigen::Dynamic> result = this->values_.array().arg();
                return unwrap(result);
            }

            /**
             * @brief Creates a vector of phases in degree
             * 
             * Calculates the phases from the complex magnitudes
             */
            Eigen::Vector<T, Eigen::Dynamic> phases_deg() const {
                Eigen::Vector<T, Eigen::Dynamic> result = this->values_.array().arg() * static_cast<T>(180) / std::numbers::pi_v<T>;
                return unwrap_deg(result);
            }

            /**
             * @brief Sets the magnitude vector of complex magnitudes
             * @param mags The magnitude vector of complex magnitudes
             * @see set_magnitudes_and_phases
             * @see set_magnitudes_dB_and_phases_deg
             */
            void set_magnitudes(Eigen::Vector<std::complex<T>, Eigen::Dynamic>& mags) {
                this->values_ = mags;
            }

            /**
             * @brief Sets the magnitude data from a vector of magnitudes and phases
             * @param mags A vector of absolute magnitudes (not dB)
             * @param phases A vector of absolute phases in rad (not degree)
             * @see set_magnitudes_dB_and_phases_deg
             */
            void set_magnitudes_and_phases(Eigen::Vector<T, Eigen::Dynamic>& mags, Eigen::Vector<T, Eigen::Dynamic>& phases) {
                assert(mags.size() == phases.size());
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
            void set_magnitudes_dB_and_phases_deg(Eigen::Vector<T, Eigen::Dynamic>& mags_dB, Eigen::Vector<T, Eigen::Dynamic>& phases_deg) {
                assert(mags_dB.size() == phases_deg.size());
                auto mags = (mags_dB.array() / T(20)).exp10();
                auto phases = phases_deg.array() * std::numbers::pi_v<T> / static_cast<T>(180);
                auto real = mags * phases.cos();
                auto imag = mags * phases.sin();
                this->values_.real() = real;
                this->values_.imag() = imag;
            }

            size_t size() const {return this->freqs_.size();}

            bool empty() const {return this->size() == 0;}
    };

    /**
     * @brief Calculates the impulse-response of frequency data.
     * 
     * Assumes the frequency-response is only the one sided analytical form containing only positive frequencies.
     * Thus it extends the frequency response to negative frequencies by mirroring it around the y axis in its conjugate complex form.
     * 
     * Uses linear piecewise first order integration in between data points
     * 
     * If the bode does not include a DC frequency sample (at 0 Hz) one will be added with the value of the smallest frequency in the data
     * 
     * @tparam T The value type used to represent numbers. Usually `float`, `double` or a custom fixpoint type
     * @param bode The frequency response data (assumed to only contain the one sided positive frequencies)
     * @param time_step 
     * @param number_of_samples 
     * @return 
     */
    template<class T>
    TimeSeries<T> to_time_series(const FrequencyResponse<T>& bode, const T time_step, const int number_of_samples){
        assert(bode.empty() == false);
        assert(number_of_samples > 0);

        const T pi = std::numbers::pi_v<T>;
        const std::complex<T> j(0, 1);

        const T start_time = 0.0;
        const T simulation_time = time_step * (number_of_samples-1);
        Eigen::Vector<T, Eigen::Dynamic> times = Eigen::Vector<T, Eigen::Dynamic>::LinSpaced(number_of_samples, start_time, simulation_time);
        Eigen::Vector<T, Eigen::Dynamic> values(number_of_samples);
        values.setZero();

        const auto f1 = bode.frequencies().head(bode.frequencies().size() - 1);
        const auto f2 = bode.frequencies().tail(bode.frequencies().size() - 1);

        Eigen::Vector<T, Eigen::Dynamic> delta_f = f2.array() - f1.array();
        T df_max = delta_f.maxCoeff();

        const auto X1 = bode.values().head(bode.values().size() - 1);
        const auto X2 = bode.values().tail(bode.values().size() - 1);

        // estimate where to switch from the small t or t=0 solution to the large t solution
        const T t_switch = 0.001 / (2 * pi * df_max);

        // complex phase of the current iteration
        Eigen::Vector<std::complex<T>, Eigen::Dynamic> e_j_2_pi_f_ti(bode.frequencies().size());
        e_j_2_pi_f_ti.setOnes();

        // complex phase turner (is also the value of the first iteration)
        Eigen::Vector<std::complex<T>, Eigen::Dynamic> e_j_2_pi_f_t1 = (j * T(2) * pi * time_step * bode.frequencies().array()).exp();

        const bool includes_dc = bode.frequencies()[0] == 0;
        const std::complex<T> dc_gain = bode.values()[0];
        const T delta_f_dc = bode.frequencies()[0];

        // special case t=0
        std::complex<T> v0_one_sided = 0;
        {
            // the one sided result for positive frequencies
            v0_one_sided = ((X1.array() + X2.array()) * delta_f.array()).sum() * T(0.5);

            // add artificial sample at 0Hz
            if(includes_dc == false){
                v0_one_sided += ((dc_gain + X1(0)) * delta_f_dc) * T(0.5);
            }

            values(0) = std::real(v0_one_sided) * 2; // assume mirrored conjugated negative frequencies
        }
        

        // linear approximation (for numerical stability)
        size_t t_itr = 1;
        std::complex<T> LinF = j * (pi/T(3)) * (delta_f.array() * (X1.array() * (T(2) * f1.array() + f2.array()) + X2.array() * (T(2) * f2.array() + f1.array()))).sum();
        
        // add artificial sample at 0Hz
        if(includes_dc == false){
            LinF += j * (pi/T(3)) * (delta_f_dc * (dc_gain * f1(0) + X1(0) * (T(2) * f1(0))));
        }

        for(; (t_itr < static_cast<size_t>(number_of_samples)) && (times(t_itr) <= t_switch); ++ t_itr){
            const T t = times(t_itr);
            std::complex<T> v = v0_one_sided + LinF * t;
            values(t_itr) = std::real(v) * 2; // assume mirrored conjugated negative frequencies (imaginary part cancels, real part adds twice)
            e_j_2_pi_f_ti.array() *= e_j_2_pi_f_t1.array();
        }

        // exact integration
        for(; t_itr < static_cast<size_t>(number_of_samples); ++ t_itr){
            const T t = times(t_itr);
            const T w = 2 * pi * t;
            
            e_j_2_pi_f_ti.array() *= e_j_2_pi_f_t1.array();

            const auto E_f1 = e_j_2_pi_f_ti.head(e_j_2_pi_f_ti.size()-1);
            const auto E_f2 = e_j_2_pi_f_ti.tail(e_j_2_pi_f_ti.size()-1);
            
            const Eigen::Vector<std::complex<T>, Eigen::Dynamic> Ia = E_f1.array() * (T(1) + j * w * delta_f.array()) - E_f2.array();
            const Eigen::Vector<std::complex<T>, Eigen::Dynamic> Ib = E_f2.array() * (T(1) - j * w * delta_f.array()) - E_f1.array();
            const Eigen::Vector<std::complex<T>, Eigen::Dynamic> I = (X1.array() * Ia.array() + X2.array() * Ib.array()) / delta_f.array();

            std::complex<T> sum = I.sum();

            // add artificial dc sample
            if(includes_dc == false){
                const std::complex<T> Ia_ = (T(1) + j * w * delta_f_dc) - E_f1(0);
                const std::complex<T> Ib_ = E_f1(0) * (T(1) - j * w * delta_f_dc) - T(1);
                const std::complex<T> I = (dc_gain * Ia_ + X1(0) * Ib_) / delta_f_dc;
                sum += I;
            }

            const std::complex<T> v = sum / (w * w);

            values(t_itr) = std::real(v) * 2; // assume mirrored conjugated negative frequencies
        }

        return TimeSeries<T>(std::move(times), std::move(values));
    }

    /**
     * \brief Calculates the time-series of a frequency response
     * 
     * Allows for arbitrary number, spacing and density-changes of samples
     * 
     * \param bode A frequency response or bode plot/measurement
     * 
     * \returns a tuple of [times, values]
     */
    template<class T>
    TimeSeries<T> to_time_series(const FrequencyResponse<T>& bode){
        const T max_freq = bode.frequencies()[bode.frequencies().size()-1];
        const T min_freq = bode.frequencies()[0];

        const T time_step = static_cast<T>(1) / (static_cast<T>(2)*max_freq);
        const T total_time = static_cast<T>(1) / min_freq;

        const int samples = static_cast<int>(total_time / time_step + static_cast<T>(0.5));

        return to_time_series(bode, time_step, samples);
    }

    // operator +
    //-----------------

    template<class T>
    FrequencyResponse<T> operator+ (const FrequencyResponse<T>& l, const FrequencyResponse<T>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies() && "Both operands need to have the same frequency vector");
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() + r.values().array();
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    FrequencyResponse<T> operator+ (const FrequencyResponse<T>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() + static_cast<std::complex<T>>(l)(r);
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    FrequencyResponse<T> operator+ (const T2& l, const FrequencyResponse<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = static_cast<std::complex<T>>(l) + r.values().array();
        return FrequencyResponse<T>(r.frequencies(), result);
    }

    template<class T>
    FrequencyResponse<T> operator+ (const FrequencyResponse<T>& b){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = +b.values().array();
        return FrequencyResponse<T>(b.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> operator+ (const FrequencyResponse<T>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> r_mags = r.eval_frequencies_Hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() + r_mags.array();
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> operator+ (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const FrequencyResponse<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> l_mags = l.eval_frequencies_Hz(r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l_mags.array() + l.values().array();
        return FrequencyResponse<T>(r.frequencies(), result);
    }

    // operator -
    //-----------------

    template<class T>
    FrequencyResponse<T> operator- (const FrequencyResponse<T>& l, const FrequencyResponse<T>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies() && "Both operands need to have the same frequency vector");
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() - r.values().array();
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    FrequencyResponse<T> operator- (const FrequencyResponse<T>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() - static_cast<std::complex<T>>(r);
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    FrequencyResponse<T> operator- (const T2& l, const FrequencyResponse<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = static_cast<std::complex<T>>(l) - r.values().array();
        return FrequencyResponse<T>(r.frequencies(), result);
    }

    template<class T>
    FrequencyResponse<T> operator- (const FrequencyResponse<T>& b){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = -b.values().array();
        return FrequencyResponse<T>(b.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> operator- (const FrequencyResponse<T>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> r_mags = r.eval_frequencies_Hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() - r_mags.array();
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> operator- (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const FrequencyResponse<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> l_mags = l.eval_frequencies_Hz(r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l_mags.array() - l.values();
        return FrequencyResponse<T>(r.frequencies(), result);
    }

    // operator *
    //-----------------

    template<class T>
    FrequencyResponse<T> operator* (const FrequencyResponse<T>& l, const FrequencyResponse<T>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies() && "Both operands need to have the same frequency vector");
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values() * r.values();
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    FrequencyResponse<T> operator* (const FrequencyResponse<T>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() * static_cast<std::complex<T>>(r);
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    FrequencyResponse<T> operator* (const T2& l, const FrequencyResponse<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = static_cast<std::complex<T>>(l) * r.values().array();
        return FrequencyResponse<T>(r.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> operator* (const FrequencyResponse<T>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> r_mags = r.eval_frequencies_Hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() * r_mags.array();
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> operator* (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const FrequencyResponse<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> l_mags = l.eval_frequencies_Hz(r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l_mags.array() * l.values().array();
        return FrequencyResponse<T>(r.frequencies(), result);
    }

    // operator /
    //-----------------

    template<class T>
    FrequencyResponse<T> operator/ (const FrequencyResponse<T>& l, const FrequencyResponse<T>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies() && "Both operands need to have the same frequency vector");
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() / r.values().array();
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    FrequencyResponse<T> operator/ (const FrequencyResponse<T>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() / std::complex<T>(r);
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    FrequencyResponse<T> operator/ (const T2& l, const FrequencyResponse<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = std::complex<T>(l) / r.values().array();
        return FrequencyResponse<T>(r.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> operator/ (const FrequencyResponse<T>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> r_mags = r.eval_frequencies_Hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() / r_mags.array();
        return FrequencyResponse<T>(l.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> operator/ (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const FrequencyResponse<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.eval_frequencies_Hz(r.frequencies()) / l.values();
        return FrequencyResponse<T>(r.frequencies(), result);
    }

    /**
     * \brief Prints a bode plot to an output stream as a `.csv` file
     * \param stream The stream to be printed to
     * \param bode The bode container with the frequencies, magnitudes and phases
     * \returns A reference to the stream object for operation chaining.
     * \see FrequencyResponse
     */
    template<class T>
    std::ostream& operator<< (std::ostream& stream, const FrequencyResponse<T>& bode){
        stream << "Frequencies (Hz), Magnitudes (dB), Phases (deg), Real, Imag" << std::endl;
        const auto frequs = bode.frequencies();
        const auto real = bode.real();
        const auto imag = bode.imag();
        const auto mags = bode.magnitudes_dB();
        const auto phases = bode.phases_deg();
        const int n = std::min({frequs.size(), mags.size(), phases.size()});
        for(int i = 0; i < n; ++i){
            stream << frequs(i) << ", " << mags(i) << ", " << phases(i) << ", " << real(i) << ", " << imag(i);
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
     * @returns The bode result as a `FrequencyResponse` struct
     * @see ContinuousTransferFunction
     * @see FrequencyResponse
     */
    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const Eigen::Vector<T, Eigen::Dynamic>& freqs_Hz
    ){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> complex_magnitudes = tf.eval_frequencies_Hz(freqs_Hz);
        const FrequencyResponse<T> result(freqs_Hz, complex_magnitudes);
        return result;
    }

    

    /**
     * \brief Calculates the bode response of a transfer function
     * \param slowest_freq_Hz The slowest/lowest frequency in Hz from which to calculate frequency responses
     * \param fastest_freq_Hz The fastest/highest frequency in Hz to which to calculate the frequency response
     * \param samples_per_decade The number of samples per decade of frequencies to be calculated
     * \returns A FrequencyResponse containing the frequencies (Hz), magnitues (dB) and phases (deg)
     * \see FrequencyResponse
     * \see ContinuousTransferFunction::eval_frequencies
     * \see bode(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, const int samples_per_decade)
     */
    template<class T, int NumOrder, int DenOrder, std::convertible_to<T> T1, std::convertible_to<T> T2>
    FrequencyResponse<T> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const T1& slowest_freq_Hz, 
            const T2& fastest_freq_Hz, 
            const int samples_per_decade=100
    ){
        const T decades = std::log10(static_cast<T>(fastest_freq_Hz)) - std::log10(static_cast<T>(slowest_freq_Hz));
        const T samples = samples_per_decade * decades;

        const Eigen::Vector<T, Eigen::Dynamic> freqs_Hz = Eigen::Vector<T, Eigen::Dynamic>::LinSpaced(samples, std::log(slowest_freq_Hz), std::log(fastest_freq_Hz)).array().exp();
        return bode<T, NumOrder, DenOrder>(tf, freqs_Hz);
    }

    /**
     * \brief Calculates the bode response of a transfer function. 
     * 
     * Infers the frequency range of the bode plot from the transfer function.
     * 
     * \param tf The transfer function
     * \param samples_per_decade The number of samples per decade of frequencies to be calculated
     * \returns A FrequencyResponse containing the frequencies (Hz), magnitues (dB) and phases (deg)
     * \see FrequencyResponse
     * \see bode(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, const T& slowest_freq_Hz, const T& fastest_freq_Hz, const int samples_per_decade)
     * \see ContinuousTransferFunction::eval_frequencies
     */
    template<class T, int NumOrder, int DenOrder>
    FrequencyResponse<T> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const int samples_per_decade=100
    ){
        const auto [slowest_freq_rad, fastest_freq_rad] = slowest_fastest_frequencies(tf);
        const T frequency_from_Hz = slowest_freq_rad / (static_cast<T>(10 * 2) * std::numbers::pi_v<T>);
        const T frequency_to_Hz = fastest_freq_rad * static_cast<T>(10 / 2) / std::numbers::pi_v<T>;
        return bode(tf, frequency_from_Hz, frequency_to_Hz, samples_per_decade);
    }

}