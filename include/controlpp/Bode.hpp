#pragma once

//std
#include <numbers>
#include <cassert>
#include <concepts>
#include <complex>
#include <algorithm>
#include <variant>

// eigen
#include <Eigen/Dense>

// csvd a csv reader
#include <csvd/csvd.hpp>

// controlpp
#include <controlpp/ContinuousTransferFunction.hpp>
#include <controlpp/TimeSeries.hpp>
#include <controlpp/algorithm.hpp>
#include <controlpp/conversion.hpp>

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
    class Bode{
        private:
            Eigen::Vector<T, Eigen::Dynamic> omegas_; // frequencies in rad
            Eigen::Vector<std::complex<T>, Eigen::Dynamic> values_; // complex magnitudes
    
        public:
            Bode() = default;

            /**
             * @brief Constructs a bode from frequencies and complex magnitudes
             * @param freqs_rad Ordered (ascending) frequency vector in rad
             * @param mags Complex Magnitudes Vector that corresponds to the frequencies
             */
            Bode(const Eigen::Vector<T, Eigen::Dynamic>& freqs_rad, const Eigen::Vector<std::complex<T>, Eigen::Dynamic>& values)
                : omegas_(freqs_rad)
                , values_(values)
            {
                // assert same size
                assert(this->omegas_.size() == this->values_.size());

                // assert only positive frequencies
                for(size_t i = 0; i < static_cast<size_t>(omegas_.size()); ++i){
                    assert(omegas_(i) >= T(0));
                }

                // assert frequencies are ascending
                for(size_t i = 1; i < static_cast<size_t>(omegas_.size()); ++i){
                    assert(omegas_(i) > omegas_(i-1));
                }
            }

            Bode(Eigen::Vector<T, Eigen::Dynamic>&& freqs_rad, Eigen::Vector<std::complex<T>, Eigen::Dynamic>&& values)
                : omegas_(std::move(freqs_rad))
                , values_(std::move(values))
            {
                // assert same size
                assert(this->omegas_.size() == this->values_.size());
                
                // assert only positive frequencies
                for(size_t i = 0; i < static_cast<size_t>(omegas_.size()); ++i){
                    assert(omegas_(i) >= T(0));
                }

                // assert frequencies are ascending
                for(size_t i = 1; i < static_cast<size_t>(omegas_.size()); ++i){
                    assert(omegas_(i) > omegas_(i-1));
                }
            }

            Bode(const Eigen::Vector<T, Eigen::Dynamic>& freqs_Hz, Eigen::Vector<std::complex<T>, Eigen::Dynamic>&& values)
                : omegas_(freqs_Hz)
                , values_(std::move(values))
            {
                // assert same size
                assert(this->omegas_.size() == this->values_.size());

                // assert only positive frequencies
                for(size_t i = 0; i < static_cast<size_t>(omegas_.size()); ++i){
                    assert(omegas_(i) >= T(0));
                }

                // assert frequencies are ascending
                for(size_t i = 1; i < static_cast<size_t>(omegas_.size()); ++i){
                    assert(omegas_(i) > omegas_(i-1));
                }
            }

            Bode(Eigen::Vector<T, Eigen::Dynamic>&& freqs_Hz, const Eigen::Vector<std::complex<T>, Eigen::Dynamic>& values)
                : omegas_(std::move(freqs_Hz))
                , values_(values)
            {
                // assert same size
                assert(this->omegas_.size() == this->values_.size());

                // assert only positive frequencies
                for(size_t i = 0; i < static_cast<size_t>(omegas_.size()); ++i){
                    assert(omegas_(i) >= T(0));
                }

                // assert frequencies are ascending
                for(size_t i = 1; i < static_cast<size_t>(omegas_.size()); ++i){
                    assert(omegas_(i) > omegas_(i-1));
                }
            }

            /**
             * \brief Returns a const-reference to the frequency vector in rad
             * \returns const-reference to an eigen vector
             */
            const Eigen::Vector<T, Eigen::Dynamic>& frequencies() const {
                return this->omegas_;
            }

            /**
             * \brief Returns a reference to the frequency vector in rad
             * \returns const-reference to an eigen vector
             */
            Eigen::Vector<T, Eigen::Dynamic>& frequencies() {
                return this->omegas_;
            }

            /**
             * \brief Returns a reference to the frequency in rad at the index n
             * \param n The index at which to return the frequency
             * \returns reference to the frequency at the n-th position
             */
            T& frequency(std::size_t n) {
                return this->omegas_(n);
            }

            /**
             * \brief Returns a reference to the frequency in rad at the index n
             * \param n The index at which to return the frequency
             * \returns reference to the frequency at the n-th position
             */
            const T& frequency(std::size_t n) const {
                return this->omegas_(n);
            }

            /**
             * \brief Returns a const-reference to the complex magnitued vector
             * \returns const-reference to an complex magnitude vector
             */
            const Eigen::Vector<std::complex<T>, Eigen::Dynamic>& values() const {
                return this->values_;
            }

            /**
             * \brief Returns a reference to the complex magnitued vector
             * \returns reference to an complex magnitude vector
             */
            Eigen::Vector<std::complex<T>, Eigen::Dynamic>& values() {
                return this->values_;
            }

            /**
             * \brief Returns a reference to the complex magnitue at the n-th position
             * \param n The n-th position to get from the magnitude vector
             * \returns reference to an complex magnitude
             */
            const std::complex<T>& value(std::size_t n) const {
                return this->values_(n);
            }

            /**
             * \brief Returns a reference to the complex magnitue at the n-th position
             * \param n The n-th position to get from the magnitude vector
             * \returns reference to an complex magnitude
             */
            std::complex<T>& value(std::size_t n) {
                return this->values_(n);
            }

            size_t size() const {return this->omegas_.size();}

            bool empty() const {return this->size() == 0;}


            /**
             * \brief returns the complex value at the given frequency using interpolation
             * 
             * If the passed frequency is not within the bode plot, data will not be extrapolated.
             * Instead the first or last available value will be returned, depending weather or not the 
             * value is below or above the spectrum.
             * 
             * Uses linear interpolation in polar coordinates
             * 
             * \param frequency The frequency (in rad/s) at which to read the complex amplitudes.
             * \returns Interpolated complex value for the passed frequency value
             */
            std::complex<T> value_at(const T& frequency) const {
                std::optional<std::pair<const T*, const T*>> result = find_enclosing(this->omegas_, frequency);
                if(result.has_value()){
                    // extract result values (iterators)
                    const T* low = result.value().first;
                    const T* high = result.value().second;

                    // calculate integral ierators
                    const size_t i_low = low - this->omegas_.data();
                    const size_t i_high = high - this->omegas_.data();

                    // get frequencies and values
                    const T f_low = *low;
                    const T f_high = *high;
                    const std::complex<T> v_low = this->values_(i_low);
                    const std::complex<T> v_high = this->values_(i_high);

                    // calculate percentage for interpolation
                    const T p = (frequency - f_low) / (f_high - f_low);

                    // linear magnitude interpolation
                    const T mag_interp = std::abs(v_low) * (static_cast<T>(1) - p) + (p) * std::abs(v_high);

                    // linear phase interpolation
                    const T phase_interp = std::arg(v_low) * (static_cast<T>(1) - p) + (p) * std::arg(v_high);

                    // calculate resulting complex value
                    const std::complex<T> result = std::polar(mag_interp, phase_interp);
                    return result;
                }else{
                    // assume same value for out of bound values | do not extrapolate
                    if(frequency <= this->omegas_(0)){
                        return this->values_(0);
                    }else{
                        return this->values_(this->values_.size()-1);
                    }
                }
            }

            /**
             * \brief Returns the phase at the passed frequency
             * 
             * Uses linear interpolation.
             * 
             * Internally uses `value_at(const T& frequency)`.
             * 
             * \param frequency The frequency (in rad/s) at which to get the phase at
             * \returns The interpolated phase (in rad) at the given frequency
             * \see value_at(const T& frequency)
             */
            T phase_at(const T& frequency) const {
                const std::complex<T> value = this->value_at(frequency);
                const T phase = std::arg(value);
                return phase;
            }

            /**
             * \brief return the phase at the given frequency in degree
             * 
             * Uses linear interpolation
             * 
             * \param frequency The frequency (in hz) at which to get the phase at
             * \returns The interpolated phase (in deg) at the given frequency
             * \see phase_at()
             * \see value_at()
             */
            T phase_deg_at(const T& frequency) const {
                const T phase_rad = this->phase_at(frequency);
                const T phase_deg = phase_rad * static_cast<T>(180) / std::numbers::pi_v<T>;
                return phase_deg;
            }

            /**
             * \brief Returns the magnitude at the passed frequency 
             * 
             * Uses linear interpolation
             * \param frequency The frequency (in hz) at which to get the magnitude at
             * \returns The interpolated magnitude (in absolute) at the given frequency
             */
            T magnitude_at(const T& frequency) const {
                const std::complex<T> value = this->value_at(frequency);
                const T mag = std::abs(value);
                return mag;
            }

            /**
             * \brief Returns the magnitude at the passed frequency 
             * 
             * Uses linear interpolation
             * \param frequency The frequency (in hz) at which to get the magnitude at
             * \returns The interpolated magnitude (in absolute) at the given frequency
             */
            T magnitude_dB_at(const T& frequency) const {
                const T mag = this->magnitude_at(frequency);
                const T mag_dB = static_cast<T>(20) * std::log10(mag);
                return mag_dB;
            }
    };

    /**
     * \brief Converts and returns the frequency vector in Hz
     * \returns const-reference to an eigen vector
     */
    template<class T>
    const Eigen::Vector<T, Eigen::Dynamic>& frequencies(const Bode<T>& bode) {
        return bode.frequencies();
    }

    /**
     * \brief Converts and returns the frequency vector in Hz
     * \returns const-reference to an eigen vector
     */
    template<class T>
    Eigen::Vector<T, Eigen::Dynamic> frequencies_hz(const Bode<T>& bode) {
        return (bode.frequencies() / 2) * std::numbers::inv_pi_v<T>;
    }

    template<class T>
    Eigen::Vector<T, Eigen::Dynamic> real(const Bode<T>& bode){
        return bode.values().real();
    }

    template<class T>
    Eigen::Vector<T, Eigen::Dynamic> imag(const Bode<T>& bode){
        return bode.values().imag();
    }

    template<class T>
    const Eigen::Vector<T, Eigen::Dynamic>& values(const Bode<T>& bode) {
        return bode.values();
    }

    /**
     * @brief Creates a vector containing the absolute magnitudes
     * 
     * Calculates the absolute magnitudes from the complex magnitudes 
     */
    template<class T>
    Eigen::Vector<T, Eigen::Dynamic> magnitudes(const Bode<T>& bode) {
        return bode.values().array().abs();
    }

    /**
     * @brief Creates a vector of magnitudes in dB
     * 
     * Calculates the magnitudes from the complex magnitudes 
     */
    template<class T>
    Eigen::Vector<T, Eigen::Dynamic> magnitudes_dB(const Bode<T>& bode) {
        return static_cast<T>(20) * bode.values().array().abs().log10();
    }

    /**
     * @brief Creates a vector of phases in rad
     * 
     * Calculates the phases from the complex magnitudes
     */
    template<class T>
    Eigen::Vector<T, Eigen::Dynamic> phases(const Bode<T>& bode) {
        Eigen::Vector<T, Eigen::Dynamic> result = bode.values().array().arg();
        return unwrap_rad(result);
    }

    /**
     * @brief Creates a vector of phases in degree
     * 
     * Calculates the phases from the complex magnitudes
     */
    template<class T>
    Eigen::Vector<T, Eigen::Dynamic> phases_deg(const Bode<T>& bode) {
        Eigen::Vector<T, Eigen::Dynamic> result = bode.values().array().arg() * static_cast<T>(180) / std::numbers::pi_v<T>;
        return unwrap_deg(result);
    }

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
     * @param time_step The timestep for the integration
     * @param simulation_time The simulation time over which to simulate. This is the minimal simulation time and may be overstepped by one timestep.
     * @return A timeseries containing the impulse frequency response of the bode data
     */
    template<class T>
    TimeSeries<T> impulse(const Bode<T>& bode, const T& time_step, const T& simulation_time){
        assert(bode.empty() == false);
        assert(simulation_time > T(0));
        assert(time_step > T(0));

        const T pi = std::numbers::pi_v<T>;
        const std::complex<T> j(0, 1);

        const T start_time = 0.0;
        const size_t number_of_samples = static_cast<size_t>(simulation_time / time_step + T(0.5));
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
     * Allows for arbitrary number, spacing and density-changes of samples.
     * 
     * Automatically estimates time-steps and the simulation time from the bode data.
     * 
     * Tip: For step responses do not do:
     * \code{.cpp}
     * stp = impulse(bode * 1/s);
     * \endcode
     * instead use the step function
     * \code{.cpp}
     * stp = step(bode * 1/s);
     * \endcode
     * 
     * \param bode A frequency response or bode plot/measurement
     * 
     * \returns a time series
     * 
     * \see template<class T> TimeSeries<T> step(const Bode<T>& bode)
     */
    template<class T>
    TimeSeries<T> impulse(const Bode<T>& bode){
        const T max_freq = bode.frequencies()[bode.frequencies().size()-1];
        const T min_freq = bode.frequencies()[0];

        const T time_step = static_cast<T>(1) / (T(2 * 4) * max_freq); // oversample 4 times
        const T total_time = static_cast<T>(1) / (T(4) * min_freq); // only use a quarter of the slowest frequency

        return impulse(bode, time_step, total_time);
    }

    /**
     * @brief Integrates the time series and writes it to out
     * 
     * Resizes `out` accordingly
     * 
     * `out` and `in` may be the same time-series object
     * 
     * @tparam T The value type
     * @param out The output time-series where the result will be written to
     * @param in The input time-series that will be integrated over
     * @param v0 The start of the integration. Integration constant C in other nomiclatures.
     */
    template<class T>
    void integrate(TimeSeries<T>& out, const TimeSeries<T>& in, const T& v0 = T(0)){
        T prev_value = in.values(0);
        T sum = 0;
        out.resize(in.size());
        out.values(0) = v0;
        for(size_t i = 1; i < in.size(); ++i){
            // first order integration
            const T sum_i = (in.values(i) + prev_value) * (in.times(i) - in.times(i-1)) * T(0.5);
            prev_value = in.values(i);
            sum += sum_i;
            out.values(i) = sum;
        }
    }

    /**
     * @brief Integrates the time series and writes it to out
     * @tparam T The value type
     * @param in The input time-series that will be integrated over
     * @param v0 The start of the integration. Integration constant C in other nomiclatures.
     * @returns The integated time-series
     */
    template<class T>
    TimeSeries<T> integrate(const TimeSeries<T>& in, const T& v0 = T(0)){
        TimeSeries<T> out;
        return integrate(out, in, v0);
    }

    /**
     * @brief Calculates the step response time-value pairs from frequency-value data
     * 
     * Internally this will calculate the impulse response first and then integrates in 
     * the time domain. This way is numerically more stable than integrating in the frequency
     * domain by `1/s` multiplication and transforming that.
     * 
     * Automatically estimates the time-step and simulation time from the bode/frequency-response data.
     * 
     * @tparam T The data type. Typically `float` or `double`.
     * @param bode The step response of the system
     * @return A TimeSeries containing time-value pairs
     */
    template<class T>
    TimeSeries<T> step(const Bode<T>& bode){
        TimeSeries<T> imp = impulse(bode);
        integrate<T>(imp, imp);
        return imp;
    }

    /**
     * @brief Calculates the step response time-value pairs from frequency-value data
     * @tparam T The data type. Typically `float` or `double`.
     * @param bode The step response of the system
     * @param time_step The time step used for the time series data
     * @param simulation_time The time until the step response will be calculated (+1 for rounding).
     * @return The TimeSeries step response of the bode data
     */
    template<class T>
    TimeSeries<T> step(const Bode<T>& bode, const T& time_step, const T& simulation_time){
        TimeSeries<T> imp = impulse(bode, time_step, simulation_time);
        integrate<T>(imp, imp);
        return imp;
    }

    // operator +
    //-----------------

    template<class T>
    Bode<T> operator+ (const Bode<T>& l, const Bode<T>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() + r.values().array();
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    Bode<T> operator+ (const Bode<T>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() + static_cast<std::complex<T>>(l)(r);
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    Bode<T> operator+ (const T2& l, const Bode<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = static_cast<std::complex<T>>(l) + r.values().array();
        return Bode<T>(r.frequencies(), result);
    }

    template<class T>
    Bode<T> operator+ (const Bode<T>& b){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = +b.values().array();
        return Bode<T>(b.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    Bode<T> operator+ (const Bode<T>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> r_mags = r.eval_frequencies_hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() + r_mags.array();
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    Bode<T> operator+ (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const Bode<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> l_mags = l.eval_frequencies_hz(r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l_mags.array() + l.values().array();
        return Bode<T>(r.frequencies(), result);
    }

    // operator -
    //-----------------

    template<class T>
    Bode<T> operator- (const Bode<T>& l, const Bode<T>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() - r.values().array();
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    Bode<T> operator- (const Bode<T>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() - static_cast<std::complex<T>>(r);
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    Bode<T> operator- (const T2& l, const Bode<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = static_cast<std::complex<T>>(l) - r.values().array();
        return Bode<T>(r.frequencies(), result);
    }

    template<class T>
    Bode<T> operator- (const Bode<T>& b){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = -b.values().array();
        return Bode<T>(b.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    Bode<T> operator- (const Bode<T>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> r_mags = r.eval_frequencies_hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() - r_mags.array();
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    Bode<T> operator- (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const Bode<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> l_mags = l.eval_frequencies_hz(r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l_mags.array() - l.values().array();
        return Bode<T>(r.frequencies(), result);
    }

    // operator *
    //-----------------

    template<class T>
    Bode<T> operator* (const Bode<T>& l, const Bode<T>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() * r.values().array();
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    Bode<T> operator* (const Bode<T>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() * static_cast<std::complex<T>>(r);
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    Bode<T> operator* (const T2& l, const Bode<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = static_cast<std::complex<T>>(l) * r.values().array();
        return Bode<T>(r.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    Bode<T> operator* (const Bode<T>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> r_mags = r.eval_frequencies_hz(l.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() * r_mags.array();
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    Bode<T> operator* (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const Bode<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> l_mags = l.eval_frequencies_hz(r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l_mags.array() * l.values().array();
        return Bode<T>(r.frequencies(), result);
    }

    // operator /
    //-----------------

    /**
     * @brief 
     * 
     * TODO: make it also work for bode that have different frequency vectors
     * 
     * @tparam T 
     * @param l 
     * @param r 
     * @return 
     */
    template<class T>
    Bode<T> operator/ (const Bode<T>& l, const Bode<T>& r){
        assert(l.size() == r.size());
        assert(l.frequencies() == r.frequencies());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() / r.values().array();
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    Bode<T> operator/ (const Bode<T>& l, const T2& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() / std::complex<T>(r);
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, std::convertible_to<T> T2>
    Bode<T> operator/ (const T2& l, const Bode<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = std::complex<T>(l) / r.values().array();
        return Bode<T>(r.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    Bode<T> operator/ (const Bode<T>& l, const ContinuousTransferFunction<T, NumOrder, DenOrder>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> r_mags = r.eval_frequencies_hz(l.frequencies().array());
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.values().array() / r_mags.array();
        return Bode<T>(l.frequencies(), result);
    }

    template<class T, int NumOrder, int DenOrder>
    Bode<T> operator/ (const ContinuousTransferFunction<T, NumOrder, DenOrder>& l, const Bode<T>& r){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> result = l.eval_frequencies_hz(r.frequencies().array()) / l.values().array();
        return Bode<T>(r.frequencies(), result);
    }

    /**
     * @brief Calculates the bode response for a pre defined frequency (rad/s) vector
     * @tparam T The value type of the transfer function
     * @tparam NumOrder The numerator order
     * @tparam DenOrder The denominator order
     * @tparam NSize The number of elements in the frequency vector (may also be `Eigen::Dynamic`)
     * @param tf The continuous transfer function to analyse
     * @param freqs The frequency vector at which to evaluate the transfer function in rad/s
     * @returns The bode result as a `Bode` struct
     * @see ContinuousTransferFunction
     * @see Bode
     */
    template<class T, int NumOrder, int DenOrder>
    Bode<T> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const Eigen::Vector<T, Eigen::Dynamic>& freqs
    ){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> complex_magnitudes = tf.eval_frequencies(freqs);
        const Bode<T> result(freqs, complex_magnitudes);
        return result;
    }

    /**
     * @overload 
     */
    template<class T, int NumOrder, int DenOrder>
    Bode<T> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            Eigen::Vector<T, Eigen::Dynamic>&& freqs
    ){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> complex_magnitudes = tf.eval_frequencies(freqs);
        const Bode<T> result(std::move(freqs), complex_magnitudes);
        return result;
    }
    
    /**
     * @brief Calculates the bode response for a pre defined frequency (Hz) vector
     * @tparam T The value type of the transfer function
     * @tparam NumOrder The numerator order
     * @tparam DenOrder The denominator order
     * @tparam NSize The number of elements in the frequency vector (may also be `Eigen::Dynamic`)
     * @param tf The continuous transfer function to analyse
     * @param freqs_Hz The frequency vector at which to evaluate the transfer function in hz
     * @returns The bode result as a `Bode` struct
     * @see ContinuousTransferFunction
     * @see Bode
     */
    template<class T, int NumOrder, int DenOrder>
    Bode<T> bode_hz(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const Eigen::Vector<T, Eigen::Dynamic>& freqs_Hz
    ){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> complex_magnitudes = tf.eval_frequencies_hz(freqs_Hz);
        Eigen::Vector<T, Eigen::Dynamic> freqs_rad = freqs_Hz.array() * 2 * std::numbers::pi_v<T>;
        const Bode<T> result(std::move(freqs_rad), complex_magnitudes);
        return result;
    }

    /**
     * @overload 
     */
    template<class T, int NumOrder, int DenOrder>
    Bode<T> bode_hz(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            Eigen::Vector<T, Eigen::Dynamic>&& freqs_Hz
    ){
        const Eigen::Vector<std::complex<T>, Eigen::Dynamic> complex_magnitudes = tf.eval_frequencies_hz(freqs_Hz);
        freqs_Hz *= 2 * std::numbers::pi_v<T>;
        const Bode<T> result(std::move(freqs_Hz), complex_magnitudes);
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
    Bode<T> bode_hz(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const T1& slowest_freq_Hz, 
            const T2& fastest_freq_Hz, 
            const int samples_per_decade=100
    ){
        const T decades = std::log10(static_cast<T>(fastest_freq_Hz)) - std::log10(static_cast<T>(slowest_freq_Hz));
        const T samples = samples_per_decade * decades;
        const Eigen::Vector<T, Eigen::Dynamic> freqs_Hz = Eigen::Vector<T, Eigen::Dynamic>::LinSpaced(samples, std::log(slowest_freq_Hz), std::log(fastest_freq_Hz)).array().exp();
        return bode_hz<T, NumOrder, DenOrder>(tf, freqs_Hz);
    }

    /**
     * \brief Calculates the bode response of a transfer function
     * \param slowest_freq_rad The slowest/lowest frequency in rad/s from which to calculate frequency responses
     * \param fastest_freq_rad The fastest/highest frequency in rad/s to which to calculate the frequency response
     * \param samples_per_decade The number of samples per decade of frequencies to be calculated
     * \returns A Bode containing the frequencies (Hz), magnitues (dB) and phases (deg)
     * \see Bode
     * \see ContinuousTransferFunction::eval_frequencies
     * \see bode(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, const int samples_per_decade)
     */
    template<class T, int NumOrder, int DenOrder, std::convertible_to<T> T1, std::convertible_to<T> T2>
    Bode<T> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const T1& slowest_freq_rad, 
            const T2& fastest_freq_rad, 
            const int samples_per_decade=100
    ){
        const T decades = std::log10(static_cast<T>(fastest_freq_rad)) - std::log10(static_cast<T>(slowest_freq_rad));
        const T samples = samples_per_decade * decades;

        const Eigen::Vector<T, Eigen::Dynamic> freqs_rad = Eigen::Vector<T, Eigen::Dynamic>::LinSpaced(samples, std::log(slowest_freq_rad), std::log(fastest_freq_rad)).array().exp();
        return bode<T, NumOrder, DenOrder>(tf, freqs_rad);
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
    Bode<T> bode_hz(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const int samples_per_decade=100
    ){
        const auto [slowest_freq_rad, fastest_freq_rad] = slowest_fastest_frequencies(tf);
        const T frequency_from_Hz = slowest_freq_rad / (static_cast<T>(10 * 2) * std::numbers::pi_v<T>);
        const T frequency_to_Hz = fastest_freq_rad * static_cast<T>(10 / 2) / std::numbers::pi_v<T>;
        return bode_hz(tf, frequency_from_Hz, frequency_to_Hz, samples_per_decade);
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
    Bode<T> bode(
            const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf, 
            const int samples_per_decade=100
    ){
        const auto [slowest_freq_rad, fastest_freq_rad] = slowest_fastest_frequencies(tf);
        const T frequency_from_rad = slowest_freq_rad / static_cast<T>(10);
        const T frequency_to_rad = fastest_freq_rad * static_cast<T>(10);
        return bode(tf, frequency_from_rad, frequency_to_rad, samples_per_decade);
    }

    /**
     * \brief Prints a bode plot to an output stream as a `.csv` file
     * \param stream The stream to be printed to
     * \param bode The bode container with the frequencies, magnitudes and phases
     * \returns A reference to the stream object for operation chaining.
     * \see Bode
     */
    template<class T>
    void write_csv (std::ostream& stream, const Bode<T>& bode){
        stream << "Frequencies (Hz), Magnitudes (dB), Phases (deg), Real, Imag" << std::endl;

        const int n = bode.size();
        for(int i = 0; i < n; ++i){
            const T omega = bode.frequency(i);
            const T f = to_hz(omega);

            const std::complex<T> value = bode.value(i);
            const T mag = to_dB(std::abs(value));
            const T phase = to_deg(std::arg(value));
            const T real = std::real(value);
            const T imag = std::imag(value);

            stream << f << ", " << mag << ", " << phase << ", " << real << ", " << imag;
            if(i < n-1) stream << "\n";
        }
    }

    /**
     * \brief Error cases for reading/parsing bode data from CSV formated data
     */
    enum class EBodeCsvReadError{
        CouldNotFindFrequencyVector, ///< Could not find the frequency axis in the csv data. Searched for a header that starts with `"f"` (case-insensitive).
        CouldNotFindAmplitudeVectors, ///< Could not find the amplitude data. Needed either: 1) Two columns that start with `"re"` and `"im"` (case-insensitive). 2) Two columns that start with `"mag"` and `"ph"` (case-insensitive).
    };

    std::ostream& operator<<(std::ostream& stream, EBodeCsvReadError val);

    /**
     * \brief Enum that determines how frequency data will be interpreted when reading CSV data.
     */
    enum class EFrequencyInterpretation{
        AutoHz,     ///< Automatically infers the unit from the header or interprets the data in **Hz** if no unit is found in the header name
        AutoRad,    ///< Automatically infers the unit from the header or interprets the data in **rad** if no unit is found in the header name
        ForceHz,    ///< When reading frequency data the data is always interpreted in **Hz** regardless of the header name
        ForceRad,   ///< When reading frequency data the data is always interpreted in **rad** regardless of the header name
    };

    /**
     * \brief Enum that determines how frequency data will be interpreted when reading CSV data.
     */
    enum class EMagnitudeInterpretation{
        Auto,    ///< Automatically infers the unit from the header and interprets as **deci-Bell** if `"dB"` has been found or interprets the data as **absolute values** if no unit is found in the header name
        ForceAbs,   ///< Forces the magnitude data to be interpreted in **absolute values** regardless of the header name
        ForceDB,    ///< Forces the magnitude data to be interpreted in **dB** regardless of the header name
    };

    /**
     * \brief Enum that determines how phase data will be interpreted when reading CSV data.
     */
    enum class EPhaseInterpretation{
        AutoRad,    ///< Automatically infers the unit from the header or interprets the data as **rad** if no unit is found in the header name
        AutoDeg,    ///< Automatically infers the unit from the header or interprets the data as **deg** if no unit is found in the header name
        ForceRad,   ///< Forces the magnitude data to be interpreted in **rad** regardless of the header name
        ForceDeg    ///< Forces the magnitude data to be interpreted in **deg** regardless of the header name
    };

    /**
     * @brief Loads bode data from csv data
     * 
     * TODO
     * 
     * Automatic data discovery:
     * ---
     * 
     * First: searches for the frequency data range. Then for the real and imaginary data ranges. 
     * If real and imaginary data could not be identified magnitudes and phases will be searched for.
     * 
     * All following textual comparisons are case insensitive.
     * 
     * - Assumes colums starting with `"f"` to be the frequency data.
     *   - If the frequency data name contains `"hz"` the frequency is interpreted as Hertz.
     *   - If the frequency data name contains `"rad"` the frequency is interpreted as Radiant.
     * - Assumes columns starting with `"re"` or `"im"` to be the real and imaginary amplitudes.
     * - If complex amplitudes could not be found assumes that columns starting with `"mag"` or `"ph"` to be absolute magnitudes and phase
     *   - If the magnitude column header contains `"dB"` the magnitude data will be interpreted in deci-Bell.
     *   - If the phase column header contains `"deg"` the phase data is interpreted in degree.
     *   - If the phase column header contains `"rad"` the phase data is interpreted in radiants.
     * 
     * @see EFrequencyInterpretation
     * @see EMagnitudeInterpretation
     * @see EPhaseInterpretation
     * 
     * @param stream 
     * @param settings 
     * @returns If successful: A bode structure. On failure/error: a variant with an error of the csv reader (`csvd::ReadError`) or an error
     * of assembing a bode from the csv `EBodeCsvReadError`.
     */
    tl::expected<
        Bode<double>, // success
        std::variant<EBodeCsvReadError, csvd::ReadError>> //error
    read_bode_from_csv(
        std::istream& stream, 
        const csvd::Settings& csv_settings = csvd::Settings(),
        EFrequencyInterpretation freq_interp = EFrequencyInterpretation::AutoHz,
        EMagnitudeInterpretation mag_interp = EMagnitudeInterpretation::Auto,
        EPhaseInterpretation phase_interp = EPhaseInterpretation::AutoDeg
    );

}