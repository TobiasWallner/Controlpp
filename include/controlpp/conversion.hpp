#pragma once

#include <numbers>
#include <cmath>
#include <Eigen/Dense>

namespace controlpp{

    /**
     * @brief Converts a number from radiants per second to herz
     * @tparam T The data type
     * @param radps The input value in radiants per second
     * @returns A value in herz
     */
    template<class T>
    T to_hz(const T& radps){
        constexpr T factor = std::numbers::inv_pi_v<T> / 2;
        return radps * factor;
    }

    template<class T, int N>
    Eigen::Vector<T, N> to_hz(const Eigen::Vector<T, N>& radps){
        constexpr T factor = std::numbers::inv_pi_v<T> / 2;
        return radps * factor;
    }

    /**
     * @brief Converts a number from herz to radiants per second
     * @tparam T The data type
     * @param radps The input value in herz
     * @returns A value in herz
     */
    template<class T>
    T to_radps(const T& hz){
        constexpr T factor = 2 * std::numbers::pi_v<T>;
        return hz * factor;
    }

    template<class T, int N>
    Eigen::Vector<T, N> to_radps(const Eigen::Vector<T, N>& hz){
        constexpr T factor = 2 * std::numbers::pi_v<T>;
        return hz * factor;
    }

    template<class T>
    T to_deg(const T& rad){
        constexpr T factor = static_cast<T>(180) / std::numbers::pi_v<T>;
        return rad * factor;
    }

    template<class T, int N>
    Eigen::Vector<T, N> to_deg(const Eigen::Vector<T, N>& rad){
        constexpr T factor = static_cast<T>(180) / std::numbers::pi_v<T>;
        return rad * factor;
    }

    template<class T>
    T to_rad(const T& deg){
        constexpr T factor = std::numbers::pi_v<T> / static_cast<T>(180);
        return deg * factor;
    }

    template<class T, int N>
    Eigen::Vector<T, N> to_rad(const Eigen::Vector<T, N>& deg){
        constexpr T factor = std::numbers::pi_v<T> / static_cast<T>(180);
        return deg * factor;
    }

    template<class T>
    T to_dB(const T& value){
        return 20 * std::log10(value);
    }

    template<class T>
    T from_dB(const T& value){
        return std::pow(static_cast<T>(10), value / static_cast<T>(20));
    }

    /**
     * @brief Prewarps frequencies for the tustin transformation
     * 
     * The tustin transformation stretches frequencies with the scaling:
     * \f[
     * \omega_d = \frac{2}{T_s} \text{atan}\left( \frac{\omega_c T_s}{2} \right)
     * \f]
     * 
     * where:
     * - \f$ \omega_c \f$ is the frequency set in continuous time
     * - \f$ \omega_d \f$ is where the frequency will be placed at after discretisation (tustin)
     * - \f$ T_s \f$ is the sample time of the discretisation
     * 
     * So for example if you place a continuous time notch filter at 300 Hz and then discretise it
     * with the tustin transformation at a sample frequency of Fs=1kHz the discretised notch will land at 240.6 Hz.
     * 
     * To have the frequencies after the discretisation exactly at the given frequency we 'pre-warp' them with the
     * inverse function which this function provides:
     * 
     * \f[
     * \omega_\text{pre} = \frac{2}{T_s} \text{tan}\left( \frac{\omega_\text{target} T_s}{2} \right)
     * \f]
     * 
     * Note that this function returns rad/s and not rad/sample
     * 
     * @tparam T The value type
     * @param omega The frequency before the pre-warping (in rad/s)
     * @param Ts The sample time that is also used for the tustin discretisation
     * @returns The pre-warped/pre-scaled frequency (in rad/s) so that after applying the tustin transformation the set frequency is exactly at `omega`.
     * 
     * @see prewarp_tustin(const T& omega, const T& Ts)
     * @see prewarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts)
     * @see unwarp_tustin(const T& omega, const T& Ts)
     * @see unwarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts)
     */
    template<class T>
    T prewarp_tustin(const T& omega, const T& Ts){
        const T result = (static_cast<T>(2) / Ts) * std::tan(omega * Ts / static_cast<T>(2));
        return result;
    }

    /**
     * @overload
     * 
     * @see prewarp_tustin(const T& omega, const T& Ts)
     * @see prewarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts)
     * @see unwarp_tustin(const T& omega, const T& Ts)
     * @see unwarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts)
     */
    template<class T, int N>
    Eigen::Vector<T, N> prewarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts){
        const Eigen::Vector<T, N> result = ((static_cast<T>(2) / Ts) * (omegas.array() * Ts / static_cast<T>(2)).tan()).matrix();
        return result;
    }

    /**
     * @brief 
     * 
     * provides the inverse to `prewarp_tustin`:
     * 
     * \f[
     * \omega_d = \frac{2}{T_s} \text{atan}\left( \frac{\omega_c T_s}{2} \right)
     * \f]
     * 
     * Calculates where the tustin transformation would actually place the frequency
     * after discretisation.
     * 
     * @tparam T The value type
     * @param omega The input frequency
     * @param Ts The sample time
     * @return The frequency where the tustin transformation would actually place the frequency after discretisation
     * 
     * @see prewarp_tustin(const T& omega, const T& Ts)
     * @see prewarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts)
     * @see unwarp_tustin(const T& omega, const T& Ts)
     * @see unwarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts)
     */
    template<class T>
    T unwarp_tustin(const T& omega, const T& Ts){
        const T result = (static_cast<T>(2) / Ts) * std::atan(omega * Ts / static_cast<T>(2));
        return result;
    }

    /**
     * @overload
     * 
     * @see prewarp_tustin(const T& omega, const T& Ts)
     * @see prewarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts)
     * @see unwarp_tustin(const T& omega, const T& Ts)
     * @see unwarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts)
     */
    template<class T, int N>
    Eigen::Vector<T, N> unwarp_tustin(const Eigen::Vector<T, N>& omegas, const T& Ts){
        const Eigen::Vector<T, N> result = ((static_cast<T>(2) / Ts) * (omegas.array() * Ts / static_cast<T>(2)).atan()).matrix();
        return result;
    }

}