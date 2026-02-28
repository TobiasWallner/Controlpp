#pragma once

#include <numbers>

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

}