#pragma once

#include "Polynom.hpp"
#include "TransferFunction.hpp"

namespace controlpp
{

    /**
     * \brief Continuous transfer functions in the s lapace plain
     * 
     * Transfer functions of the shape:
     * 
     * \f[
     * G(s) = \frac{b_0 + b_1 s \cdots b_m s^m}{a_0 + a_1 s \cdots a_n s^n}
     * \f]
     * 
     * \tparam T The value type that the transfer function should use
     * \tparam NumberOrder The order of the numerator. m in the equation.
     * \tparam DenOrder The order of the denominator. n in the equation.
     */
    template<class T, int NumOrder, int DenOrder>
    class ContinuousTransferFunction{
        public:
            using value_type = T;
            using transfer_function_type = TransferFunction<T, NumOrder, DenOrder>;
            using num_type = typename transfer_function_type::num_type;
            using den_type = typename transfer_function_type::den_type;
            using num_vector_type = typename transfer_function_type::num_vector_type;
            using den_vector_type = typename transfer_function_type::den_vector_type;

        private:
            TransferFunction<T, NumOrder, DenOrder> tf_;
        public:

            constexpr ContinuousTransferFunction() = default;
            constexpr ContinuousTransferFunction(const ContinuousTransferFunction&) = default;
            constexpr ContinuousTransferFunction& operator=(const ContinuousTransferFunction&) = default;

            constexpr explicit ContinuousTransferFunction(const Polynom<T, NumOrder>& num, const Polynom<T, DenOrder>& den)
                : tf_(num, den){}

            constexpr explicit ContinuousTransferFunction(const TransferFunction<T, NumOrder, DenOrder>& transfer_function)
                : tf_(transfer_function){}

            constexpr explicit ContinuousTransferFunction(const Eigen::Vector<T, NumOrder+1>& num, const Eigen::Vector<T, DenOrder+1>& den)
                : tf_(num, den){}

            constexpr explicit ContinuousTransferFunction(const T(&num)[NumOrder+1], const T(&den)[DenOrder+1])
                : tf_(num, den){}

            constexpr num_type& num() {return this->tf_.num();}
            constexpr const num_type& num() const {return this->tf_.num();}

            constexpr T& num(size_t i) {return this->tf_.num(i);}
            constexpr const T& num(size_t i) const {return this->tf_.num(i);}

            constexpr den_type& den() {return this->tf_.den();}
            constexpr const den_type& den() const {return this->tf_.den();}

            constexpr T& den(size_t i) {return this->tf_.den(i);}
            constexpr const T& den(size_t i) const {return this->tf_.den(i);}

            constexpr transfer_function_type& transfer_function() {return this->tf_;}
            constexpr const transfer_function_type& transfer_function() const {return this->tf_;}

            /**
             * \brief Evaluates the rational polynomial at `x`
             * \param x The variable used to evaluate the polynomial
             * \returns The result of the polynomial
             */
            constexpr T eval(const T& x) const {
                return this->tf_.eval(x);
            }

            /**
             * \brief Evaluates the rational polynomial elementwise at every element of the input vector
             * \param x_vec The parameter vector at which the rational polynomial is being evaluated at
             * \returns An eigen vector with the results of the rational polynomial
             */
            template<int M>
            constexpr Eigen::Vector<T, M> eval(const Eigen::Vector<T, M>& x_vec) const {
                return this->tf_.eval(x_vec);
            }

            /**
             * \brief Evaluates the rational polynomial at a complex `x`
             * \param x The complex variable used to evaluate the polynomial
             * \returns The result of the polynomial as a complex value
             */
            constexpr std::complex<T> eval(const std::complex<T>& x) const {
                return this->tf_.eval(x);
            }

            /**
             * \brief Evaluates the rational polynomial elementwise at every complex element of the input vector
             * \param x_vec The complex parameter vector at which the rational polynomial is being evaluated at
             * \returns An eigen vector with the complex results of the rational polynomial
             */
            template<int M>
            constexpr Eigen::Vector<std::complex<T>, M> eval(const Eigen::Vector<std::complex<T>, M>& x_vec) const {
                return this->tf_.eval(x_vec);
            }

            /**
             * \brief Evaluates the transfer function at the given frequency
             * 
             * equivalent to calling `.eval(std::complex<T>(0, f))`.
             * 
             * \param frequency The frequency (in radiants per second) at which to evaluate the transfer function at
             * \returns The complex result of the frequency evaluation/analysis.
             */
            constexpr std::complex<T> eval_frequency(const T& frequency) const {
                return this->tf_.eval_frequency(frequency);
            }

            constexpr std::complex<T> eval_frequency_Hz(const T& frequency) const {
                return this->tf_.eval_frequency_Hz(frequency);
            }

            /**
             * \brief Evaluates the transfer function (rad/s) at the given frequencies
             * \param frequencies The frequencies (rad/s) at which to evaluate the transfer function at
             * \returns The complex result of the frequency evaluation/analysis.
             */
            template<int M>
            constexpr Eigen::Vector<std::complex<T>, M> eval_frequencies(const Eigen::Vector<T, M>& frequencies) const {
                return this->tf_.eval_frequencies(frequencies);
            }

            /**
             * \brief Evaluates the transfer function (Hz) at the given frequencies
             * \param frequencies The frequencies (Hz) at which to evaluate the transfer function at
             * \returns The complex result of the frequency evaluation/analysis.
             */
            template<int M>
            constexpr Eigen::Vector<std::complex<T>, M> eval_frequencies_Hz(const Eigen::Vector<T, M>& frequencies) const {
                return this->tf_.eval_frequencies_Hz(frequencies);
            }

            friend std::ostream& operator<<(std::ostream& stream, const ContinuousTransferFunction& ctf){
                ctf.transfer_function().print(stream, "s");
                return stream;
            }
    };

    // operator +
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator+(const ContinuousTransferFunction<T, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.transfer_function() + rhs.transfer_function());
    }

    template<class Tpoly, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const Tscalar& lhs, const ContinuousTransferFunction<Tpoly, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs + rhs.transfer_function());
    }

    template<class Tpoly, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const ContinuousTransferFunction<Tpoly, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.transfer_function() + rhs);
    }

    // operator -
    // -----------

    template<class T, int NumOrder, int DenOrder>
    constexpr ContinuousTransferFunction<T, NumOrder, DenOrder> operator-(const ContinuousTransferFunction<T, NumOrder, DenOrder>& a){
        return ContinuousTransferFunction(-a.transfer_function());
    }

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator-(const ContinuousTransferFunction<T, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.transfer_function() - rhs.transfer_function());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const Tscalar& lhs, const ContinuousTransferFunction<T, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs - rhs.transfer_function());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const ContinuousTransferFunction<T, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.transfer_function() - rhs);
    }

    // operator *
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator*(const ContinuousTransferFunction<T, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.transfer_function() * rhs.transfer_function());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const Tscalar& lhs, const ContinuousTransferFunction<T, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs * rhs.transfer_function());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const ContinuousTransferFunction<T, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.transfer_function() * rhs);
    }

    // operator /
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator/(const ContinuousTransferFunction<T, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.transfer_function() / rhs.transfer_function());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const Tscalar& lhs, const ContinuousTransferFunction<T, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs / rhs.transfer_function());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const ContinuousTransferFunction<T, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        ContinuousTransferFunction result(lhs.transfer_function() / rhs);
        return result;
    }

    namespace tf{
        template<class T=double>
        static inline const ContinuousTransferFunction<T, 1, 0> s({T(0), T(1)}, {T(1)});
    }

// analysis

    template<class T, int NumOrder, int DenOrder>
    Eigen::Vector<std::complex<T>, NumOrder> zeros(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf){
        const Eigen::Vector<std::complex<T>, NumOrder> result = zeros(tf.num());
        return result;
    }   

    template<class T, int NumOrder, int DenOrder>
    Eigen::Vector<std::complex<T>, DenOrder> poles(const ContinuousTransferFunction<T, NumOrder, DenOrder>& tf){
        const Eigen::Vector<std::complex<T>, DenOrder> result = zeros(tf.den());
        return result;
    }

} // namespace controlpp