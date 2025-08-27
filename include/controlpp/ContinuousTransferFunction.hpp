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
            using ratpoly_type = TransferFunction<T, NumOrder, DenOrder>;
            using num_type = typename ratpoly_type::num_type;
            using den_type = typename ratpoly_type::den_type;
            using num_vector_type = typename ratpoly_type::num_vector_type;
            using den_vector_type = typename ratpoly_type::den_vector_type;

        private:
            TransferFunction<T, NumOrder, DenOrder> tf_;
        public:

            constexpr ContinuousTransferFunction() = default;
            constexpr ContinuousTransferFunction(const ContinuousTransferFunction&) = default;
            constexpr ContinuousTransferFunction& operator=(const ContinuousTransferFunction&) = default;

            constexpr explicit ContinuousTransferFunction(const Polynom<T, NumOrder>& num, const Polynom<T, DenOrder>& den)
                : tf_(num, den){}

            constexpr explicit ContinuousTransferFunction(const TransferFunction<T, NumOrder, DenOrder>& ratpoly)
                : tf_(ratpoly){}

            constexpr explicit ContinuousTransferFunction(const num_vector_type& num, const den_vector_type& den)
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

            constexpr ratpoly_type& ratpoly() {return this->tf_;}
            constexpr const ratpoly_type& ratpoly() const {return this->tf_;}

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
            constexpr std::complex<T> eval_frequency(const T& frequency){
                return this->tf_.eval_frequency(frequency);
            }

            /**
             * \brief Evaluates the transfer function at the given frequencies
             * 
             * \param frequencies The frequencies (in radiants per second) at which to evaluate the transfer function at
             * \returns The complex result of the frequency evaluation/analysis.
             */
            template<int M>
            constexpr Eigen::Vector<std::complex<T>, M> eval_frequencies(const Eigen::Vector<T, M>& frequencies){
                return this->tf_.eval_frequencies(frequencies);
            }

            friend std::ostream& operator<<(std::ostream& stream, const ContinuousTransferFunction& ctf){
                ctf.ratpoly().print(stream, "s");
                return stream;
            }
    };

    // operator +
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator+(const ContinuousTransferFunction<T, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const Tscalar& lhs, const ContinuousTransferFunction<Tpoly, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const ContinuousTransferFunction<Tpoly, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() + rhs);
    }

    // operator -
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator-(const ContinuousTransferFunction<T, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() - rhs.ratpoly());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const Tscalar& lhs, const ContinuousTransferFunction<T, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs - rhs.ratpoly());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const ContinuousTransferFunction<T, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() - rhs);
    }

    // operator *
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator*(const ContinuousTransferFunction<T, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() * rhs.ratpoly());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const Tscalar& lhs, const ContinuousTransferFunction<T, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs * rhs.ratpoly());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const ContinuousTransferFunction<T, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() * rhs);
    }

    // operator /
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator/(const ContinuousTransferFunction<T, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() / rhs.ratpoly());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const Tscalar& lhs, const ContinuousTransferFunction<T, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs / rhs.ratpoly());
    }

    template<class T, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const ContinuousTransferFunction<T, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        ContinuousTransferFunction result(lhs.ratpoly() / rhs);
        return result;
    }

    namespace tf{
        template<class T=double>
        static inline const ContinuousTransferFunction<T, 1, 0> s({T(0), T(1)}, {T(1)});
    }

} // namespace controlpp