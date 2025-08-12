#pragma once

// std
#include <initializer_list>

#include "Polynom.hpp"

namespace controlpp
{
    template<class T, size_t NumSize, size_t DenSize>
    class TransferFunction{
        public:
            using value_type = T;
            using num_type = Polynom<T, NumSize>;
            using den_type = Polynom<T, DenSize>;
            using num_vector_type = typename num_type::vector_type;
            using den_vector_type = typename den_type::vector_type;

        private:
            num_type _num;
            den_type _den;

        public:
            constexpr TransferFunction() = default;
            constexpr TransferFunction(const TransferFunction&) = default;
            constexpr TransferFunction& operator=(const TransferFunction&) = default;

            constexpr TransferFunction(const Polynom<T, NumSize>& num, const Polynom<T, DenSize>& den)
                : _num(num)
                , _den(den){}

            constexpr TransferFunction(const num_vector_type& num, const den_vector_type& den)
                : _num(num)
                , _den(den){}

            constexpr TransferFunction(const T(&num)[NumSize], const T(&den)[DenSize])
                : _num(num)
                , _den(den){}

            constexpr num_type& num() {return this->_num;}
            constexpr const num_type& num() const {return this->_num;}

            constexpr T& num(size_t i) {return this->_num.at(i);}
            constexpr const T& num(size_t i) const {return this->_num.at(i);}

            constexpr den_type& den() {return this->_den;}
            constexpr const den_type& den() const {return this->_den;}

            constexpr T& den(size_t i) {return this->_den.at(i);}
            constexpr const T& den(size_t i) const {return this->_den.at(i);}

            void print(std::ostream& stream, std::string_view var="x") const {
                stream << "num: "; this->num().prsize_t(stream, var);
                stream << "\nden: "; this->den().prsize_t(stream, var); stream << '\n';
            }

            friend std::ostream& operator<<(std::ostream& stream, const TransferFunction& rpoly){
                rpoly.prsize_t(stream);
                return stream;
            }
    };
    // -----------------------------------------------------------------------------------------------
    //                                      Comparison Operators
    // -----------------------------------------------------------------------------------------------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator==(const TransferFunction<T, NumSize1, DenSize1>& lhs, const TransferFunction<T, NumSize1, DenSize1>& rhs){
        return ((lhs.num() == rhs.num()) && (lhs.den() == rhs.den()));
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr bool operator!=(const TransferFunction<T, NumSize1, DenSize1>& lhs, const TransferFunction<T, NumSize1, DenSize1>& rhs){
        return ((lhs.num() != rhs.num()) || (lhs.den() != rhs.den()));
    }
    
    // -----------------------------------------------------------------------------------------------
    //                                      Arithmetic Operators
    // -----------------------------------------------------------------------------------------------

    // operator +
    // -----------

    /**
     * \brief rational addition of polynomials
     * 
     * calculates:
     * \f[
     *      \frac{num_1(x)}{den_1(x)} \frac{num_2(x)}{den_2(x)}
     * \f]
     * 
     * \tparam T The value type of the polynomial parameters
     * \tparam NumSize1 The size of the numerator of the left-hand-side addition argument
     * \tparam DenSize1 The size of the denominator of the left-hand-side addition argument
     * \tparam NumSize2 The size of the numberator of the right-hand-side addition argument
     * \tparam DenSize1 The size of the denominator of the right-hand-side addition argument
     * 
     * \param lhs The left-hand-side additino argument as a rational polynom
     * \param rhs The right-hand-side addition argument as a rational polynom
     */
    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator+(const TransferFunction<T, NumSize1, DenSize1>& lhs, const TransferFunction<T, NumSize2, DenSize2>& rhs){
        return TransferFunction(lhs.num() * rhs.den() + rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const Tscalar& lhs, const TransferFunction<T, NumSize, DenSize>& rhs){
        return TransferFunction(static_cast<T>(lhs) * rhs.den() + rhs.num(), rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const TransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return TransferFunction(lhs.num() + lhs.den() * static_cast<T>(rhs), lhs.den());
    }

    // operator -
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1>
    constexpr TransferFunction<T, NumSize1, DenSize1> operator-(const TransferFunction<T, NumSize1, DenSize1>& poly){
        return TransferFunction<T, NumSize1, DenSize1>(-poly.num(), poly.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator-(const TransferFunction<T, NumSize1, DenSize1>& lhs, const TransferFunction<T, NumSize2, DenSize2>& rhs){
        return TransferFunction(lhs.num() * rhs.den() - rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    
    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const Tscalar& lhs, const TransferFunction<T, NumSize, DenSize>& rhs){
        return TransferFunction(static_cast<T>(lhs) * rhs.den() - rhs.num(), rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const TransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return TransferFunction(lhs.num() - lhs.den() * static_cast<T>(rhs), lhs.den());
    }

    // operator *
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator*(const TransferFunction<T, NumSize1, DenSize1>& lhs, const TransferFunction<T, NumSize2, DenSize2>& rhs){
        return TransferFunction(lhs.num() * rhs.num(), lhs.den() * rhs.den());
    }

    
    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const Tscalar& lhs, const TransferFunction<T, NumSize, DenSize>& rhs){
        return TransferFunction(static_cast<T>(lhs) * rhs.num(), rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const TransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return TransferFunction(lhs.num() * static_cast<T>(rhs), lhs.den());
    }

    // operator /
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator/(const TransferFunction<T, NumSize1, DenSize1>& lhs, const TransferFunction<T, NumSize2, DenSize2>& rhs){
        return TransferFunction(lhs.num() * rhs.den(), lhs.den() * rhs.num());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const Tscalar& lhs, const TransferFunction<T, NumSize, DenSize>& rhs){
        return TransferFunction(static_cast<T>(lhs) * rhs.den(), rhs.num());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const TransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return TransferFunction(lhs.num(), lhs.den() * static_cast<T>(rhs));
    }



    ///--------------------------------------------------------------------------------------------------------------------
    ///--------------------------------------------------------------------------------------------------------------------
    ///--------------------------------------------------------------------------------------------------------------------

    /**
     * \brief Fixed sized polynomial
     * 
     * This rational polynomial does not automatically grow as a result of operations like its sibling: `controlpp::TransferFunction`
     * 
     * \see controlpp::TransferFunction
     */
    template<class T, size_t N>
    class FixedRationalPolynom{
        public:
            using value_type = T;
            using num_type = FixedPolynom<T, N>;
            using den_type = FixedPolynom<T, N>;
            using num_vector_type = typename num_type::vector_type;
            using den_vector_type = typename den_type::vector_type;

        private:
            num_type _num;
            den_type _den;

        public:
            constexpr FixedRationalPolynom() = default;
            constexpr FixedRationalPolynom(const FixedRationalPolynom&) = default;
            constexpr FixedRationalPolynom& operator=(const FixedRationalPolynom&) = default;

            template<std::convertible_to<T> U>
            constexpr FixedRationalPolynom(const U& scalar){
                this->_num.setZero();
                this->_num.at(0) = static_cast<T>(scalar);
                this->_den.setZero();
                this->_den.at(0) = static_cast<T>(0);
            }

            constexpr FixedRationalPolynom(const Polynom<T, N>& num, const Polynom<T, N>& den)
                : _num(num)
                , _den(den){}

            constexpr FixedRationalPolynom(const FixedPolynom<T, N>& num, const FixedPolynom<T, N>& den)
                : _num(num)
                , _den(den){}

            constexpr FixedRationalPolynom(const num_vector_type& num, const den_vector_type& den)
                : _num(num)
                , _den(den){}

            constexpr FixedRationalPolynom(const T(&num)[N], const T(&den)[N])
                : _num(num)
                , _den(den){}

            constexpr operator TransferFunction<T, N, N> () const {
                return TransferFunction<T, N, N>(this->_num, this->_den);
            }

            constexpr TransferFunction<T, N, N> toRationalPolynom() const {
                return TransferFunction<T, N, N>(this->_num, this->_den);
            }

            constexpr num_type& num() {return this->_num;}
            constexpr const num_type& num() const {return this->_num;}

            constexpr den_type& den() {return this->_den;}
            constexpr const den_type& den() const {return this->_den;}

            void prsize_t(std::ostream& stream, std::string_view var="x") const {
                stream << "num: "; this->num().prsize_t(stream, var);
                stream << "\nden: "; this->den().prsize_t(stream, var); stream << '\n';
            }

            friend std::ostream& operator<<(std::ostream& stream, const FixedRationalPolynom& rpoly){
                rpoly.prsize_t(stream);
                return stream;
            }
    };
    // -----------------------------------------------------------------------------------------------
    //                                      Comparison Operators
    // -----------------------------------------------------------------------------------------------

    template<class T, size_t N>
    constexpr bool operator==(const FixedRationalPolynom<T, N>& lhs, const FixedRationalPolynom<T, N>& rhs){
        return ((lhs.num() == rhs.num()) && (lhs.den() == rhs.den()));
    }

    template<class T, size_t N>
    constexpr bool operator!=(const FixedRationalPolynom<T, N>& lhs, const FixedRationalPolynom<T, N>& rhs){
        return ((lhs.num() != rhs.num()) || (lhs.den() != rhs.den()));
    }
    
    // -----------------------------------------------------------------------------------------------
    //                                      Arithmetic Operators
    // -----------------------------------------------------------------------------------------------

    // operator +
    // -----------

    /**
     * \brief rational addition of polynomials
     * 
     * calculates:
     * \f[
     *      \frac{num_1(x)}{den_1(x)} \frac{num_2(x)}{den_2(x)}
     * \f]
     * 
     * \tparam T The value type of the polynomial parameters
     * \tparam NumSize1 The size of the numerator of the left-hand-side addition argument
     * \tparam DenSize1 The size of the denominator of the left-hand-side addition argument
     * \tparam NumSize2 The size of the numberator of the right-hand-side addition argument
     * \tparam DenSize1 The size of the denominator of the right-hand-side addition argument
     * 
     * \param lhs The left-hand-side additino argument as a rational polynom
     * \param rhs The right-hand-side addition argument as a rational polynom
     */
    template<class T, size_t N>
    constexpr auto operator+(const FixedRationalPolynom<T, N>& lhs, const FixedRationalPolynom<T, N>& rhs){
        return FixedRationalPolynom(lhs.num() * rhs.den() + rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t N>
    constexpr auto operator+(const Tscalar& lhs, const FixedRationalPolynom<T, N>& rhs){
        return FixedRationalPolynom(static_cast<T>(lhs) * rhs.den() + rhs.num(), rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t N>
    constexpr auto operator+(const FixedRationalPolynom<T, N>& lhs, const Tscalar& rhs){
        return FixedRationalPolynom(lhs.num() + lhs.den() * static_cast<T>(rhs), lhs.den());
    }

    // operator -
    // -----------

    template<class T, size_t N>
    constexpr FixedRationalPolynom<T, N> operator-(const FixedRationalPolynom<T, N>& poly){
        return FixedRationalPolynom<T, N>(-poly.num(), poly.den());
    }

    template<class T, size_t N>
    constexpr auto operator-(const FixedRationalPolynom<T, N>& lhs, const FixedRationalPolynom<T, N>& rhs){
        return FixedRationalPolynom(lhs.num() * rhs.den() - rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    
    template<class T, std::convertible_to<T> Tscalar, size_t N>
    constexpr auto operator-(const Tscalar& lhs, const FixedRationalPolynom<T, N>& rhs){
        return FixedRationalPolynom(static_cast<T>(lhs) * rhs.den() - rhs.num(), rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t N>
    constexpr auto operator-(const FixedRationalPolynom<T, N>& lhs, const Tscalar& rhs){
        return FixedRationalPolynom(lhs.num() - lhs.den() * static_cast<T>(rhs), lhs.den());
    }

    // operator *
    // -----------

    template<class T, size_t N>
    constexpr FixedRationalPolynom<T, N> operator*(const FixedRationalPolynom<T, N>& lhs, const FixedRationalPolynom<T, N>& rhs){
        return FixedRationalPolynom(lhs.num() * rhs.num(), lhs.den() * rhs.den());
    }

    
    template<class T, std::convertible_to<T> Tscalar, size_t N>
    constexpr FixedRationalPolynom<T, N> operator*(const Tscalar& lhs, const FixedRationalPolynom<T, N>& rhs){
        return FixedRationalPolynom(static_cast<T>(lhs) * rhs.num(), rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t N>
    constexpr FixedRationalPolynom<T, N> operator*(const FixedRationalPolynom<T, N>& lhs, const Tscalar& rhs){
        return FixedRationalPolynom(lhs.num() * static_cast<T>(rhs), lhs.den());
    }

    // operator /
    // -----------

    template<class T, size_t N>
    constexpr auto operator/(const FixedRationalPolynom<T, N>& lhs, const FixedRationalPolynom<T, N>& rhs){
        return FixedRationalPolynom(lhs.num() * rhs.den(), rhs.den() * rhs.num());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t N>
    constexpr auto operator/(const Tscalar& lhs, const FixedRationalPolynom<T, N>& rhs){
        return FixedRationalPolynom(static_cast<T>(lhs) * rhs.den(), rhs.num());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t N>
    constexpr auto operator/(const FixedRationalPolynom<T, N>& lhs, const Tscalar& rhs){
        return FixedRationalPolynom(lhs.num(), rhs.den() * static_cast<T>(rhs));
    }

} // namespace controlpp


namespace Eigen {
    template<typename T, size_t N>
    struct NumTraits<controlpp::FixedRationalPolynom<T, N>> {
        using Self = controlpp::FixedRationalPolynom<T, N>;
        using Real = Self; // or T if you want to treat real parts differently
        using NonInteger = Self;
        using Literal = Self;
        using Nested = Self;

        enum {
            IsComplex = 0,
            IsInteger = 0,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost = 5,
            MulCost = 10
        };

        static inline Self epsilon() { return Self(); }
        static inline Self dummy_precision() { return Self(); }
        static inline Self highest() { return Self(); }
        static inline Self lowest() { return Self(); }
    };

    // ResultType for operator +
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedRationalPolynom<Tpoly, N>,
        internal::scalar_sum_op<Tscalar, controlpp::FixedRationalPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() + std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedRationalPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_sum_op<controlpp::FixedRationalPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>() + std::declval<Tscalar>());
    };

    template <class Tpoly, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedRationalPolynom<Tpoly, N>,
        controlpp::FixedRationalPolynom<Tpoly, N>,
        internal::scalar_sum_op<controlpp::FixedRationalPolynom<Tpoly, N>, controlpp::FixedRationalPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>() + std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>());
    };

    // ResultType for operator -
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedRationalPolynom<Tpoly, N>,
        internal::scalar_difference_op<Tscalar, controlpp::FixedRationalPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() - std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedRationalPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_difference_op<controlpp::FixedRationalPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>() - std::declval<Tscalar>());
    };

    template <class Tpoly, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedRationalPolynom<Tpoly, N>,
        controlpp::FixedRationalPolynom<Tpoly, N>,
        internal::scalar_difference_op<controlpp::FixedRationalPolynom<Tpoly, N>, controlpp::FixedRationalPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>() - std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>());
    };

    // ResultType for operator *
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedRationalPolynom<Tpoly, N>,
        internal::scalar_product_op<Tscalar, controlpp::FixedRationalPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() * std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedRationalPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_product_op<controlpp::FixedRationalPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>() * std::declval<Tscalar>());
    };

    template <class Tpoly, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedRationalPolynom<Tpoly, N>,
        controlpp::FixedRationalPolynom<Tpoly, N>,
        internal::scalar_product_op<controlpp::FixedRationalPolynom<Tpoly, N>, controlpp::FixedRationalPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>() * std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>());
    };

    // ResultType for operator /
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedRationalPolynom<Tpoly, N>,
        internal::scalar_quotient_op<Tscalar, controlpp::FixedRationalPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() / std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedRationalPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_quotient_op<controlpp::FixedRationalPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>() / std::declval<Tscalar>());
    };

    template <class Tpoly, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedRationalPolynom<Tpoly, N>,
        controlpp::FixedRationalPolynom<Tpoly, N>,
        internal::scalar_quotient_op<controlpp::FixedRationalPolynom<Tpoly, N>, controlpp::FixedRationalPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>() / std::declval<controlpp::FixedRationalPolynom<Tpoly, N>>());
    };
}// Eigen