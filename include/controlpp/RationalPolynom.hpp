#pragma once

// std
#include <initializer_list>

#include "Polynom.hpp"

namespace controlpp
{
    template<class T, size_t NumSize, size_t DenSize>
    class RationalPolynom{
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
            constexpr RationalPolynom() = default;
            constexpr RationalPolynom(const RationalPolynom&) = default;
            constexpr RationalPolynom& operator=(const RationalPolynom&) = default;

            constexpr RationalPolynom(const Polynom<T, NumSize>& num, const Polynom<T, DenSize>& den)
                : _num(num)
                , _den(den){}

            constexpr RationalPolynom(const num_vector_type& num, const den_vector_type& den)
                : _num(num)
                , _den(den){}

            constexpr RationalPolynom(const T(&num)[NumSize], const T(&den)[DenSize])
                : _num(num)
                , _den(den){}

            constexpr num_type& num() {return this->_num;}
            constexpr const num_type& num() const {return this->_num;}

            constexpr den_type& den() {return this->_den;}
            constexpr const den_type& den() const {return this->_den;}

            void print(std::ostream& stream, std::string_view var="x") const {
                stream << "num: "; this->num().print(stream, var);
                stream << "\nden: "; this->den().print(stream, var); stream << '\n';
            }

            friend std::ostream& operator<<(std::ostream& stream, const RationalPolynom& rpoly){
                rpoly.print(stream);
                return stream;
            }
    };
    // -----------------------------------------------------------------------------------------------
    //                                      Comparison Operators
    // -----------------------------------------------------------------------------------------------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator==(const RationalPolynom<T, NumSize1, DenSize1>& lhs, const RationalPolynom<T, NumSize1, DenSize1>& rhs){
        return ((lhs.num() == rhs.num()) && (lhs.den() == rhs.den()));
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr bool operator!=(const RationalPolynom<T, NumSize1, DenSize1>& lhs, const RationalPolynom<T, NumSize1, DenSize1>& rhs){
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
    constexpr auto operator+(const RationalPolynom<T, NumSize1, DenSize1>& lhs, const RationalPolynom<T, NumSize2, DenSize2>& rhs){
        return RationalPolynom(lhs.num() * rhs.den() + rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const Tscalar& lhs, const RationalPolynom<T, NumSize, DenSize>& rhs){
        return RationalPolynom(static_cast<T>(lhs) * rhs.den() + rhs.num(), rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const RationalPolynom<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return RationalPolynom(lhs.num() + lhs.den() * static_cast<T>(rhs), lhs.den());
    }

    // operator -
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator-(const RationalPolynom<T, NumSize1, DenSize1>& lhs, const RationalPolynom<T, NumSize2, DenSize2>& rhs){
        return RationalPolynom(lhs.num() * rhs.den() - rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    
    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const Tscalar& lhs, const RationalPolynom<T, NumSize, DenSize>& rhs){
        return RationalPolynom(static_cast<T>(lhs) * rhs.den() - rhs.num(), rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const RationalPolynom<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return RationalPolynom(lhs.num() - lhs.den() * static_cast<T>(rhs), lhs.den());
    }

    // operator *
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator*(const RationalPolynom<T, NumSize1, DenSize1>& lhs, const RationalPolynom<T, NumSize2, DenSize2>& rhs){
        return RationalPolynom(lhs.num() * rhs.num(), lhs.den() * rhs.den());
    }

    
    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const Tscalar& lhs, const RationalPolynom<T, NumSize, DenSize>& rhs){
        return RationalPolynom(static_cast<T>(lhs) * rhs.num(), rhs.den());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const RationalPolynom<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return RationalPolynom(lhs.num() * static_cast<T>(rhs), lhs.den());
    }

    // operator /
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator/(const RationalPolynom<T, NumSize1, DenSize1>& lhs, const RationalPolynom<T, NumSize2, DenSize2>& rhs){
        return RationalPolynom(lhs.num() * rhs.den(), rhs.den() * rhs.num());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const Tscalar& lhs, const RationalPolynom<T, NumSize, DenSize>& rhs){
        return RationalPolynom(static_cast<T>(lhs) * rhs.den(), rhs.num());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const RationalPolynom<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return RationalPolynom(lhs.num(), rhs.den() * static_cast<T>(rhs));
    }




} // namespace controlpp
