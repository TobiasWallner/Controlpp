#pragma once

#include "Polynom.hpp"

namespace control
{
    template<class T, size_t NumSize, size_t DenSize>
    class RationalPolynomial{
        public:
            using num_type = Polynom<T, NumeratorSize>;
            using den_type = Polynom<T, DenominatorSize>;

        private:
            num_type _num;
            den_type _den;

        public:
            constexpr RationalPolynomial() = default;
            constexpr RationalPolynomial(const RationalPolynomial&) = default;
            constexpr RationalPolynomial& operator=(const RationalPolynomial&) = default;

            constexpr RationalPolynomial(const Polynom<T, NumeratorSize>& num, const Polynom<T, DenominatorSize>& den)
                : _num(num)
                , _den(den){}

            constexpr num_type& num() {return this->_num;}
            constexpr const num_type& num() const {return this->_num;}

            constexpr num_type& den() {return this->_den;}
            constexpr const num_type& den() const {return this->_den;}
    };

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator+(const RationalPolynomial<T, NumSize1, DenSize1>& lhs, const RationalPolynomial<T, NumSize1, DenSize1>& rhs){
        return RationalPolynomial(lhs.num() * rhs.den() + rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator-(const RationalPolynomial<T, NumSize1, DenSize1>& lhs, const RationalPolynomial<T, NumSize1, DenSize1>& rhs){
        return RationalPolynomial(lhs.num() * rhs.den() - rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator*(const RationalPolynomial<T, NumSize1, DenSize1>& lhs, const RationalPolynomial<T, NumSize1, DenSize1>& rhs){
        return RationalPolynomial(lhs.num() * rhs.num(), rhs.den() * rhs.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator/(const RationalPolynomial<T, NumSize1, DenSize1>& lhs, const RationalPolynomial<T, NumSize1, DenSize1>& rhs){
        return RationalPolynomial(lhs.num() * rhs.den(), rhs.den() * rhs.num());
    }


} // namespace control
