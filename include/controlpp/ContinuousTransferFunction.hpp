#pragma once

#include "Polynom.hpp"
#include "RationalPolynom.hpp"

namespace control
{
    template<class T, size_t NumSize, size_t DenSize>
    class ContinuousTransferFunction{
        public:
            using value_type = T;
            using ratpoly_type = RationalPolynomial<T, NumSize, DenSize>;
            using num_type = typename ratpoly_type::num_type;
            using den_type = typename ratpoly_type::den_type;

        private:
            ratpoly_type ratpoly;

        public:

            constexpr ContinuousTransferFunction() = default;
            constexpr ContinuousTransferFunction(const ContinuousTransferFunction&) = default;
            constexpr ContinuousTransferFunction& operator=(const ContinuousTransferFunction&) = default;

            constexpr ContinuousTransferFunction(const Polynom<T, NumeratorSize>& num, const Polynom<T, DenominatorSize>& den)
                : ratpoly_type(num, den){}

            constexpr ContinuousTransferFunction(const RationalPolynomial<T, NumSize, DenSize>& ratpoly)
                : ratpoly_type(ratpoly){}

            constexpr num_type& num() {return this->ratpoly.num();}
            constexpr const num_type& num() const {return this->ratpoly.num();}

            constexpr num_type& den() {return this->ratpoly.den();}
            constexpr const num_type& den() const {return this->ratpoly.den();}
    };

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator+(const ContinuousTransferFunction<T, NumSize1, DenSize1>& lhs, const ContinuousTransferFunction<T, NumSize1, DenSize1>& rhs){
        return ContinuousTransferFunction(lhs.num() * rhs.den() + rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator-(const ContinuousTransferFunction<T, NumSize1, DenSize1>& lhs, const ContinuousTransferFunction<T, NumSize1, DenSize1>& rhs){
        return ContinuousTransferFunction(lhs.num() * rhs.den() - rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator*(const ContinuousTransferFunction<T, NumSize1, DenSize1>& lhs, const ContinuousTransferFunction<T, NumSize1, DenSize1>& rhs){
        return ContinuousTransferFunction(lhs.num() * rhs.num(), rhs.den() * rhs.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator/(const ContinuousTransferFunction<T, NumSize1, DenSize1>& lhs, const ContinuousTransferFunction<T, NumSize1, DenSize1>& rhs){
        return ContinuousTransferFunction(lhs.num() * rhs.den(), rhs.den() * rhs.num());
    }


} // namespace control