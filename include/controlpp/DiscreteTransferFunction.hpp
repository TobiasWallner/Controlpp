#pragma once

namespace control
{
    template<class ValueType, class TimeType, size_t NumSize, size_t DenSize>
    class DiscreteTransferFunction{
        public:
            using value_type = ValueType;
            using time_type = TimeType;
            using ratpoly_type = RationalPolynomial<T, NumSize, DenSize>;
            using num_type = typename ratpoly_type::num_type;
            using den_type = typename ratpoly_type::den_type;

        private:
            ratpoly_type _ratpoly;
            time_type _time;

        public:

            constexpr DiscreteTransferFunction() = default;
            constexpr DiscreteTransferFunction(const DiscreteTransferFunction&) = default;
            constexpr DiscreteTransferFunction& operator=(const DiscreteTransferFunction&) = default;

            constexpr DiscreteTransferFunction(const Polynom<T, NumeratorSize>& num, const Polynom<T, DenominatorSize>& den)
                : _ratpoly(num, den){}

            constexpr DiscreteTransferFunction(const RationalPolynomial<T, NumSize, DenSize>& ratpoly)
                : _ratpoly(ratpoly){}

            constexpr num_type& num() {return this->_ratpoly.num();}
            constexpr const num_type& num() const {return this->_ratpoly.num();}

            constexpr num_type& den() {return this->_ratpoly.den();}
            constexpr const num_type& den() const {return this->_ratpoly.den();}

            constexpr time_type& sampleTime() {return this->_ratpoly._time();}
            constexpr const time_type& sampleTime() const {return this->_ratpoly._time();}
    };

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator+(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize1, DenSize1>& rhs){
        return DiscreteTransferFunction(lhs.num() * rhs.den() + rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator-(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize1, DenSize1>& rhs){
        return DiscreteTransferFunction(lhs.num() * rhs.den() - rhs.num() * lhs.den(), lhs.den() * rhs.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator*(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize1, DenSize1>& rhs){
        return DiscreteTransferFunction(lhs.num() * rhs.num(), rhs.den() * rhs.den());
    }

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator/(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize1, DenSize1>& rhs){
        return DiscreteTransferFunction(lhs.num() * rhs.den(), rhs.den() * rhs.num());
    }
} // namespace control