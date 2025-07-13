#pragma once

namespace controlpp
{
    template<class ValueType, size_t NumSize, size_t DenSize>
    class DiscreteTransferFunction{
        public:
            using value_type = ValueType;
            using ratpoly_type = RationalPolynom<ValueType, NumSize, DenSize>;
            using num_type = typename ratpoly_type::num_type;
            using den_type = typename ratpoly_type::den_type;
            using num_vector_type = typename ratpoly_type::num_vector_type;
            using den_vector_type = typename ratpoly_type::den_vector_type;

        private:
            ratpoly_type _ratpoly;
            value_type _sample_time;

        public:

            constexpr DiscreteTransferFunction() = default;
            constexpr DiscreteTransferFunction(const DiscreteTransferFunction&) = default;
            constexpr DiscreteTransferFunction& operator=(const DiscreteTransferFunction&) = default;

            constexpr DiscreteTransferFunction(
                const Polynom<ValueType, NumSize>& num, 
                const Polynom<ValueType, DenSize>& den)
                : _ratpoly(num, den){}

            constexpr DiscreteTransferFunction(
                const RationalPolynom<ValueType, NumSize, DenSize>& ratpoly)
                : _ratpoly(ratpoly){}

            constexpr explicit DiscreteTransferFunction(
                const num_vector_type& num, 
                const den_vector_type& den)
                : _ratpoly(num, den){}

            constexpr explicit DiscreteTransferFunction(
                const value_type(&num)[NumSize], 
                const value_type(&den)[DenSize])
                : _ratpoly(num, den){}

            constexpr ratpoly_type& ratpoly() {return this->_ratpoly;}
            constexpr const ratpoly_type& ratpoly() const {return this->_ratpoly;}

            constexpr num_type& num() {return this->_ratpoly.num();}
            constexpr const num_type& num() const {return this->_ratpoly.num();}

            constexpr ValueType& num(size_t i) {return this->_ratpoly.num(i);}
            constexpr const ValueType& num(size_t i) const {return this->_ratpoly.num(i);}

            constexpr den_type& den() {return this->_ratpoly.den();}
            constexpr const den_type& den() const {return this->_ratpoly.den();}

            constexpr ValueType& den(size_t i) {return this->_ratpoly.den(i);}
            constexpr const ValueType& den(size_t i) const {return this->_ratpoly.den(i);}

            friend inline std::ostream& operator<<(std::ostream& stream, const DiscreteTransferFunction& dtf){
                dtf.ratpoly().print(stream, "z");
                return stream;
            }
    };

    // operator +
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator+(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize2, DenSize2>& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() + rhs.ratpoly());
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const Tscalar& lhs, const DiscreteTransferFunction<Tpoly, NumSize, DenSize>& rhs){
        return DiscreteTransferFunction(lhs + rhs.ratpoly());
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const DiscreteTransferFunction<Tpoly, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() + rhs);
    }

    // operator -
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator-(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize2, DenSize2>& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() - rhs.ratpoly());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const Tscalar& lhs, const DiscreteTransferFunction<T, NumSize, DenSize>& rhs){
        return DiscreteTransferFunction(lhs + rhs.ratpoly());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const DiscreteTransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() - rhs);
    }

    // operator *
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator*(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize2, DenSize2>& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() * rhs.ratpoly());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const Tscalar& lhs, const DiscreteTransferFunction<T, NumSize, DenSize>& rhs){
        return DiscreteTransferFunction(lhs * rhs.ratpoly());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const DiscreteTransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() * rhs);
    }

    // operator /
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator/(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize2, DenSize2>& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() / rhs.ratpoly());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const Tscalar& lhs, const DiscreteTransferFunction<T, NumSize, DenSize>& rhs){
        return DiscreteTransferFunction(lhs / rhs.ratpoly());
    }

    template<class T, std::convertible_to<T> Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const DiscreteTransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() / rhs);
    }

    namespace tf{
        template<class ValueType=double>
        static inline const DiscreteTransferFunction<ValueType, 2, 1> z = DiscreteTransferFunction<ValueType, 2, 1>({ValueType(0), ValueType(1)}, {ValueType(1)});
    }

} // namespace controlpp



