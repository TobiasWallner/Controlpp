#pragma once

#include "Polynom.hpp"
#include "RationalPolynom.hpp"

namespace controlpp
{
    template<class ValueType, size_t NumSize, size_t DenSize>
    class ContinuousTransferFunction{
        public:
            using value_type = ValueType;
            using ratpoly_type = RationalPolynom<ValueType, NumSize, DenSize>;
            using num_type = typename ratpoly_type::num_type;
            using den_type = typename ratpoly_type::den_type;
            using num_vector_type = typename ratpoly_type::num_vector_type;
            using den_vector_type = typename ratpoly_type::den_vector_type;

        private:
            ratpoly_type _ratpoly;

        public:

            constexpr ContinuousTransferFunction() = default;
            constexpr ContinuousTransferFunction(const ContinuousTransferFunction&) = default;
            constexpr ContinuousTransferFunction& operator=(const ContinuousTransferFunction&) = default;

            constexpr explicit ContinuousTransferFunction(const Polynom<ValueType, NumSize>& num, const Polynom<ValueType, DenSize>& den)
                : _ratpoly(num, den){}

            constexpr explicit ContinuousTransferFunction(const RationalPolynom<ValueType, NumSize, DenSize>& ratpoly)
                : _ratpoly(ratpoly){}

            constexpr explicit ContinuousTransferFunction(const num_vector_type& num, const den_vector_type& den)
                : _ratpoly(num, den){}

            constexpr explicit ContinuousTransferFunction(const ValueType(&num)[NumSize], const ValueType(&den)[DenSize])
                : _ratpoly(num, den){}

            constexpr num_type& num() {return this->_ratpoly.num();}
            constexpr const num_type& num() const {return this->_ratpoly.num();}

            constexpr num_type& den() {return this->_ratpoly.den();}
            constexpr const num_type& den() const {return this->_ratpoly.den();}

            constexpr ratpoly_type& ratpoly() {return this->_ratpoly;}
            constexpr const ratpoly_type& ratpoly() const {return this->_ratpoly;}

            friend std::ostream& operator<<(std::ostream& stream, const ContinuousTransferFunction& ctf){
                ctf.ratpoly().print(stream, "s");
                return stream;
            }
    };

    // operator +
    // -----------

    template<class ValueType, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator+(const ContinuousTransferFunction<ValueType, NumSize1, DenSize1>& lhs, const ContinuousTransferFunction<ValueType, NumSize2, DenSize2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const Tscalar& lhs, const ContinuousTransferFunction<Tpoly, NumSize, DenSize>& rhs){
        return ContinuousTransferFunction(lhs + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const ContinuousTransferFunction<Tpoly, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() + rhs);
    }

    // operator -
    // -----------

    template<class ValueType, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator-(const ContinuousTransferFunction<ValueType, NumSize1, DenSize1>& lhs, const ContinuousTransferFunction<ValueType, NumSize2, DenSize2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() - rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const Tscalar& lhs, const ContinuousTransferFunction<ValueType, NumSize, DenSize>& rhs){
        return ContinuousTransferFunction(lhs + rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const ContinuousTransferFunction<ValueType, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() - rhs);
    }

    // operator *
    // -----------

    template<class ValueType, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator*(const ContinuousTransferFunction<ValueType, NumSize1, DenSize1>& lhs, const ContinuousTransferFunction<ValueType, NumSize2, DenSize2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() * rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const Tscalar& lhs, const ContinuousTransferFunction<ValueType, NumSize, DenSize>& rhs){
        return ContinuousTransferFunction(lhs * rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const ContinuousTransferFunction<ValueType, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() * rhs);
    }

    // operator /
    // -----------

    template<class ValueType, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator/(const ContinuousTransferFunction<ValueType, NumSize1, DenSize1>& lhs, const ContinuousTransferFunction<ValueType, NumSize2, DenSize2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() / rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const Tscalar& lhs, const ContinuousTransferFunction<ValueType, NumSize, DenSize>& rhs){
        return ContinuousTransferFunction(lhs / rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const ContinuousTransferFunction<ValueType, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() / rhs);
    }

    namespace tf{
        template<class ValueType=double>
        static inline const ContinuousTransferFunction<ValueType, 2, 1> s({ValueType(0), ValueType(1)}, {ValueType(1)});
    }

} // namespace controlpp