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
     * \tparam ValueType The value type that the transfer function should use
     * \tparam NumberOrder The order of the numerator. m in the equation.
     * \tparam DenOrder The order of the denominator. n in the equation.
     */
    template<class ValueType, int NumOrder, int DenOrder>
    class ContinuousTransferFunction{
        public:
            using value_type = ValueType;
            using ratpoly_type = TransferFunction<ValueType, NumOrder, DenOrder>;
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

            constexpr explicit ContinuousTransferFunction(const Polynom<ValueType, NumOrder>& num, const Polynom<ValueType, DenOrder>& den)
                : _ratpoly(num, den){}

            constexpr explicit ContinuousTransferFunction(const TransferFunction<ValueType, NumOrder, DenOrder>& ratpoly)
                : _ratpoly(ratpoly){}

            constexpr explicit ContinuousTransferFunction(const num_vector_type& num, const den_vector_type& den)
                : _ratpoly(num, den){}

            constexpr explicit ContinuousTransferFunction(const ValueType(&num)[NumOrder+1], const ValueType(&den)[DenOrder+1])
                : _ratpoly(num, den){}

            constexpr num_type& num() {return this->_ratpoly.num();}
            constexpr const num_type& num() const {return this->_ratpoly.num();}

            constexpr ValueType& num(size_t i) {return this->_ratpoly.num(i);}
            constexpr const ValueType& num(size_t i) const {return this->_ratpoly.num(i);}

            constexpr den_type& den() {return this->_ratpoly.den();}
            constexpr const den_type& den() const {return this->_ratpoly.den();}

            constexpr ValueType& den(size_t i) {return this->_ratpoly.den(i);}
            constexpr const ValueType& den(size_t i) const {return this->_ratpoly.den(i);}

            constexpr ratpoly_type& ratpoly() {return this->_ratpoly;}
            constexpr const ratpoly_type& ratpoly() const {return this->_ratpoly;}

            friend std::ostream& operator<<(std::ostream& stream, const ContinuousTransferFunction& ctf){
                ctf.ratpoly().print(stream, "s");
                return stream;
            }
    };

    // operator +
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator+(const ContinuousTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
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

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator-(const ContinuousTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() - rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const Tscalar& lhs, const ContinuousTransferFunction<ValueType, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs + rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const ContinuousTransferFunction<ValueType, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() - rhs);
    }

    // operator *
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator*(const ContinuousTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() * rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const Tscalar& lhs, const ContinuousTransferFunction<ValueType, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs * rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const ContinuousTransferFunction<ValueType, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() * rhs);
    }

    // operator /
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator/(const ContinuousTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const ContinuousTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return ContinuousTransferFunction(lhs.ratpoly() / rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const Tscalar& lhs, const ContinuousTransferFunction<ValueType, NumOrder, DenOrder>& rhs){
        return ContinuousTransferFunction(lhs / rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const ContinuousTransferFunction<ValueType, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        ContinuousTransferFunction result(lhs.ratpoly() / rhs);
        return result;
    }

    namespace tf{
        template<class ValueType=double>
        static inline const ContinuousTransferFunction<ValueType, 1, 0> s({ValueType(0), ValueType(1)}, {ValueType(1)});
    }

} // namespace controlpp