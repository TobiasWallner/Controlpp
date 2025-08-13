#pragma once

#include "Polynom.hpp"
#include "TransferFunction.hpp"

namespace controlpp
{
    template<class ValueType, int NumOrder, int DenOrder>
    class BilinearTransferFunction{
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

            constexpr BilinearTransferFunction() = default;
            constexpr BilinearTransferFunction(const BilinearTransferFunction&) = default;
            constexpr BilinearTransferFunction& operator=(const BilinearTransferFunction&) = default;

            constexpr explicit BilinearTransferFunction(const Polynom<ValueType, NumOrder>& num, const Polynom<ValueType, DenOrder>& den)
                : _ratpoly(num, den){}

            constexpr explicit BilinearTransferFunction(const TransferFunction<ValueType, NumOrder, DenOrder>& ratpoly)
                : _ratpoly(ratpoly){}

            constexpr explicit BilinearTransferFunction(const num_vector_type& num, const den_vector_type& den)
                : _ratpoly(num, den){}

            constexpr explicit BilinearTransferFunction(const ValueType(&num)[NumOrder+1], const ValueType(&den)[DenOrder+1])
                : _ratpoly(num, den){}

            constexpr num_type& num() {return this->_ratpoly.num();}
            constexpr const num_type& num() const {return this->_ratpoly.num();}

            constexpr den_type& den() {return this->_ratpoly.den();}
            constexpr const den_type& den() const {return this->_ratpoly.den();}

            constexpr ratpoly_type& ratpoly() {return this->_ratpoly;}
            constexpr const ratpoly_type& ratpoly() const {return this->_ratpoly;}

            friend std::ostream& operator<<(std::ostream& stream, const BilinearTransferFunction& ctf){
                ctf.ratpoly().print(stream, "s");
                return stream;
            }
    };

    // operator +
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator+(const BilinearTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const BilinearTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return BilinearTransferFunction(lhs.ratpoly() + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const Tscalar& lhs, const BilinearTransferFunction<Tpoly, NumOrder, DenOrder>& rhs){
        return BilinearTransferFunction(lhs + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const BilinearTransferFunction<Tpoly, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.ratpoly() + rhs);
    }

    // operator -
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator-(const BilinearTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const BilinearTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return BilinearTransferFunction(lhs.ratpoly() - rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const Tscalar& lhs, const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& rhs){
        return BilinearTransferFunction(lhs + rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.ratpoly() - rhs);
    }

    // operator *
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator*(const BilinearTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const BilinearTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return BilinearTransferFunction(lhs.ratpoly() * rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const Tscalar& lhs, const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& rhs){
        return BilinearTransferFunction(lhs * rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.ratpoly() * rhs);
    }

    // operator /
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator/(const BilinearTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const BilinearTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return BilinearTransferFunction(lhs.ratpoly() / rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const Tscalar& lhs, const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& rhs){
        return BilinearTransferFunction(lhs / rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.ratpoly() / rhs);
    }

    namespace tf{
        template<class ValueType=double>
        static inline const BilinearTransferFunction<ValueType, 1, 0> q({ValueType(0), ValueType(1)}, {ValueType(1)});
    }

} // namespace controlpp