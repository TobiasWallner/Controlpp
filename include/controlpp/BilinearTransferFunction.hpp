#pragma once

#include "Polynom.hpp"
#include "TransferFunction.hpp"

namespace controlpp
{
    template<class ValueType, size_t NumSize, size_t DenSize>
    class BilinearTransferFunction{
        public:
            using value_type = ValueType;
            using ratpoly_type = TransferFunction<ValueType, NumSize, DenSize>;
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

            constexpr explicit BilinearTransferFunction(const Polynom<ValueType, NumSize>& num, const Polynom<ValueType, DenSize>& den)
                : _ratpoly(num, den){}

            constexpr explicit BilinearTransferFunction(const TransferFunction<ValueType, NumSize, DenSize>& ratpoly)
                : _ratpoly(ratpoly){}

            constexpr explicit BilinearTransferFunction(const num_vector_type& num, const den_vector_type& den)
                : _ratpoly(num, den){}

            constexpr explicit BilinearTransferFunction(const ValueType(&num)[NumSize], const ValueType(&den)[DenSize])
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

    template<class ValueType, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator+(const BilinearTransferFunction<ValueType, NumSize1, DenSize1>& lhs, const BilinearTransferFunction<ValueType, NumSize2, DenSize2>& rhs){
        return BilinearTransferFunction(lhs.ratpoly() + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const Tscalar& lhs, const BilinearTransferFunction<Tpoly, NumSize, DenSize>& rhs){
        return BilinearTransferFunction(lhs + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const BilinearTransferFunction<Tpoly, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.ratpoly() + rhs);
    }

    // operator -
    // -----------

    template<class ValueType, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator-(const BilinearTransferFunction<ValueType, NumSize1, DenSize1>& lhs, const BilinearTransferFunction<ValueType, NumSize2, DenSize2>& rhs){
        return BilinearTransferFunction(lhs.ratpoly() - rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const Tscalar& lhs, const BilinearTransferFunction<ValueType, NumSize, DenSize>& rhs){
        return BilinearTransferFunction(lhs + rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const BilinearTransferFunction<ValueType, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.ratpoly() - rhs);
    }

    // operator *
    // -----------

    template<class ValueType, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator*(const BilinearTransferFunction<ValueType, NumSize1, DenSize1>& lhs, const BilinearTransferFunction<ValueType, NumSize2, DenSize2>& rhs){
        return BilinearTransferFunction(lhs.ratpoly() * rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const Tscalar& lhs, const BilinearTransferFunction<ValueType, NumSize, DenSize>& rhs){
        return BilinearTransferFunction(lhs * rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const BilinearTransferFunction<ValueType, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.ratpoly() * rhs);
    }

    // operator /
    // -----------

    template<class ValueType, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator/(const BilinearTransferFunction<ValueType, NumSize1, DenSize1>& lhs, const BilinearTransferFunction<ValueType, NumSize2, DenSize2>& rhs){
        return BilinearTransferFunction(lhs.ratpoly() / rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const Tscalar& lhs, const BilinearTransferFunction<ValueType, NumSize, DenSize>& rhs){
        return BilinearTransferFunction(lhs / rhs.ratpoly());
    }

    template<class ValueType, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const BilinearTransferFunction<ValueType, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.ratpoly() / rhs);
    }

    namespace tf{
        template<class ValueType=double>
        static inline const BilinearTransferFunction<ValueType, 2, 1> q({ValueType(0), ValueType(1)}, {ValueType(1)});
    }

} // namespace controlpp