#pragma once

#include "Polynom.hpp"
#include "TransferFunction.hpp"

namespace controlpp
{
    template<class ValueType, int NumOrder, int DenOrder>
    class BilinearTransferFunction{
        public:
            using value_type = ValueType;
            using transfer_function_type = TransferFunction<ValueType, NumOrder, DenOrder>;
            using num_type = typename transfer_function_type::num_type;
            using den_type = typename transfer_function_type::den_type;
            using num_vector_type = typename transfer_function_type::num_vector_type;
            using den_vector_type = typename transfer_function_type::den_vector_type;

        private:
            transfer_function_type tf_;

        public:

            constexpr BilinearTransferFunction() = default;
            constexpr BilinearTransferFunction(const BilinearTransferFunction&) = default;
            constexpr BilinearTransferFunction& operator=(const BilinearTransferFunction&) = default;

            constexpr explicit BilinearTransferFunction(const Polynom<ValueType, NumOrder>& num, const Polynom<ValueType, DenOrder>& den)
                : tf_(num, den){}

            constexpr explicit BilinearTransferFunction(const TransferFunction<ValueType, NumOrder, DenOrder>& transfer_function)
                : tf_(transfer_function){}

            constexpr explicit BilinearTransferFunction(const num_vector_type& num, const den_vector_type& den)
                : tf_(num, den){}

            constexpr explicit BilinearTransferFunction(const ValueType(&num)[NumOrder+1], const ValueType(&den)[DenOrder+1])
                : tf_(num, den){}

            constexpr num_type& num() {return this->tf_.num();}
            constexpr const num_type& num() const {return this->tf_.num();}

            constexpr den_type& den() {return this->tf_.den();}
            constexpr const den_type& den() const {return this->tf_.den();}

            constexpr transfer_function_type& transfer_function() {return this->tf_;}
            constexpr const transfer_function_type& transfer_function() const {return this->tf_;}

            friend std::ostream& operator<<(std::ostream& stream, const BilinearTransferFunction& ctf){
                ctf.transfer_function().print(stream, "s");
                return stream;
            }
    };

    // operator +
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator+(const BilinearTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const BilinearTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return BilinearTransferFunction(lhs.transfer_function() + rhs.transfer_function());
    }

    template<class Tpoly, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const Tscalar& lhs, const BilinearTransferFunction<Tpoly, NumOrder, DenOrder>& rhs){
        return BilinearTransferFunction(lhs + rhs.transfer_function());
    }

    template<class Tpoly, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const BilinearTransferFunction<Tpoly, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.transfer_function() + rhs);
    }

    // operator -
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator-(const BilinearTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const BilinearTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return BilinearTransferFunction(lhs.transfer_function() - rhs.transfer_function());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const Tscalar& lhs, const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& rhs){
        return BilinearTransferFunction(lhs - rhs.transfer_function());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.transfer_function() - rhs);
    }

    // operator *
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator*(const BilinearTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const BilinearTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return BilinearTransferFunction(lhs.transfer_function() * rhs.transfer_function());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const Tscalar& lhs, const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& rhs){
        return BilinearTransferFunction(lhs * rhs.transfer_function());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.transfer_function() * rhs);
    }

    // operator /
    // -----------

    template<class ValueType, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator/(const BilinearTransferFunction<ValueType, NumOrder1, DenOrder1>& lhs, const BilinearTransferFunction<ValueType, NumOrder2, DenOrder2>& rhs){
        return BilinearTransferFunction(lhs.transfer_function() / rhs.transfer_function());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const Tscalar& lhs, const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& rhs){
        return BilinearTransferFunction(lhs / rhs.transfer_function());
    }

    template<class ValueType, class Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return BilinearTransferFunction(lhs.transfer_function() / rhs);
    }

    namespace tf{
        template<class ValueType=double>
        static inline const BilinearTransferFunction<ValueType, 1, 0> q({ValueType(0), ValueType(1)}, {ValueType(1)});
    }

} // namespace controlpp