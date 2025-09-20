#pragma once

#include <controlpp/TransferFunction.hpp>

namespace controlpp
{
    /**
     * \brief Continuous transfer functions in the s lapace plain
     * 
     * Transfer functions of the shape:
     * 
     * \f[
     * G(s) = \frac{b_0 + b_1 d \cdots b_m d^m}{a_0 + a_1 d \cdots a_n d^n}
     * \f]
     * 
     * Note that the transfer function uses the delay operator: \f$d = z^{-1}\f$.
     * 
     * \tparam T The value type that the transfer function should use. Like `double` or `float`.
     * \tparam NumberOrder The order of the numerator. m in the equation.
     * \tparam DenOrder The order of the denominator. n in the equation.
     */
    template<class T, int NumOrder, int DenOrder>
    class DiscreteTransferFunction{
        public:
            using value_type = T;
            using transfer_function_type = TransferFunction<T, NumOrder, DenOrder>;
            using num_type = typename transfer_function_type::num_type;
            using den_type = typename transfer_function_type::den_type;
            using num_vector_type = typename transfer_function_type::num_vector_type;
            using den_vector_type = typename transfer_function_type::den_vector_type;

        private:
            transfer_function_type tf_;

        public:

            constexpr DiscreteTransferFunction() = default;
            constexpr DiscreteTransferFunction(const DiscreteTransferFunction&) = default;

            constexpr DiscreteTransferFunction& operator=(const DiscreteTransferFunction&) = default;

            
            template<int NumOrder2, int DenOrder2>
            requires(((NumOrder2 < NumOrder) || (DenOrder2 < DenOrder)) && (NumOrder2 <= NumOrder) && (DenOrder2 <= DenOrder))
            DiscreteTransferFunction(const DiscreteTransferFunction<T, NumOrder2, DenOrder2>& other)
                : tf_(other.transfer_function()){}

            template<int NumOrder2, int DenOrder2>
            requires(((NumOrder2 < NumOrder) || (DenOrder2 < DenOrder)) && (NumOrder2 <= NumOrder) && (DenOrder2 <= DenOrder))
            DiscreteTransferFunction& operator=(const DiscreteTransferFunction<T, NumOrder2, DenOrder2>& other){
                this->tf_ = other.transfer_function();
                return *this;
            }

            constexpr DiscreteTransferFunction(
                const Polynom<T, NumOrder>& num, 
                const Polynom<T, DenOrder>& den)
                : tf_(num, den){}

            constexpr DiscreteTransferFunction(const TransferFunction<T, NumOrder, DenOrder>& transfer_function)
                : tf_(transfer_function){}

            constexpr explicit DiscreteTransferFunction(
                const num_vector_type& num, 
                const den_vector_type& den
            )
                : tf_(num, den){}

            constexpr explicit DiscreteTransferFunction(
                const value_type(&num)[NumOrder+1], 
                const value_type(&den)[DenOrder+1])
                : tf_(num, den){}
            

            constexpr T eval(const Eigen::Vector<T, NumOrder+1>& input_series, const Eigen::Vector<T, DenOrder>& output_series){
                const T B = this->num().vector().dot(input_series);
                const T A = this->den().vector().tail(DenOrder).dot(output_series);
                const T y = (B - A) / this->den(0);
                return y;
            }
            
            constexpr transfer_function_type& transfer_function() {return this->tf_;}
            constexpr const transfer_function_type& transfer_function() const {return this->tf_;}
            
            constexpr num_type& num() {return this->tf_.num();}
            constexpr const num_type& num() const {return this->tf_.num();}

            constexpr T& num(int i) {return this->tf_.num(i);}
            constexpr const T& num(int i) const {return this->tf_.num(i);}

            constexpr den_type& den() {return this->tf_.den();}
            constexpr const den_type& den() const {return this->tf_.den();}

            constexpr T& den(int i) {return this->tf_.den(i);}
            constexpr const T& den(int i) const {return this->tf_.den(i);}

            friend inline std::ostream& operator<<(std::ostream& stream, const DiscreteTransferFunction& dtf){
                dtf.transfer_function().print(stream, "z^{-1}");
                return stream;
            }
    };

    // operator +
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator+(const DiscreteTransferFunction<T, NumOrder1, DenOrder1>& lhs, const DiscreteTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return DiscreteTransferFunction(lhs.transfer_function() + rhs.transfer_function());
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const Tscalar& lhs, const DiscreteTransferFunction<Tpoly, NumOrder, DenOrder>& rhs){
        return DiscreteTransferFunction(lhs + rhs.transfer_function());
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator+(const DiscreteTransferFunction<Tpoly, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.transfer_function() + rhs);
    }

    // operator -
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator-(const DiscreteTransferFunction<T, NumOrder1, DenOrder1>& lhs, const DiscreteTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return DiscreteTransferFunction(lhs.transfer_function() - rhs.transfer_function());
    }

    template<class T, std::convertible_to<T> Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const Tscalar& lhs, const DiscreteTransferFunction<T, NumOrder, DenOrder>& rhs){
        return DiscreteTransferFunction(lhs - rhs.transfer_function());
    }

    template<class T, std::convertible_to<T> Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator-(const DiscreteTransferFunction<T, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.transfer_function() - rhs);
    }

    // operator *
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator*(const DiscreteTransferFunction<T, NumOrder1, DenOrder1>& lhs, const DiscreteTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return DiscreteTransferFunction(lhs.transfer_function() * rhs.transfer_function());
    }

    template<class T, std::convertible_to<T> Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const Tscalar& lhs, const DiscreteTransferFunction<T, NumOrder, DenOrder>& rhs){
        return DiscreteTransferFunction(lhs * rhs.transfer_function());
    }

    template<class T, std::convertible_to<T> Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator*(const DiscreteTransferFunction<T, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.transfer_function() * rhs);
    }

    // operator /
    // -----------

    template<class T, int NumOrder1, int DenOrder1, int NumOrder2, int DenOrder2>
    constexpr auto operator/(const DiscreteTransferFunction<T, NumOrder1, DenOrder1>& lhs, const DiscreteTransferFunction<T, NumOrder2, DenOrder2>& rhs){
        return DiscreteTransferFunction(lhs.transfer_function() / rhs.transfer_function());
    }

    template<class T, std::convertible_to<T> Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const Tscalar& lhs, const DiscreteTransferFunction<T, NumOrder, DenOrder>& rhs){
        return DiscreteTransferFunction(lhs / rhs.transfer_function());
    }

    template<class T, std::convertible_to<T> Tscalar, int NumOrder, int DenOrder>
    constexpr auto operator/(const DiscreteTransferFunction<T, NumOrder, DenOrder>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.transfer_function() / rhs);
    }

    namespace tf{

        /**
         * \brief The delay operator \f$z^{-1}\f$
         * 
         * \tparam T The value type of the delay operator. Usually `double`, `float` or a custom fixpoint.
         */
        template<class T=double>
        static inline const DiscreteTransferFunction<T, 1, 0> z_1({T(0), T(1)}, {T(1)});
    }

} // namespace controlpp



