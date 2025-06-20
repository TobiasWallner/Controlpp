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
                const Polynom<ValueType, DenSize>& den, 
                const ValueType& sample_time)
                : _ratpoly(num, den)
                , _sample_time(sample_time){}

            constexpr DiscreteTransferFunction(
                const RationalPolynom<ValueType, NumSize, DenSize>& ratpoly,
                const ValueType& sample_time)
                : _ratpoly(ratpoly)
                , _sample_time(sample_time){}

            constexpr explicit DiscreteTransferFunction(
                const num_vector_type& num, 
                const den_vector_type& den, 
                const ValueType& sample_time)
                : _ratpoly(num, den)
                , _sample_time(sample_time){}

            constexpr explicit DiscreteTransferFunction(
                const value_type(&num)[NumSize], 
                const value_type(&den)[DenSize], 
                const ValueType& sample_time)
                : _ratpoly(num, den)
                , _sample_time(sample_time){}

            constexpr ratpoly_type& ratpoly() {return this->_ratpoly;}
            constexpr const ratpoly_type& ratpoly() const {return this->_ratpoly;}

            constexpr num_type& num() {return this->_ratpoly.num();}
            constexpr const num_type& num() const {return this->_ratpoly.num();}

            constexpr num_type& den() {return this->_ratpoly.den();}
            constexpr const num_type& den() const {return this->_ratpoly.den();}

            constexpr value_type& sample_time() {return this->_ratpoly._time();}
            constexpr const value_type& sample_time() const {return this->_ratpoly._time();}

            friend std::ostream& operator<<(std::ostream& stream, const DiscreteTransferFunction& dtf){
                dtf.ratpoly().print(stream, "z");
                stream << "Sample time: " << dtf.sample_time();
                return stream;
            }
    };

    // operator +
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator+(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize2, DenSize2>& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const Tscalar& lhs, const DiscreteTransferFunction<Tpoly, NumSize, DenSize>& rhs){
        return DiscreteTransferFunction(lhs + rhs.ratpoly());
    }

    template<class Tpoly, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator+(const DiscreteTransferFunction<Tpoly, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() + rhs);
    }

    // operator -
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator-(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize2, DenSize2>& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() - rhs.ratpoly());
    }

    template<class T, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const Tscalar& lhs, const DiscreteTransferFunction<T, NumSize, DenSize>& rhs){
        return DiscreteTransferFunction(lhs + rhs.ratpoly());
    }

    template<class T, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator-(const DiscreteTransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() - rhs);
    }

    // operator *
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator*(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize2, DenSize2>& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() * rhs.ratpoly());
    }

    template<class T, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const Tscalar& lhs, const DiscreteTransferFunction<T, NumSize, DenSize>& rhs){
        return DiscreteTransferFunction(lhs * rhs.ratpoly());
    }

    template<class T, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator*(const DiscreteTransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() * rhs);
    }

    // operator /
    // -----------

    template<class T, size_t NumSize1, size_t DenSize1, size_t NumSize2, size_t DenSize2>
    constexpr auto operator/(const DiscreteTransferFunction<T, NumSize1, DenSize1>& lhs, const DiscreteTransferFunction<T, NumSize2, DenSize2>& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() / rhs.ratpoly());
    }

    template<class T, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const Tscalar& lhs, const DiscreteTransferFunction<T, NumSize, DenSize>& rhs){
        return DiscreteTransferFunction(lhs / rhs.ratpoly());
    }

    template<class T, class Tscalar, size_t NumSize, size_t DenSize>
    constexpr auto operator/(const DiscreteTransferFunction<T, NumSize, DenSize>& lhs, const Tscalar& rhs){
        return DiscreteTransferFunction(lhs.ratpoly() / rhs);
    }

    namespace tf{
        template<class ValueType=double>
        constexpr DiscreteTransferFunction<ValueType, 2, 1> z(const ValueType& sample_time){
            return DiscreteTransferFunction<ValueType, 2, 1>({ValueType(0), ValueType(1)}, {ValueType(1)}, sample_time);
        }
    }

} // namespace controlpp

namespace Eigen {
    template<typename T, size_t N, size_t D>
    struct NumTraits<controlpp::DiscreteTransferFunction<T, N, D>> {
        using Self = controlpp::DiscreteTransferFunction<T, N, D>;
        using Real = Self; // or T if you want to treat real parts differently
        using NonInteger = Self;
        using Literal = Self;
        using Nested = Self;

        enum {
            IsComplex = 0,
            IsInteger = 0,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost = 5,
            MulCost = 10
        };

        static inline Self epsilon() { return Self(); }
        static inline Self dummy_precision() { return Self(); }
        static inline Self highest() { return Self(); }
        static inline Self lowest() { return Self(); }
    };
}