#pragma once

#include <cstdint>

#include <Eigen/Core>

namespace control
{

    /**
     * \brief Describes a mathematical polynomial like $a_1 + a_2 x+ a_3 x^2 ..$
     */
    template<class T, size_t N>
    class Polynomial{
        public:
            using vector_type = Eigen::Matrix<T, N, 1>;

        private:
            vector_type _vector;

        public:

            constexpr Polynomial() = default;
            explicit constexpr Polynomial(const T (&values)[N]) : _vector(values){}
            explicit constexpr Polynomial(const vector_type& vector) : _vector(vector){}
            constexpr Polynomial(const Polynomial&) = default;
            constexpr Polynomial& operator=(const Polynomial&) = default;
            constexpr Polynomial& operator=(T (&values)[N]) {
                this->_vector = values;
                return *this;
            }
            constexpr Polynomial& operator=(const vector_type& vector){
                this->_vector = vector;
                return *this;
            }

            constexpr const T& operator[](size_t i) const {return this->_vector[i];}
            constexpr T& operator[](size_t i) {return this->_vector[i];}

            constexpr const T& at(size_t i) const {return this->_vector[i];}
            constexpr T& at(size_t i) {return this->_vector[i];}

            constexpr vector_type& vector(){return this->_vector;}
            constexpr const vector_type& vector() const {return this->_vector;}

            constexpr size_t order() const {
                size_t result = N-1;
                const T zero(0);
                for(size_t i = 0; i < N; ++i){
                    if(this->at(N-i-1) != zero){
                        break;
                    }
                    --result;
                }
                return result;
            }

            constexpr void setZero() {this->_vector.setZero();}
    };

    // -----------------------------------------------------------------------------------------------
    //                                      Comparison Operators
    // -----------------------------------------------------------------------------------------------

    template<class T, size_t N>
    constexpr bool operator==(const Polynomial<T, N>& lhs, const Polynomial<T, N>& rhs){
        return lhs.vector() == rhs.vector();
    }

    template<class T, size_t N>
    constexpr bool operator!=(const Polynomial<T, N>& lhs, const Polynomial<T, N>& rhs){
        return lhs.vector() != rhs.vector();
    }

    // -----------------------------------------------------------------------------------------------
    //                                      Arithmetic Operators
    // -----------------------------------------------------------------------------------------------

    template<class T, size_t N>
    constexpr Polynomial<T, N> operator+(const Polynomial<T, N>& lhs, const Polynomial<T, N>& rhs){
        return Polynomial<T, N>(lhs.vector() + rhs.vector());
    }

    template<class T, size_t N>
    constexpr Polynomial<T, N> operator-(const Polynomial<T, N>& lhs, const Polynomial<T, N>& rhs){
        return Polynomial<T, N>(lhs.vector() - rhs.vector());
    }

    template<class T, size_t N>
    constexpr Polynomial<T, N> operator*(const Polynomial<T, N>& lhs, const T& rhs){
        const Eigen::Vector<T, N> result = lhs.vector() * rhs;
        return Polynomial<T, N>(result);
    }

    template<class T, size_t N>
    constexpr Polynomial<T, N> operator*(const T& lhs, const Polynomial<T, N>& rhs){
        return (rhs * lhs);
    }

    template<class T, size_t N1, size_t N2>
    constexpr Polynomial<T, N1+N2-1> operator*(const Polynomial<T, N1>& lhs, const Polynomial<T, N2>& rhs){
        Polynomial<T, N1+N2-1> result;
        result.setZero();

        for(size_t i = 0; i < N2; ++i){
            result.vector().segment(i, N1) += lhs.vector() * rhs[i];
        }

        return result;
    }

} // namespace control
