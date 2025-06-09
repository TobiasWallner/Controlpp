#pragma once

#include <cstdint>

#include <Eigen/Core>

namespace control
{

    /**
     * \brief Describes a mathematical polynomial like \f$ a_1 + a_2 x + a_3 x^2 .. \f$
     * 
     * Indices correspond to the parameter number and the potence
     */
    template<class T, size_t N>
    class Polynom{
        public:
            using vector_type = Eigen::Matrix<T, N, 1>;

        private:
            vector_type _vector;

        public:

            constexpr Polynom() = default;
            explicit constexpr Polynom(const T (&values)[N]) : _vector(values){}
            explicit constexpr Polynom(const vector_type& vector) : _vector(vector){}
            constexpr Polynom(const Polynom&) = default;
            constexpr Polynom& operator=(const Polynom&) = default;
            constexpr Polynom& operator=(T (&values)[N]) {
                this->_vector = values;
                return *this;
            }
            constexpr Polynom& operator=(const vector_type& vector){
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
    constexpr bool operator==(const Polynom<T, N>& lhs, const Polynom<T, N>& rhs){
        return lhs.vector() == rhs.vector();
    }

    template<class T, size_t N>
    constexpr bool operator!=(const Polynom<T, N>& lhs, const Polynom<T, N>& rhs){
        return lhs.vector() != rhs.vector();
    }

    // -----------------------------------------------------------------------------------------------
    //                                      Arithmetic Operators
    // -----------------------------------------------------------------------------------------------

    template<class T, size_t N>
    constexpr Polynom<T, N> operator+(const Polynom<T, N>& lhs, const Polynom<T, N>& rhs){
        return Polynom<T, N>(lhs.vector() + rhs.vector());
    }

    template<class T, size_t N>
    constexpr Polynom<T, N> operator-(const Polynom<T, N>& lhs, const Polynom<T, N>& rhs){
        return Polynom<T, N>(lhs.vector() - rhs.vector());
    }

    template<class T, size_t N>
    constexpr Polynom<T, N> operator*(const Polynom<T, N>& lhs, const T& rhs){
        const Eigen::Vector<T, N> result = lhs.vector() * rhs;
        return Polynom<T, N>(result);
    }

    template<class T, size_t N>
    constexpr Polynom<T, N> operator*(const T& lhs, const Polynom<T, N>& rhs){
        return (rhs * lhs);
    }

    template<class T, size_t N1, size_t N2>
    constexpr Polynom<T, N1+N2-1> operator*(const Polynom<T, N1>& lhs, const Polynom<T, N2>& rhs){
        Polynom<T, N1+N2-1> result;
        result.setZero();

        for(size_t i = 0; i < N2; ++i){
            result.vector().segment(i, N1) += lhs.vector() * rhs[i];
        }

        return result;
    }

} // namespace control
