#pragma once

// std
#include <cstdint>
#include <ostream>
#include <initializer_list>
#include <iterator>
#include <concepts>

// eigen
#include <Eigen/Core>

namespace controlpp
{

    /**
     * \brief Describes a mathematical polynomial like \f$ a_1 + a_2 x + a_3 x^2 .. \f$
     * 
     * Indices correspond to the parameter number and the potence
     */
    template<class T, size_t N>
    class Polynom{
        public:
            using value_type = T;
            using vector_type = Eigen::Matrix<T, N, 1, Eigen::ColMajor | Eigen::AutoAlign>;

        private:
            vector_type _vector;

        public:

            constexpr Polynom() = default;
            constexpr Polynom(const Polynom&) = default;
            constexpr Polynom& operator=(const Polynom&) = default;

            constexpr explicit Polynom(const T (&values)[N]){
                for(size_t i = 0; i < N; ++i) this->_vector[i] = values[i];
            }

            constexpr explicit Polynom(const vector_type& vector) : _vector(vector){}
            
            constexpr Polynom& operator=(const T (&values)[N]) {
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


            constexpr size_t size() const {return N;}

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

            constexpr void set_zero() {this->_vector.setZero();}
            
            void print (std::ostream& stream, std::string_view var="x") const {
                std::string_view plus = "";
                for(size_t i = 0; i < this->size(); ++i){
                    if(this->at(i) != static_cast<T>(0)){
                        stream << plus << this->at(i) << ' ' << var << '^' << i;
                        plus = " + ";
                    }
                }
            }

            friend std::ostream& operator<< (std::ostream& stream, const Polynom& poly){
                poly.print(stream, "x");
                return stream;
            }
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

    // operator +
    // ----------

    template<class T, size_t Nl, size_t Nr>
    constexpr Polynom<T, (Nl > Nr) ? Nl : Nr> operator+(const Polynom<T, Nl>& lhs, const Polynom<T, Nr>& rhs){
        if constexpr (Nl > Nr){
            Polynom<T, Nl> result;
            result.vector().head(Nr) = lhs.vector().head(Nr) + rhs.vector();
            result.vector().tail(Nl - Nr) = lhs.vector.tail(Nl - Nr);
            return result;
        }else if constexpr (Nl == Nr){
            return Polynom<T, Nl>(lhs.vector() + rhs.vector());
        }else{
            Polynom<T, Nr> result;
            result.vector().head(Nl) = lhs.vector() + rhs.vector().head(Nl);
            result.vector().tail(Nr - Nl) = rhs.vector().tail(Nr - Nl);
            return result;
        }
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    requires (N >= 1)
    constexpr Polynom<Tpoly, N> operator+(const Tscalar& lhs, const Polynom<Tpoly, N>& rhs){
        Polynom<Tpoly, N> result(rhs);
        result[0] += static_cast<Tpoly>(lhs);
        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    requires (N >= 1)
    constexpr Polynom<Tpoly, N> operator+(const Polynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return (static_cast<Tpoly>(rhs) + lhs);
    }

    // operator -
    // ----------

    template<class T, size_t Nl, size_t Nr>
    constexpr Polynom<T, (Nl > Nr) ? Nl : Nr> operator-(const Polynom<T, Nl>& lhs, const Polynom<T, Nr>& rhs){
        if constexpr (Nl > Nr){
            Polynom<T, Nl> result;
            result.vector().head(Nr) = lhs.vector().head(Nr) - rhs.vector();
            result.vector().tail(Nl - Nr) = lhs.vector.tail(Nl - Nr);
            return result;
        }else if constexpr (Nl == Nr){
            return Polynom<T, Nl>(lhs.vector() - rhs.vector());
        }else{
            Polynom<T, Nr> result;
            result.vector().head(Nl) = lhs.vector() + rhs.vector().head(Nl);
            result.vector().tail(Nr - Nl) = -rhs.vector.tail(Nr - Nl);
            return result;
        }
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    requires (N >= 1)
    constexpr Polynom<Tpoly, N> operator-(const Tscalar& lhs, const Polynom<Tpoly, N>& rhs){
        Polynom<Tpoly, N> result(-rhs);
        result[0] += static_cast<Tpoly>(lhs);
        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    requires (N >= 1)
    constexpr Polynom<Tpoly, N> operator-(const Polynom<Tpoly, N>& lhs, const Tscalar& rhs){
        Polynom<Tpoly, N> result(lhs);
        result[0] -= static_cast<Tpoly>(rhs);
        return result;
    }

    template<class T, size_t N>
    constexpr Polynom<T, N> operator-(const Polynom<T, N>& values){
        return Polynom<T, N>(-values.vector());
    }

    // operator *
    // ----------

    template<class T, size_t N1, size_t N2>
    constexpr Polynom<T, N1+N2-1> operator*(const Polynom<T, N1>& lhs, const Polynom<T, N2>& rhs){
        Polynom<T, N1+N2-1> result;
        result.set_zero();

        for(size_t i = 0; i < N2; ++i){
            result.vector().segment(i, N1) += lhs.vector() * rhs[i];
        }

        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    constexpr Polynom<Tpoly, N> operator*(const Polynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return Polynom<Tpoly, N>(lhs.vector() * static_cast<Tpoly>(rhs));
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    constexpr Polynom<Tpoly, N> operator*(const Tscalar& lhs, const Polynom<Tpoly, N>& rhs){
        return (rhs * static_cast<Tpoly>(lhs));
    }

    // operator /
    // ----------

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    constexpr Polynom<Tpoly, N> operator/(const Polynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return Polynom<Tpoly, N>(lhs.vector() / static_cast<Tscalar>(rhs));
    }

    

} // namespace controlpp
