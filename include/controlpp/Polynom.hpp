#pragma once

// std
#include <cstdint>
#include <ostream>
#include <initializer_list>
#include <iterator>
#include <concepts>
#include <complex>

// eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include "math.hpp"

namespace controlpp
{

    /**
     * \brief Describes a mathematical polynomial
     * 
     * Describes polynomials in the form of:
     *  
     * \f[
     * a_1 + a_2 x + a_3 x^2 .. 
     * \f]
     * 
     * , where indices correspond to the parameter number and the potence
     * 
     * \tparam T The datatype of the polynomial parameters (e.g: `float`, `double`, `complex`)
     * \tparam N The size of the polynomial, aka. number of parameters, aka. max order plus one
     * \tparam AutoGrow If true for at least one polynomial of an operation the resulting polynomial will 
     * grow automatically in size to contain all values. Usually set `false` if one needs to sum over polynomials in a loop.
     */
    template<class T, size_t N>
    class Polynom{
        public:
            using value_type = T;
            using vector_type = Eigen::Vector<T, N>;

        private:
            vector_type _vector;

        public:

            /// @brief decault constructs the polynomials with the default values of their data type
            constexpr Polynom() = default;

            /// @brief Copy constructor
            constexpr Polynom(const Polynom&) = default;

            template<size_t M>
            requires(M < N)
            constexpr Polynom(const Polynom<T, M>& other){
                size_t i = 0;
                for(; i < M; ++i){
                    this->at(i) = other.at(i);
                }
                for(; i < N; ++i){
                    this->at(i) = static_cast<T>(0);
                }
            }

            /// @brief Copy assignment operator
            constexpr Polynom& operator=(const Polynom&) = default;

            template<size_t M>
            requires(M < N)
            constexpr Polynom& operator=(const Polynom<T, M>& other){
                size_t i = 0;
                for(; i < M; ++i){
                    this->at(i) = other.at(i);
                }
                for(; i < N; ++i){
                    this->at(i) = static_cast<T>(0);
                }
                return *this;
            }

            /// @brief Constructs a polynomial from an array
            /// @param values the values to be assigned to the polynomial
            constexpr explicit Polynom(const T (&values)[N]){
                for(size_t i = 0; i < N; ++i){
                    this->_vector[i] = values[i];
                }
            }

            template<size_t M>
            requires(M < N)
            constexpr explicit Polynom(const T (&values)[M]){
                size_t i = 0;
                for(; i < M; ++i){
                    this->_vector[i] = values[i];
                } 
                for(; i < N; ++i){
                    this->_vector[i] = static_cast<T>(0);
                } 
            }
            

            /// @brief Constructs a polynomial from a vector
            /// @param vector Eigen::Vector object
            constexpr explicit Polynom(const Eigen::Vector<T, N>& vector) : _vector(vector){}

            template<size_t M>
            requires(M < N)
            constexpr explicit Polynom(const Eigen::Vector<T, M>& vector){
                size_t i = 0;
                for(; i < M; ++i){
                    this->_vector[i] = vector(i);
                } 
                for(; i < N; ++i){
                    this->_vector[i] = static_cast<T>(0);
                } 
            }
            
            constexpr Polynom& operator=(const T (&values)[N]) {
                this->_vector = values;
                return *this;
            }

            template<size_t M>
            requires(M < N)
            constexpr Polynom& operator=(const T (&values)[M]){
                size_t i = 0;
                for(; i < M; ++i){
                    this->_vector[i] = values[i];
                } 
                for(; i < N; ++i){
                    this->_vector[i] = static_cast<T>(0);
                } 
                return *this;
            }
            
            /// @brief Assigns the values of a vector to this polynomial
            /// @param vector An Eigen::Vector holding the values to be assigned to this polynomial
            constexpr Polynom& operator=(const Eigen::Vector<T, N>& vector){
                this->_vector = vector;
                return *this;
            }

            template<size_t M>
            constexpr Polynom& operator=(const Eigen::Vector<T, M>& vector){
                size_t i = 0;
                for(; i < M; ++i){
                    this->_vector[i] = vector(i);
                } 
                for(; i < N; ++i){
                    this->_vector[i] = static_cast<T>(0);
                } 
                return *this;
            }

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return an immutable reference to the element
            constexpr const T& operator[](size_t i) const {return this->_vector[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return a mutable reference to the element
            constexpr T& operator[](size_t i) {return this->_vector[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return an immutable reference to the element
            constexpr const T& at(size_t i) const {return this->_vector[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return a mutable reference to the element
            constexpr T& at(size_t i) {return this->_vector[i];}

            /// @brief Returns the underlying vector that holds the values
            /// @return A mutable reference to an Eigen::Vector object holding the values
            constexpr vector_type& vector(){return this->_vector;}

            /// @brief Returns the underlying vector that holds the values
            /// @return An immutable reference to an Eigen::Vector object holding the values
            constexpr const vector_type& vector() const {return this->_vector;}

            /// @brief returns the size of the polynomial
            /// @details Note that the size is one larger than the order of the polynomial (assuming non-zero entries)
            /// @return an size_teger value holding the size
            constexpr size_t size() const {return N;}

            /// @brief returns the order of the polynomial
            /// @details the order can maxially be one smaller then the size of the polynomial. Takes zero data entries size_to account.
            /// @return an size_teger value holding the order of the polynomial
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

            /// @brief sets all entries to zero
            constexpr void setZero() {this->_vector.setZero();}
            
            /// @brief prsize_ts polynomial to an output character stream with a variale
            /// @param stream the stream to be prsize_ted to
            /// @param var the symbol name of the variable
            void prsize_t (std::ostream& stream, std::string_view var="x") const {
                std::string_view plus = "";
                for(size_t i = 0; i < this->size(); ++i){
                    if(this->at(i) != static_cast<T>(0)){
                        stream << plus << this->at(i) << ' ' << var << '^' << i;
                        plus = " + ";
                    }
                }
            }

            /// @brief Prsize_ts the polynomial (pretty) to an output stream with `"x"` as the symbol name of the variable
            /// @param stream a reference to the output stream prsize_ted to
            /// @param poly the polynomial to be prsize_ted
            /// @return retuns the reference of the output stream
            friend std::ostream& operator<< (std::ostream& stream, const Polynom& poly){
                poly.prsize_t(stream, "x");
                return stream;
            }

            template<size_t M>
            requires(M <= N)
            Polynom& operator+=(const Polynom<T, M>& other){
                for(size_t i = 0; i < M; ++i){
                    this->at(i) += other.at(i);
                }
                return *this;
            }

            template<size_t M>
            requires(M <= N)
            Polynom& operator-=(const Polynom<T, M>& other){
                for(size_t i = 0; i < M; ++i){
                    this->at(i) -= other.at(i);
                }
                return *this;
            }

            template<std::convertible_to<T> U>
            Polynom& operator*=(const U& num){
                this->vector() *= num;
                return *this;
            }

            template<std::convertible_to<T> U>
            Polynom& operator/=(const U& num){
                this->vector() /= num;
                return *this;
            }
    };

    // -----------------------------------------------------------------------------------------------
    //                                        Roots, Zeros
    // -----------------------------------------------------------------------------------------------
    
    /**
     * \brief calculates the zeros/roots of a polynomial
     * 
     * Specialisation for polynomials of degree 1.
     * 
     * Solves the following equation:
     * 
     * \f[
     * 0 = a_0 + a_1 x
     * \f]
     * 
     * \returns an eigen vector of complex numbers which are the roots of the polynomial
     */
    template<class T>
    constexpr Eigen::Vector<std::complex<T>, 1> zeros(const Polynom<T, 2>& polynom){
        const std::complex x0 = - polynom[0] / polynom[1];
        Eigen::Vector<std::complex<T>, 1> result(x0);
        return result;
    }

    /**
     * \brief calculates the zeros/roots of a polynomial
     * 
     * Specialisation for polynomials of degree 1.
     * 
     * \f[
     * 0 = a_0 + a_1 x + a_2 x^2
     * \f]
     * 
     * \returns an eigen vector of complex numbers which are the roots of the polynomial
     */
    template<class T>
    constexpr Eigen::Vector<std::complex<T>, 2> zeros(const Polynom<T, 3>& polynom){
        const std::complex c = polynom[0];
        const std::complex b = polynom[1];
        const std::complex a = polynom[2];
        const std::complex x1 = (-b - std::sqrt(b * b - static_cast<T>(4) * a * c)) / (static_cast<T>(2) * a);
        const std::complex x2 = (-b + std::sqrt(b * b - static_cast<T>(4) * a * c)) / (static_cast<T>(2) * a);
        const Eigen::Vector<std::complex<T>, 2> result(x1, x2);
        return result;
    }

    /**
     * \brief Calculates the zeros of a polynomial
     * 
     * Uses the companion matrix to solve for the zeros of the polynomial
     * 
     */
    template<class T, size_t N>
    requires(N > 1)
    constexpr Eigen::Vector<std::complex<T>, N-1> zeros(const Polynom<T, N>& polynom){
        const Eigen::Matrix<T, N-1, N-1> C = controlpp::companion(polynom.vector());
        const Eigen::Vector<std::complex<T>, N-1> result = C.eigenvalues();
        return result;
    }

    // -----------------------------------------------------------------------------------------------
    //                                      Comparison Operators
    // -----------------------------------------------------------------------------------------------

    /// @brief Compares two polynomials for equality
    /// @tparam T The datatype of the polynomials
    /// @tparam N The size of the polynomials
    /// @param lhs the left-hand-side polynomial of the comparison
    /// @param rhs the right-hand-sode polynomial of the comparison
    /// @return a boolean value that is true if both polynomials have the same parameters
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
            result.vector().tail(Nl - Nr) = lhs.vector().tail(Nl - Nr);
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

    template<class T, size_t N>
    constexpr Polynom<T, N> operator-(const Polynom<T, N>& poly){
        return Polynom<T, N>(-poly.vector());
    }

    template<class T, size_t Nl, size_t Nr>
    constexpr Polynom<T, (Nl > Nr) ? Nl : Nr> operator-(const Polynom<T, Nl>& lhs, const Polynom<T, Nr>& rhs){
        if constexpr (Nl > Nr){
            Polynom<T, Nl> result;
            result.vector().head(Nr) = lhs.vector().head(Nr) - rhs.vector();
            result.vector().tail(Nl - Nr) = lhs.vector().tail(Nl - Nr);
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

    // operator *
    // ----------

    template<class T, size_t N1, size_t N2>
    constexpr Polynom<T, N1+N2-1> operator*(const Polynom<T, N1>& lhs, const Polynom<T, N2>& rhs){
        Polynom<T, N1+N2-1> result;
        result.setZero();

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

    
// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------

    /**
     * \brief Describes a mathematical polynomial of fixed size
     * 
     * Does not grow in size like its sibling `controlpp::Polynom` as a result of operations.
     * 
     * Describes polynomials in the form of:
     *  
     * \f[
     * a_1 + a_2 x + a_3 x^2 .. 
     * \f]
     * 
     * , where indices correspond to the parameter number and the potence
     * 
     * \tparam T The datatype of the polynomial parameters (e.g: `float`, `double`, `complex`)
     * \tparam N The size of the polynomial, aka. number of parameters, aka. max order plus one
     * \tparam AutoGrow If true for at least one polynomial of an operation the resulting polynomial will 
     * grow automatically in size to contain all values. Usually set `false` if one needs to sum over polynomials in a loop.
     */
    template<class T, size_t N>
    class FixedPolynom{
        public:
            using value_type = T;
            using vector_type = Eigen::Vector<T, N>;

        private:
            vector_type _vector;

        public:

            /// @brief decault constructs the polynomials with the default values of their data type
            constexpr FixedPolynom() = default;

            /// @brief Copy constructor
            constexpr FixedPolynom(const FixedPolynom&) = default;

            template<std::convertible_to<T> U>
            constexpr FixedPolynom(const U& value) : _vector(Eigen::Vector<T, N>::Zero()){
                _vector(0) = static_cast<T>(value);
            }

            template<size_t M>
            requires(M < N)
            constexpr FixedPolynom(const FixedPolynom<T, M>& other){
                size_t i = 0;
                for(; i < M; ++i){
                    this->at(i) = other.at(i);
                }
                for(; i < N; ++i){
                    this->at(i) = static_cast<T>(0);
                }
            };

            /// @brief Copy assignment operator
            constexpr FixedPolynom& operator=(const FixedPolynom&) = default;

            template<size_t M>
            requires(M < N)
            constexpr FixedPolynom& operator=(const FixedPolynom<T, M>& other){
                size_t i = 0;
                for(; i < M; ++i){
                    this->at(i) = other.at(i);
                }
                for(; i < N; ++i){
                    this->at(i) = static_cast<T>(0);
                }
                return *this;
            };


            /// @brief Constructs a polynomial from an array
            /// @param values the values to be assigned to the polynomial
            constexpr explicit FixedPolynom(const T (&values)[N]){
                for(size_t i = 0; i < N; ++i){
                    this->_vector[i] = values[i];
                }
            }

            template<size_t M>
            requires(M < N)
            constexpr explicit FixedPolynom(const T (&values)[M]){
                size_t i = 0;
                for(; i < M; ++i){
                    this->_vector[i] = values[i];
                }
                for(; i < N; ++i){
                    this->_vector[i] = static_cast<T>(0);
                }
            }

            /// @brief Constructs a polynomial from a vector
            /// @param vector Eigen::Vector object
            constexpr explicit FixedPolynom(const Eigen::Vector<T, N>& vector) : _vector(vector){}

            template<size_t M>
            requires(M < N)
            constexpr explicit FixedPolynom(const Eigen::Vector<T, M>& vector){
                size_t i = 0;
                for(; i < M; ++i){
                    this->_vector[i] = vector[i];
                }
                for(; i < N; ++i){
                    this->_vector[i] = static_cast<T>(0);
                }
            }
            
            constexpr FixedPolynom& operator=(const T (&values)[N]) {
                this->_vector = values;
                return *this;
            }

            template<size_t M>
            requires(M < N)
            constexpr FixedPolynom& operator=(const T (&values)[M]) {
                size_t i = 0;
                for(; i < M; ++i){
                    this->_vector[i] = values[i];
                }
                for(; i < N; ++i){
                    this->_vector[i] = static_cast<T>(0);
                }
                return *this;
            }
            
            constexpr explicit FixedPolynom(const Polynom<T, N>& poly) : FixedPolynom(poly.vector()){}

            template<size_t M>
            requires(M < N)
            constexpr explicit FixedPolynom(const Polynom<T, N>& poly){
                size_t i = 0;
                for(; i < M; ++i){
                    this->_vector[i] = poly[i];
                }
                for(; i < N; ++i){
                    this->_vector[i] = static_cast<T>(0);
                }
            }

            constexpr operator Polynom<T, N>() const {return Polynom<T, N>(this->_vector);}

            /// @brief Assigns the values of a vector to this polynomial
            /// @param vector An Eigen::Vector holding the values to be assigned to this polynomial
            constexpr FixedPolynom& operator=(const Eigen::Vector<T, N>& vector){
                this->_vector = vector;
                return *this;
            }

            template<size_t M>
            requires(M < N)
            constexpr FixedPolynom& operator=(const Eigen::Vector<T, M>& vector){
                size_t i = 0;
                for(; i < M; ++i){
                    this->_vector[i] = vector[i];
                }
                for(; i < N; ++i){
                    this->_vector[i] = static_cast<T>(0);
                }
                return *this;
            }

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return an immutable reference to the element
            constexpr const T& operator[](size_t i) const {return this->_vector[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return a mutable reference to the element
            constexpr T& operator[](size_t i) {return this->_vector[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return an immutable reference to the element
            constexpr const T& at(size_t i) const {return this->_vector[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return a mutable reference to the element
            constexpr T& at(size_t i) {return this->_vector[i];}

            /// @brief Returns the underlying vector that holds the values
            /// @return A mutable reference to an Eigen::Vector object holding the values
            constexpr vector_type& vector(){return this->_vector;}

            /// @brief Returns the underlying vector that holds the values
            /// @return An immutable reference to an Eigen::Vector object holding the values
            constexpr const vector_type& vector() const {return this->_vector;}

            /// @brief returns the size of the polynomial
            /// @details Note that the size is one larger than the order of the polynomial (assuming non-zero entries)
            /// @return an size_teger value holding the size
            constexpr size_t size() const {return N;}

            /// @brief returns the order of the polynomial
            /// @details the order can maxially be one smaller then the size of the polynomial. Takes zero data entries size_to account.
            /// @return an size_teger value holding the order of the polynomial
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

            /// @brief sets all entries to zero
            constexpr void setZero() {this->_vector.setZero();}
            
            /// @brief prsize_ts polynomial to an output character stream with a variale
            /// @param stream the stream to be prsize_ted to
            /// @param var the symbol name of the variable
            void prsize_t (std::ostream& stream, std::string_view var="x") const {
                std::string_view plus = "";
                for(size_t i = 0; i < this->size(); ++i){
                    if(this->at(i) != static_cast<T>(0)){
                        stream << plus << this->at(i) << ' ' << var << '^' << i;
                        plus = " + ";
                    }
                }
            }

            /// @brief Prsize_ts the polynomial (pretty) to an output stream with `"x"` as the symbol name of the variable
            /// @param stream a reference to the output stream prsize_ted to
            /// @param poly the polynomial to be prsize_ted
            /// @return retuns the reference of the output stream
            friend std::ostream& operator<< (std::ostream& stream, const FixedPolynom& poly){
                poly.prsize_t(stream, "x");
                return stream;
            }

            template<size_t M>
            requires(M <= N)
            FixedPolynom& operator+=(const FixedPolynom<T, M>& other){
                for(size_t i = 0; i < M; ++i){
                    this->at(i) += other.at(i);
                }
                return *this;
            }

            template<size_t M>
            requires(M <= N)
            FixedPolynom& operator-=(const FixedPolynom<T, M>& other){
                for(size_t i = 0; i < M; ++i){
                    this->at(i) -= other.at(i);
                }
                return *this;
            }
    };

    // -----------------------------------------------------------------------------------------------
    //                                      Comparison Operators
    // -----------------------------------------------------------------------------------------------

    /// @brief Compares two polynomials for equality
    /// @tparam T The datatype of the polynomials
    /// @tparam N The size of the polynomials
    /// @param lhs the left-hand-side polynomial of the comparison
    /// @param rhs the right-hand-sode polynomial of the comparison
    /// @return a boolean value that is true if both polynomials have the same parameters
    template<class T, size_t N>
    constexpr bool operator==(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        return lhs.vector() == rhs.vector();
    }

    template<class T, size_t N>
    constexpr bool operator!=(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        return lhs.vector() != rhs.vector();
    }

    // -----------------------------------------------------------------------------------------------
    //                                      Arithmetic Operators
    // -----------------------------------------------------------------------------------------------

    // operator +
    // ----------

    template<class T, size_t N>
    constexpr FixedPolynom<T, N> operator+(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        return FixedPolynom<T, N>(lhs.vector() + rhs.vector());
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    requires (N >= 1)
    constexpr FixedPolynom<Tpoly, N> operator+(const Tscalar& lhs, const FixedPolynom<Tpoly, N>& rhs){
        FixedPolynom<Tpoly, N> result(rhs);
        result[0] += static_cast<Tpoly>(lhs);
        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    requires (N >= 1)
    constexpr FixedPolynom<Tpoly, N> operator+(const FixedPolynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return (static_cast<Tpoly>(rhs) + lhs);
    }

    // operator -
    // ----------

    template<class T, size_t N>
    constexpr FixedPolynom<T, N> operator-(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        return FixedPolynom<T, N>(lhs.vector() - rhs.vector());
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    requires (N >= 1)
    constexpr FixedPolynom<Tpoly, N> operator-(const Tscalar& lhs, const FixedPolynom<Tpoly, N>& rhs){
        FixedPolynom<Tpoly, N> result(-rhs);
        result[0] += static_cast<Tpoly>(lhs);
        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    requires (N >= 1)
    constexpr FixedPolynom<Tpoly, N> operator-(const FixedPolynom<Tpoly, N>& lhs, const Tscalar& rhs){
        FixedPolynom<Tpoly, N> result(lhs);
        result[0] -= static_cast<Tpoly>(rhs);
        return result;
    }

    template<class T, size_t N>
    constexpr FixedPolynom<T, N> operator-(const FixedPolynom<T, N>& values){
        return FixedPolynom<T, N>(-values.vector());
    }

    // operator *
    // ----------

    /**
     * \brief Non growing multiplication
     * 
     * Note that this multiplication does not grow the result to fit all posible values.
     * The correct size has to be allocated beforehand by the user.
     */
    template<class T, size_t N>
    constexpr FixedPolynom<T, N> operator*(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        FixedPolynom<T, N> result;
        result.setZero();

        for(size_t i = 0; i < N; ++i){
            result.vector().tail(N-i) += lhs.vector().head(N-i) * rhs[i];
        }

        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    constexpr FixedPolynom<Tpoly, N> operator*(const FixedPolynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return FixedPolynom<Tpoly, N>(lhs.vector() * static_cast<Tpoly>(rhs));
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    constexpr FixedPolynom<Tpoly, N> operator*(const Tscalar& lhs, const FixedPolynom<Tpoly, N>& rhs){
        return (rhs * static_cast<Tpoly>(lhs));
    }

    // operator /
    // ----------

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    constexpr FixedPolynom<Tpoly, N> operator/(const FixedPolynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return FixedPolynom<Tpoly, N>(lhs.vector() / static_cast<Tscalar>(rhs));
    }
} // namespace controlpp

namespace Eigen {
    template<typename T, size_t N>
    struct NumTraits<controlpp::FixedPolynom<T, N>> {
        using Self = controlpp::FixedPolynom<T, N>;
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

    // ResultType for operator +
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_sum_op<Tscalar, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() + std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_sum_op<controlpp::FixedPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() + std::declval<Tscalar>());
    };

    template <class Tpoly, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_sum_op<controlpp::FixedPolynom<Tpoly, N>, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() + std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    // ResultType for operator -
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_difference_op<Tscalar, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() - std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_difference_op<controlpp::FixedPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() - std::declval<Tscalar>());
    };

    template <class Tpoly, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_difference_op<controlpp::FixedPolynom<Tpoly, N>, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() - std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    // ResultType for operator *
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_product_op<Tscalar, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() * std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_product_op<controlpp::FixedPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() * std::declval<Tscalar>());
    };

    template <class Tpoly, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_product_op<controlpp::FixedPolynom<Tpoly, N>, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() * std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    // ResultType for operator /
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_quotient_op<Tscalar, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() / std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_quotient_op<controlpp::FixedPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() / std::declval<Tscalar>());
    };

    template <class Tpoly, size_t N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_quotient_op<controlpp::FixedPolynom<Tpoly, N>, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() / std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };
}// Eigen
