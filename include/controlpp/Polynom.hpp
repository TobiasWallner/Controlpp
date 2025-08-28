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
     * \tparam Order The Order of the polynomial
     * \tparam AutoGrow If true for at least one polynomial of an operation the resulting polynomial will 
     * grow automatically in size to contain all values. Usually set `false` if one needs to sum over polynomials in a loop.
     */
    template<class T, int Order>
    class Polynom{
        public:
            using value_type = T;
            using vector_type = Eigen::Vector<T, Order+1>;

        private:
            vector_type vector_;

        public:

            /// @brief decault constructs the polynomials with the default values of their data type
            constexpr Polynom() = default;

            /// @brief Copy constructor
            constexpr Polynom(const Polynom&) = default;

            template<int OtherOrder>
            requires(OtherOrder < Order)
            constexpr Polynom(const Polynom<T, OtherOrder>& other){
                size_t i = 0;
                for(; i < other.size(); ++i){
                    this->at(i) = other.at(i);
                }
                for(; i < this->size(); ++i){
                    this->at(i) = static_cast<T>(0);
                }
            }

            /// @brief Copy assignment operator
            constexpr Polynom& operator=(const Polynom&) = default;

            template<int M>
            requires(M < Order)
            constexpr Polynom& operator=(const Polynom<T, M>& other){
                size_t i = 0;
                for(; i < other.size(); ++i){
                    this->at(i) = other.at(i);
                }
                for(; i < this->size(); ++i){
                    this->at(i) = static_cast<T>(0);
                }
                return *this;
            }



            /// @brief Constructs a polynomial from an array
            /// @param values the values to be assigned to the polynomial
            constexpr explicit Polynom(const T (&values)[Order+1])
                : vector_(values){}

            /// @brief Constructs a polynomial from an array
            /// @param values the values to be assigned to the polynomial
            template<int M>
            requires(M < Order+1)
            constexpr explicit Polynom(const T (&values)[M]){
                size_t i = 0;
                for(; i < M; ++i){
                    this->vector_[i] = values[i];
                } 
                for(; i < this->size(); ++i){
                    this->vector_[i] = static_cast<T>(0);
                } 
            }

            /// @brief Constructs a polynomial from a vector
            /// @param vector Eigen::Vector object
            constexpr explicit Polynom(const Eigen::Vector<T, Order+1>& vector) : vector_(vector){}

            template<class U>
            requires(std::constructible_from<Eigen::Vector<T, Order+1>, U>)
            constexpr explicit Polynom(const U& vector_expression) : vector_(vector_expression){}

            template<int M>
            requires(M < Order+1)
            constexpr explicit Polynom(const Eigen::Vector<T, M>& vector){
                size_t i = 0;
                for(; i < vector.size(); ++i){
                    this->vector_[i] = vector(i);
                } 
                for(; i < this->size(); ++i){
                    this->vector_[i] = static_cast<T>(0);
                } 
            }
            
            constexpr Polynom& operator=(const T (&values)[Order+1]) {
                this->vector_ = values;
                return *this;
            }

            template<int M>
            requires(M < Order+1)
            constexpr Polynom& operator=(const T (&values)[M]){
                size_t i = 0;
                for(; i < M; ++i){
                    this->vector_[i] = values[i];
                } 
                for(; i < this->size(); ++i){
                    this->vector_[i] = static_cast<T>(0);
                } 
                return *this;
            }
            
            /// @brief Assigns the values of a vector to this polynomial
            /// @param vector An Eigen::Vector holding the values to be assigned to this polynomial
            constexpr Polynom& operator=(const Eigen::Vector<T, Order+1>& vector){
                this->vector_ = vector;
                return *this;
            }

            /// @brief Assigns the values of a vector expression to this polynomial
            /// @param vector An Eigen::Vector holding the values to be assigned to this polynomial
            template<class U>
            requires(std::constructible_from<Eigen::Vector<T, Order+1>, U>)
            constexpr Polynom& operator=(const U& vector){
                this->vector_ = vector;
                return *this;
            }

            template<int M>
            requires(M < Order+1)
            constexpr Polynom& operator=(const Eigen::Vector<T, M>& vector){
                size_t i = 0;
                for(; i < vector.size(); ++i){
                    this->vector_[i] = vector(i);
                } 
                for(; i < this->size(); ++i){
                    this->vector_[i] = static_cast<T>(0);
                } 
                return *this;
            }

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return an immutable reference to the element
            constexpr const T& operator[](size_t i) const {return this->vector_[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return a mutable reference to the element
            constexpr T& operator[](size_t i) {return this->vector_[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return an immutable reference to the element
            constexpr const T& at(size_t i) const {return this->vector_[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return a mutable reference to the element
            constexpr T& at(size_t i) {return this->vector_[i];}

            /// @brief Returns the underlying vector that holds the values
            /// @return A mutable reference to an Eigen::Vector object holding the values
            constexpr vector_type& vector(){return this->vector_;}

            /// @brief Returns the underlying vector that holds the values
            /// @return An immutable reference to an Eigen::Vector object holding the values
            constexpr const vector_type& vector() const {return this->vector_;}

            /// @brief returns the size of the polynomial
            /// @details Note that the size is one larger than the order of the polynomial (assuming non-zero entries)
            /// @return an size_teger value holding the size
            constexpr size_t size() const {return this->vector_.size();}

            /// @brief returns the order of the polynomial
            /// @details the order can maxially be one smaller then the size of the polynomial. Takes zero data entries size_to account.
            /// @return an size_teger value holding the order of the polynomial
            constexpr size_t order() const {
                size_t result = Order;
                const T zero(0);
                for(size_t i = 0; i < this->size(); ++i){
                    if(this->at(this->size()-1-i) != zero){
                        break;
                    }
                    --result;
                }
                return result;
            }

            /// @brief sets all entries to zero
            constexpr void setZero() {this->vector_.setZero();}
            
            /**
             * \brief Evaluates the polynomial at position x
             * \returns The result of the polynomial evaluated at the position x
             */
            constexpr T eval(const T& x) const {
                Eigen::Vector<T, Order+1> X;
                T p = static_cast<T>(1);
                for(int i = 0; i < X.size(); ++i){
                    X(i) = p;
                    p *= x;
                }
                const T result = this->vector_.dot(X);
                return result;
            }

            /**
             * \brief Evaluates the polynomial at position `x`
             * \param x the input variable of the polynomial
             * \returns The result of the polynomial evaluated at the position `x`
             */
            constexpr T operator() (const T& x) const {
                return this->eval(x);
            }

            /**
             * \brief Evaluates the polynomial at positions of the vector `x_vec`
             * \param x_vec A vector of input variable at which to evaluate the polynomial
             * \returns The results of the polynomial evaluated at each of the vector elements
             */
            template<int M>
            constexpr Eigen::Vector<T, M> eval(const Eigen::Vector<T, M>& x_vec) const {
                Eigen::Matrix<T, M, Order+1> X;
                Eigen::Vector<T, M> p = Eigen::Vector<T, M>::Ones();
                for(int i = 0; i < X.cols(); ++i){
                    X.col(i) = p;
                    p.array() *= x_vec.array(); // element wise multiplication
                }
                const Eigen::Vector<T, M> result = X * this->vector_;
                return result;
            }

            /**
             * \brief Evaluates the polynomial at positions of the vector `x_vec`
             * \param x_vec A vector of input variable at which to evaluate the polynomial
             * \returns The results of the polynomial evaluated at each of the vector elements
             */
            template<int M>
            constexpr Eigen::Vector<T, M> operator() (const Eigen::Vector<T, M>& x_vec) const {
                return this->eval(x_vec);
            }

            /**
             * \brief Evaluates the polynomial at position x
             * \returns The result of the polynomial evaluated at the position x
             */
            constexpr std::complex<T> eval(const std::complex<T>& x) const {
                Eigen::Vector<std::complex<T>, Order+1> X;
                std::complex<T> p = static_cast<std::complex<T>>(1);
                for(int i = 0; i < X.size(); ++i){
                    X(i) = p;
                    p *= x;
                }
                const std::complex<T> result = this->vector_.dot(X);
                return result;
            }

            /**
             * \brief Evaluates the polynomial at position `x`
             * \param x the input variable of the polynomial
             * \returns The result of the polynomial evaluated at the position `x`
             */
            constexpr std::complex<T> operator() (const std::complex<T>& x) const {
                return this->eval(x);
            }

            /**
             * \brief Evaluates the polynomial at positions of the vector `x_vec`
             * \param x_vec A vector of input variable at which to evaluate the polynomial
             * \returns The results of the polynomial evaluated at each of the vector elements
             */
            template<int M>
            constexpr Eigen::Vector<std::complex<T>, M> eval(const Eigen::Vector<std::complex<T>, M>& x_vec) const {
                Eigen::Matrix<std::complex<T>, M, Order+1> X;
                Eigen::Vector<std::complex<T>, M> p = Eigen::Vector<T, M>::Ones();
                for(int i = 0; i < X.cols(); ++i){
                    X.col(i) = p;
                    p.array() *= x_vec.array(); // element wise multiplication
                }
                const Eigen::Vector<std::complex<T>, M> result = X * this->vector_;
                return result;
            }

            /**
             * \brief Evaluates the polynomial at positions of the vector `x_vec`
             * \param x_vec A vector of input variable at which to evaluate the polynomial
             * \returns The results of the polynomial evaluated at each of the vector elements
             */
            template<int M>
            constexpr Eigen::Vector<std::complex<T>, M> operator() (const Eigen::Vector<std::complex<T>, M>& x_vec) const {
                return this->eval(x_vec);
            }

            /// @brief prints polynomial to an output character stream with a variale
            /// @param stream the stream to be printed to
            /// @param var the symbol name of the variable
            void print (std::ostream& stream, std::string_view var="x") const {
                for(size_t i = 0; i < this->size(); ++i){
                    if(this->at(i) != static_cast<T>(0)){
                        const auto value = this->at(i);
                        if(value > 0){
                            if(i > 0) stream << " + ";
                            
                            stream << value;

                            if(i == 1) stream << ' ' << var;
                            else if(i > 1) stream << ' ' << var << '^' << i;
                        }else if(value < 0){
                            if(i == 0) stream << "- ";
                            else if(i > 0) stream << " - ";

                            stream << (-value);

                            if(i == 1) stream << ' ' << var;
                            else if(i > 1) stream << ' ' << var << '^' << i;
                        }
                    }
                }
            }

            /// @brief prints the polynomial (pretty) to an output stream with `"x"` as the symbol name of the variable
            /// @param stream a reference to the output stream printed to
            /// @param poly the polynomial to be printed
            /// @return retuns the reference of the output stream
            friend std::ostream& operator<< (std::ostream& stream, const Polynom& poly){
                poly.print(stream, "x");
                return stream;
            }

            template<int M>
            requires(M <= Order)
            Polynom& operator+=(const Polynom<T, M>& other){
                for(size_t i = 0; i < other.size(); ++i){
                    this->at(i) += other.at(i);
                }
                return *this;
            }

            template<int M>
            requires(M <= Order)
            Polynom& operator-=(const Polynom<T, M>& other){
                for(size_t i = 0; i < other.size(); ++i){
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
    constexpr Eigen::Vector<std::complex<T>, 1> zeros(const Polynom<T, 1>& polynom){
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
    constexpr Eigen::Vector<std::complex<T>, 2> zeros(const Polynom<T, 2>& polynom){
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
    template<class T, int N>
    requires(N > 1)
    constexpr Eigen::Vector<std::complex<T>, N> zeros(const Polynom<T, N>& polynom){
        const Eigen::Matrix<T, N, N> C = controlpp::companion(polynom.vector());
        const Eigen::Vector<std::complex<T>, N> result = C.eigenvalues();
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
    template<class T, int N>
    constexpr bool operator==(const Polynom<T, N>& lhs, const Polynom<T, N>& rhs){
        return lhs.vector() == rhs.vector();
    }

    template<class T, int N>
    constexpr bool operator!=(const Polynom<T, N>& lhs, const Polynom<T, N>& rhs){
        return lhs.vector() != rhs.vector();
    }

    // -----------------------------------------------------------------------------------------------
    //                                      Arithmetic Operators
    // -----------------------------------------------------------------------------------------------

    // operator +
    // ----------

    template<class T, int Nl, int Nr>
    constexpr Polynom<T, (Nl > Nr) ? Nl : Nr> operator+(const Polynom<T, Nl>& lhs, const Polynom<T, Nr>& rhs){
        if constexpr (Nl > Nr){
            Polynom<T, Nl> result;
            result.vector().head(rhs.size()) = lhs.vector().head(rhs.size()) + rhs.vector();
            result.vector().tail(lhs.size() - rhs.size()) = lhs.vector().tail(lhs.size() - rhs.size());
            return result;
        }else if constexpr (Nl == Nr){
            return Polynom<T, Nl>(lhs.vector() + rhs.vector());
        }else{
            Polynom<T, Nr> result;
            result.vector().head(lhs.size()) = lhs.vector() + rhs.vector().head(lhs.size());
            result.vector().tail(rhs.size() - lhs.size()) = rhs.vector().tail(rhs.size() - lhs.size());
            return result;
        }
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    constexpr Polynom<Tpoly, N> operator+(const Tscalar& lhs, const Polynom<Tpoly, N>& rhs){
        Polynom<Tpoly, N> result(rhs);
        result[0] += static_cast<Tpoly>(lhs);
        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    constexpr Polynom<Tpoly, N> operator+(const Polynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return (static_cast<Tpoly>(rhs) + lhs);
    }

    // operator -
    // ----------

    template<class T, int N>
    constexpr Polynom<T, N> operator-(const Polynom<T, N>& poly){
        return Polynom<T, N>(-poly.vector());
    }

    template<class T, int Nl, int Nr>
    constexpr Polynom<T, (Nl > Nr) ? Nl : Nr> operator-(const Polynom<T, Nl>& lhs, const Polynom<T, Nr>& rhs){
        if constexpr (Nl > Nr){
            Polynom<T, Nl> result;
            result.vector().head(rhs.size()) = lhs.vector().head(rhs.size()) - rhs.vector();
            result.vector().tail(lhs.size() - rhs.size()) = lhs.vector().tail(lhs.size() - rhs.size());
            return result;
        }else if constexpr (Nl == Nr){
            return Polynom<T, Nl>(lhs.vector() - rhs.vector());
        }else{
            Polynom<T, Nr> result;
            result.vector().head(lhs.size()) = lhs.vector() + rhs.vector().head(lhs.size());
            result.vector().tail(rhs.size() - lhs.size()) = -rhs.vector().tail(rhs.size() - lhs.size());
            return result;
        }
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    constexpr Polynom<Tpoly, N> operator-(const Tscalar& lhs, const Polynom<Tpoly, N>& rhs){
        Polynom<Tpoly, N> result(-rhs);
        result[0] += static_cast<Tpoly>(lhs);
        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    constexpr Polynom<Tpoly, N> operator-(const Polynom<Tpoly, N>& lhs, const Tscalar& rhs){
        Polynom<Tpoly, N> result(lhs);
        result[0] -= static_cast<Tpoly>(rhs);
        return result;
    }

    // operator *
    // ----------

    template<class T, int lOrder, int rOrder>
    constexpr Polynom<T, lOrder + rOrder> operator*(const Polynom<T, lOrder>& lhs, const Polynom<T, rOrder>& rhs){
        Polynom<T, lOrder + rOrder> result;
        result.setZero();

        for(size_t i = 0; i < rhs.size(); ++i){
            result.vector().segment(i, lhs.size()) += lhs.vector() * rhs[i];
        }

        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int Order>
    constexpr Polynom<Tpoly, Order> operator*(const Polynom<Tpoly, Order>& lhs, const Tscalar& rhs){
        Polynom<Tpoly, Order> result(lhs.vector() * static_cast<Tpoly>(rhs));
        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int Order>
    constexpr Polynom<Tpoly, Order> operator*(const Tscalar& lhs, const Polynom<Tpoly, Order>& rhs){
        return (rhs * static_cast<Tpoly>(lhs));
    }

    // operator /
    // ----------

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int Order>
    constexpr Polynom<Tpoly, Order> operator/(const Polynom<Tpoly, Order>& lhs, const Tscalar& rhs){
        return Polynom<Tpoly, Order>(lhs.vector() / static_cast<Tscalar>(rhs));
    }

    
    namespace polynom{

        template<class T>
        inline const Polynom<T, 1> x({0, 1});
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
    template<class T, int Order>
    class FixedPolynom{
        public:
            using value_type = T;
            using vector_type = Eigen::Vector<T, Order+1>;

        private:
            vector_type vector_;

        public:

            /// @brief decault constructs the polynomials with the default values of their data type
            constexpr FixedPolynom() = default;

            /// @brief Copy constructor
            constexpr FixedPolynom(const FixedPolynom&) = default;

            template<std::convertible_to<T> U>
            constexpr FixedPolynom(const U& value) : vector_(Eigen::Vector<T, Order+1>::Zero()){
                vector_(0) = static_cast<T>(value);
            }

            template<int M>
            requires(M < Order)
            constexpr FixedPolynom(const FixedPolynom<T, M>& other){
                size_t i = 0;
                for(; i < other.size(); ++i){
                    this->at(i) = other.at(i);
                }
                for(; i < this->size(); ++i){
                    this->at(i) = static_cast<T>(0);
                }
            };

            /// @brief Copy assignment operator
            constexpr FixedPolynom& operator=(const FixedPolynom&) = default;

            template<int M>
            requires(M < Order)
            constexpr FixedPolynom& operator=(const FixedPolynom<T, M>& other){
                size_t i = 0;
                for(; i < other.size(); ++i){
                    this->at(i) = other.at(i);
                }
                for(; i < this->size(); ++i){
                    this->at(i) = static_cast<T>(0);
                }
                return *this;
            };


            /// @brief Constructs a polynomial from an array
            /// @param values the values to be assigned to the polynomial
            constexpr explicit FixedPolynom(const T (&values)[Order+1]){
                for(size_t i = 0; i < this->size(); ++i){
                    this->vector_[i] = values[i];
                }
            }

            template<int M>
            requires(M < Order+1)
            constexpr explicit FixedPolynom(const T (&values)[M]){
                size_t i = 0;
                for(; i < M; ++i){
                    this->vector_[i] = values[i];
                }
                for(; i < this->size(); ++i){
                    this->vector_[i] = static_cast<T>(0);
                }
            }

            /// @brief Constructs a polynomial from a vector
            /// @param vector Eigen::Vector object
            constexpr explicit FixedPolynom(const Eigen::Vector<T, Order+1>& vector) : vector_(vector){}

            template<int M>
            requires(M < Order+1)
            constexpr explicit FixedPolynom(const Eigen::Vector<T, M>& vector){
                size_t i = 0;
                for(; i < vector.size(); ++i){
                    this->vector_[i] = vector[i];
                }
                for(; i < this->size(); ++i){
                    this->vector_[i] = static_cast<T>(0);
                }
            }
            
            constexpr FixedPolynom& operator=(const T (&values)[Order+1]) {
                this->vector_ = values;
                return *this;
            }

            template<int M>
            requires(M < Order+1)
            constexpr FixedPolynom& operator=(const T (&values)[M]) {
                size_t i = 0;
                for(; i < M; ++i){
                    this->vector_[i] = values[i];
                }
                for(; i < this->size(); ++i){
                    this->vector_[i] = static_cast<T>(0);
                }
                return *this;
            }
            
            constexpr explicit FixedPolynom(const Polynom<T, Order>& poly) : FixedPolynom(poly.vector()){}

            template<int M>
            requires(M < Order)
            constexpr explicit FixedPolynom(const Polynom<T, M>& poly){
                size_t i = 0;
                for(; i < M; ++i){
                    this->vector_[i] = poly[i];
                }
                for(; i < this->size(); ++i){
                    this->vector_[i] = static_cast<T>(0);
                }
            }

            constexpr operator Polynom<T, Order>() const {return Polynom<T, Order>(this->vector_);}

            /// @brief Assigns the values of a vector to this polynomial
            /// @param vector An Eigen::Vector holding the values to be assigned to this polynomia
            constexpr FixedPolynom& operator=(const Eigen::Vector<T, Order+1>& vector){
                this->vector_ = vector;
                return *this;
            }

            template<int M>
            requires(M < Order+1)
            constexpr FixedPolynom& operator=(const Eigen::Vector<T, M>& vector){
                size_t i = 0;
                for(; i < M; ++i){
                    this->vector_[i] = vector[i];
                }
                for(; i < this->size(); ++i){
                    this->vector_[i] = static_cast<T>(0);
                }
                return *this;
            }

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return an immutable reference to the element
            constexpr const T& operator[](size_t i) const {return this->vector_[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return a mutable reference to the element
            constexpr T& operator[](size_t i) {return this->vector_[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return an immutable reference to the element
            constexpr const T& at(size_t i) const {return this->vector_[i];}

            /// @brief Access elements at the i-th position
            /// @param i the value index to be accessed. Indices correspond to the order of the parameter
            /// @return a mutable reference to the element
            constexpr T& at(size_t i) {return this->vector_[i];}

            /// @brief Returns the underlying vector that holds the values
            /// @return A mutable reference to an Eigen::Vector object holding the values
            constexpr vector_type& vector(){return this->vector_;}

            /// @brief Returns the underlying vector that holds the values
            /// @return An immutable reference to an Eigen::Vector object holding the values
            constexpr const vector_type& vector() const {return this->vector_;}

            /// @brief returns the size of the polynomial
            /// @details Note that the size is one larger than the order of the polynomial (assuming non-zero entries)
            /// @return an size_teger value holding the size
            constexpr size_t size() const {return Order+1;}

            /// @brief returns the order of the polynomial
            /// @details the order can maxially be one smaller then the size of the polynomial. Takes zero data entries size_to account.
            /// @return an size_teger value holding the order of the polynomial
            constexpr size_t order() const {
                size_t result = Order;
                const T zero(0);
                for(size_t i = 0; i < this->size(); ++i){
                    if(this->at(this->size() - i - 1) != zero){
                        break;
                    }
                    --result;
                }
                return result;
            }

            /// @brief sets all entries to zero
            constexpr void setZero() {this->vector_.setZero();}
            
            /// @brief prints polynomial to an output character stream with a variale
            /// @param stream the stream to be printed to
            /// @param var the symbol name of the variable
            void print (std::ostream& stream, std::string_view var="x") const {
                std::string_view plus = "";
                for(size_t i = 0; i < this->size(); ++i){
                    if(this->at(i) != static_cast<T>(0)){
                        stream << plus << this->at(i) << ' ' << var << '^' << i;
                        plus = " + ";
                    }
                }
            }

            /// @brief prints the polynomial (pretty) to an output stream with `"x"` as the symbol name of the variable
            /// @param stream a reference to the output stream printed to
            /// @param poly the polynomial to be printed
            /// @return retuns the reference of the output stream
            friend std::ostream& operator<< (std::ostream& stream, const FixedPolynom& poly){
                poly.print(stream, "x");
                return stream;
            }

            template<int M>
            requires(M <= Order)
            FixedPolynom& operator+=(const FixedPolynom<T, M>& other){
                for(size_t i = 0; i < other->size(); ++i){
                    this->at(i) += other.at(i);
                }
                return *this;
            }

            template<int M>
            requires(M <= Order)
            FixedPolynom& operator-=(const FixedPolynom<T, M>& other){
                for(size_t i = 0; i < other->size(); ++i){
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
    template<class T, int N>
    constexpr bool operator==(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        return lhs.vector() == rhs.vector();
    }

    template<class T, int N>
    constexpr bool operator!=(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        return lhs.vector() != rhs.vector();
    }

    // -----------------------------------------------------------------------------------------------
    //                                      Arithmetic Operators
    // -----------------------------------------------------------------------------------------------

    // operator +
    // ----------

    template<class T, int N>
    constexpr FixedPolynom<T, N> operator+(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        return FixedPolynom<T, N>(lhs.vector() + rhs.vector());
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    requires (N >= 1)
    constexpr FixedPolynom<Tpoly, N> operator+(const Tscalar& lhs, const FixedPolynom<Tpoly, N>& rhs){
        FixedPolynom<Tpoly, N> result(rhs);
        result[0] += static_cast<Tpoly>(lhs);
        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    requires (N >= 1)
    constexpr FixedPolynom<Tpoly, N> operator+(const FixedPolynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return (static_cast<Tpoly>(rhs) + lhs);
    }

    // operator -
    // ----------

    template<class T, int N>
    constexpr FixedPolynom<T, N> operator-(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        return FixedPolynom<T, N>(lhs.vector() - rhs.vector());
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    requires (N >= 1)
    constexpr FixedPolynom<Tpoly, N> operator-(const Tscalar& lhs, const FixedPolynom<Tpoly, N>& rhs){
        FixedPolynom<Tpoly, N> result(-rhs);
        result[0] += static_cast<Tpoly>(lhs);
        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    requires (N >= 1)
    constexpr FixedPolynom<Tpoly, N> operator-(const FixedPolynom<Tpoly, N>& lhs, const Tscalar& rhs){
        FixedPolynom<Tpoly, N> result(lhs);
        result[0] -= static_cast<Tpoly>(rhs);
        return result;
    }

    template<class T, int N>
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
    template<class T, int N>
    constexpr FixedPolynom<T, N> operator*(const FixedPolynom<T, N>& lhs, const FixedPolynom<T, N>& rhs){
        FixedPolynom<T, N> result;
        result.setZero();

        for(size_t i = 0; i < lhs.size(); ++i){
            result.vector().tail(lhs.size()-i) += lhs.vector().head(lhs.size()-i) * rhs[i];
        }

        return result;
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    constexpr FixedPolynom<Tpoly, N> operator*(const FixedPolynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return FixedPolynom<Tpoly, N>(lhs.vector() * static_cast<Tpoly>(rhs));
    }

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    constexpr FixedPolynom<Tpoly, N> operator*(const Tscalar& lhs, const FixedPolynom<Tpoly, N>& rhs){
        return (rhs * static_cast<Tpoly>(lhs));
    }

    // operator /
    // ----------

    template<class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    constexpr FixedPolynom<Tpoly, N> operator/(const FixedPolynom<Tpoly, N>& lhs, const Tscalar& rhs){
        return FixedPolynom<Tpoly, N>(lhs.vector() / static_cast<Tscalar>(rhs));
    }
} // namespace controlpp

namespace Eigen {
    template<typename T, int N>
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

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_sum_op<Tscalar, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() + std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_sum_op<controlpp::FixedPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() + std::declval<Tscalar>());
    };

    template <class Tpoly, int N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_sum_op<controlpp::FixedPolynom<Tpoly, N>, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() + std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    // ResultType for operator -
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_difference_op<Tscalar, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() - std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_difference_op<controlpp::FixedPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() - std::declval<Tscalar>());
    };

    template <class Tpoly, int N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_difference_op<controlpp::FixedPolynom<Tpoly, N>, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() - std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    // ResultType for operator *
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_product_op<Tscalar, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() * std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_product_op<controlpp::FixedPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() * std::declval<Tscalar>());
    };

    template <class Tpoly, int N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_product_op<controlpp::FixedPolynom<Tpoly, N>, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() * std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    // ResultType for operator /
    // -------------------------

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    struct ScalarBinaryOpTraits<
        Tscalar,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_quotient_op<Tscalar, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<Tscalar>() / std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };

    template <class Tpoly, std::convertible_to<Tpoly> Tscalar, int N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        Tscalar,
        internal::scalar_quotient_op<controlpp::FixedPolynom<Tpoly, N>, Tscalar>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() / std::declval<Tscalar>());
    };

    template <class Tpoly, int N>
    struct ScalarBinaryOpTraits<
        controlpp::FixedPolynom<Tpoly, N>,
        controlpp::FixedPolynom<Tpoly, N>,
        internal::scalar_quotient_op<controlpp::FixedPolynom<Tpoly, N>, controlpp::FixedPolynom<Tpoly, N>>> 
    {
        using ReturnType = decltype(std::declval<controlpp::FixedPolynom<Tpoly, N>>() / std::declval<controlpp::FixedPolynom<Tpoly, N>>());
    };
}// Eigen
