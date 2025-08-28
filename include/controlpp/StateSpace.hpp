#pragma once

// std
#include <tuple>
#include <concepts>

// eigen
#include <Eigen/Core>

// controlpp
#include <controlpp/math.hpp>
#include <controlpp/TransferFunction.hpp>

namespace controlpp
{

    enum class Domain{s, q, z};
    
    constexpr std::string_view to_string(Domain domain){
        switch(domain){
            case Domain::s : return "s";
            case Domain::q : return "q";
            case Domain::z : return "z";
            default : return "None";
        }
    }

    /**
     * \brief Base class for the state space representation of a linear time invariant system
     * 
     * System describeing the matrices (A, B, C, D) of the following linear system:
     * 
     * \f[
     *      x = \mathbf{A} x + \mathbf{B} x
     *      y = \mathbf{C} x + \mathbf{D} u
     * \f]
     * 
     * Note that this class only stores the A, B, C, and D matrices and **not** the state x.
     * 
     */
    template<class T, int internal_states, int inputs, int outputs>
    class StateSpace{
        public:
            using value_type = T;

            using A_matrix_type = Eigen::Matrix<T, internal_states, internal_states>;
            using B_matrix_type = Eigen::Matrix<T, internal_states, inputs>;
            using C_matrix_type = Eigen::Matrix<T, outputs, internal_states>;
            using D_matrix_type = Eigen::Matrix<T, outputs, inputs>;

            //constexpr static int number_of_states = internal_states;
            //constexpr static int number_of_inputs = inputs;
            //constexpr static int number_of_outputs = outputs;

        private:
            A_matrix_type _A;
            B_matrix_type _B;
            C_matrix_type _C;
            D_matrix_type _D;

        public:
            constexpr StateSpace() = default;
            constexpr StateSpace(const StateSpace&) = default;
            constexpr StateSpace& operator=(const StateSpace&) = default;

            constexpr StateSpace(
                const Eigen::Matrix<T, internal_states, internal_states>& A,
                const Eigen::Matrix<T, internal_states, inputs>& B,
                const Eigen::Matrix<T, outputs, internal_states>& C,
                const Eigen::Matrix<T, outputs, inputs>& D
            )
                : _A(A)
                , _B(B)
                , _C(C)
                , _D(D){}

            /**
             * \brief calculates the new system states and outupts
             * 
             * Usage Example:
             * ```
             * StateSpace ss = some_calculation();
             * Eigen::Vector<float, 2> = some_measurement();
             * auto [new_internal_states, new_output] = ss(internal_states, input);
             * ```
             */
            constexpr std::tuple<Eigen::Vector<T, internal_states>, Eigen::Vector<T, outputs>> eval(const Eigen::Vector<T, internal_states>& x, const Eigen::Vector<T, inputs>& u) const {
                const Eigen::Vector<T, internal_states> result_x = this->A() * x + this->B() * u;
                const Eigen::Vector<T, outputs> result_y = this->C() * x + this->D() * u;
                return std::tuple(result_x, result_y);
            }

            /**
             * \brief calculates the new system states and outupts for SISO (single input, single output) systems
             * 
             * Overload for SI systems that accepts scalar values as inputs.
             * 
             * Usage Example:
             * ```
             * StateSpace ss = some_calculation();
             * float input = some_measurement();
             * auto [new_internal_states, new_output] = ss(internal_states, input);
             * ```
             */
            template<std::same_as<T> U>
                requires(inputs == 1 && outputs != 1)
            constexpr std::tuple<Eigen::Vector<U, internal_states>, Eigen::Vector<U, outputs>> eval(const Eigen::Vector<U, internal_states>& x, const U& u_scalar) const {
                const Eigen::Vector<T, 1> u(u_scalar);
                return this->eval(x, u);
            }

            /**
             * \brief calculates the new system states and outupts for SISO (single input, single output) systems
             * 
             * Overload for SISO systems that accepts scalar values as inputs.
             * 
             * Usage Example:
             * ```
             * StateSpace ss = some_calculation();
             * float input = some_measurement();
             * auto [new_internal_states, new_output] = ss(internal_states, input);
             * ```
             */
            template<std::same_as<T> U>
                requires(inputs == 1 && outputs == 1)
            constexpr std::tuple<Eigen::Vector<U, internal_states>, U> eval(const Eigen::Vector<U, internal_states>& x, const U& u_scalar) const {
                const Eigen::Vector<T, 1> u(u_scalar);
                const auto [new_x, y] = this->eval(x, u);
                return {new_x, y(0)};
            }

            constexpr const A_matrix_type& A() const {return this->_A;}
            constexpr const B_matrix_type& B() const {return this->_B;}
            constexpr const C_matrix_type& C() const {return this->_C;}
            constexpr const D_matrix_type& D() const {return this->_D;}

            constexpr A_matrix_type& A() {return this->_A;}
            constexpr B_matrix_type& B() {return this->_B;}
            constexpr C_matrix_type& C() {return this->_C;}
            constexpr D_matrix_type& D() {return this->_D;}

            friend std::ostream& operator<<(std::ostream& stream, const StateSpace& state_space){
                stream << "A:\n" << state_space.A() << '\n';
                stream << "B:\n" << state_space.B() << '\n';
                stream << "C:\n" << state_space.C() << '\n';
                stream << "D:\n" << state_space.D() << '\n';
                return stream;
            }
    };

    /**
     * \brief calculates the state space representation from a rational polynomial
     */
    template<class T, int NumOrder, int DenOrder>
    constexpr StateSpace<T, DenOrder, 1, 1> to_state_space(const TransferFunction<T, NumOrder, DenOrder>& rp){
        static constexpr int number_of_states = DenOrder;
        StateSpace<T, number_of_states, 1, 1> result;
        const T a_n = rp.den()[rp.den().order()];

        // normalise
        const Polynom<T, DenOrder> a = -(rp.den() / a_n);
        const Polynom<T, NumOrder> b = rp.num() / a_n;
        
        // write A matrix
        if constexpr (number_of_states > 0){
            const auto I = Eigen::Matrix<T, number_of_states-1, number_of_states-1>::Identity();
            result.A().template block<number_of_states-1, number_of_states-1>(0, 1) = I; 
            result.A().col(0).head(number_of_states-1).setZero(); 
            result.A().row(number_of_states-1) = a.vector().head(a.vector().size()-1);
        }

        // write B matrix
        if constexpr (number_of_states > 0){
            result.B().col(0).head(number_of_states-1).setZero();
            result.B()(number_of_states-1, 0) = T(1);
        }
        
        // write C matrix
        if constexpr (number_of_states > 0){
            result.C().row(0).head(b.size()) = b.vector();
            result.C().row(0).tail(number_of_states - b.size()).setZero();
        }
        
        // write D matrix
        result.D()(0, 0) = (b.size() > (number_of_states)) ? b[number_of_states] : T(0);
        return result;
    }


    /**
     * \brief transforms a continuous state space system representation to a transfer function representation
     * 
     * Uses the formular:
     * 
     * \f[
     * G(s) = \mathbf{C} \left( s \mathbf{I} - \mathbf{A} \right)^{-1} \mathbf{B} + D
     * \f]
     */
    template<class T, int states>
    constexpr TransferFunction<T, states+1, states+1> to_transfer_function(const StateSpace<T, states, 1, 1>& css){
        const FixedPolynom<T, states+1> s({0, 1});
        const auto I = Eigen::Matrix<T, states, states>::Identity();
        const Eigen::Matrix<FixedPolynom<T, states+1>, states, states> sI_min_A = s * I - css.A();
        const Eigen::Matrix<FixedPolynom<T, states+1>, states, states> adj_sI_min_A = controlpp::adj(sI_min_A);
        const FixedPolynom<T, states+1> num = (css.C() * adj_sI_min_A * css.B() + css.D())(0, 0);
        const FixedPolynom<T, states+1> den = sI_min_A.determinant();
        return TransferFunction<T, states+1, states+1>(num.vector(), den.vector());
    }

/*
    template<class T, int LStates, int Linputs, int Loutputs, int RStates, int Rinputs, int Routputs>
    StateSpace<...> operator+(const StateSpace<>& lhs, const StateSpace<>& rhs){
        // TODO:    
    }

    template<class T, int LStates, int Linputs, int Loutputs, int RStates, int Rinputs, int Routputs>
    StateSpace<...> operator-(const StateSpace<>& lhs, const StateSpace<>& rhs){
        // TODO:
    }

    template<class T, int LStates, int Linputs, int Loutputs, int RStates, int Rinputs, int Routputs>
    StateSpace<...> operator*(const StateSpace<>& lhs, const StateSpace<>& rhs){
        // TODO:
    }

    template<class T, int LStates, int Linputs, int Loutputs, int RStates, int Rinputs, int Routputs>
    StateSpace<...> operator/(const StateSpace<>& lhs, const StateSpace<>& rhs){
        // TODO:
    }

    template<class T, int LStates, int Linputs, int Loutputs, int RStates, int Rinputs, int Routputs>
    StateSpace<...> inverse(const StateSpace<>& Sys){
        // TODO:
    }
*/
} // namespace controlpp
