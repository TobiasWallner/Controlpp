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
    template<class T, int NStates, int NInputs, int NOutputs>
    class StateSpace{
        public:
            using value_type = T;

            using A_matrix_type = Eigen::Matrix<T, NStates, NStates>;
            using B_matrix_type = Eigen::Matrix<T, NStates, NInputs>;
            using C_matrix_type = Eigen::Matrix<T, NOutputs, NStates>;
            using D_matrix_type = Eigen::Matrix<T, NOutputs, NInputs>;

            //constexpr static int number_of_states = NStates;
            //constexpr static int number_of_inputs = NInputs;
            //constexpr static int number_of_outputs = NOutputs;

        private:
            A_matrix_type _A;
            B_matrix_type _B;
            C_matrix_type _C;
            D_matrix_type _D;

        public:
            StateSpace() = default;
            StateSpace(const StateSpace&) = default;
            StateSpace& operator=(const StateSpace&) = default;

            StateSpace(
                const Eigen::Matrix<T, NStates, NStates>& A,
                const Eigen::Matrix<T, NStates, NInputs>& B,
                const Eigen::Matrix<T, NOutputs, NStates>& C,
                const Eigen::Matrix<T, NOutputs, NInputs>& D
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
             * auto [new_internal_states, new_output] = ss(NStates, input);
             * ```
             */
            std::tuple<Eigen::Vector<T, NStates>, Eigen::Vector<T, NOutputs>> eval(const Eigen::Vector<T, NStates>& x, const Eigen::Vector<T, NInputs>& u) const {
                const Eigen::Vector<T, NStates> result_x = this->A() * x + this->B() * u;
                const Eigen::Vector<T, NOutputs> result_y = this->C() * x + this->D() * u;
                return std::tuple(result_x, result_y);
            }

            /**
             * \brief calculates the new system states and outupts for SISO (single input, single output) systems
             * 
             * Overload for SI systems that accepts scalar values as NInputs.
             * 
             * Usage Example:
             * ```
             * StateSpace ss = some_calculation();
             * float input = some_measurement();
             * auto [new_internal_states, new_output] = ss(NStates, input);
             * ```
             */
            template<std::same_as<T> U>
                requires(NInputs == 1 && NOutputs != 1)
            std::tuple<Eigen::Vector<U, NStates>, Eigen::Vector<U, NOutputs>> eval(const Eigen::Vector<U, NStates>& x, const U& u_scalar) const {
                const Eigen::Vector<T, 1> u(u_scalar);
                return this->eval(x, u);
            }

            /**
             * \brief calculates the new system states and outupts for SISO (single input, single output) systems
             * 
             * Overload for SISO systems that accepts scalar values as NInputs.
             * 
             * Usage Example:
             * ```
             * StateSpace ss = some_calculation();
             * float input = some_measurement();
             * auto [new_internal_states, new_output] = ss(NStates, input);
             * ```
             */
            template<std::same_as<T> U>
                requires(NInputs == 1 && NOutputs == 1)
            std::tuple<Eigen::Vector<U, NStates>, U> eval(const Eigen::Vector<U, NStates>& x, const U& u_scalar) const {
                const Eigen::Vector<T, 1> u(u_scalar);
                const auto [new_x, y] = this->eval(x, u);
                return {new_x, y(0)};
            }

            const A_matrix_type& A() const {return this->_A;}
            const B_matrix_type& B() const {return this->_B;}
            const C_matrix_type& C() const {return this->_C;}
            const D_matrix_type& D() const {return this->_D;}

            A_matrix_type& A() {return this->_A;}
            B_matrix_type& B() {return this->_B;}
            C_matrix_type& C() {return this->_C;}
            D_matrix_type& D() {return this->_D;}

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
     * 
     * The transfer function:
     * 
     * \f[
     * Tf(s) = \frac{b_0 + b_1 s + \cdots + b_n s^{n}}{a_0 + a_1 s + \cdots + a_n s^{n}}
     * \f]
     * 
     * will be turned into the state space system:
     * 
     * \f[
     * \dot{x} = A x + B u
     * y = C x + D u
     * \f]
     * 
     * with:
     * 
     * \f[
     * \hat{a}_j = a_j / a_n
     * \f]
     * 
     * \f[
     * \hat{b}_j = b_j / a_n
     * \f]
     * 
     * 
     * \f[
     * A = 
     * \begin{bmatrix}
     *  0           & 1             & 0             & 0             & \cdots & 0                \\
     *  0           & 0             & 1             & 0             & \cdots & 0                \\
     *  0           & 0             & 0             & 1             & \cdots & 0                \\
     *  \vdots      & \vdots        & \vdots        & \vdots        & \ddots & \vdots           \\
     *  0           & 0             & 0             & 0             & \cdots & 1                \\
     *  -\hat{a}_0  & -\hat{a}_1    & -\hat{a}_2    & -\hat{a}_3    & \cdots & -\hat{a}_{n-1}
     * \end{bmatrix}
     * \f]
     * 
     * \f[
     * B = 
     * \begin{bmatrix}
     *  0      \\
     *  \vdots \\
     *  1
     * \end{bmatrix}
     * \f]
     * 
     * \f[
     * C = \begin{bmatrix}
     * \hat{b}_0 - \hat{a}_0 \hat{b}_n & b_1 - \hat{a}_1 \hat{b}_n & \cdots & \hat{b}_{n-1} - \hat{a}_{n-1} \hat{b}_n
     * \end{bmatrix}
     * \f]
     * 
     * \f[
     * D = \hat{b}_n
     * \f]
     * 
     */
    template<class T, int NumOrder, int DenOrder>
    StateSpace<T, DenOrder, 1, 1> to_state_space(const TransferFunction<T, NumOrder, DenOrder>& rp){
        constexpr int number_of_states = DenOrder;
        StateSpace<T, number_of_states, 1, 1> result;
        const T a_n = rp.den()[rp.den().order()];

        // normalise
        const Polynom<T, DenOrder> a = -(rp.den() / a_n);
        const Polynom<T, NumOrder> b = rp.num() / a_n;
        
        const T bn = b.at(b.size()-1);

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
            const int l = (b.size() < number_of_states) ? b.size() : number_of_states;
            if(b.size() > (number_of_states)){
                result.C().row(0).head(l) = b.vector().head(l) + a.vector().head(l) * bn;
            }else{
                result.C().row(0).head(l) = b.vector().head(l);
                result.C().row(0).tail(number_of_states - b.size()).setZero();
            }
        }
        
        // write D matrix
        result.D()(0, 0) = (b.size() > (number_of_states)) ? bn : T(0);
        return result;
    }


    /**
     * @brief Generates the block diagonal state space representation of the system of transfer functions
     * 
     * The matrix of transfer functions represensts a multiple input and multiple output system where:  
     *  - the number of rows corresponds to the systems outputs
     *  - and the number of columns corresponds to the systems inputs
     * 
     * @tparam T The data type of the matrix elements and transfer function parameters. (Usually `float` or `double`)
     * @param Mtf A Matrix of transfer functions
     * @return The state space representation of the matrix of transfer functions
     */
    template<class T, int NRows, int NCols, int Opts, int NMaxRows, int NMaxCols>
    StateSpace<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic> 
    to_state_space(const Eigen::Matrix<TransferFunction<T, Eigen::Dynamic, Eigen::Dynamic>, NRows, NCols, Opts, NMaxRows, NMaxCols>& Mtf){
        // count the number of states for matrix pre-allocation
        int states = 0;
        for(int row = 0; row < Mtf.rows(); ++row){
            for(int col = 0; col < Mtf.cols(); ++col){
                states += Mtf.at(row, col).den().order();
            }
        }

        int inputs = Mtf.cols();
        int outputs = Mtf.rows();

        // allocate matrices
        Eigen::Matrix<T, Eigen::Dynamic,  Eigen::Dynamic> A(states, states); A.setZero();
        Eigen::Matrix<T, Eigen::Dynamic,  Eigen::Dynamic> B(states, inputs); B.setZero();
        Eigen::Matrix<T, Eigen::Dynamic,  Eigen::Dynamic> C(outputs, states); C.setZero();
        Eigen::Matrix<T, Eigen::Dynamic,  Eigen::Dynamic> D(outputs, inputs); D.setZero();

        // turn each transfer function into its state space representation and build the total system
        int state_itr = 0;
        for(int row = 0; row < Mtf.rows(); ++row){
            for(int col = 0; col < Mtf.cols(); ++col){
                const StateSpace<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic> ss = to_state_space(Mtf(row, col));
                A.block(state_itr, state_itr, ss.states(), ss.states()) = ss.A();
                B.block(state_itr, col, ss.states(), 1) = ss.B();
                C.block(row, state_itr, 1, ss.states()) = ss.C();
                D.block(row, col, 1, 1) = ss.D();
                state_itr += ss.states();
            }
        }
        return StateSpace<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic>(A, B, C, D);
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
    TransferFunction<T, states+1, states+1> to_transfer_function(const StateSpace<T, states, 1, 1>& css){
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
} // nam