#pragma once

#include <controlpp/StateSpace.hpp>
#include <controlpp/DiscreteTransferFunction.hpp>

namespace controlpp
{
    /**
     * \brief Matrix (A, B, C, D) representation of a linear time invariant system
     * 
     * \f[
     * \hat{x} = \mathbf{A} x + \mathbf{B} u
     * \f]
     * 
     * \f[
     * y = \mathbf{C} x + \mathbf{D} u
     * \f]
     * 
     * - u is the input
     * - x are the internal states
     * - y is the output
     * 
     * This is only a representation using the A, B, C, D matrices.
     * This does **not** contain an internal state. 
     * You may use the `.eval()` method, but have to provide the state manually.
     * 
     * If you want a system with internal state have a look at: `DssFilter`.
     * 
     * This is a type wrapper around `controlpp::StateSpace`
     * 
     * \see controlpp::StateSpace
     * \see controlpp::DiscreteTransferFunction
     * \see controlpp::DssFilter
     * 
     * \tparam ValueType The value type of the matrix elements and arithmetic calculations (usually `double` or `float`)
     * \tparam internal_states The number of internal states, this influences the size of the A, B and C matrices
     * \tparam inputs The number of input to the system, this influences the size of the B and D matrix
     * \tparam outputs The number of outputs of the system, this influences the size of the C and D matrix
     */
    template<class ValueType, int internal_states, int inputs, int outputs>
    class DiscreteStateSpace{
        public:
            using value_type = ValueType;

            using state_space_type = StateSpace<ValueType, internal_states, inputs, outputs>;

            using A_matrix_type = typename state_space_type::A_matrix_type;
            using B_matrix_type = typename state_space_type::B_matrix_type;
            using C_matrix_type = typename state_space_type::C_matrix_type;
            using D_matrix_type = typename state_space_type::D_matrix_type;

            constexpr static int number_of_states = internal_states;
            constexpr static int number_of_inputs = inputs;
            constexpr static int number_of_outputs = outputs;

        private:
            state_space_type _state_space;

        public:
            constexpr DiscreteStateSpace() = default;
            constexpr DiscreteStateSpace(const DiscreteStateSpace&) = default;
            constexpr DiscreteStateSpace& operator=(const DiscreteStateSpace&) = default;

            constexpr DiscreteStateSpace(
                const Eigen::Matrix<ValueType, internal_states, internal_states>& A,
                const Eigen::Matrix<ValueType, internal_states, inputs>& B,
                const Eigen::Matrix<ValueType, outputs, internal_states>& C,
                const Eigen::Matrix<ValueType, outputs, inputs>& D)
                : _state_space(A, B, C, D){}

            constexpr DiscreteStateSpace(const state_space_type& state_space) 
                : _state_space(state_space){}

            constexpr state_space_type& state_space(){return this->_state_space;}
            constexpr const state_space_type& state_space() const {return this->_state_space;}

            /**
             * \brief calculates the next states and calculates the output from the previous states and new inputs
             * \returns a tuple of `[nest_states, output]`
             */
            constexpr std::tuple<Eigen::Vector<ValueType, internal_states>, Eigen::Vector<ValueType, outputs>> eval(const Eigen::Vector<ValueType, internal_states>& x, const Eigen::Vector<ValueType, inputs>& u) const {
                return this->_state_space.eval(x, u);
            }

            /**
             * \brief calculates the next states and calculates the output from the previous states and new inputs
             * \returns a tuple of `[nest_states, output]`
             */
            template<std::same_as<ValueType> U>
                requires(inputs == 1 && outputs != 1)
            constexpr std::tuple<Eigen::Vector<U, internal_states>, Eigen::Vector<U, outputs>> eval(const Eigen::Vector<U, internal_states>& x, const U& u_scalar) const {
                return this->_state_space.eval(x, u_scalar);
            }

            /**
             * \brief calculates the next states and calculates the output from the previous states and new inputs
             * \returns a tuple of `[nest_states, output]`
             */
            template<std::same_as<ValueType> U>
                requires(inputs == 1 && outputs == 1)
            constexpr std::tuple<Eigen::Vector<U, internal_states>, U> eval(const Eigen::Vector<U, internal_states>& x, const U& u_scalar) const {
                return this->_state_space.eval(x, u_scalar);
            }
            
            
            constexpr A_matrix_type& A() {return this->_state_space.A();}
            constexpr B_matrix_type& B() {return this->_state_space.B();}
            constexpr C_matrix_type& C() {return this->_state_space.C();}
            constexpr D_matrix_type& D() {return this->_state_space.D();}

            constexpr const A_matrix_type& A() const {return this->_state_space.A();}
            constexpr const B_matrix_type& B() const {return this->_state_space.B();}
            constexpr const C_matrix_type& C() const {return this->_state_space.C();}
            constexpr const D_matrix_type& D() const {return this->_state_space.D();}

            friend std::ostream& operator<<(std::ostream& stream, const DiscreteStateSpace& dss){
                stream << dss.state_space();
                return stream;
            }
    };

    /**
     * \brief constructs a discrete state space function from a rational polynom
     */
    template<class T, int NumOrder, int DenOrder>
    requires(NumOrder <= DenOrder)
    constexpr DiscreteStateSpace<T, DenOrder, 1, 1> to_discrete_state_space(const TransferFunction<T, NumOrder, DenOrder>& rp){
        static constexpr int number_of_states = DenOrder;
        StateSpace<T, number_of_states, 1, 1> result;
        const T a_0 = rp.den(0);
        const T b_0 = rp.num(0);

        // normalise
        const Polynom<T, DenOrder> a = -(rp.den() / a_0);
        const Polynom<T, NumOrder> b = rp.num() / a_0;
        
        // write A matrix
        if constexpr (number_of_states > 0){
            const auto I = Eigen::Matrix<T, number_of_states-1, number_of_states-1>::Identity();
            result.A().setOnes();
            result.A().template block<number_of_states-1, number_of_states-1>(1, 0) = I; 
            result.A().col(number_of_states-1).tail(number_of_states-1).setZero(); 
            result.A().row(0) = a.vector().tail(a.vector().size()-1);
        }

        // write B matrix
        if constexpr (number_of_states > 0){
            result.B()(0, 0) = T(1);
            result.B().col(0).tail(number_of_states-1).setZero();
        }
        
        // write C matrix
        if constexpr (number_of_states > 0){
            Eigen::Vector<double, DenOrder> aa = a.vector().tail(DenOrder) * b_0;
            Eigen::Vector<double, DenOrder> bb = Eigen::Vector<double, DenOrder>::Zero();
            bb.head(NumOrder) = b.vector().tail(b.size()-1);
            result.C().row(0) = bb + aa;
        }
        
        // write D matrix
        result.D()(0, 0) = b_0;
        return result;
    }

    /**
     * \brief constructs a continuous state space function from a continuous transfer function
     */
    template<class ValueType, int NumOrder, int DenOrder>
    constexpr DiscreteStateSpace<ValueType, DenOrder, 1, 1> to_discrete_state_space(const DiscreteTransferFunction<ValueType, NumOrder, DenOrder>& dtf){
        return DiscreteStateSpace<ValueType, DenOrder, 1, 1>(to_discrete_state_space(dtf.transfer_function()));
    }

    /**
     * \brief constructs a continuous state space function from a continuous transfer function
     */
    template<class ValueType, int NumOrder, int DenOrder>
    constexpr DiscreteStateSpace<ValueType, DenOrder, 1, 1> to_state_space(const DiscreteTransferFunction<ValueType, NumOrder, DenOrder>& dtf){
        return DiscreteStateSpace<ValueType, DenOrder, 1, 1>(to_discrete_state_space(dtf.transfer_function()));
    }

    /**
     * \brief Transforms a discrete state space system to a discrete transfer function
     * \returns a discrete transfer function `controlpp::DiscreteTransferFunction`
     * \see controlpp::DiscreteTransferFunction
     */
    // redo:
    //template<class T, int states>
    //constexpr DiscreteTransferFunction<T, states+1, states+1> to_transfer_function(const DiscreteStateSpace<T, states, 1, 1>& dss){
    //    return DiscreteTransferFunction<T, states+1, states+1>(to_transfer_function(dss.state_space()));
    //}

} // namespace control