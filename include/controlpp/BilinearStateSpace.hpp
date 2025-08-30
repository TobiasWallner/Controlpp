#pragma once

#include "StateSpace.hpp"
#include "BilinearTransferFunction.hpp"

namespace controlpp
{
    template<class ValueType, int internal_states, int inputs, int outputs>
    class BilinearStateSpace{
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
            constexpr BilinearStateSpace() = default;
            constexpr BilinearStateSpace(const BilinearStateSpace&) = default;
            constexpr BilinearStateSpace& operator=(const BilinearStateSpace&) = default;

            constexpr BilinearStateSpace(
                const Eigen::Matrix<ValueType, internal_states, internal_states>& A,
                const Eigen::Matrix<ValueType, internal_states, inputs>& B,
                const Eigen::Matrix<ValueType, outputs, internal_states>& C,
                const Eigen::Matrix<ValueType, outputs, inputs>& D
            )
                : _state_space(A, B, C, D){}

            constexpr BilinearStateSpace(const state_space_type& state_space) 
                : _state_space(state_space){}

            constexpr state_space_type& state_space(){return this->_state_space;}
            constexpr const state_space_type& state_space() const {return this->_state_space;}

            constexpr std::tuple<Eigen::Vector<ValueType, internal_states>, Eigen::Vector<ValueType, outputs>> eval(const Eigen::Vector<ValueType, internal_states>& x, const Eigen::Vector<ValueType, inputs>& u) const {
                this->_state_space.eval(x, u);
            }

            template<std::same_as<ValueType> U>
                requires(inputs == 1)
            constexpr std::tuple<Eigen::Vector<U, internal_states>, Eigen::Vector<U, outputs>> eval(const Eigen::Vector<U, internal_states>& x, const U& u_scalar) const {
                this->_state_space.eval(x, u_scalar);
            }
            
            
            constexpr A_matrix_type& A() {return this->_state_space.A();}
            constexpr B_matrix_type& B() {return this->_state_space.B();}
            constexpr C_matrix_type& C() {return this->_state_space.C();}
            constexpr D_matrix_type& D() {return this->_state_space.D();}

            constexpr const A_matrix_type& A() const {return this->_state_space.A();}
            constexpr const B_matrix_type& B() const {return this->_state_space.B();}
            constexpr const C_matrix_type& C() const {return this->_state_space.C();}
            constexpr const D_matrix_type& D() const {return this->_state_space.D();}

            friend std::ostream& operator<<(std::ostream& stream, const BilinearStateSpace& dss){
                stream << dss.state_space();
                return stream;
            }
    };

    /**
     * \brief constructs a continuous state space function from a rational polynom
     */
    template<class ValueType, int NumOrder, int DenOrder>
    constexpr BilinearStateSpace<ValueType, DenOrder, 1, 1> to_bilinear_state_space(const TransferFunction<ValueType, NumOrder, DenOrder>& rp){
        return BilinearStateSpace<ValueType, DenOrder, 1, 1>(to_state_space(rp));
    }

    /**
     * \brief constructs a continuous state space function from a continuous transfer function
     */
    template<class ValueType, int NumOrder, int DenOrder>
    constexpr BilinearStateSpace<ValueType, DenOrder, 1, 1> to_bilinear_state_space(const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& dtf){
        return BilinearStateSpace<ValueType, DenOrder, 1, 1>(to_state_space(dtf.transfer_function()));
    }

    template<class ValueType, int NumOrder, int DenOrder>
    constexpr BilinearStateSpace<ValueType, DenOrder, 1, 1> to_state_space(const BilinearTransferFunction<ValueType, NumOrder, DenOrder>& dtf){
        return BilinearStateSpace<ValueType, DenOrder, 1, 1>(to_state_space(dtf.transfer_function()));
    }

    /**
     * \brief Transforms a discrete state space system to a discrete transfer function
     * \returns a discrete transfer function `controlpp::DiscreteTransferFunction`
     * \see controlpp::DiscreteTransferFunction
     */
    template<class T, int states>
    constexpr BilinearTransferFunction<T, states+1, states+1> to_transfer_function(const BilinearStateSpace<T, states, 1, 1>& dss){
        return BilinearTransferFunction<T, states+1, states+1>(to_transfer_function(dss.state_space()));
    }


} // namespace controlpp