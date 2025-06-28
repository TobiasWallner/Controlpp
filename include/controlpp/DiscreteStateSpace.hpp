#pragma once

#include "StateSpace.hpp"

namespace controlpp
{
    template<class ValueType, size_t internal_states, size_t inputs, size_t outputs>
    class DiscreteStateSpace{
        public:
            using value_type = ValueType;

            using state_space_type = StateSpace<ValueType, internal_states, inputs, outputs>;

            using A_matrix_type = typename state_space_type::A_matrix_type;
            using B_matrix_type = typename state_space_type::B_matrix_type;
            using C_matrix_type = typename state_space_type::C_matrix_type;
            using D_matrix_type = typename state_space_type::D_matrix_type;

            constexpr static size_t number_of_states = internal_states;
            constexpr static size_t number_of_inputs = inputs;
            constexpr static size_t number_of_outputs = outputs;

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
     * \brief constructs a continuous state space function from a rational polynom
     */
    template<class ValueType, size_t num_size, size_t den_size>
    constexpr DiscreteStateSpace<ValueType, den_size-1, 1, 1> to_DiscreteStateSpace(const RationalPolynom<ValueType, num_size, den_size>& rp){
        return DiscreteStateSpace<ValueType, den_size-1, 1, 1>(to_StateSpace(rp));
    }

    /**
     * \brief constructs a continuous state space function from a continuous transfer function
     */
    template<class ValueType, size_t num_size, size_t den_size>
    constexpr DiscreteStateSpace<ValueType, den_size-1, 1, 1> to_DiscreteStateSpace(const DiscreteTransferFunction<ValueType, num_size, den_size>& dtf){
        return DiscreteStateSpace<ValueType, den_size-1, 1, 1>(to_StateSpace(dtf.ratpoly()));
    }

    /**
     * \brief constructs a continuous state space function from a continuous transfer function
     */
    template<class ValueType, size_t num_size, size_t den_size>
    constexpr DiscreteStateSpace<ValueType, den_size-1, 1, 1> to_StateSpace(const DiscreteTransferFunction<ValueType, num_size, den_size>& dtf){
        return DiscreteStateSpace<ValueType, den_size-1, 1, 1>(to_StateSpace(dtf.ratpoly()));
    }

} // namespace controlpp