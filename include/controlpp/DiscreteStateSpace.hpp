#pragma once

#include "StateSpace.hpp"

namespace control
{
    template<class ValueType, class TimeType, size_t internal_states, size_t inputs, size_t outputs>
    class DiscreteStateSpace{
        public:
            using value_type = ValueType;
            using time_type = TimeType;

            using state_space_type = StateSpace<T, internal_states, inputs, outputs>

            using A_matrix_type = typename state_space_type::A_matrix_type;
            using B_matrix_type = typename state_space_type::B_matrix_type;
            using C_matrix_type = typename state_space_type::C_matrix_type;
            using D_matrix_type = typename state_space_type::D_matrix_type;

            constexpr static size_t number_of_states = internal_states;
            constexpr static size_t number_of_inputs = inputs;
            constexpr static size_t number_of_outputs = outputs;

        private:
            state_space_type _state_space;
            time_type _sample_time;

        public:
            constexpr DiscreteStateSpace() = default;
            constexpr DiscreteStateSpace(const DiscreteStateSpace&) = default;
            constexpr DiscreteStateSpace& operator=(const DiscreteStateSpace&) = default;

            constexpr DiscreteStateSpace(
                const A_matrix_type& A,
                const B_matrix_type& B,
                const C_matrix_type& C,
                const D_matrix_type& D,
                const time_type& sample_time
            )
                : _state_space(A, B, C, D)
                , _sample_time(sample_time){}

            constexpr DiscreteStateSpace(const state_space_type& state_space, const time_type& sample_time) 
                : _state_space(state_space)
                , _sample_time(sample_time){}
            
            constexpr const state_space_type& state_space() const {return this->_state_space;}
            constexpr const time_type& sample_time() const {return this->_sample_time;}
            constexpr const A_matrix_type& A() const {return this->_state_space.A();}
            constexpr const B_matrix_type& B() const {return this->_state_space.B();}
            constexpr const C_matrix_type& C() const {return this->_state_space.C();}
            constexpr const D_matrix_type& D() const {return this->_state_space.D();}
    };
} // namespace control