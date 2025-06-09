#pragma once

#include "StateSpace.hpp"

namespace control
{
    template<class T, size_t internal_states, size_t inputs, size_t outputs>
    class ContinuousStateSpace{
        public:
            using value_type = T;

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

        public:
            constexpr ContinuousStateSpace() = default;
            constexpr ContinuousStateSpace(const ContinuousStateSpace&) = default;
            constexpr ContinuousStateSpace& operator=(const ContinuousStateSpace&) = default;

            constexpr ContinuousStateSpace(
                const A_matrix_type& A,
                const B_matrix_type& B,
                const C_matrix_type& C,
                const D_matrix_type& D
            ) : _state_space(A, B, C, D){}

            constexpr ContinuousStateSpace(const state_space_type& state_space) : _state_space(state_space){}
            
            constexpr state_space_type& state_space(){return this->_state_space;}
            constexpr const state_space_type& state_space() const {return this->_state_space;}

            constexpr A_matrix_type& A(){return this->_state_space.A();}
            constexpr const A_matrix_type& A()const{return this->_state_space.A();}
            
            constexpr B_matrix_type& B(){return this->_state_space.B();}
            constexpr const B_matrix_type& B() const {return this->_state_space.B();}
            
            constexpr C_matrix_type& C(){return this->_state_space.C();}
            constexpr const C_matrix_type& C() const {return this->_state_space.C();}
            
            constexpr D_matrix_type D() {return this->_state_space.D();}
            constexpr const D_matrix_type D() const {return this->_state_space.D();}
    };
} // namespace control