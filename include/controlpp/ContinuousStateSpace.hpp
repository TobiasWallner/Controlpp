#pragma once

#include "StateSpace.hpp"
#include "RationalPolynom.hpp"
#include "ContinuousTransferFunction.hpp"

namespace controlpp
{
    template<class T, size_t internal_states, size_t inputs, size_t outputs>
    class ContinuousStateSpace{
        public:
            using value_type = T;

            using state_space_type = StateSpace<T, internal_states, inputs, outputs>;

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
                const Eigen::Matrix<T, internal_states, internal_states>& A,
                const Eigen::Matrix<T, internal_states, inputs>& B,
                const Eigen::Matrix<T, outputs, internal_states>& C,
                const Eigen::Matrix<T, outputs, inputs>& D
            ) : _state_space(A, B, C, D){}

            constexpr ContinuousStateSpace(const StateSpace<T, internal_states, inputs, outputs>& state_space) : _state_space(state_space){}
            
            constexpr StateSpace<T, internal_states, inputs, outputs>& state_space(){return this->_state_space;}
            constexpr const StateSpace<T, internal_states, inputs, outputs>& state_space() const {return this->_state_space;}

            constexpr std::tuple<Eigen::Vector<T, internal_states>, Eigen::Vector<T, outputs>> eval(const Eigen::Vector<T, internal_states>& x, const Eigen::Vector<T, inputs>& u) const {
                this->_state_space(x, u);
            }

            template<std::same_as<T> U>
                requires(inputs == 1)
            constexpr std::tuple<Eigen::Vector<U, internal_states>, Eigen::Vector<U, outputs>> eval(const Eigen::Vector<U, internal_states>& x, const U& u_scalar) const {
                this->_state_space(x, u_scalar);
            }

            constexpr A_matrix_type& A() {return this->_state_space.A();}
            constexpr B_matrix_type& B() {return this->_state_space.B();}
            constexpr C_matrix_type& C() {return this->_state_space.C();}
            constexpr D_matrix_type& D() {return this->_state_space.D();}

            constexpr const A_matrix_type& A() const {return this->_state_space.A();}
            constexpr const B_matrix_type& B() const {return this->_state_space.B();}
            constexpr const C_matrix_type& C() const {return this->_state_space.C();}
            constexpr const D_matrix_type& D() const {return this->_state_space.D();}

            friend std::ostream& operator<<(std::ostream& stream, const ContinuousStateSpace& css){
                stream << css.state_space();
                return stream;
            }
    };

    /**
     * \brief constructs a continuous state space function from a rational polynom
     */
    template<class T, size_t num_size, size_t den_size>
    constexpr ContinuousStateSpace<T, den_size-1, 1, 1> to_ContinuousStateSpace(const RationalPolynom<T, num_size, den_size>& rp){
        return ContinuousStateSpace<T, den_size-1, 1, 1>(to_StateSpace(rp));
    }

    /**
     * \brief constructs a continuous state space function from a continuous transfer function
     */
    template<class T, size_t num_size, size_t den_size>
    constexpr ContinuousStateSpace<T, den_size-1, 1, 1> to_ContinuousStateSpace(const ContinuousTransferFunction<T, num_size, den_size>& ctf){
        return ContinuousStateSpace<T, den_size-1, 1, 1>(to_StateSpace(ctf.ratpoly()));
    }

    /**
     * \brief constructs a continuous state space function from a continuous transfer function
     */
    template<class T, size_t num_size, size_t den_size>
    constexpr ContinuousStateSpace<T, den_size-1, 1, 1> to_StateSpace(const ContinuousTransferFunction<T, num_size, den_size>& ctf){
        return ContinuousStateSpace<T, den_size-1, 1, 1>(to_StateSpace(ctf.ratpoly()));
    }
    
} // namespace controlpp