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
            value_type _sample_time;

        public:
            constexpr DiscreteStateSpace() = default;
            constexpr DiscreteStateSpace(const DiscreteStateSpace&) = default;
            constexpr DiscreteStateSpace& operator=(const DiscreteStateSpace&) = default;

            constexpr DiscreteStateSpace(
                const A_matrix_type& A,
                const B_matrix_type& B,
                const C_matrix_type& C,
                const D_matrix_type& D,
                const value_type& sample_time
            )
                : _state_space(A, B, C, D)
                , _sample_time(sample_time){}

            constexpr DiscreteStateSpace(const state_space_type& state_space, const value_type& sample_time) 
                : _state_space(state_space)
                , _sample_time(sample_time){}

            constexpr state_space_type& state_space(){return this->_state_space;}
            constexpr const state_space_type& state_space() const {return this->_state_space;}

            constexpr const value_type& sample_time() const {return this->_sample_time;}
            constexpr value_type& sample_time() {return this->_sample_time;}

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

            friend std::ostream& operator<<(std::ostream& stream, const DiscreteStateSpace& dss){
                stream << dss.state_space();
                stream << "Sample time: " << dss.sample_time() << '\n';
                return stream;
            }
    };

    /**
     * \brief constructs a continuous state space function from a rational polynom
     */
    template<class ValueType, size_t num_size, size_t den_size>
    constexpr DiscreteStateSpace<ValueType, den_size-1, 1, 1> to_DiscreteStateSpace(const RationalPolynom<ValueType, num_size, den_size>& rp, const ValueType& sample_time){
        return DiscreteStateSpace<ValueType, den_size-1, 1, 1>(to_state_space(rp), sample_time);
    }

    /**
     * \brief constructs a continuous state space function from a continuous transfer function
     */
    template<class ValueType, size_t num_size, size_t den_size>
    constexpr DiscreteStateSpace<ValueType, den_size-1, 1, 1> to_DiscreteStateSpace(const DiscreteTransferFunction<ValueType, num_size, den_size>& dtf){
        return DiscreteStateSpace<ValueType, den_size-1, 1, 1>(to_state_space(dtf.ratpoly()));
    }

    // template<class ValueType, size_t num_size, size_t den_size>
    // constexpr DiscreteTransferFunction<ValueType, den_size-1, 1, 1> to_DiscreteStateSpace(const DiscreteStateSpace<ValueType, num_size, den_size>& dtf){
    //     return DiscreteTransferFunction<ValueType, den_size-1, 1, 1>(to_state_space(dtf.ratpoly()));
    // }

} // namespace controlpp