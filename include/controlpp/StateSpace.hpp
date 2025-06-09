#pragma once

#include <Eigen/Core>

namespace control
{
    template<class T, size_t internal_states, size_t inputs = 1, size_t outputs = 1>
    class StateSpace{
        public:
            template<size_t rows, size_t cols>
            using matrix_type = Eigen::Matrix<T, rows, cols>;

            using A_matrix_type = matrix_type<internal_states, internal_states>;
            using B_matrix_type = matrix_type<internal_states, inputs>;
            using C_matrix_type = matrix_type<outputs, internal_states>;
            using D_matrix_type = matrix_type<outputs, inputs>;

            constexpr static size_t number_of_states = internal_states;
            constexpr static size_t number_of_inputs = inputs;
            constexpr static size_t number_of_outputs = outputs;

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
                const A_matrix_type& A,
                const B_matrix_type& B,
                const C_matrix_type& C,
                const D_matrix_type& D
            )
                : _A(A)
                , _B(B)
                , _C(C)
                , _D(D){}

            constexpr const state_space_type& state_space() const {return this->_state_space;}
            constexpr const A_matrix_type& A() const {return this->_A;}
            constexpr const B_matrix_type& B() const {return this->_B;}
            constexpr const C_matrix_type& C() const {return this->_C;}
            constexpr const D_matrix_type& D() const {return this->_D;}

            constexpr A_matrix_type& A() {return this->_A;}
            constexpr B_matrix_type& B() {return this->_B;}
            constexpr C_matrix_type& C() {return this->_C;}
            constexpr D_matrix_type& D() {return this->_D;}
    };

    template<class T, size_t num_size, size_t den_size>
    constexpr StateSpace<T, den_size-1, 1, 1> to_state_space(const RationalPolynom<T, num_size, den_size>& rp){
        StateSpace<T, den_size-1, 1, 1> result;
        const T a_n = rp.den[rp.den().order()];
        const Polynom<T, den_size> a = -(rp.den/a_n);
        const Polynom<T, num_size> b = rp.num / a_n;
        
        result.A().block()

        result.A().template sub_matrix<result.A().rows()-1, result.A().columns()-1>(0, 1).set_unity();
        result.A().template column<result.A().rows()-1>(0).set_zero();
        result.A().row(result.A().rows()-1) = a.vec.template sub_vector<a.vec.size()-1>();
        
        result.B().template column<result.B().size()-1>(0).set_zero();
        result.B().at(-1, 0) = T(1);
        
        result.C().template row<b.size()>(0) = b.vec;
        result.C().template row<result.C().columns() - b.size()>(0, b.size()).set_zero();
        
        result.D().at(0, 0) = (b.size() > (den_size-1)) ? b[den_size-1] : T(0);
        
        return result;
    }


} // namespace control
