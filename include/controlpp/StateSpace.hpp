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

    /**
     * \brief calculates the state space representation from a rational polynomial
     * 
     * TODO: Testing
     */
    template<class T, size_t num_size, size_t den_size>
    constexpr StateSpace<T, den_size-1, 1, 1> to_state_space(const RationalPolynom<T, num_size, den_size>& rp){
        StateSpace<T, den_size-1, 1, 1> result;
        const T a_n = rp.den[rp.den().order()];
        const Polynom<T, den_size> a = -(rp.den/a_n);
        const Polynom<T, num_size> b = rp.num / a_n;
        
        result.A().template block<result.A().rows()-1, result.A().columns()-1>(0, 1) = Eigen::Matrix<result.A().rows()-1, result.A().columns()-1>::Identity();
        result.A().col(result.A().rows()-1).setZero();
        result.A().row(result.A().rows()-1) = a.vector().head(a.vector().size()-1);
        
        result.B().col(0).head(result.B().size()-1).setZero();
        result.B()(-1, 0) = T(1);
        
        result.C().row(0).head(b.size()) = b.vector();
        result.C().row(0).tail(result.C().columns() - b.size()).setZero();
        
        result.D()(0, 0) = (b.size() > (den_size-1)) ? b[den_size-1] : T(0);
        
        return result;
    }


} // namespace control
