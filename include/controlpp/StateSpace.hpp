#pragma once

#include <Eigen/Core>

namespace control
{
    template<class T, size_t internal_states, size_t inputs, size_t outputs>
    class StateSpace{
        public:
            template<size_t rows, size_t cols>
            using matrix_type = Eigen::Matrix<T, rows, cols>;

            using A_matrix_type = matrix_type<internal_states, internal_states>;
            using B_matrix_type = matrix_type<internal_states, inputs>;
            

        private:
            A_matrix_type _A;
            Matrix<T, internal_states, inputs> _B;
            Matrix<T, outputs, internal_states> _C;
            Matrix<T, outputs, inputs> _D;

        public:
            constexpr StateSpace() = default;
            
            constexpr static size_t number_of_states = internal_states;
            constexpr static size_t number_of_inputs = inputs;
            constexpr static size_t number_of_outputs = outputs;
            
            constexpr Matrix<T, internal_states, internal_states>& A(){return this->_A;}
            constexpr const Matrix<T, internal_states, internal_states>& A()const{return this->_A;}
            
            constexpr Matrix<T, internal_states, inputs>& B(){return this->_B;}
            constexpr const Matrix<T, internal_states, inputs>& B() const {return this->_B;}
            
            constexpr Matrix<T, outputs, internal_states>& C(){return this->_C;}
            constexpr const Matrix<T, outputs, internal_states>& C() const {return this->_C;}
            
            constexpr Matrix<T, outputs, inputs> D() {return this->_D;}
            constexpr const Matrix<T, outputs, inputs> D() const {return this->_D;}
    };
} // namespace control
