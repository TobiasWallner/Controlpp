#pragma once

#include <Eigen/Dense>

#include <controlpp/math.hpp>
#include <controlpp/ContinuousStateSpace.hpp>


namespace controlpp{

    /**
     * @brief Synthesizes the continuous LQR gain from plant matrices (A, B) and weights (Q, R)
     * @tparam T The valuetype of the parameters (usually `float` or `double`)
     * @tparam NStates The number of states in the plant
     * @tparam NInputs The number of inputs of the plant
     * @param A Plant state transition matrix
     * @param B Plant input matrix
     * @param R Control wight matrix (Penalises control effort). Expected to be symetric positive definite
     * @param Q State weight matrix (Penalises large plant states)
     * @returns An LQR controller in contunuous state space form
     * @see ContinuousStateSpace
     */
    template<class T, int NStates, int NInputs>
    Eigen::Matrix<T, NInputs, NStates> lqr_continuous(
        const Eigen::Matrix<T, NStates, NStates>& A,
        const Eigen::Matrix<T, NStates, NInputs>& B,
        const Eigen::Matrix<T, NStates, NStates>& Q = Eigen::Matrix<T, NStates, NStates>::Identity(),
        const Eigen::Matrix<T, NInputs, NInputs>& R = Eigen::Matrix<T, NInputs, NInputs>::Identity()
    ){
        const Eigen::Matrix<T, NStates, NStates> P = care_solver(A, B, Q, R);
        const Eigen::Matrix<T, NInputs, NStates> K = R.llt().solve(B.transpose() * P);
        return K;
    }

    /**
     * @brief Synthesizes the continuous LQR gain from plant matrices (A, B) and weights (Q, R)
     * @tparam T The valuetype of the parameters (usually `float` or `double`)
     * @tparam NStates The number of states in the plant
     * @tparam NInputs The number of inputs of the plant
     * @param A Plant state transition matrix
     * @param B Plant input matrix
     * @param Q State weight matrix (Penalises large plant states)
     * @param R Control wight matrix (Penalises control effort). Expected to be symetric positive definite
     * @returns An LQR controller in contunuous state space form
     * @see ContinuousStateSpace
     */
    template<class T, int NStates, int NInputs>
    Eigen::Matrix<T, NInputs, NStates> lqr_discrete(
        const Eigen::Matrix<T, NStates, NStates>& A,
        const Eigen::Matrix<T, NStates, NInputs>& B,
        const Eigen::Matrix<T, NStates, NStates>& Q = Eigen::Matrix<T, NStates, NStates>::Identity(),
        const Eigen::Matrix<T, NInputs, NInputs>& R = Eigen::Matrix<T, NInputs, NInputs>::Identity()
    ){
        
        std::cout << "A:\n  " << A << std::endl;
        std::cout << "B:\n  " << B << std::endl;
        const Eigen::Matrix<T, NStates, NStates> P = dare_solver(A, B, Q, R);
        std::cout << "P:\n  " << P << std::endl;

        const Eigen::Matrix<T, NInputs, NInputs> M1 = R + B.transpose() * P * B;
        const Eigen::Matrix<T, NInputs, NStates> M2 = B.transpose() * P * A;
        const Eigen::Matrix<T, NInputs, NStates> K = M1.ldlt().solve(M2);

        std::cout << "K:\n  " << K << std::endl;
        return K;
    }
}