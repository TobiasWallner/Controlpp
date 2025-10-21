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
        const Eigen::Matrix<T, NStates, NStates> P = dare_solver(A, B, Q, R);
        const Eigen::Matrix<T, NInputs, NInputs> M1 = R + B.transpose() * P * B;
        const Eigen::Matrix<T, NInputs, NStates> M2 = B.transpose() * P * A;
        const Eigen::Matrix<T, NInputs, NStates> K = M1.ldlt().solve(M2);
        return K;
    }

    /**
     * @brief Computes a discrete LQR controller from a discrete state space plant model
     * @tparam T The value types of the parameters
     * @tparam NStates The number of states of the plant
     * @tparam NInputs The number of inputs of the plant
     * @tparam NOutputs Then number of outputs of the plant
     * @param Q The state weight matrix (penalizes large states)
     * @param R The control weight matrix (penalizes control effort)
     * @return The gain matrix of the LQR controller as an `Eigen::Matrix`
     * 
     * @see lqr_feed_forward
     */
    template<class T, int NStates, int NInputs, int NOutputs>
    Eigen::Matrix<T, NInputs, NStates> lqr(
        const DiscreteStateSpace<T, NStates, NInputs, NOutputs> Gss,
        const Eigen::Matrix<T, NStates, NStates>& Q,
        const Eigen::Matrix<T, NInputs, NInputs>& R
    ){
        return lqr_discrete(Gss.A(), Gss.B(), Q, R);
    }

    /**
     * @brief Construcs a discrete LQR controller with automatic weights Q and R
     * 
     * The weights Q (state penalty) and R (control penalty) are chosen automatically.
     * 
     * TODO: explain how Q and R are chosen
     * 
     * 
     * @tparam T The value type of the parameters
     * @tparam NStates The number of states of the plant
     * @tparam NInputs The number of inputs of the plant
     * @tparam NOutputs The number of outputs of the plant
     * @param Gss The state space description of the plant
     * @param p
     * @returns The LQR gain matrix as an `Eigen::Matrix`
     * 
     * @see lqr_feed_forward
     */
    template<class T, int NStates, int NInputs, int NOutputs, std::convertible_to<T> U1 = T, std::convertible_to<T> U2 = T>
    Eigen::Matrix<T, NInputs, NStates> lqr(
        const DiscreteStateSpace<T, NStates, NInputs, NOutputs> Gss
    ){
        const Eigen::Matrix<T, NStates, NStates> Iq = Eigen::Matrix<T, NStates, NStates>::Identity(); 
        const Eigen::Matrix<T, NStates, NStates> Q1 = Gss.C().transpose() * Gss.C();
        const Eigen::Matrix<T, NStates, NStates> Q = Q1 + Iq * (Q1.norm() * 0.000001);

        const Eigen::Matrix<T, NInputs, NInputs> R = Eigen::Matrix<T, NInputs, NInputs>::Ones()
        
        return lqr_discrete(Gss.A(), Gss.B(), Q, R);
    }

    /**
     * @brief Construcs a discrete LQR controller with automatic weights Q and R
     * 
     * The weights Q (state penalty) and R (control penalty) are chosen automatically.
     * 
     * TODO: explain how Q and R are chosen
     * 
     * 
     * @tparam T The value type of the parameters
     * @tparam NStates The number of states of the plant
     * @tparam NInputs The number of inputs of the plant
     * @tparam NOutputs The number of outputs of the plant
     * @param Gss The state space description of the plant
     * @param p
     * @returns The LQR gain matrix as an `Eigen::Matrix`
     * 
     * @see lqr_feed_forward
     */
    template<class T, int NStates, int NInputs, int NOutputs, std::convertible_to<T> U1 = T, std::convertible_to<T> U2 = T>
    Eigen::Matrix<T, NInputs, NStates> lqr(
        const DiscreteStateSpace<T, NStates, NInputs, NOutputs> Gss,
        const Eigen::Vector<T, NOutputs>& y_max
        const Eigen::Vector<T, NInputs>& u_max
    ){
        
        const Eigen::Vector<T, NOutputs> q = static_cast<T>(1) / y_max.array().square();
        const Eigen::Matrix<T, NStates, NStates> Iq = Eigen::Matrix<T, NStates, NStates>::Identity(); 
        const Eigen::Matrix<T, NStates, NStates> Q1 = Gss.C().transpose() * q.asDiagonal() * Gss.C();
        const Eigen::Matrix<T, NStates, NStates> Q = Q1 + Iq * (Q1.norm() * 0.000001);

        const Eigen::Vector<T, NOutputs> r = static_cast<T>(1) / u_max.array().square();
        const Eigen::Matrix<T, NInputs, NInputs> R = r.asDiagonal()
        
        return lqr_discrete(Gss.A(), Gss.B(), Q, R);
    }

    template<class T, int NStates, int NInputs, int NOutputs>
    Eigen::Matrix<T, NInputs, NOutputs> lqr_feed_forward(
        const DiscreteStateSpace<T, NStates, NInputs, NOutputs> Gss,
        Eigen::Matrix<T, NInputs, NStates> LQR
    ){
        // LQR feed forward
        const Eigen::Matrix<double, NStates, NStates> I = controlpp::identity_like(Gss.A());
        const Eigen::Matrix<double, NStates, NStates> M1 = I - Gss.A() + Gss.B() * LQR;
        const Eigen::Matrix<double, 1, 1> M = Gss.C() * M1.partialPivLu().solve(Gss.B());
        const Eigen::Matrix<T, NInputs, NOutputs> F = M.inverse().eval();
        return F;
    }
}