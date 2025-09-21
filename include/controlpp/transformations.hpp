#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>


// controlpp
#include <controlpp/math.hpp>

#include <controlpp/ContinuousStateSpace.hpp>
#include <controlpp/ContinuousTransferFunction.hpp>

#include <controlpp/DiscreteStateSpace.hpp>
#include <controlpp/DiscreteTransferFunction.hpp>


namespace controlpp
{
    /**
     * \brief transform s-domain into z-domain using zero-order-hold
     * 
     * TODO: Testing
     * 
     * Use this function to discretize a plant.
     * 
     * To correctly design controllers, transform the plant first with zoh into the z-domain
     * and then from the z-domain into the q-domain using the tustin transformation.
     * Perform the actual controller design in q-domain and then transform the controller back into the z-domain.
     * 
     * Note: This method uses exact discretisation - but not quite - it uses an approximated matrix exponent `control::exp_taylor_scale(M)`
     * to calculate:
     * 
     * \f[
     *  \exp{\mathbf{M}}
     * \f]
     * 
     * Controller design chain:
     * -------------------------
     * - G ... Plant
     * - R ... Controller
     * 
     * 1. mathematical model --> G(s)
     * 2. discretise_zoh (zero-order-hold) --> G(z) 
     * 3. discrete_to_bilinear (tustin) --> controller design --> R(q) 
     * 4. bilinear_to_discrete (tustin) --> R(z)
     * 
     * \tparam ValueType The value type of the matrix entries (e.g.: float)
     * \tparam TimePoint The representation of the time type (e.g.: float)
     * \tparam states The number of internal states of the system
     * 
     * \param sys A continuous state space sytem
     * \param sample_time The sample time in seconds
     * 
     * \returns The discretised version of `sys`
     */
    template<class ValueType, int states>
    DiscreteStateSpace<ValueType, states, 1, 1> discretise_zoh(
            const ContinuousStateSpace<ValueType, states, 1, 1>& sys, 
            ValueType sample_time,
            int approximation_order=8
    ){
        // allocation
        DiscreteStateSpace<ValueType, states, 1, 1> result;
        Eigen::Matrix<ValueType, states+1, states+1> M;

        // scale the energy of the states

        // preparation
        M.template block<states, states>(0, 0) = sys.A();
        M.col(states).head(states) = sys.B().col(0);
        M.row(states).setZero();
        M *= sample_time;
        
        // calculation
        Eigen::Matrix<ValueType, states+1, states+1> Md = controlpp::mexp(M, approximation_order);

        // re-assignment
        result.A() = Md.template block<states, states>(0, 0);
        result.B().col(0) = Md.col(states).head(states);
        result.C() = sys.C();
        result.D() = sys.D();

        return result;
    }

    /**
     * \brief discretises a continuous state space to a discrete one with the Tustin transformation
     * 
     * Applies the following transformation:
     * 
     * \f[
     * A_d = \left( I - \frac{Ts}{2} A \right)^{-1} \left( I + \frac{Ts}{2} A \right)\\
     * B_d = \left( I - \frac{Ts}{2} A \right)^{-1} Ts B
     * C_d = C \left( I - \frac{Ts}{2} A \right)^{-1}
     * D_d = D + C * \left( I - \frac{Ts}{2} A \right)^{-1} * Ts * B / 2
     * \f]
     * 
     * where:
     *  - A, B, C, D are the continuous time system matrices
     *  - A_d, B_d, C_d, D_d are the discrete time system matrices
     *  - Ts is the sample time
     * 
     * \tparam T The value type of the matrices/systems. Usually `float` or `double`.
     * \tparam NStates The number of states of the systems
     * \tparam NInputs The number of inputs of the systems
     * \tparam NOutputs The number of outputs of the systems
     * 
     * \param sys The continuous time state space system about to be discretised
     * \param sample_time The sample time used for the discretisation
     * 
     * \returns A discrete state space system
     */
    template<class T, int NStates, int NInputs, int NOutputs>
    DiscreteStateSpace<T, NStates, NInputs, NOutputs> discretise_tustin(
            const ContinuousStateSpace<T, NStates, NInputs, NOutputs>& sys, 
            const T& sample_time
    ){
        const Eigen::Matrix<T, NStates, NStates> I = Eigen::Matrix<T, NStates, NStates>::Identity();
        const Eigen::Matrix<T, NStates, NStates> M = (I - (sample_time / static_cast<T>(2)) * sys.A()).eval();
        const Eigen::Matrix<T, NStates, NStates> P = (I + (sample_time / static_cast<T>(2)) * sys.A()).eval();

        Eigen::PartialPivLU<Eigen::Matrix<T, NStates, NStates>> Mfactor(M);

        const Eigen::Matrix<T, NStates, NStates> A_d = Mfactor.solve(P);
        const Eigen::Matrix<T, NStates, NInputs> B_d = Mfactor.solve(sample_time * sys.B());
        const Eigen::Matrix<T, NOutputs, NStates> C_d = M.transpose().partialPivLu().solve(sys.C().transpose()).transpose();
        const Eigen::Matrix<T, NOutputs, NInputs> X = (sys.C() * Mfactor.solve(sys.B() * (sample_time / static_cast<T>(2))));
        const Eigen::Matrix<T, NOutputs, NInputs> D_d = (sys.D() + X);

        const DiscreteStateSpace<T, NStates, NInputs, NOutputs> result(A_d, B_d, C_d, D_d);
        return result;
    }

    
} // namespace controlpp

