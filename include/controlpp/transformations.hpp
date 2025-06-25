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

#include <controlpp/BilinearStateSpace.hpp>


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
     * 2. continuous_to_discrete (zero-order-hold) --> G(z) 
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
    template<class ValueType, size_t states>
    constexpr DiscreteStateSpace<ValueType, states, 1, 1> continuous_to_discrete(const ContinuousStateSpace<ValueType, states, 1, 1>& sys, ValueType sample_time){
        // allocation
        DiscreteStateSpace<ValueType, states, 1, 1> result;
        Eigen::Matrix<ValueType, states+1, states+1> M;

        // preparation
        M.template block<states, states>(0, 0) = sys.A();
        M.col(states).head(states) = sys.B().col(0);
        M.row(states).setZero();
        M *= sample_time;
        
        // calculation
        M = controlpp::mexp(M);

        // re-assignment
        result.A() = M.template block<states, states>(0, 0);
        result.B().col(0) = M.col(states).head(states);
        result.C() = sys.C();
        result.D() = sys.D();

        return result;
    }

    /**
     * \brief transform s-domain into z-domain using zero-order-hold
     * 
     * alias for the function: `controlpp::continuous_to_discrete()`.
     */
    template<class ValueType, size_t states>
    constexpr DiscreteStateSpace<ValueType, states, 1, 1> s_to_z(const ContinuousStateSpace<ValueType, states, 1, 1>& sys, ValueType sample_time){
        return continuous_to_discrete<ValueType, states>(sys, sample_time);
    }

    /**
     * \brief transform from z-domain into q-domain using the tustin (bilinear) transformation
     * 
     * Controller design chain:
     * -------------------------
     * - G ... Plant
     * - R ... Controller
     * 
     * 1. mathematical model --> G(s)
     * 2. continuous_to_discrete (zero-order-hold) --> G(z) 
     * 3. discrete_to_bilinear (tustin) --> controller design --> R(q) 
     * 4. bilinear_to_discrete (tustin) --> R(z)
     */
    template<class ValueType, size_t internal_states, size_t inputs, size_t outputs>
    BilinearStateSpace<ValueType, internal_states, inputs, outputs> discrete_to_bilinear(const DiscreteStateSpace<ValueType, internal_states, inputs, outputs>& dss){
        const auto I = identity_like(dss.A());
        const auto ApI = (dss.A() + I).eval();

        const auto inv_ApI = ApI.inverse().eval();

        const auto sqrt_2 = std::sqrt(static_cast<ValueType>(2));

        const auto A_q = ((dss.A() - I) * inv_ApI).eval();
        const auto B_q = (sqrt_2 * inv_ApI * dss.B()).eval();
        const auto C_q = (sqrt_2 * dss.C() * inv_ApI).eval();
        const auto D_q = (dss.D() - dss.C() * inv_ApI * dss.B()).eval();

        BilinearStateSpace<ValueType, internal_states, inputs, outputs> result(A_q, B_q, C_q, D_q);
        return result;
    }

    /**
     * \brief transform from z-domain into q-domain using the tustin (bilinear) transformation
     * 
     * alias for the function `controlpp::discrete_to_bilinear()`.
     */
    template<class ValueType, size_t internal_states, size_t inputs, size_t outputs>
    BilinearStateSpace<ValueType, internal_states, inputs, outputs> z_to_q(const DiscreteStateSpace<ValueType, internal_states, inputs, outputs>& dss){
        return discrete_to_bilinear(dss);
    }

    /**
     * \brief transorms a system in the continuous s-domain into the bilinear q-domain
     * 
     * Internally performs:
     * 
     * 1. a zero-order-hold to go from the s-domain into the z-domain
     * 2. a forward bilinear tustin to go from the z-domain into the q-domain
     */
    template<class ValueType, size_t internal_states>
    BilinearStateSpace<ValueType, internal_states, 1, 1> continuous_to_bilinear(const ContinuousStateSpace<ValueType, internal_states, 1, 1>& css){
        return discrete_to_bilinear(continuous_to_discrete(css));
    }

    /**
     * \brief transorms a system in the continuous s-domain into the bilinear q-domain
     * 
     * alias for the function `controlpp::continuous_to_bilinear()`.
     */
    template<class ValueType, size_t internal_states>
    BilinearStateSpace<ValueType, internal_states, 1, 1> s_to_q(const ContinuousStateSpace<ValueType, internal_states, 1, 1>& css){
        return z_to_q(s_to_z(css));
    }

    /**
     * \brief transform from q-domain into z-domain using the tustin (inverse bilinear) transformation
     * 
     * Use this function to convert a controller designed in the q-domain back into the z-domain (discrete sampled world).
     * 
     * Controller design chain:
     * -------------------------
     * - G ... Plant
     * - R ... Controller
     * 
     * 1. mathematical model --> G(s)
     * 2. continuous_to_discrete (zero-order-hold) --> G(z) 
     * 3. discrete_to_bilinear (tustin) --> controller design --> R(q) 
     * 4. bilinear_to_discrete (tustin) --> R(z)
     * 
     * ---
     * 
     * Tip: For controll systems that are oversampled (\f$f_s >> f_{-3dB}\f$) the q- and s-domains are approximatelly similar.
     * Thus many desing their controller directly in the q-domain without converting to the z-domain first and then only perform
     * one bilinear transformation to go from the q-domain to the z-domain. (They may also not write q but s in their equations).
     */
    //DiscreteStateSpace bilinear_to_discrete(const BilinearTransferFunction& css, float sample_time){//TODO}
    template<class ValueType, size_t internal_states, size_t inputs, size_t outputs>
    DiscreteStateSpace<ValueType, internal_states, inputs, outputs> bilinear_to_discrete(const BilinearStateSpace<ValueType, internal_states, inputs, outputs>& bss){
        const auto I = identity_like(bss.A());
        const auto AmI = (bss.A() - I).eval();
        const auto inv_AmI = AmI.inverse().eval();
        
        const auto sqrt_2 = std::sqrt(static_cast<ValueType>(2));

        const auto A_z = (I + bss.A()) * inv_AmI;
        const auto B_z = sqrt_2 * inv_AmI * bss.B();
        const auto C_z = sqrt_2 * bss.C() * inv_AmI;
        const auto D_z = bss.D() + bss.C() * inv_AmI * bss.B();

        DiscreteStateSpace<ValueType, internal_states, inputs, outputs> result(A_z, B_z, C_z, D_z);
        return result;
    }

    /**
     * \brief transorms a system in the continuous s-domain into the bilinear q-domain
     * 
     * alias for the function `controlpp::continuous_to_bilinear()`.
     */
    template<class ValueType, size_t internal_states, size_t inputs, size_t outputs>
    DiscreteStateSpace<ValueType, internal_states, inputs, outputs> q_to_z(const BilinearStateSpace<ValueType, internal_states, inputs, outputs>& qss){
        return bilinear_to_discrete(qss);
    }

    
} // namespace controlpp

