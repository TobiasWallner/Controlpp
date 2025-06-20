#pragma once

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
        result.sample_time() = sample_time;

        return result;
    }

    /**
     * \brief convenience overload that transfrom from s-domain into z-domain using zero-order-hold
     * 
     * Use this function to discretize a plant
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
    //template<class T, size_t NumSize, size_t DenSize>
    //DiscreteTransferFunction<T, NumSize, DenSize> continuous_to_discrete(const ContinuousTransferFunction<T, NumSize, DenSize>& ctf, float sample_time){
        
    //}

    /**
     * \brief transform from z-domain into q-domain using the tustin (bilinear) transformation
     * 
     * Use this function to continuize a plant and perform the digital controller desing in the q-domain.
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
    //BilinearTransferFunction discrete_to_bilinear(const DiscreteTransferFunction& css){//TODO}

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
     */
    //DiscreteStateSpace bilinear_to_discrete(const BilinearTransferFunction& css, float sample_time){//TODO}

    /**
     * \brief convenience transform from s-domain into q-domain
     * 
     * Convenience transformation that will simplify the following processing chain:
     * 
     * Will internally do `continuous_to_discrete()` followed by `discrete_to_bilinear()`
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
    //BilinearTransferFunction continuous_to_bilinear(const ContinuousTransferFuntion& css, float sample_time){//TODO}
    
} // namespace controlpp

