#pragma once

#include "math.hpp"
#include "ContinuousTransferFunction.hpp"

namespace controlpp
{
    /**
     * \brief Creates a continuous transfer function approximating a delay using the pade approximation
     * 
     * 
     * > The standard Pade approximation where the order of the numerator is equal to the
     * > order of the denominator exhibits a jump at t=0 in its step response. To avoid this
     * > the use of pade approximations where the numerator order is one less than that of the
     * > denominator is reccomended. This gives a better step response
     * >
     * > *SOME REMARKS ON PADÉ-APPROXIMATION* by *M.Vajta*
     * 
     * \param delay The delay that should be approximated using a 
     * 
     * \tparam T The value type of the ContinuousTransferFunction
     * \tparam NumOrder The order of the numerator
     * \tparam DenOrder The order of the denominator
     * 
     * \see https://ris.utwente.nl/ws/portalfiles/portal/134422804/Some_remarks_on_Pade-approximations.pdf
     */
    template<class T, int NumOrder, int DenOrder=NumOrder>
    requires(NumOrder <= DenOrder)
    ContinuousTransferFunction<T, NumOrder, DenOrder> pade_delay(T delay){
        ContinuousTransferFunction<T, NumOrder, DenOrder> result;
        
        for(int k = 0; k < (NumOrder+1); ++k){
            result.num().at(k) = pade_num_param<T>(NumOrder, DenOrder, k) * controlpp::pow(-delay, k);
        }

        for(int k = 0; k < (DenOrder+1); ++k){
            result.den().at(k) = pade_den_param<T>(NumOrder, DenOrder, k) * controlpp::pow(delay, k);
        }

        return result;
    }


} // namespace controlpp
