#pragma once

#include "math.hpp"
#include "ContinuousTransferFunction.hpp"

namespace controlpp
{

	/**
	 * \brief Calculates the numerator parameters of the pade approximation
     * 
     * Calculates the parameter:
     * 
     * \f[
     * P_k = \frac{(m + n - k)! m!}{(m + n)! k! (m - k)!}
     * \f]
     * 
     * For a Pade approximation like:
     * 
     * \f[
     * A = \frac{\sum_{k = 0}^{m} P_k s^k}{\sum_{k = 0}^{n} Q_k s^k} 
     * \f]
     * 
     * \param m Is the order of the numerator
     * \param n Is the order of the denominator
     * \param k Is the order of the parameter
     * 
     * \tparam The result type of the division
     * 
     * \returns The pade parameter of the numerator at the k-th order
     * 
     * \see pade_den_param
     * \see https://ris.utwente.nl/ws/portalfiles/portal/134422804/Some_remarks_on_Pade-approximations.pdf
	 */
    template<class T = double>
	constexpr T pade_num_param(unsigned long m, unsigned long n, unsigned long k){
		const unsigned long num1 = product_over(k+1, m+1);
		const unsigned long den1 = product_over(m + n - k + 1, m + n + 1);
		const unsigned long den2 = product_over(1, m-k + 1);
		const T result = static_cast<T>(num1) / static_cast<T>(den1 * den2);
		return result;
	}

    /**
	 * \brief Calculates the numerator parameters of the pade approximation
     * 
     * Calculates the parameter:
     * 
     * \f[
     * Q_k = \frac{(n + m - k)! n!}{(n + m)! k! (n - k)!}
     * \f]
     * 
     * For a Pade approximation like:
     * 
     * \f[
     * A = \frac{\sum_{k = 0}^{m} P_k s^k}{\sum_{k = 0}^{n} Q_k s^k} 
     * \f]
     * 
     * \param m Is the order of the numerator
     * \param n Is the order of the denominator
     * \param k Is the order of the parameter
     * 
     * \tparam T the result type of the function (used for the final division)
     * 
     * \returns The pade parameter of the numerator at the k-th order
     * 
     * \see pade_num_param
     * \see https://ris.utwente.nl/ws/portalfiles/portal/134422804/Some_remarks_on_Pade-approximations.pdf
	 */
    template<class T = double>
	constexpr T pade_den_param(unsigned long m, unsigned long n, unsigned long k){
		return pade_num_param(n, m, k);
	}

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
    ContinuousTransferFunction<T, NumOrder, DenOrder> pade(T delay){
        ContinuousTransferFunction<T, NumOrder, DenOrder> result;
        
        for(int k = 0; k < NumOrder; ++k){
            result.num().at(k) = pade_num_param<T>(NumOrder, DenOrder, k) * controlpp::pow(-delay, k);
        }

        for(int k = 0; k < DenOrder; ++k){
            result.den().at(k) = pade_num_param<T>(NumOrder, DenOrder, k) * controlpp::pow(delay, k);
        }

        return result;
    }


} // namespace controlpp
