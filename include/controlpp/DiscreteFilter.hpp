#pragma once

#include <ostream>

#include "ContinuousTransferFunction.hpp"
#include "ContinuousStateSpace.hpp"
#include "DiscreteTransferFunction.hpp"
#include "DiscreteStateSpace.hpp"
#include "transformations.hpp"

namespace controlpp{

    

    /**
     * \brief Controller from a discrete state space
     */
    template<class T, int NStates_, int NInputs_=1, int NOutputs_=1>
    class DssFilter{
        private:
        
        Eigen::Vector<T, NStates_> states_ = Eigen::Vector<T, NStates_>::Zero();
        DiscreteStateSpace<T, NStates_, NInputs_, NOutputs_> dss_;

        public:

        static constexpr int n_states = NStates_;
        static constexpr int n_inputs = NInputs_;
        static constexpr int n_outputs = NOutputs_;

        /**
         * \brief construct a discrete filter from a state space
         * 
         * A filter will also internally track 
         */
        DssFilter(const DiscreteStateSpace<T, NStates_, NInputs_, NOutputs_>& dss)
            : dss_(dss){}

        template<int N>
        DssFilter(const DiscreteTransferFunction<T, N, NStates_>& dtf)
            : DssFilter(to_state_space(dtf)){}

        DssFilter(
            const ContinuousStateSpace<T, NStates_, NInputs_, NOutputs_>& css, 
            double Ts, 
            EDiscretisation method = EDiscretisation::tustin
        )
            : DssFilter(discretise(css, Ts, method)){}

        template<int N>
        DssFilter(
            const ContinuousTransferFunction<T, N, NStates_>& ctf,
            double Ts,
            EDiscretisation method = EDiscretisation::tustin
        )
            : DssFilter(to_state_space(ctf), Ts, method){}

        Eigen::Vector<T, NOutputs_> input(const Eigen::Vector<T, NInputs_>& u){
            const auto [x, y] = this->dss_.eval(this->states_, u);
            this->states_ = x;
            return y;
        }

        /**
         * \brief input a new value to advance the internal state and calculate the new output
         * \returns  
         */
        template<std::same_as<T> U>
        requires(NInputs_ == 1 && NOutputs_ > 1)
        Eigen::Vector<T, NOutputs_> input(const U& u){
            const auto [x, y] = this->dss_.eval(this->states_, u);
            this->states_ = x;
            return y;
        }

        /**
         * \brief input a new value to advance the internal state and calculate the new output
         * \returns  
         */
        template<std::same_as<T> U>
        requires(NInputs_ == 1 && NOutputs_ == 1)
        T input(const U& u){
            const auto [x, y] = this->dss_.eval(this->states_, u);
            this->states_ = x;
            return y;
        }

        template<std::same_as<T> U>
        requires(NInputs_ == 1 && NOutputs_ == 1)
        T operator() (const U& u){
            return this->input(u);
        }

        /**
         * \brief sets the interanal state vector to zero
         */
        void clear(){
            this->states_.setZero();
        }

        /**
         * \brief returns the contained state space object
         * \returns a const reference of a discrete state space object
         */
        const DiscreteStateSpace<T, NStates_, NInputs_, NOutputs_>& dss() const {
            return this->dss_;
        }

        /**
         * \brief returns the contained state space object
         * \returns a mutable reference of a discrete state space object
         */
        DiscreteStateSpace<T, NStates_, NInputs_, NOutputs_>& dss() {
            return this->dss_;
        }

        friend std::ostream& operator<< (std::ostream& stream, const DssFilter& dssf){
            stream << dssf.dss() << "\n";
            stream << "states: " << dssf.states_.transpose() << "\n";
            return stream;
        }
    };



    /**
     * \brief Controller from a discrete transfer function
     */
    template<class T, int NumOrder, int DenOrder>
    class DtfFilter{
        private:
        
        Eigen::Vector<T, NumOrder> uk_ = Eigen::Vector<T, NumOrder>::Zero();
        Eigen::Vector<T, DenOrder> yk_ = Eigen::Vector<T, DenOrder>::Zero();
        DiscreteTransferFunction<T, NumOrder, DenOrder> dtf_;
        public:

        /**
         * \brief construct a discrete filter from a Transfer Function
         * 
         * A filter will also internally track 
         */
        DtfFilter(const DiscreteTransferFunction<T, NumOrder, DenOrder>& dtf)
            : dtf_(dtf)
        {}

        T input(const T& u){
            const Eigen::Vector<T, NumOrder+1> uk;
            uk(0) = u;
            uk.tail(NumOrder) = this->uk_;

            const auto y = this->dtf_.eval(uk, this->yk_);

            std::copy_backward(this->uk_.data(), this->uk_.data()+this->uk_.size(), this->uk_.data()+1);
            this->uk_(0) = u;

            std::copy_backward(this->yk_.data(), this->yk_.data()+this->yk_.size(), this->yk_.data()+1);
            this->yk_(0) = y;
            return y;
        }

        /**
         * \brief sets the interanal state vector to zero
         */
        void clear(){
            this->uk_.setZero();
            this->yk_.setZero();
        }

        const DiscreteTransferFunction<T, NumOrder, DenOrder>& dft() const {
            return this->dtf_;
        }

        DiscreteTransferFunction<T, NumOrder, DenOrder>& dft() {
            return this->dtf_;
        }

    };

}