#pragma once

#include "DiscreteStateSpace.hpp"

namespace controlpp{

    template<class T, int NStates_, int NInputs_, int NOutputs_>
    class DiscreteStateSpaceFilter{
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
        DiscreteStateSpaceFilter(const DiscreteStateSpace<T, NStates_, NInputs_, NOutputs_>& dss)
            : dss_(dss)
        {}

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

    };

}