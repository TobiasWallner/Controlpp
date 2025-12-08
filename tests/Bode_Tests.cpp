#include <fstream>
#include <iostream>

#include <controlpp.hpp>

#include <gtest/gtest.h>


/*
    Tests if one can perform loop shaping from measured data
*/
TEST(FrequencyResponse, Calculation){
    const controlpp::ContinuousTransferFunction s = controlpp::tf::s<double>;

    controlpp::ContinuousTransferFunction Gs = 1 / (1 + s + s * s);
    controlpp::FrequencyResponse Gb = controlpp::bode(Gs);
    std::ofstream("bode_G.csv") << Gb;

    std::ofstream("step_G.csv") << controlpp::step(Gs);

    controlpp::ContinuousTransferFunction Rs = 0.5 + 2 / s;
    
    controlpp::FrequencyResponse L = Gb * Rs;
    std::ofstream("bode_L.csv") << L; // test by comparison

    controlpp::FrequencyResponse L2 = (1 + L);
    std::ofstream("bode_L2.csv") << L2; // test by comparison

    controlpp::FrequencyResponse Try = L / (1 + L);
    std::ofstream("bode_Try.csv") << Try; // test by comparison

    controlpp::FrequencyResponse Tdy = 1 / (1 + L);
    std::ofstream("bode_Tdy.csv") << Tdy; // test by comparison

    controlpp::TimeSeries ImpulseResponse = to_time_series(Gb);
    std::ofstream("ImpulseResponse.csv") << ImpulseResponse;

    controlpp::TimeSeries StepResponse = to_time_series(Gb * (1/s));
    std::ofstream("StepResponse.csv") << StepResponse;
}