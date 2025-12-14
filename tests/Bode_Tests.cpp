#include <fstream>
#include <iostream>

#include <controlpp.hpp>

#include <gtest/gtest.h>


/*
    Tests if one can perform loop shaping from measured data
*/
TEST(Bode, Calculation){
    const controlpp::ContinuousTransferFunction s = controlpp::tf::s<double>;

    controlpp::ContinuousTransferFunction Gs = 1 / (1 + s + s * s);
    controlpp::Bode Gb = controlpp::bode(Gs);
    std::ofstream("bode_G.csv") << Gb;

    std::ofstream("step_G.csv") << controlpp::step(Gs);

    controlpp::ContinuousTransferFunction Rs = 0.5 + 2 / s;
    
    controlpp::Bode L = Gb * Rs;
    std::ofstream("bode_L.csv") << L; // test by comparison

    controlpp::Bode L2 = (1 + L);
    std::ofstream("bode_L2.csv") << L2; // test by comparison

    controlpp::Bode Try = L / (1 + L);
    std::ofstream("bode_Try.csv") << Try; // test by comparison

    controlpp::Bode Tdy = 1 / (1 + L);
    std::ofstream("bode_Tdy.csv") << Tdy; // test by comparison

    controlpp::TimeSeries ImpulseResponse = controlpp::impulse(Gb);
    std::ofstream("ImpulseResponse.csv") << ImpulseResponse;

    std::ofstream("G_1_s.csv") << (Gb * (1/s));

    controlpp::TimeSeries StepResponse = controlpp::step(Gb);
    std::ofstream("StepResponse.csv") << StepResponse;
}

TEST(Bode, set_magnitudes_dB_and_phases_deg){

    Eigen::Vector<double, Eigen::Dynamic> mags_dB(3);
    mags_dB << 1.0, 1.0, 1.0;

    Eigen::Vector<double, Eigen::Dynamic> phases_deg(3);
    phases_deg << 0.0, 1.0, 2.0;

    controlpp::Bode bode;
    bode.set_magnitudes_dB_and_phases_deg(mags_dB, phases_deg);

}