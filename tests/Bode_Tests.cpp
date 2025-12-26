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
    controlpp::write_csv(std::ofstream("bode_G.csv"), Gb);

    std::ofstream("step_G.csv") << controlpp::step(Gs);

    controlpp::ContinuousTransferFunction Rs = 0.5 + 2 / s;
    
    controlpp::Bode L = Gb * Rs;
    controlpp::write_csv(std::ofstream("bode_L.csv"), L); // test by comparison

    controlpp::Bode L2 = (1 + L);
    controlpp::write_csv(std::ofstream("bode_L2.csv"), L2); // test by comparison

    controlpp::Bode Try = L / (1 + L);
    controlpp::write_csv(std::ofstream("bode_Try.csv"),Try); // test by comparison

    controlpp::Bode Tdy = 1 / (1 + L);
    controlpp::write_csv(std::ofstream("bode_Tdy.csv"), Tdy); // test by comparison

    controlpp::TimeSeries ImpulseResponse = controlpp::impulse(Gb);
    controlpp::write_csv(std::ofstream("ImpulseResponse.csv"), ImpulseResponse);

    controlpp::write_csv(std::ofstream("G_1_s.csv"), (Gb * (1/s)));

    controlpp::TimeSeries StepResponse = controlpp::step(Gb);
    controlpp::write_csv(std::ofstream("StepResponse.csv"), StepResponse);
}

TEST(Bode, set_magnitudes_dB_and_phases_deg){

    Eigen::Vector<double, Eigen::Dynamic> mags_dB(3);
    mags_dB << 1.0, 1.0, 1.0;

    Eigen::Vector<double, Eigen::Dynamic> phases_deg(3);
    phases_deg << 0.0, 1.0, 2.0;

    controlpp::Bode bode;
    bode.set_magnitudes_dB_and_phases_deg(mags_dB, phases_deg);

}

TEST(Bode, read_csv){
    using namespace controlpp;

    std::stringstream csv;
    
    const controlpp::ContinuousTransferFunction Gs = 1 / (1 + s + s * s);
    const controlpp::Bode Gb = controlpp::bode(Gs);

    std::stringstream stream;

    write_csv(stream, Gb);

    controlpp Gb = read_csv(stream);
}