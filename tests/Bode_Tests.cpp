#include <controlpp.hpp>

#include <gtest/gtest.h>

#include <fstream>

/*
    Tests if one can perform loop shaping from measured data
*/
TEST(Bode, Calculation){
    const controlpp::ContinuousTransferFunction s = controlpp::tf::s<double>;

    controlpp::ContinuousTransferFunction Gs = 1 / (1 + s + s * s);
    controlpp::Bode Gb = controlpp::bode(Gs);
    std::ofstream("bode_G.csv") << Gb;

    controlpp::ContinuousTransferFunction Rs = 0.5 + 2 / s;
    
    controlpp::Bode L = Gb * Rs;
    std::ofstream("bode_L.csv") << L; // test by comparison

    controlpp::Bode Try = L / (1 + L);
    std::ofstream("bode_Try.csv") << Try; // test by comparison

    controlpp::Bode Tdy = 1 / (1 + L);
    std::ofstream("bode_Tdy.csv") << Tdy; // test by comparison


}