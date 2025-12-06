#include <controlpp.hpp>

#include <gtest/gtest.h>

#include <fstream>

/*
    Tests if one can perform loop shaping from measured data
*/
TEST(Bode, Calculation){
    const auto s = controlpp::tf::s<double>;

    auto Gs = 1 / (1 + s + s * s);
    auto Gb = controlpp::bode(Gs);
    std::ofstream("bode_G.csv") << Gb;

    auto Rs = 0.5 + 2 / s;
    
    auto L = Gb * Rs;
    std::ofstream("bode_L.csv") << L; // test by comparison

    auto Try = L / (1 + L);
    std::ofstream("bode_Try.csv") << Try; // test by comparison

    auto Tdy = 1 / (1 + L);
    std::ofstream("bode_Tdy.csv") << Tdy; // test by comparison

    
}