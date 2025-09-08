#include <cmath>
#include <random>
#include <fstream>
// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/Pade.hpp>
#include <controlpp/transformations.hpp>
#include <controlpp/analysis.hpp>

TEST(Pade, pade_coefficients_1_1){
    const double NumOrder = 1;
    const double DenOrder = 1;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);

    ASSERT_NEAR(a0/a1, 2, 1e-9);

    ASSERT_NEAR(b0/a1, 2, 1e-9);
    ASSERT_NEAR(b1/a1, 1, 1e-9);
}

TEST(Pade, pade_coefficients_2_2){
    const double NumOrder = 2;
    const double DenOrder = 2;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);

    ASSERT_NEAR(a0/a2, 12, 1e-9);
    ASSERT_NEAR(a1/a2, 6, 1e-9);

    ASSERT_NEAR(b0/a2, 12, 1e-9);
    ASSERT_NEAR(b1/a2, 6, 1e-9);
    ASSERT_NEAR(b2/a2, 1, 1e-9);
}

TEST(Pade, pade_coefficients_3_3){
    const double NumOrder = 3;
    const double DenOrder = 3;
    
    const double a0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);
    const double a3 = controlpp::pade_num_param(NumOrder, DenOrder, 3);

    const double b0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);
    const double b3 = controlpp::pade_den_param(NumOrder, DenOrder, 3);

    ASSERT_NEAR(a0/a3, 120, 1e-9);
    ASSERT_NEAR(a1/a3, 60, 1e-9);
    ASSERT_NEAR(a2/a3, 12, 1e-9);

    ASSERT_NEAR(b0/a3, 120, 1e-9);
    ASSERT_NEAR(b1/a3, 60, 1e-9);
    ASSERT_NEAR(b2/a3, 12, 1e-9);
    ASSERT_NEAR(b3/a3, 1, 1e-9);
}

TEST(Pade, pade_coefficients_4_4){
    const double NumOrder = 4;
    const double DenOrder = 4;
    
    const double a0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);
    const double a3 = controlpp::pade_num_param(NumOrder, DenOrder, 3);
    const double a4 = controlpp::pade_num_param(NumOrder, DenOrder, 4);

    const double b0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);
    const double b3 = controlpp::pade_den_param(NumOrder, DenOrder, 3);
    const double b4 = controlpp::pade_den_param(NumOrder, DenOrder, 4);

    ASSERT_NEAR(a0/a4, 1680, 1e-9);
    ASSERT_NEAR(a1/a4, 840, 1e-9);
    ASSERT_NEAR(a2/a4, 180, 1e-9);
    ASSERT_NEAR(a3/a4, 20, 1e-9);

    ASSERT_NEAR(b0/a4, 1680, 1e-9);
    ASSERT_NEAR(b1/a4, 840, 1e-9);
    ASSERT_NEAR(b2/a4, 180, 1e-9);
    ASSERT_NEAR(b3/a4, 20, 1e-9);
    ASSERT_NEAR(b4/a4, 1, 1e-9);
}

// ----------------------------------------------------------------------------

TEST(Pade, pade_coefficients_0_1){
    const double NumOrder = 0;
    const double DenOrder = 1;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);

    ASSERT_NEAR(a0/a1, 1, 1e-9);

    ASSERT_NEAR(b0/a1, 1, 1e-9);
}

TEST(Pade, pade_coefficients_1_2){
    const double NumOrder = 1;
    const double DenOrder = 2;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);

    ASSERT_NEAR(a0/a2, 6, 1e-9);
    ASSERT_NEAR(a1/a2, 4, 1e-9);

    ASSERT_NEAR(b0/a2, 6, 1e-9);
    ASSERT_NEAR(b1/a2, 2, 1e-9);
}

TEST(Pade, pade_coefficients_2_3){
    const double NumOrder = 2;
    const double DenOrder = 3;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);
    const double a3 = controlpp::pade_den_param(NumOrder, DenOrder, 3);

    ASSERT_NEAR(a0/a3, 60, 1e-9);
    ASSERT_NEAR(a1/a3, 36, 1e-9);
    ASSERT_NEAR(a2/a3, 9, 1e-9);

    ASSERT_NEAR(b0/a3, 60, 1e-9);
    ASSERT_NEAR(b1/a3, 24, 1e-9);
    ASSERT_NEAR(b2/a3, 3, 1e-9);
}

TEST(Pade, pade_coefficients_3_4){
    const double NumOrder = 3;
    const double DenOrder = 4;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);
    const double b3 = controlpp::pade_num_param(NumOrder, DenOrder, 3);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);
    const double a3 = controlpp::pade_den_param(NumOrder, DenOrder, 3);
    const double a4 = controlpp::pade_den_param(NumOrder, DenOrder, 4);

    ASSERT_NEAR(b0/a4, 840, 1e-9);
    ASSERT_NEAR(b1/a4, 360, 1e-9);
    ASSERT_NEAR(b2/a4, 60, 1e-9);
    ASSERT_NEAR(b3/a4, 4, 1e-9);

    ASSERT_NEAR(a0/a4, 840, 1e-9);
    ASSERT_NEAR(a1/a4, 480, 1e-9);
    ASSERT_NEAR(a2/a4, 120, 1e-9);
    ASSERT_NEAR(a3/a4, 16, 1e-9);

}

TEST(Pade, pade_polynom_1_1){
    const double Td = 0.1;
    controlpp::ContinuousTransferFunction<double, 1, 1> ctf = controlpp::pade<double, 1, 1>(Td);
    controlpp::ContinuousTransferFunction<double, 1, 1> expected({2, -Td}, {2, Td});

    // norm
    {
        const auto n = ctf.den(ctf.den().size()-1);
        ctf.num() /= n;
        ctf.den() /= n;
    }

    {
        const auto n = expected.den(expected.den().size()-1);
        expected.num() /= n;
        expected.den() /= n;
    }

    for(size_t i = 0; i < ctf.num().size(); ++i){
        ASSERT_NEAR(ctf.num(i), expected.num(i), 1e-9);
    }
    
    for(size_t i = 0; i < ctf.den().size(); ++i){
        ASSERT_NEAR(ctf.den(i), expected.den(i), 1e-9);
    }
}

TEST(Pade, pade_polynom_2_2){
    const double Td = 0.1;
    controlpp::ContinuousTransferFunction<double, 2, 2> ctf = controlpp::pade<double, 2, 2>(Td);
    controlpp::ContinuousTransferFunction<double, 2, 2> expected({12, -6*Td, Td*Td}, {12, 6*Td, Td*Td});

    // norm
    {
        const auto n = ctf.den(ctf.den().size()-1);
        ctf.num() /= n;
        ctf.den() /= n;
    }

    {
        const auto n = expected.den(expected.den().size()-1);
        expected.num() /= n;
        expected.den() /= n;
    }

    for(size_t i = 0; i < ctf.num().size(); ++i){
        ASSERT_NEAR(ctf.num(i), expected.num(i), 1e-9);
    }
    
    for(size_t i = 0; i < ctf.den().size() ; ++i){
        ASSERT_NEAR(ctf.den(i), expected.den(i), 1e-9);
    }
}

TEST(Pade, pade_polynom_3_4){

    const double Td = 0.1;
    controlpp::ContinuousTransferFunction<double, 3, 4> ctf = controlpp::pade<double, 3, 4>(Td);
    controlpp::ContinuousTransferFunction<double, 3, 4> expected({840, -360*Td, 60*Td*Td, -4*Td*Td*Td}, {840, 480*Td, 120*Td*Td, 16*Td*Td*Td, Td*Td*Td*Td});

    // norm
    {
        const auto n = ctf.den(ctf.den().size()-1);
        ctf.num() /= n;
        ctf.den() /= n;
    }

    {
        const auto n = expected.den(expected.den().size()-1);
        expected.num() /= n;
        expected.den() /= n;
    }

    for(size_t i = 0; i < ctf.num().size(); ++i){
        ASSERT_NEAR(ctf.num(i), expected.num(i), 1e-8);
    }
    
    for(size_t i = 0; i < ctf.den().size(); ++i){
        ASSERT_NEAR(ctf.den(i), expected.den(i), 1e-8);
    }
}

TEST(Pade, pade_3_4_step_response){
    const auto P = controlpp::pade<double, 3, 4>(600e-6);

    const auto [slowest_freq, fastest_freq] = slowest_fastest_frequencies(P);
    const double sample_time = (0.05) / fastest_freq;
    const double simulation_time = (20) / slowest_freq;
    std::cout << "sample_time: " << sample_time << std::endl;

    const auto Pss = controlpp::to_state_space(P);
    std::cout << "Pss:\n" << Pss << std::endl;

    const auto Pz = controlpp::discretise_zoh(Pss, sample_time);
    std::cout << "Pz:\n" << Pz << std::endl;

    const auto step = controlpp::step(P, sample_time, simulation_time);
    
    std::ofstream file("step.csv");
    file << step << std::endl;
}