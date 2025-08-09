#include <iostream>
#include <cmath>
#include <random>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/Pade.hpp>

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

    ASSERT_NEAR(a0/a4, 840, 1e-9);
    ASSERT_NEAR(a1/a4, 480, 1e-9);
    ASSERT_NEAR(a2/a4, 120, 1e-9);
    ASSERT_NEAR(a3/a4, 16, 1e-9);

    ASSERT_NEAR(b0/a4, 840, 1e-9);
    ASSERT_NEAR(b1/a4, 360, 1e-9);
    ASSERT_NEAR(b2/a4, 60, 1e-9);
    ASSERT_NEAR(b3/a4, 4, 1e-9);
}