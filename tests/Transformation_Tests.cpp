#include <iostream>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp.hpp>

TEST(Transformation, continuous_ss_to_discrete_ss){
    // preparation
    const auto s = controlpp::tf::s<double>;
    const auto G_s = 1 / (1 + s + s*s);
    const auto Sys_s = controlpp::to_ContinuousStateSpace(G_s);
    const double sample_time = 0.01; // seconds
    
    // expected values
    const Eigen::Matrix<double, 2, 2> expected_A({
        {0.99995, 0.00995},
        {-0.00995, 0.99}
    });

    const Eigen::Matrix<double, 2, 1> expected_B({
        {4.98333e-05},
        {0.00995}
    });

    const Eigen::Matrix<double, 1, 2> expected_C({
        {1, 0},
    });

    const Eigen::Matrix<double, 1, 1> expected_D({
        {0},
    });

    // to test
    const auto Sys_z = controlpp::continuous_to_discrete(Sys_s, sample_time);
   
    // check
    Sys_z.A().isApprox(expected_A, 0.01);
    Sys_z.B().isApprox(expected_B, 0.01);
    Sys_z.C().isApprox(expected_C, 0.01);
    Sys_z.D().isApprox(expected_D, 0.01);
}

TEST(Transformation, continuous_tf_to_discrete_tf){
    // preparation
    const auto s = controlpp::tf::s<double>;
    const auto G_s = 1 / (1 + s + s*s);
    const auto Sys_s = controlpp::to_ContinuousStateSpace(G_s);
    const double sample_time = 0.01; // seconds
    const auto Sys_z = controlpp::continuous_to_discrete(Sys_s, sample_time);
    
    // transform from state space to transfer function
    const auto z = controlpp::tf::z<double>;
    const auto I = controlpp::identity_like(Sys_s.A());

    const auto temp = Sys_s.C() * Eigen::Inverse(I * z - Sys_z.A()) * Sys_s.B() + Sys_s.D();
    
    
    // calculate: C * (I * z - A)^{-1} * B + D
    const auto X = I * z - Sys_z.A();
    const auto Y = Eigen::Inverse(X).eval();
    std::cout << "Y(0,0): " << Y(0,0) << std::endl;
    // const auto Z = (Sys_z.C() * Y).eval();

    //const auto G_z = Sys_s.C() * Eigen::inverse(Iz - Sys_s.A()) * Sys_s.B() + Sys_s.D();

}