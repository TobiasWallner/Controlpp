// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/transformations.hpp>

TEST(Transformation, continuous_ss_to_discrete_ss){
    // preparation
    const controlpp::ContinuousTransferFunction<double, 1, 0> s = controlpp::tf::s<double>;
    const controlpp::ContinuousTransferFunction<double, 0, 2> G_s = 1 / (1 + s + s*s);
    const auto Sys_s = controlpp::to_state_space(G_s);
    const double sample_time = 0.01; // seconds
    
    // expected values
    const Eigen::Matrix<double, 2, 2> expected_A({
        {0.99995, -0.00995},
        {0.00995, 0.99}
    });

    const Eigen::Matrix<double, 2, 1> expected_B({
        {0.009999},
        {4.983333e-05}
    });

    const Eigen::Matrix<double, 1, 2> expected_C({
        {0, 1},
    });

    const Eigen::Matrix<double, 1, 1> expected_D({
        {0},
    });

    // to test
    const auto Sys_z = controlpp::discretise_zoh(Sys_s, sample_time);
   
    // check
    ASSERT_TRUE(Sys_z.A().isApprox(expected_A, 0.001)) << "Sys_z.A():\n" << Sys_z.A() << "\n" << "expected_A\n" << expected_A;
    ASSERT_TRUE(Sys_z.B().isApprox(expected_B, 0.001)) << "Sys_z.B():\n" << Sys_z.B() << "\n" << "expected_B\n" << expected_B;
    ASSERT_TRUE(Sys_z.C().isApprox(expected_C, 0.001)) << "Sys_z.C():\n" << Sys_z.C() << "\n" << "expected_C\n" << expected_C;;
    ASSERT_TRUE(Sys_z.D().isApprox(expected_D, 0.001)) << "Sys_z.D():\n" << Sys_z.D() << "\n" << "expected_D\n" << expected_D;;
}

TEST(Transformation, state_space_to_transfer_function){
    const Eigen::Matrix<double, 2, 2> A({
        {0, 1},
        {-1, -3}
    });
    const Eigen::Matrix<double, 2, 1> B(
        0,
        1
    );
    const Eigen::Matrix<double, 1, 2> C(1, 1);
    const Eigen::Matrix<double, 1, 1> D(0);

    const controlpp::StateSpace sys(A, B, C, D);

    const auto G = controlpp::to_transfer_function(sys);

    ASSERT_NEAR(G.num(0), 1.0, 1e-6);
    ASSERT_NEAR(G.num(1), 1.0, 1e-6);
    ASSERT_NEAR(G.num(2), 0.0, 1e-6);
    
    ASSERT_NEAR(G.den(0), 1.0, 1e-6);
    ASSERT_NEAR(G.den(1), 3.0, 1e-6);
    ASSERT_NEAR(G.den(2), 1.0, 1e-6);
}