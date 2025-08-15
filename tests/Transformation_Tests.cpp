// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp.hpp>

TEST(Transformation, continuous_ss_to_discrete_ss){
    // preparation
    const controlpp::ContinuousTransferFunction<double, 1, 0> s = controlpp::tf::s<double>;
    const controlpp::ContinuousTransferFunction<double, 0, 2> G_s = 1 / (1 + s + s*s);
    const auto Sys_s = controlpp::to_state_space(G_s);
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
    ASSERT_TRUE(Sys_z.A().isApprox(expected_A, 0.01));
    ASSERT_TRUE(Sys_z.B().isApprox(expected_B, 0.01));
    ASSERT_TRUE(Sys_z.C().isApprox(expected_C, 0.01));
    ASSERT_TRUE(Sys_z.D().isApprox(expected_D, 0.01));
}

TEST(Transformation, discrete_ss_to_bilinear_ss){
    // preparation
    const auto z = controlpp::tf::z<double>;
    const auto G_z = (1 + z) / (1 + 3*z + z*z);
    const auto Sys_z = controlpp::to_state_space(G_z);

    // function under test
    const auto Sys_q = controlpp::discrete_to_bilinear(Sys_z);

    // expected values
    const Eigen::Matrix<double, 2, 2> expected_A({
        {-3, -2},
        {2, 3}
    });

    const Eigen::Matrix<double, 2, 1> expected_B({
        {std::sqrt(2.0)},
        {-std::sqrt(2.0)}
    });

    const Eigen::Matrix<double, 1, 2> expected_C({
        {std::sqrt(2.0), 0},
    });

    const Eigen::Matrix<double, 1, 1> expected_D({
        {0},
    });

    // check
    ASSERT_TRUE(Sys_q.A().isApprox(expected_A, 0.01));
    ASSERT_TRUE(Sys_q.B().isApprox(expected_B, 0.01));
    ASSERT_TRUE(Sys_q.C().isApprox(expected_C, 0.01));
    ASSERT_NEAR(Sys_q.D()(0, 0), Sys_q.D()(0, 0), 0.01);
}

TEST(Transformation, continuous_ss_to_bilinear_ss){
    // preparation
    const auto z = controlpp::tf::z<double>;
    const auto G_z = (1 + z) / (1 + 3*z + z*z);
    const auto Sys_z = controlpp::to_state_space(G_z);
    const auto Sys_q = controlpp::discrete_to_bilinear(Sys_z);

    // function under test
    const auto Sys_z2 = controlpp::bilinear_to_discrete(Sys_q);

    // check
    if(Sys_z2.A().isApprox(Sys_z.A(), 0.01)){
        ASSERT_TRUE(Sys_z2.A().isApprox(Sys_z.A(), 0.01));
        ASSERT_TRUE(Sys_z2.B().isApprox(Sys_z.B(), 0.01));
        ASSERT_TRUE(Sys_z2.C().isApprox(Sys_z.C(), 0.01));
        ASSERT_NEAR(Sys_z2.D()(0, 0), Sys_z.D()(0, 0), 0.01);
    }else{
        // maybe signs have flipped
        ASSERT_TRUE(Sys_z2.A().isApprox((-Sys_z.A()).eval(), 0.01));
        ASSERT_TRUE(Sys_z2.B().isApprox((-Sys_z.B()).eval(), 0.01));
        ASSERT_TRUE(Sys_z2.C().isApprox((-Sys_z.C()).eval(), 0.01));
        ASSERT_NEAR(Sys_z2.D()(0, 0), -Sys_z.D()(0, 0), 0.01);
    }

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