#include <cmath>
#include <random>
#include <numbers>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/math.hpp>

TEST(math, solve_riccati){
    const Eigen::Matrix<double, 2, 2> A({
        {0, 1},
        {1, 1}
    });

    const Eigen::Matrix<double, 2, 1> B({
        {0},
        {1}
    });

    const Eigen::Matrix<double, 2, 2> Q({
        {1, 0},
        {0, 1}
    });

    const Eigen::Matrix<double, 1, 1> R({
        {1}
    });

    Eigen::Matrix<double, 2, 2> X = controlpp::solve_continuous_lqr_riccati(
        A, B, R, Q);
    Eigen::Matrix<double, 2, 2> care = A.transpose() * X + X * A - X * B * R.inverse() * B.transpose() * X + Q;

    ASSERT_NEAR(care(0, 0), 0, 1e-10);
    ASSERT_NEAR(care(1, 0), 0, 1e-10);
    ASSERT_NEAR(care(0, 1), 0, 1e-10);
    ASSERT_NEAR(care(1, 1), 0, 1e-10);
}

TEST(math, phase_unwrap_rad_positive){
    using namespace std::numbers;
    Eigen::VectorXd phases(26);
    phases 
        << 0.0, (0.2*pi), (0.4*pi), (0.6*pi), (0.8*pi), (pi)
        , (-pi * 0.8), (-pi * 0.6), (-pi * 0.4), (-pi * 0.2), 0.0, (0.2*pi), (0.4*pi), (0.6*pi), (0.8*pi), (pi)
        , (-pi * 0.8), (-pi * 0.6), (-pi * 0.4), (-pi * 0.2), 0.0, (0.2*pi), (0.4*pi), (0.6*pi), (0.8*pi), (pi);

    Eigen::VectorXd expected_result(26);
    expected_result 
        << 0.0, (0.2*pi), (0.4*pi), (0.6*pi), (0.8*pi), (1.0 * pi)
        , (1.2*pi), (1.4*pi), (1.6*pi), (1.8*pi), (2.0 * pi), (2.2*pi), (2.4*pi), (2.6*pi), (2.8*pi), (3.0 * pi)
        , (3.2*pi), (3.4*pi), (3.6*pi), (3.8*pi), (4.0 * pi), (4.2*pi), (4.4*pi), (4.6*pi), (4.8*pi), (5.0 * pi);

    Eigen::VectorXd unwrapped_phases = controlpp::phase_unwrap_rad(phases);

    for(int i = 0; i < phases.size(); ++i){
        ASSERT_NEAR(unwrapped_phases(i), expected_result(i), 1e-9) << "at index: " << i;
    }
}