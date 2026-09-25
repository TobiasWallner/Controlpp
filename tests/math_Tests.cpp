#include <cmath>
#include <random>
#include <numbers>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/math.hpp>

TEST(math, join_to_diagonal){
    Eigen::Vector2d left;
    left << 2.0, 3.0;
    Eigen::Vector<double, 1> right;
    right << 5.0;

    const auto result = controlpp::join_to_diagonal(left, right);
    Eigen::Matrix3d expected = Eigen::Matrix3d::Zero();
    expected.diagonal() << 2.0, 3.0, 5.0;
    EXPECT_EQ(result, expected);
}

TEST(math, solve_riccati){
    // TODO: this
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

    Eigen::VectorXd unwrapped_phases = controlpp::unwrap_rad(phases);

    for(int i = 0; i < phases.size(); ++i){
        ASSERT_NEAR(unwrapped_phases(i), expected_result(i), 1e-9) << "at index: " << i;
    }
}
