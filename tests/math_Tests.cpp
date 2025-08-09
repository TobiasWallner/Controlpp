#include <iostream>
#include <cmath>
#include <random>

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

    Eigen::Matrix<double, 2, 2> X = controlpp::solve_continuous_riccati(A, B, Q, R);
    Eigen::Matrix<double, 2, 2> care = A.transpose() * X + X * A - X * B * R.inverse() * B.transpose() * X + Q;

    ASSERT_NEAR(care(0, 0), 0, 1e-10);
    ASSERT_NEAR(care(1, 0), 0, 1e-10);
    ASSERT_NEAR(care(0, 1), 0, 1e-10);
    ASSERT_NEAR(care(1, 1), 0, 1e-10);
}