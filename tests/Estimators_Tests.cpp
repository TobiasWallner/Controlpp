#include <iostream>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp.hpp>

TEST(Estimators, least_squares){

    const double a0 = 1;
    const double a1 = 3;
    const double a2 = -5;

    const auto x = Eigen::VectorXd::LinSpaced(100, -10, 10);
    const auto xpow2 = x.cwiseProduct(x);
    const auto y = (a0 + x.array() * a1 + xpow2.array() * a2).matrix().eval();
    const auto data = (y + Eigen::VectorXd::Random(y.size())).eval();

    Eigen::Matrix<double, Eigen::Dynamic, 3> X(data.size(), 3);
    X.col(0).setOnes();
    X.col(1) = x;
    X.col(2) = xpow2;

    const Eigen::Vector<double, 3> param_a = controlpp::least_squares(X, data);

    ASSERT_NEAR(param_a(0), a0, 0.1);
    ASSERT_NEAR(param_a(1), a1, 0.1);
    ASSERT_NEAR(param_a(2), a2, 0.1);
}