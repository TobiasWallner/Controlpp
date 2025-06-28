#include <iostream>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp.hpp>

TEST(Estimators, least_squares){
    const double a0 = 1;
    const double a1 = 3;
    const double a2 = -5;

    // construct noisy data (measurement)
    const auto x = Eigen::VectorXd::LinSpaced(100, -10, 10);
    const auto xpow2 = x.cwiseProduct(x);
    const auto y = (a0 + x.array() * a1 + xpow2.array() * a2).matrix().eval();
    const auto data = (y + Eigen::VectorXd::Random(y.size())).eval();

    // construct solution matrix 
    Eigen::Matrix<double, Eigen::Dynamic, 3> X(data.size(), 3);
    X.col(0).setOnes();
    X.col(1) = x;
    X.col(2) = xpow2;

    // function under test
    // estimate parameters from noisy data
    const Eigen::Vector<double, 3> param_a = controlpp::least_squares(X, data);

    // check
    ASSERT_NEAR(param_a(0), a0, 0.1);
    ASSERT_NEAR(param_a(1), a1, 0.1);
    ASSERT_NEAR(param_a(2), a2, 0.1);
}

TEST(Estimators, reccursive_least_squares_1memory){
    const double a0 = 1;
    const double a1 = 3;
    const double a2 = -5;
  
    // construct noisy data (measurement)
    const auto x = Eigen::VectorXd::LinSpaced(100, -10, 10);
    const auto xpow2 = x.cwiseProduct(x);
    const auto y = (a0 + x.array() * a1 + xpow2.array() * a2).matrix().eval();
    const auto data = (y + Eigen::VectorXd::Random(y.size())).eval();

    // construct solution matrix
    Eigen::Matrix<double, Eigen::Dynamic, 3> X(data.size(), 3);
    X.col(0).setOnes();
    X.col(1) = x;
    X.col(2) = xpow2;

    controlpp::ReccursiveLeastSquares<double, 3> estimator;

    for(int i = 0; i < data.size(); ++i){
        Eigen::Vector<double, 3> x_ = X.row(i).eval();
        double y_ = data(i);
        estimator.add(y_, x_);
    }

    
    // check
    ASSERT_NEAR(estimator.estimate()(0), a0, 0.1);
    ASSERT_NEAR(estimator.estimate()(1), a1, 0.1);
    ASSERT_NEAR(estimator.estimate()(2), a2, 0.1);

}

TEST(Estimators, DTFEstimator){
    const auto s = controlpp::tf::s<double>;

    const auto Gs = (1) / (1 + s);

    const auto Sz = controlpp::s_to_z(controlpp::to_StateSpace(Gs), 0.1);

    controlpp::DiscreteStateSpaceFilter dssf(Sz);
    controlpp::DTFEstimator<double, 1, 2> dtf_est;

    for(int i = 0; i < 100; ++i){
        const double u = (i == 0) ? 0.0 : 1.0;
        const double y = dssf.input(u);
        dtf_est.add(y, u);
    }

    const auto Gz_est = dtf_est.estimate();
    std::cout << "Gz_est:\n" << Gz_est << std::endl;

    // ASSERT_NEAR(Gz_est.num()[0], 0.09529, 0.005);
    // ASSERT_NEAR(Gz_est.den()[0], 0.9048, 0.005);
    // ASSERT_NEAR(Gz_est.den()[1], 1.0, 0.005);

    const auto Sz_est = to_StateSpace(Gz_est);
    controlpp::DiscreteStateSpaceFilter dssf_est(Sz_est);
    dssf.clear();
    for(int i = 0; i < 100; ++i){
        const double u = (i == 0) ? 0.0 : 1.0;
        const double y = dssf.input(u);
        const double y_est = dssf_est.input(u);
        std::cout << "u: " << u << ", y: " << y << ", y_est: " << y_est << std::endl;
    }

}