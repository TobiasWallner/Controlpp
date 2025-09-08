// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/ContinuousTransferFunction.hpp>
#include <controlpp/DiscreteTransferFunction.hpp>
#include <controlpp/DiscreteStateSpace.hpp>
#include <controlpp/DiscreteFilter.hpp>
#include <controlpp/transformations.hpp>
#include <controlpp/Estimators.hpp>

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
    ASSERT_NEAR(estimator.estimate()(0), a0, 0.2);
    ASSERT_NEAR(estimator.estimate()(1), a1, 0.2);
    ASSERT_NEAR(estimator.estimate()(2), a2, 0.2);

}

TEST(Estimators, DtfEstimator){
    const auto s = controlpp::tf::s<double>;

    const auto Gs = (8) / (1 + 3*s);

    const auto Sz = controlpp::discretise_zoh(controlpp::to_state_space(Gs), 0.1);

    controlpp::DssFilter dssf(Sz);
    controlpp::DtfEstimator<double, 1, 1> dtf_est;

    for(int i = 0; i < 100; ++i){
        const double u = (i == 0) ? 0.0 : 1.0;
        const double y = dssf.input(u);
        dtf_est.add(y, u);
    }

    const auto Gz_est = dtf_est.estimate();

    const auto Sz_est = to_state_space(Gz_est);
    controlpp::DssFilter dssf_est(Sz_est);
    dssf.clear();

    // check by comparing step responses
    for(int i = 0; i < 100; ++i){
        const double u = 1.0;
        const double y = dssf.input(u);
        const double y_est = dssf_est.input(u);
        //std::cout << "y: " << y << ", y_est: " << y_est << std::endl;
        ASSERT_NEAR(y, y_est, 0.01);
    }

}