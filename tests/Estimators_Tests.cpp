// google test
#include <gtest/gtest.h>

#include <fstream>

// controlpp
#include <controlpp/ContinuousTransferFunction.hpp>
#include <controlpp/DiscreteTransferFunction.hpp>
#include <controlpp/DiscreteStateSpace.hpp>
#include <controlpp/DiscreteFilter.hpp>
#include <controlpp/transformations.hpp>
#include <controlpp/Estimators.hpp>
#include <controlpp/TimeSeries.hpp>
#include <controlpp/analysis.hpp>

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
        estimator.input(y_, x_);
    }

    
    // check
    ASSERT_NEAR(estimator.estimate()(0), a0, 0.2);
    ASSERT_NEAR(estimator.estimate()(1), a1, 0.2);
    ASSERT_NEAR(estimator.estimate()(2), a2, 0.2);

}

TEST(Estimators, DtfEstimator){
    const auto s = controlpp::tf::s<double>;

    const auto Gs = (8) / (1 + 3*s);

//    const auto Sz = controlpp::discretise_zoh(controlpp::to_state_space(Gs), 0.1);

    controlpp::DssFilter dssf(Gs, 0.1, controlpp::EDiscretisation::zero_order_hold);
    controlpp::DtfEstimator<double, 1, 1> dtf_est;

    for(int i = 0; i < 100; ++i){
        const double u = (i == 0) ? 0.0 : 1.0;
        const double y = dssf.input(u);
        dtf_est.input(y, u);
        // std::cout << "iteration: " << i << std::endl;
        // std::cout << "- estimage:" << std::endl;
        // std::cout << dtf_est.estimate() << std::endl;
        // std::cout << "- covariance:" << std::endl;
        // std::cout << dtf_est.cov() << std::endl;
        // std::cout << "- gain:" << std::endl;
        // std::cout << dtf_est.gain().transpose() << '\n' << std::endl;
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
        // std::cout << "iteration: " << i << ", y: " << y << ", y_est: " << y_est << std::endl;
        ASSERT_NEAR(y, y_est, 0.05);
    }

}

TEST(Estimators, DtfEstimatorNaNInitialization){
    const Eigen::Vector2d num_uncertainty = Eigen::Vector2d::Constant(1000.0);
    const Eigen::Vector<double, 1> den_uncertainty = Eigen::Vector<double, 1>::Constant(1000.0);
    const auto covariance = controlpp::join_to_diagonal(num_uncertainty, den_uncertainty);

    const Eigen::Matrix3d expected = Eigen::Matrix3d::Identity() * 1000.0;
    ASSERT_EQ(covariance, expected);

    controlpp::ReccursiveLeastSquares<double, 3> estimator(
        Eigen::Vector3d::Zero(), covariance, 0.995);
    estimator.input(0.0, Eigen::Vector3d::Zero());

    EXPECT_TRUE(estimator.cov().allFinite());
    EXPECT_TRUE(estimator.gain().allFinite());
    EXPECT_TRUE(estimator.estimate().allFinite());
}

TEST(Estimators, DtfEstimatorSecondOrderHistory){
    const controlpp::DiscreteTransferFunction<double, 2, 2> hint({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
    controlpp::DtfEstimator<double, 2, 2> estimator(hint, 1000.0, 1.0);
    controlpp::ReccursiveLeastSquares<double, 5> reference(
        Eigen::Vector<double, 5>::Zero(),
        Eigen::Matrix<double, 5, 5>::Identity() * 1000.0,
        1.0);

    Eigen::Vector<double, 6> inputs;
    inputs << 0.0, 1.0, 0.5, -0.25, 0.75, 0.0;
    Eigen::Vector<double, 6> outputs;
    outputs << 0.0, 0.2, -0.4, 0.7, -0.1, 1.1;

    for(int i = 0; i < inputs.size(); ++i){
        Eigen::Vector<double, 5> regressor;
        regressor << inputs(i),
            i > 0 ? inputs(i - 1) : 0.0,
            i > 1 ? inputs(i - 2) : 0.0,
            i > 0 ? -outputs(i - 1) : 0.0,
            i > 1 ? -outputs(i - 2) : 0.0;

        reference.input(outputs(i), regressor);
        estimator.input(outputs(i), inputs(i));

        const auto estimated = estimator.estimate();
        for(int j = 0; j < 3; ++j){
            ASSERT_NEAR(estimated.num(j), reference.estimate()(j), 1e-10) << "sample " << i;
        }
        for(int j = 0; j < 2; ++j){
            ASSERT_NEAR(estimated.den(j + 1), reference.estimate()(j + 3), 1e-10) << "sample " << i;
        }
    }
}

TEST(Estimators, dft_estimate){
    const auto s = controlpp::tf::s<double>;

    const auto Gs = (8 + 1*s) / (1 + 3*s + 9*s*s);

    controlpp::TimeSeries ts = controlpp::step(Gs);

    Eigen::Vector<double, 10> init = Eigen::Vector<double, 10>::Zero();
    Eigen::VectorXd step(ts.size());
    step.setOnes();

    Eigen::VectorXd u(ts.size() + init.size());
    int mod = 2;
    for(int i = 0; i < u.size(); ++i) {
        u(i) = (i < init.size()) ? 0.0 : ((i % mod == 0) ? 1.1 : 0.9);
        ++mod;
    }

    Eigen::VectorXd y = controlpp::join_to_vector(init, ts.values());
    ASSERT_EQ(y.size(), init.size() + ts.size());

    const auto dest_err = controlpp::dft_estimate<double, 1, 2>(u, y, 1e-6);

    ASSERT_TRUE(dest_err.has_value());

    const auto dest = dest_err.value();

    const auto estimated_simulation = dest.eval(u);

    for(size_t i = 0; (i < ts.size()) && (i < static_cast<size_t>(estimated_simulation.size())); ++i){
        const double y_est = estimated_simulation(init.size() + i);
        const double y_true = ts.values(i);
        // std::cout << "iteration: " << i << ", y_true: " << y_true << ", y_est: " << y_est << std::endl;
        ASSERT_NEAR(y_true, y_est, 0.05);
    }
}
