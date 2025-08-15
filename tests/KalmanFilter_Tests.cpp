#include <cmath>
#include <random>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp.hpp>

TEST(KalmanFilter, noisy_speed_estimate_pos_acc){
    const double Ts = 0.01;

    const Eigen::Matrix<double, 2, 2> A({
        {1, Ts},
        {0, 1},
    });
    const Eigen::Matrix<double, 1, 2> H(1, 0);

    auto kalman = controlpp::KalmanFilter(A, H);
    std::default_random_engine noise_engine;
    std::normal_distribution noise(0.0, 0.5);

    for(double t = 0; t < 2; t += Ts){
        const double dy = 1.6;
        const double y = 3 + dy * t;
        const double noisy_y = y + noise(noise_engine);

        kalman.add(noisy_y);

        Eigen::Vector<double, 2> p = kalman.estimate();

        // give the filter some time to learn
        // then check if the estimate is useful
        if(t > 1){
            ASSERT_NEAR(y, p(0), 0.2);
            ASSERT_NEAR(dy, p(1), 0.3);
        }
    }

}