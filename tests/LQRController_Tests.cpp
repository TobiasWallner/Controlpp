#include <numbers>
#include <fstream>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/transformations.hpp>
#include <controlpp/DiscreteFilter.hpp>
#include <controlpp/LQRController.hpp>


TEST(LQRController, for_PT2){
    using namespace controlpp;
    const auto s = tf::s<double>;
    const double Ts = 1./1'000.;

    // plant description
    const double f_p = 100;
    const double omega_p = f_p * 2 * std::numbers::pi;
    const double K_p = 11;
    const double D = 0.2;
    const auto P = 
                            K_p 
    / //-------------------------------------------------- 
    (1 + 2 * D * s / omega_p + s * s / (omega_p * omega_p));
    
    // constexpr int NStates = 2;

    const auto Pss = to_state_space(P);
    const auto Pz = discretise_zoh(Pss, Ts);
    
    // LQR gain
    const Eigen::Matrix<double, 1, 2> K = controlpp::lqr(
        Pz, 
        Eigen::Matrix<double, 2, 2>::Identity().eval(), 
        Eigen::Matrix<double, 1, 1>::Identity().eval()
    );

    EXPECT_NEAR(K(0), 0.00020199, 1e-6);
    EXPECT_NEAR(K(1), -0.08006232, 1e-6);

    // LQR feed forward
    const auto F = lqr_feed_forward(Pz, K);
    EXPECT_NEAR(F(0), 0.0616135, 1e-6);
}