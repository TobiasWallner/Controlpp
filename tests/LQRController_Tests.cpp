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
    
    constexpr int NStates = 2;

    const auto Pss = to_state_space(P);
    const auto Pz = discretise_zoh(Pss, Ts);
    
    // LQR gain
    const Eigen::Matrix<double, 1, NStates> K = controlpp::lqr(Pz);

    // LQR feed forward
    const auto F = lqr_feed_forward(Pz, K);

    controlpp::DssFilter Pf(Pz);

    std::ofstream file("data.csv");
 
    file << "time, uk, uf, u, x1, x2, y, set_point" << std::endl;
    double y = 0;

    Eigen::Vector<double, NStates> x;
    x.setZero();
    for(double time = 0; time < 0.1; time+=Ts){
        const double set_point = (time < 0.01) ? 0.0 : 1.0;
        const double uk = ((-K * x)).eval()(0,0);
        const double uf = (F * set_point).eval()(0,0);
        const double u = uk + uf;
        y = Pf(u);
        x = Pf.states();
        file << time << ", " << uk << ", " << uf << ", " << u << ", " << x(0) << ", " << x(1) << ", " << y << ", " << set_point << std::endl;
    }

    // TODO: somehow verify the results

}