#include <numbers>
#include <fstream>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/transformations.hpp>
#include <controlpp/DiscreteFilter.hpp>
#include <controlpp/H2Controller.hpp>
#include <controlpp/Pade.hpp>

TEST(H2Controller, for_PT1){
    using namespace controlpp;
    const auto s = tf::s<double>;
    const double Ts = 1./4000.;

    // plant description
    const double f_p = 62;
    const double omega_p = f_p * 2 * std::numbers::pi;
    const double D_p = 0.2;
    const double K_p = 11;
    const auto P = K_p/(1 + 2 * D_p * s / omega_p + (s * s) / (omega_p * omega_p));

    // measurement description
    const double f_m = 500;
    const double omega_m = f_m * 2 * std::numbers::pi;
    const auto M = 1 / (1 + s / omega_m);

    // disturbance shape
    const auto Wd = 1 / (1 + s / (2 * std::numbers::pi * 200));

    // performance output weight
    const auto Wz = 1 / (s*(1 + s / (2 * std::numbers::pi * 200)));

    const double control_penalty = 500000.0;
    const double measurement_noise = 0.001;

    const auto H2 = continous_h2_controller(P, M, Wd, Wz, control_penalty, measurement_noise);
    
    const auto Pz = s_to_z(to_state_space(P), Ts);
    const auto Mz = s_to_z(to_state_space(M), Ts);
    const auto H2z = s_to_z(H2, Ts);

    DssFilter Pf(Pz);
    DssFilter Mf(Mz);
    DssFilter H2f(H2z);

    std::ofstream file("data.csv");

    file << "time, r, err, u, z, y" << std::endl;
    double y = 0;
    for(double time = 0; time < 0.1; time+=Ts){
        const double r = 0.;
        const double err = r - y;
        const double u = H2f(-err);
        const double z = Pf(u) + ((time > 0.01)? 1. : 0.);
        y = Mf(z);
        file << time << ", " << r << ", " << err << ", " << u << ", " << z << ", " << y << std::endl;
    }

}