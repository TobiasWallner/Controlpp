#include <numbers>
#include <fstream>

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/transformations.hpp>
#include <controlpp/DiscreteFilter.hpp>
#include <controlpp/H2Controller.hpp>

TEST(H2Controller, for_PT1){
    using namespace controlpp;

    const auto s = tf::s<double>;

    // plant description
    const double f_p = 62;
    const double omega_p = f_p * 2 * std::numbers::pi;
    const double D_p = 0.2;
    const double K_p = 11;
    const auto P_tf = K_p/(1 + 2 * D_p * s / omega_p + (s * s) / (omega_p * omega_p));

    const auto P = to_state_space(P_tf);

    // measurement description
    const double f_m = 210;
    const double omega_m = f_m * 2 * std::numbers::pi;
    const auto M = to_state_space(1 / (1 + s / omega_m));

    // Generalised system description
    ContinuousGeneralisedPlant<double, 3> Gss;
    Gss.A.setZero();
    Gss.A.topLeftCorner(P.A().rows(), P.A().cols()) = P.A();
    const auto BmCp = (M.B() * P.C()).eval();
    Gss.A.bottomLeftCorner(BmCp.rows(), BmCp.cols()) = BmCp;
    Gss.A.bottomRightCorner(M.A().rows(), M.A().cols()) = M.A();
    Gss.Bw.fill(1e-3); // some process noise is needed for it to be non triveal
    Gss.Bu.setZero();
    Gss.Bu.head(P.B().size()) = P.B();
    Gss.Cz.setZero();
    Gss.Cz.head(P.C().size()) = P.C();
    Gss.Cy.setZero();
    Gss.Cy.tail(M.C().size()) = M.C();
    Gss.Duz.setZero();
    Gss.Dwy.setZero();
    
    const auto H2 = continous_H2_controller(
            Gss.A, Gss.Bw, Gss.Bu,
            Gss.Cz, Gss.Duz,
            Gss.Cy, Gss.Dwy,
            3000,1e-3
        );
    std::cout << "H2:\n" << H2 << "\n" << std::endl;

    const double Ts = 1./2800.;
    const auto Pz = s_to_z(P, Ts);
    const auto Mz = s_to_z(M, Ts);
    const auto H2z = s_to_z(H2, Ts);

    std::cout << "H2z:\n" << H2z << "\n" << std::endl;

    DssFilter Pf(Pz);
    DssFilter Mf(Mz);
    DssFilter H2f(H2z);

    std::ofstream file("data.csv");

    file << "time, r, err, u, z, y" << std::endl;
    double y = 0;
    for(double time = 0; time < 0.1; time+=Ts){
        const double r = 0.;
        const double err = r - y;
        const double u = H2f(err);
        const double z = Pf(u) + ((time > 0.01)? 1. : 0.);
        y = Mf(z);

        file << time << ", " << r << ", " << err << ", " << u << ", " << z << ", " << y << std::endl;
    }

}