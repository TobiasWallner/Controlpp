#include <controlpp.hpp>
#include <cstddef>
#include <fstream>
#include <iostream>

/**
 * @brief Example for an LQR controller
 * 
 * Demonstrates on how to use:
 * - `controlpp::lqr()`
 * - `controlpp::lqr_feed_forward()`
 */
int main(){
    using namespace controlpp::tf;

    const auto G = ((s - 1) * (s - 3)) / ( (s - 2) * (s - 4) * (s*s + s - 1));
    std::cout << "G:\n" << G << "\n" << std::endl;

    const auto Gss = controlpp::to_state_space(G);
    std::cout << "Gss:\n" << Gss << "\n" << std::endl;

    const double Ts = 1. / 1000.;
    const auto Gz = controlpp::discretise_zoh(Gss, Ts);
    std::cout << "Gz:\n" << Gz << "\n" << std::endl;
    auto Gzf = controlpp::DssFilter(Gz);

    const auto K = controlpp::lqr(Gz);
    std::cout << "K:\n" << K << "\n" << std::endl;

    std::ofstream csv("lqr.csv");
    if(!csv){
        std::cerr << "Failed to open lqr.csv for writing\n";
        return 1;
    }
    csv << "time, reference, control_output, plant_output\n";

    const auto F = controlpp::lqr_feed_forward(Gz, K);

    for(double t = 0; t < 20; t += Ts){
        const double r = controlpp::step(t - 0.1);
        const Eigen::Vector<double, 4> x = Gzf.states();
        const double u = (-K * x + F * r)(0);
        const double y = Gzf.input(u);
        csv << t << ", " << r << ", " << u << ", " << y << "\n";
    }

    csv.flush();
    if(!csv){
        std::cerr << "Failed while writing lqr.csv\n";
        return 1;
    }
}