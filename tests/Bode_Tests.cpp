#include <fstream>
#include <iostream>

#include <controlpp.hpp>

#include <gtest/gtest.h>


/*
    Tests if one can perform loop shaping from measured data
*/
TEST(Bode, Calculation){
    const controlpp::ContinuousTransferFunction s = controlpp::tf::s<double>;

    controlpp::ContinuousTransferFunction Gs = 1 / (1 + s + s * s);
    controlpp::Bode Gb = controlpp::bode(Gs);
    {
        std::ofstream file("bode_G.csv");
        controlpp::write_csv(file, Gb);
    }
    
    {
        std::ofstream file("step_G.csv");
        controlpp::write_csv(file, controlpp::step(Gs));
    }

    controlpp::ContinuousTransferFunction Rs = 0.5 + 2 / s;
    
    controlpp::Bode L = Gb * Rs;
    {
        std::ofstream file("bode_L.csv");
        controlpp::write_csv(file, L); // test by comparison
    }

    controlpp::Bode L2 = (1 + L);

    {
        std::ofstream file("bode_L2.csv");
        controlpp::write_csv(file, L2); // test by comparison
    }

    controlpp::Bode Try = L / (1 + L);

    {
        std::ofstream file("bode_Try.csv");
        controlpp::write_csv(file,Try); // test by comparison
    }   

    controlpp::Bode Tdy = 1 / (1 + L);
    {
        std::ofstream file("bode_Tdy.csv");
        controlpp::write_csv(file, Tdy); // test by comparison
    }

    controlpp::TimeSeries ImpulseResponse = controlpp::impulse(Gb);

    {
        std::ofstream file("ImpulseResponse.csv");
        controlpp::write_csv(file, ImpulseResponse);
    }

    {
        std::ofstream file("G_1_s.csv");
        controlpp::write_csv(file, (Gb * (1/s)));
    }

    controlpp::TimeSeries StepResponse = controlpp::step(Gb);

    {
        std::ofstream file("StepResponse.csv");
        controlpp::write_csv(file, StepResponse);
    }
}

TEST(Bode, set_magnitudes_dB_and_phases_deg){

    Eigen::Vector<double, Eigen::Dynamic> mags_dB(3);
    mags_dB << 1.0, 1.0, 1.0;

    Eigen::Vector<double, Eigen::Dynamic> phases_deg(3);
    phases_deg << 0.0, 1.0, 2.0;

    controlpp::Bode bode;
    bode.set_magnitudes_dB_and_phases_deg(mags_dB, phases_deg);

}

TEST(Bode, read_csv_hz_real_image){
    using namespace controlpp;

    std::stringstream stream;
    stream << "Frequencies (Hz), Real, Imag\n";
    stream << "1, 1, 0\n";
    stream << "10, 0.5, 0.125\n";
    stream << "100, 0.25, 0.25\n";

    std::expected<Bode<double>, std::variant<EBodeCsvReadError, csvd::ReadError>> bode_result = read_bode_from_csv(stream);
    std::stringstream err_stream;
    if(bode_result.has_value() == false){
        std::visit([&](const auto& val){err_stream << val;}, bode_result.error());
    }
    ASSERT_TRUE(bode_result.has_value()) << err_stream.str();

    Bode<double>& bode = bode_result.value();

    ASSERT_EQ(bode.frequencies()[0], 1.0);
    ASSERT_EQ(bode.frequencies()[1], 10.0);
    ASSERT_EQ(bode.frequencies()[2], 100.0);

    ASSERT_EQ(bode.values()[0].real(), 1.0);
    ASSERT_EQ(bode.values()[1].real(), 0.5);
    ASSERT_EQ(bode.values()[2].real(), 0.25);

    ASSERT_EQ(bode.values()[0].imag(), 0.0);
    ASSERT_EQ(bode.values()[1].imag(), 0.125);
    ASSERT_EQ(bode.values()[2].imag(), 0.25);
}