#include <cmath>
#include <random>
#include <fstream>
// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/Pade.hpp>
#include <controlpp/transformations.hpp>
#include <controlpp/analysis.hpp>

TEST(Pade, pade_coefficients_1_1){
    const double NumOrder = 1;
    const double DenOrder = 1;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);

    ASSERT_NEAR(a0/a1, 2, 1e-9);

    ASSERT_NEAR(b0/a1, 2, 1e-9);
    ASSERT_NEAR(b1/a1, 1, 1e-9);
}

TEST(Pade, pade_coefficients_2_2){
    const double NumOrder = 2;
    const double DenOrder = 2;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);

    ASSERT_NEAR(a0/a2, 12, 1e-9);
    ASSERT_NEAR(a1/a2, 6, 1e-9);

    ASSERT_NEAR(b0/a2, 12, 1e-9);
    ASSERT_NEAR(b1/a2, 6, 1e-9);
    ASSERT_NEAR(b2/a2, 1, 1e-9);
}

TEST(Pade, pade_coefficients_3_3){
    const double NumOrder = 3;
    const double DenOrder = 3;
    
    const double a0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);
    const double a3 = controlpp::pade_num_param(NumOrder, DenOrder, 3);

    const double b0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);
    const double b3 = controlpp::pade_den_param(NumOrder, DenOrder, 3);

    ASSERT_NEAR(a0/a3, 120, 1e-9);
    ASSERT_NEAR(a1/a3, 60, 1e-9);
    ASSERT_NEAR(a2/a3, 12, 1e-9);

    ASSERT_NEAR(b0/a3, 120, 1e-9);
    ASSERT_NEAR(b1/a3, 60, 1e-9);
    ASSERT_NEAR(b2/a3, 12, 1e-9);
    ASSERT_NEAR(b3/a3, 1, 1e-9);
}

TEST(Pade, pade_coefficients_4_4){
    const double NumOrder = 4;
    const double DenOrder = 4;
    
    const double a0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);
    const double a3 = controlpp::pade_num_param(NumOrder, DenOrder, 3);
    const double a4 = controlpp::pade_num_param(NumOrder, DenOrder, 4);

    const double b0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);
    const double b3 = controlpp::pade_den_param(NumOrder, DenOrder, 3);
    const double b4 = controlpp::pade_den_param(NumOrder, DenOrder, 4);

    ASSERT_NEAR(a0/a4, 1680, 1e-9);
    ASSERT_NEAR(a1/a4, 840, 1e-9);
    ASSERT_NEAR(a2/a4, 180, 1e-9);
    ASSERT_NEAR(a3/a4, 20, 1e-9);

    ASSERT_NEAR(b0/a4, 1680, 1e-9);
    ASSERT_NEAR(b1/a4, 840, 1e-9);
    ASSERT_NEAR(b2/a4, 180, 1e-9);
    ASSERT_NEAR(b3/a4, 20, 1e-9);
    ASSERT_NEAR(b4/a4, 1, 1e-9);
}

// ----------------------------------------------------------------------------

TEST(Pade, pade_coefficients_0_1){
    const double NumOrder = 0;
    const double DenOrder = 1;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);

    ASSERT_NEAR(a0/a1, 1, 1e-9);

    ASSERT_NEAR(b0/a1, 1, 1e-9);
}

TEST(Pade, pade_coefficients_1_2){
    const double NumOrder = 1;
    const double DenOrder = 2;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);

    ASSERT_NEAR(a0/a2, 6, 1e-9);
    ASSERT_NEAR(a1/a2, 4, 1e-9);

    ASSERT_NEAR(b0/a2, 6, 1e-9);
    ASSERT_NEAR(b1/a2, 2, 1e-9);
}

TEST(Pade, pade_coefficients_2_3){
    const double NumOrder = 2;
    const double DenOrder = 3;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);
    const double a3 = controlpp::pade_den_param(NumOrder, DenOrder, 3);

    ASSERT_NEAR(a0/a3, 60, 1e-9);
    ASSERT_NEAR(a1/a3, 36, 1e-9);
    ASSERT_NEAR(a2/a3, 9, 1e-9);

    ASSERT_NEAR(b0/a3, 60, 1e-9);
    ASSERT_NEAR(b1/a3, 24, 1e-9);
    ASSERT_NEAR(b2/a3, 3, 1e-9);
}

TEST(Pade, pade_coefficients_3_4){
    const double NumOrder = 3;
    const double DenOrder = 4;
    
    const double b0 = controlpp::pade_num_param(NumOrder, DenOrder, 0);
    const double b1 = controlpp::pade_num_param(NumOrder, DenOrder, 1);
    const double b2 = controlpp::pade_num_param(NumOrder, DenOrder, 2);
    const double b3 = controlpp::pade_num_param(NumOrder, DenOrder, 3);

    const double a0 = controlpp::pade_den_param(NumOrder, DenOrder, 0);
    const double a1 = controlpp::pade_den_param(NumOrder, DenOrder, 1);
    const double a2 = controlpp::pade_den_param(NumOrder, DenOrder, 2);
    const double a3 = controlpp::pade_den_param(NumOrder, DenOrder, 3);
    const double a4 = controlpp::pade_den_param(NumOrder, DenOrder, 4);

    ASSERT_NEAR(b0/a4, 840, 1e-9);
    ASSERT_NEAR(b1/a4, 360, 1e-9);
    ASSERT_NEAR(b2/a4, 60, 1e-9);
    ASSERT_NEAR(b3/a4, 4, 1e-9);

    ASSERT_NEAR(a0/a4, 840, 1e-9);
    ASSERT_NEAR(a1/a4, 480, 1e-9);
    ASSERT_NEAR(a2/a4, 120, 1e-9);
    ASSERT_NEAR(a3/a4, 16, 1e-9);

}

TEST(Pade, pade_polynom_1_1){
    const double Td = 0.1;
    controlpp::ContinuousTransferFunction<double, 1, 1> ctf = controlpp::pade_delay<double, 1, 1>(Td);
    controlpp::ContinuousTransferFunction<double, 1, 1> expected({2, -Td}, {2, Td});

    // norm
    {
        const auto n = ctf.den(ctf.den().size()-1);
        ctf.num() /= n;
        ctf.den() /= n;
    }

    {
        const auto n = expected.den(expected.den().size()-1);
        expected.num() /= n;
        expected.den() /= n;
    }

    for(size_t i = 0; i < ctf.num().size(); ++i){
        ASSERT_NEAR(ctf.num(i), expected.num(i), 1e-9);
    }
    
    for(size_t i = 0; i < ctf.den().size(); ++i){
        ASSERT_NEAR(ctf.den(i), expected.den(i), 1e-9);
    }
}

TEST(Pade, pade_polynom_2_2){
    const double Td = 0.1;
    controlpp::ContinuousTransferFunction<double, 2, 2> ctf = controlpp::pade_delay<double, 2, 2>(Td);
    controlpp::ContinuousTransferFunction<double, 2, 2> expected({12, -6*Td, Td*Td}, {12, 6*Td, Td*Td});

    // norm
    {
        const auto n = ctf.den(ctf.den().size()-1);
        ctf.num() /= n;
        ctf.den() /= n;
    }

    {
        const auto n = expected.den(expected.den().size()-1);
        expected.num() /= n;
        expected.den() /= n;
    }

    for(size_t i = 0; i < ctf.num().size(); ++i){
        ASSERT_NEAR(ctf.num(i), expected.num(i), 1e-9);
    }
    
    for(size_t i = 0; i < ctf.den().size() ; ++i){
        ASSERT_NEAR(ctf.den(i), expected.den(i), 1e-9);
    }
}

TEST(Pade, pade_polynom_3_4){

    const double Td = 0.1;
    controlpp::ContinuousTransferFunction<double, 3, 4> ctf = controlpp::pade_delay<double, 3, 4>(Td);
    controlpp::ContinuousTransferFunction<double, 3, 4> expected({840, -360*Td, 60*Td*Td, -4*Td*Td*Td}, {840, 480*Td, 120*Td*Td, 16*Td*Td*Td, Td*Td*Td*Td});

    // norm
    {
        const auto n = ctf.den(ctf.den().size()-1);
        ctf.num() /= n;
        ctf.den() /= n;
    }

    {
        const auto n = expected.den(expected.den().size()-1);
        expected.num() /= n;
        expected.den() /= n;
    }

    for(size_t i = 0; i < ctf.num().size(); ++i){
        ASSERT_NEAR(ctf.num(i), expected.num(i), 1e-8);
    }
    
    for(size_t i = 0; i < ctf.den().size(); ++i){
        ASSERT_NEAR(ctf.den(i), expected.den(i), 1e-8);
    }
}

TEST(Pade, pade_3_4_step_response){
    const auto P = controlpp::pade_delay<double, 3, 4>(600e-6);
    const auto step = controlpp::step(P);
    
    const double expected[] = {0, -0.0287996, -0.0502435, -0.0652444, -0.0746436, -0.0792143, -0.0796658, -0.0766465, -0.0707481, -0.0625084, -0.052415, -0.0409077, -0.0283821, -0.0151921, -0.00165241, 0.0119583, 0.0253955, 0.038446, 0.0509262, 0.0626794, 0.0735741, 0.0835021, 0.0923764, 0.10013, 0.106712, 0.112091, 0.116249, 0.119179, 0.120891, 0.121401, 0.120739, 0.11894, 0.116048, 0.112114, 0.107195, 0.10135, 0.0946467, 0.0871515, 0.0789359, 0.0700725, 0.0606355, 0.0506999, 0.0403408, 0.0296334, 0.018652, 0.00747036, -0.00383938, -0.0152065, -0.0265623, -0.03784, -0.0489754, -0.0599065, -0.0705742, -0.0809223, -0.0908973, -0.100449, -0.10953, -0.118095, -0.126106, -0.133522, -0.14031, -0.146439, -0.151881, -0.15661, -0.160604, -0.163845, -0.166317, -0.168006, -0.168904, -0.169001, -0.168295, -0.166782, -0.164463, -0.161341, -0.157421, -0.152709, -0.147216, -0.140952, -0.13393, -0.126165, -0.117674, -0.108474, -0.0985838, -0.0880247, -0.0768182, -0.0649868, -0.0525543, -0.0395451, -0.0259847, -0.0118988, 0.00268602, 0.0177428, 0.0332442, 0.0491625, 0.06547, 0.0821387, 0.0991405, 0.116447, 0.134031, 0.151865, 0.16992, 0.18817, 0.206588, 0.225147, 0.243821, 0.262584, 0.281411, 0.300279, 0.319162, 0.338037, 0.356882, 0.375676, 0.394396, 0.413022, 0.431534, 0.449914, 0.468143, 0.486204, 0.50408, 0.521756, 0.539216, 0.556446, 0.573434, 0.590166, 0.606631, 0.622817, 0.638716, 0.654316, 0.669611, 0.684592, 0.699252, 0.713584, 0.727584, 0.741246, 0.754566, 0.767541, 0.780168, 0.792444, 0.804369, 0.81594, 0.827158, 0.838023, 0.848535, 0.858696, 0.868507, 0.87797, 0.887088, 0.895864, 0.904302, 0.912404, 0.920176, 0.927622, 0.934746, 0.941553, 0.94805, 0.954241, 0.960133, 0.96573, 0.971041, 0.976071, 0.980826, 0.985314, 0.989541, 0.993514, 0.99724, 1.00073, 1.00398, 1.00701, 1.00982, 1.01242, 1.01481, 1.01701, 1.01902, 1.02085, 1.0225, 1.02399, 1.02531, 1.02648, 1.0275, 1.02839, 1.02913, 1.02975, 1.03025, 1.03064, 1.03092, 1.03109, 1.03117, 1.03115, 1.03106, 1.03088, 1.03063, 1.03031, 1.02992, 1.02948, 1.02898, 1.02843, 1.02783, 1.0272, 1.02653, 1.02582, 1.02508, 1.02432, 1.02354, 1.02274, 1.02192, 1.02109, 1.02025, 1.0194, 1.01854, 1.01769, 1.01683, 1.01598, 1.01513, 1.01429, 1.01345, 1.01263, 1.01181, 1.01101, 1.01022, 1.00945, 1.0087, 1.00796, 1.00724, 1.00653, 1.00585, 1.00519, 1.00455, 1.00393, 1.00333, 1.00275, 1.0022, 1.00167, 1.00116, 1.00067, 1.00021, 0.999765, 0.999345, 0.998948, 0.998572, 0.998218, 0.997886, 0.997575, 0.997285, 0.997015, 0.996765, 0.996534, 0.996323, 0.99613, 0.995956, 0.995798, 0.995658, 0.995535, 0.995427, 0.995335, 0.995258, 0.995195, 0.995146, 0.995109, 0.995086, 0.995074, 0.995074, 0.995084, 0.995105, 0.995135, 0.995175, 0.995223, 0.995279, 0.995343, 0.995414, 0.995491, 0.995574, 0.995663, 0.995757, 0.995855, 0.995958, 0.996064, 0.996174, 0.996286, 0.996401, 0.996518, 0.996637, 0.996757, 0.996879, 0.997001, 0.997123, 0.997246, 0.997368, 0.99749, 0.997612, 0.997733, 0.997852, 0.99797, 0.998087, 0.998202, 0.998316, 0.998427, 0.998536, 0.998643, 0.998748, 0.99885, 0.998949, 0.999046, 0.99914, 0.999232, 0.99932, 0.999406, 0.999488, 0.999568, 0.999644, 0.999718, 0.999788, 0.999856, 0.99992, 0.999982, 1.00004, 1.0001, 1.00015, 1.0002, 1.00024, 1.00029, 1.00033, 1.00037, 1.0004, 1.00044, 1.00047, 1.0005, 1.00052, 1.00054, 1.00057, 1.00058, 1.0006, 1.00062, 1.00063, 1.00064, 1.00065, 1.00065, 1.00066, 1.00066, 1.00067, 1.00067, 1.00066, 1.00066, 1.00066, 1.00065, 1.00065, 1.00064, 1.00063, 1.00062, 1.00062, 1.0006, 1.00059, 1.00058, 1.00057, 1.00056, 1.00054, 1.00053, 1.00051, 1.0005, 1.00048, 1.00047, 1.00045, 1.00044, 1.00042, 1.00041, 1.00039, 1.00037, 1.00036, 1.00034, 1.00033, 1.00031, 1.0003, 1.00028, 1.00027, 1.00025, 1.00024, 1.00022, 1.00021, 1.00019, 1.00018, 1.00017, 1.00015, 1.00014, 1.00013, 1.00012, 1.00011, 1.0001, 1.00008, 1.00007, 1.00006, 1.00005, 1.00005, 1.00004, 1.00003, 1.00002, 1.00001, 1, 0.999998, 0.999991, 0.999985, 0.999979, 0.999974, 0.999969, 0.999964, 0.999959, 0.999955, 0.999951, 0.999948, 0.999945, 0.999942, 0.999939, 0.999937, 0.999934, 0.999933, 0.999931, 0.999929, 0.999928, 0.999927, 0.999927, 0.999926, 0.999926, 0.999926, 0.999926, 0.999926, 0.999926, 0.999927, 0.999927, 0.999928, 0.999929, 0.99993, 0.999931, 0.999932, 0.999933, 0.999934, 0.999936, 0.999937, 0.999939, 0.999941, 0.999942, 0.999944, 0.999946, 0.999947, 0.999949, 0.999951, 0.999953, 0.999955, 0.999956, 0.999958, 0.99996, 0.999962, 0.999964, 0.999965, 0.999967, 0.999969, 0.999971, 0.999972, 0.999974, 0.999976, 0.999977, 0.999979, 0.999981, 0.999982, 0.999984, 0.999985, 0.999986, 0.999988, 0.999989, 0.99999, 0.999991, 0.999993, 0.999994, 0.999995, 0.999996, 0.999997, 0.999998, 0.999999, 1, 1, 1, 1, 1, 1, 1, 1, 1.00001};

    // std::ofstream file("step.csv");
    // file << step << std::endl;

    for(size_t i = 0; i < step.size(); ++i){
        ASSERT_NEAR(step.values(i), expected[i], 1e-3);
    }

}