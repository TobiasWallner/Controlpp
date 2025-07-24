// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/math.hpp>
#include <controlpp/Polynom.hpp>

TEST(Polynom, default_construction){
    controlpp::Polynom<float, 3>();
    controlpp::Polynom<double, 10>();
    controlpp::Polynom<int, 1>();
}

TEST(Polynom, array_construction){
    controlpp::Polynom<int, 3> p({1,2,3});

    ASSERT_EQ(p[0], 1);
    ASSERT_EQ(p[1], 2);
    ASSERT_EQ(p[2], 3);

    ASSERT_EQ(p.order(), 2);
}

#include <iostream>
TEST(Polynom, vector_construction){
    Eigen::Vector<float, 3> vector({1.0f,2.0f,3.0f});
    controlpp::Polynom<float, 3> p(vector);

    ASSERT_EQ(p.at(0), 1.0f);
    ASSERT_EQ(p.at(1), 2.0f);
    ASSERT_EQ(p.at(2), 3.0f);
}


TEST(Polynom, addition){
    const int arr_a[] = {1, 2, 3};
    const int arr_b[] = {4, 5, 6};
    controlpp::Polynom<int, 3> a(arr_a);
    controlpp::Polynom<int, 3> b(arr_b);

    const auto result = a + b;

    for(int i = 0; i <= 2; ++i){
        ASSERT_EQ(result[i], (arr_a[i] + arr_b[i]));
    }
}

TEST(Polynom, scalar_multiplication){
    controlpp::Polynom<int, 3> a({1, 2, 3});
    controlpp::Polynom<int, 3> expected({2, 4, 6});

    const auto result1 = a * 2;
    const auto result2 = 2 * a;
    const auto result3 = a * 3;

    ASSERT_EQ(result1, expected);
    ASSERT_EQ(result2, expected);
    ASSERT_NE(result3, expected);
}

TEST(Polynom, multiplication){
    const controlpp::Polynom<int, 3> a({1, 2, 3});
    const controlpp::Polynom<int, 2> b({1, 2});

    // 1, 2, 3
    //    2, 4, 6
    // --------------
    // 1, 4, 7, 6
    const controlpp::Polynom<int, 4> expected({1, 4, 7, 6});

    const auto result = a * b;

    ASSERT_EQ(result.order(), expected.order());
    ASSERT_EQ(result, expected);
}

TEST(Polynom, zeros_order_1){
    const controlpp::Polynom<double, 2> three({3, -1});
    
    const Eigen::Vector<std::complex<double>, 1> zeros = controlpp::zeros(three);

    ASSERT_NEAR(std::real(zeros(0)), 3, 1e-6);
    ASSERT_NEAR(std::imag(zeros(0)), 0, 1e-6);
}

TEST(Polynom, zeros_order_2){
    const controlpp::Polynom<double, 2> one({1, -1});
    const controlpp::Polynom<double, 2> three({3, -1});

    const controlpp::Polynom<double, 3> p = one * three;

    const Eigen::Vector<std::complex<double>, 2> zeros = controlpp::zeros(p);

    // sort and seperate
    Eigen::Vector<double, 2> real_zeros = controlpp::real(zeros);
    Eigen::Vector<double, 2> imag_zeros = controlpp::imag(zeros);

    std::sort(real_zeros.data(), real_zeros.data() + real_zeros.size());
    
    // check the real parts
    ASSERT_NEAR(real_zeros(0), 1, 1e-6);
    ASSERT_NEAR(real_zeros(1), 3, 1e-6);

    ASSERT_NEAR(imag_zeros(0), 0, 1e-6);
    ASSERT_NEAR(imag_zeros(1), 0, 1e-6);
}

TEST(Polynom, zeros_order_n){
    // create polynomials with zeros at: 1, 2, 3 and 4
    const controlpp::Polynom<double, 2> one({1, -1});
    const controlpp::Polynom<double, 2> two({2, -1});
    const controlpp::Polynom<double, 2> three({3, -1});
    const controlpp::Polynom<double, 2> four({4, -1});

    // multiply them to get a polynomial that has all those zeros
    const controlpp::Polynom<double, 5> p = one * two * three * four;
    
    // function under test
    const auto zeros = controlpp::zeros(p);

    // sort and seperate
    Eigen::Vector<double, 4> real_zeros = controlpp::real(zeros);
    Eigen::Vector<double, 4> imag_zeros = controlpp::imag(zeros);

    std::sort(real_zeros.data(), real_zeros.data() + real_zeros.size());
    
    // check the real parts
    ASSERT_NEAR(real_zeros(0), 1, 1e-6);
    ASSERT_NEAR(real_zeros(1), 2, 1e-6);
    ASSERT_NEAR(real_zeros(2), 3, 1e-6);
    ASSERT_NEAR(real_zeros(3), 4, 1e-6);

    ASSERT_NEAR(imag_zeros(0), 0, 1e-6);
    ASSERT_NEAR(imag_zeros(1), 0, 1e-6);
    ASSERT_NEAR(imag_zeros(2), 0, 1e-6);
    ASSERT_NEAR(imag_zeros(3), 0, 1e-6);
}