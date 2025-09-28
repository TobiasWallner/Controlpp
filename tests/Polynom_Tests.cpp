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

TEST(Polynom, vector_construction){
    Eigen::Vector<float, 3> vector({1.0f,2.0f,3.0f});
    controlpp::Polynom<float, 2> p(vector);

    ASSERT_EQ(p.at(0), 1.0f);
    ASSERT_EQ(p.at(1), 2.0f);
    ASSERT_EQ(p.at(2), 3.0f);
}


TEST(Polynom, addition){
    const int arr_a[] = {1, 2, 3};
    const int arr_b[] = {4, 5, 6};
    controlpp::Polynom<int, 2> a(arr_a);
    controlpp::Polynom<int, 2> b(arr_b);

    const auto result = a + b;

    for(int i = 0; i <= 2; ++i){
        ASSERT_EQ(result[i], (arr_a[i] + arr_b[i]));
    }
}

TEST(Polynom, scalar_multiplication){
    controlpp::Polynom<int, 2> a(1, 2, 3);
    controlpp::Polynom<int, 2> expected(2, 4, 6);

    const auto result1 = a * 2;
    const auto result2 = 2 * a;
    const auto result3 = a * 3;

    ASSERT_EQ(result1, expected);
    ASSERT_EQ(result2, expected);
    ASSERT_NE(result3, expected);
}

TEST(Polynom, multiplication){
    const controlpp::Polynom<int, 2> a(1, 2, 3);
    const controlpp::Polynom<int, 1> b(1, 2);

    // 1, 2, 3
    //    2, 4, 6
    // --------------
    // 1, 4, 7, 6
    const controlpp::Polynom<int, 3> expected(1, 4, 7, 6);

    const auto result = a * b;

    ASSERT_EQ(result.order(), expected.order());
    ASSERT_EQ(result, expected);
}

TEST(Polynom, zeros_order_1){
    const controlpp::Polynom<double, 1> three(3, -1);
    
    const Eigen::Vector<std::complex<double>, 1> zeros = controlpp::zeros(three);

    ASSERT_NEAR(std::real(zeros(0)), 3, 1e-6);
    ASSERT_NEAR(std::imag(zeros(0)), 0, 1e-6);
}

TEST(Polynom, zeros_order_2){
    const controlpp::Polynom<double, 1> one(1, -1);
    const controlpp::Polynom<double, 1> three(3, -1);

    const controlpp::Polynom<double, 2> p = one * three;

    const Eigen::Vector<std::complex<double>, 2> zeros = controlpp::zeros(p);

    // sort and seperate
    Eigen::Vector<double, 2> real_zeros = zeros.real();
    Eigen::Vector<double, 2> imag_zeros = zeros.imag();

    std::sort(real_zeros.data(), real_zeros.data() + real_zeros.size());
    
    // check the real parts
    ASSERT_NEAR(real_zeros(0), 1, 1e-6);
    ASSERT_NEAR(real_zeros(1), 3, 1e-6);

    ASSERT_NEAR(imag_zeros(0), 0, 1e-6);
    ASSERT_NEAR(imag_zeros(1), 0, 1e-6);
}

TEST(Polynom, zeros_order_n){
    // create polynomials with zeros at: 1, 2, 3 and 4
    const controlpp::Polynom<double, 1> one(1, -1);
    const controlpp::Polynom<double, 1> two(2, -1);
    const controlpp::Polynom<double, 1> three(3, -1);
    const controlpp::Polynom<double, 1> four(4, -1);

    // multiply them to get a polynomial that has all those zeros
    const controlpp::Polynom<double, 4> p = one * two * three * four;
    
    // function under test
    const auto zeros = controlpp::zeros(p);

    // sort and seperate
    Eigen::Vector<double, 4> real_zeros = zeros.real();
    Eigen::Vector<double, 4> imag_zeros = zeros.imag();

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

TEST(Polynom, scalar_evaluation){
    const auto x = controlpp::polynom::x<double>;

    const auto p = (1 - 5 * x + 2 * x * x);
    auto f = [](const double x){return 1 - 5 * x + 2 * x * x;};

    for(double x=-152.548; x < 152.548; x += 5.37){
        const double xf = f(x);
        const double xp = p(x);
        ASSERT_NEAR(xf, xp, 1e-9);
    }
}

TEST(Polynom, vector_evaluation){
    const auto x = controlpp::polynom::x<double>;

    const auto p = (1 - 5 * x + 2 * x * x);
    auto f = [](const double x){return 1 - 5 * x + 2 * x * x;};

    const Eigen::Vector<double, 1000> X = Eigen::Vector<double, 1000>::LinSpaced(-15.548, 15.548);
    const Eigen::Vector<double, 1000> Xp = p(X);

    for(int i = 0; i < 1000; ++i){
        const double xf = f(X(i));
        ASSERT_NEAR(xf, Xp(i), 1e-9);
    }
}

TEST(Polynom, complex_scalar_eval){
    using namespace std::complex_literals;
    const auto x = controlpp::polynom::x<double>;

    const auto p = (1 - 5 * x + 2 * x * x);
    auto f = [](const std::complex<double> x){return 1. - 5. * x + 2. * x * x;};

    Eigen::Vector<std::complex<double>, 4> X;
    X(0) = -3.0 - 1.0i;
    X(1) = +3.0 - 1.0i;
    X(2) = +1.0 - 3.0i;
    X(3) = +1.0 + 3.0i;

    for(int i = 0; i < 4; ++i){
        const auto xp = p(X(i));
        const auto xf = f(X(i));
        ASSERT_NEAR(xp.real(), xf.real(), 1e-9);
        ASSERT_NEAR(xp.imag(), xf.imag(), 1e-9);
    }
}

TEST(Polynom, complex_vector_eval){
    using namespace std::complex_literals;
    const auto x = controlpp::polynom::x<double>;

    const auto p = (1 - 5 * x + 2 * x * x);
    auto f = [](const std::complex<double> x){return 1. - 5. * x + 2. * x * x;};

    Eigen::Vector<std::complex<double>, 4> X;
    const auto Xp = p(X);

    for(int i = 0; i < 4; ++i){
        const auto xp = Xp(i);
        const auto xf = f(X(i));
        ASSERT_NEAR(xp.real(), xf.real(), 1e-9);
        ASSERT_NEAR(xp.imag(), xf.imag(), 1e-9);
    }
}

TEST(Polynom, static_construction_from_variadic_arguments){
    const controlpp::Polynom<int, 3> p(1, 2, 3, 4);
    ASSERT_EQ(p.size(), 4);
    ASSERT_EQ(p.order(), 3);
    ASSERT_EQ(p[0], 1);
    ASSERT_EQ(p[1], 2);
    ASSERT_EQ(p[2], 3);
}

TEST(Polynom, dynamic_construction_from_array){
    const int values[] = {1, 2, 3};

    const controlpp::Polynom<int, Eigen::Dynamic> p(values);

    ASSERT_EQ(p.size(), 3);
    ASSERT_EQ(p.order(), 2);
    ASSERT_EQ(p[0], 1);
    ASSERT_EQ(p[1], 2);
    ASSERT_EQ(p[2], 3);
}

TEST(Polynom, dynamic_construction_from_vector){
    const Eigen::Vector<int, 3> values(1, 2, 3);

    const controlpp::Polynom<int, Eigen::Dynamic> p(values);

    ASSERT_EQ(p.size(), 3);
    ASSERT_EQ(p.order(), 2);
    ASSERT_EQ(p[0], 1);
    ASSERT_EQ(p[1], 2);
    ASSERT_EQ(p[2], 3);
}

TEST(Polynom, dynamic_construction_from_variadic){
    const controlpp::Polynom<int, Eigen::Dynamic> p(1, 2, 3);

    ASSERT_EQ(p.size(), 3);
    ASSERT_EQ(p.order(), 2);
    ASSERT_EQ(p[0], 1);
    ASSERT_EQ(p[1], 2);
    ASSERT_EQ(p[2], 3);
}