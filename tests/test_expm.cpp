#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <iomanip>
#include <iostream>
#include <cmath>

#include "controlpp/expm.hpp"

TEST(PadeTest, Coefficients)
{
    using namespace controlpp;

    EXPECT_NEAR(
        pade_num_param(3, 3, 0),
        1.0,
        1e-15
    );

    EXPECT_NEAR(
        pade_num_param(3, 3, 1),
        0.5,
        1e-15
    );

    EXPECT_NEAR(
        pade_num_param(3, 3, 2),
        0.1,
        1e-15
    );

    EXPECT_NEAR(
        pade_num_param(3, 3, 3),
        1.0 / 120.0,
        1e-15
    );
}

TEST(PadeTest, HighOrderCoefficients)
{
    using namespace controlpp;

    EXPECT_NEAR(pade_num_param(5, 5, 0), 1.0, 1e-15);
    EXPECT_NEAR(pade_num_param(5, 5, 1), 0.5, 1e-15);
    EXPECT_NEAR(pade_num_param(5, 5, 2), 1.0 / 9.0, 1e-15);

    EXPECT_NEAR(pade_num_param(9, 9, 0), 1.0, 1e-15);
    EXPECT_NEAR(pade_num_param(9, 9, 1), 0.5, 1e-15);

    EXPECT_NEAR(pade_num_param(13, 13, 0), 1.0, 1e-15);
    EXPECT_NEAR(pade_num_param(13, 13, 1), 0.5, 1e-15);
    EXPECT_NEAR(pade_num_param(13, 13, 2), 0.12, 1e-15);
}

template<typename DerivedA, typename DerivedB>
void expect_matrix_near(
    const Eigen::MatrixBase<DerivedA>& actual,
    const Eigen::MatrixBase<DerivedB>& expected,
    double tolerance)
{
    using namespace controlpp;
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());

    const double error =
        (actual - expected).template cast<double>().norm();

    const double scale =
        std::max(1.0, expected.template cast<double>().norm());

    EXPECT_LE(error, tolerance * scale);
}

TEST(ExpmTest, Identity){
    using namespace controlpp;
    Eigen::Matrix3d A = Eigen::Matrix3d::Identity();

    const Eigen::Matrix3d expected =
        std::exp(1.0) * Eigen::Matrix3d::Identity();

    const Eigen::Matrix3d actual = expm(A);

    std::cout << "actual:\n" << std::setprecision(17)
              << actual << "\n";

    std::cout << "expected:\n"
              << expected << "\n";

    expect_matrix_near(actual, expected, 1e-14);
}

TEST(ExpmTest, Zero){
    using namespace controlpp;
    Eigen::Matrix4d A = Eigen::Matrix4d::Zero();

    const Eigen::Matrix4d expected =
        Eigen::Matrix4d::Identity();

    const Eigen::Matrix4d actual = expm(A);

    expect_matrix_near(actual, expected, 1e-15);
}

TEST(ExpmTest, Diagonal){
    using namespace controlpp;
    Eigen::Matrix3d A;
    A << 1.0, 0.0, 0.0,
         0.0, 2.0, 0.0,
         0.0, 0.0, -1.0;

    Eigen::Matrix3d expected;
    expected << std::exp(1.0), 0.0, 0.0,
                0.0, std::exp(2.0), 0.0,
                0.0, 0.0, std::exp(-1.0);

    const Eigen::Matrix3d actual = expm(A);

    expect_matrix_near(actual, expected, 1e-14);
}

TEST(ExpmTest, Rotation){
    using namespace controlpp;
    const double theta = 0.7;

    Eigen::Matrix2d A;
    A << 0.0,  theta,
        -theta, 0.0;

    Eigen::Matrix2d expected;
    expected << std::cos(theta),  std::sin(theta),
               -std::sin(theta),  std::cos(theta);

    const Eigen::Matrix2d actual = expm(A);

    expect_matrix_near(actual, expected, 1e-14);
}

TEST(ExpmTest, NilpotentJordanBlock){
    using namespace controlpp;
    Eigen::Matrix3d A;
    A << 0.0, 1.0, 0.0,
         0.0, 0.0, 1.0,
         0.0, 0.0, 0.0;

    Eigen::Matrix3d expected;
    expected << 1.0, 1.0, 0.5,
                0.0, 1.0, 1.0,
                0.0, 0.0, 1.0;

    const Eigen::Matrix3d actual = expm(A);

    expect_matrix_near(actual, expected, 1e-14);
}

TEST(ExpmTest, DefectiveMatrix){
    using namespace controlpp;
    Eigen::Matrix2d A;
    A << -1.0, 1.0,
          0.0, -1.0;

    const double e = std::exp(-1.0);

    Eigen::Matrix2d expected;
    expected << e, e,
                0.0, e;

    const Eigen::Matrix2d actual = expm(A);

    expect_matrix_near(actual, expected, 1e-14);
}

TEST(ExpmTest, MatlabExample){
    using namespace controlpp;
    Eigen::Matrix3d A;
    A << 1.0, 1.0,  0.0,
         0.0, 0.0,  2.0,
         0.0, 0.0, -1.0;

    Eigen::Matrix3d expected;
    expected << 2.718281828459045,
                1.718281828459045,
                1.086161269630487,
                0.0,
                1.0,
                1.264241117657115,
                0.0,
                0.0,
                0.367879441171442;

    const Eigen::Matrix3d actual = expm(A);

    expect_matrix_near(actual, expected, 1e-14);
}

TEST(ExpmTest, Ward1977TestCase1){
    using namespace controlpp;
    Eigen::Matrix3d A;
    A << 4.0, 2.0, 0.0,
         1.0, 4.0, 1.0,
         1.0, 1.0, 4.0;

    Eigen::Matrix3d expected;
    expected <<
        147.86662244637000, 183.76513864636857,  71.79703239999643,
        127.78108552318250, 183.76513864636877,  91.88256932318409,
        127.78108552318204, 163.67960172318047, 111.96810624637124;

    const Eigen::Matrix3d actual = expm(A);

    expect_matrix_near(actual, expected, 1e-13);
}

TEST(ExpmTest, Scaling){
    using namespace controlpp;
    Eigen::Matrix2d A = 20.0 * Eigen::Matrix2d::Identity();

    const Eigen::Matrix2d expected =
        std::exp(20.0) * Eigen::Matrix2d::Identity();

    const Eigen::Matrix2d actual = expm(A);

    expect_matrix_near(actual, expected, 1e-13);
}

TEST(ExpmPadeParamsTest, Scaling){
    using namespace controlpp;
    Eigen::Matrix2d A = 20.0 * Eigen::Matrix2d::Identity();

    const ExpmPadeParams params = expm_pade_params(A);

    EXPECT_EQ(params.order, 13);
    EXPECT_EQ(params.scaling, 2);
}

TEST(ExpmPadeParamsTest, SelectsCorrectOrders){
    using namespace controlpp;
    constexpr double theta3  = 1.495585217958292e-2;
    constexpr double theta5  = 2.539398330063230e-1;
    constexpr double theta7  = 9.504178996162932e-1;
    constexpr double theta9  = 2.097847961257068;
    constexpr double theta13 = 5.371920351148152;

    {
        Eigen::Matrix2d A =
            (0.9 * theta3) * Eigen::Matrix2d::Identity();

        const auto p = expm_pade_params(A);

        EXPECT_EQ(p.order, 3);
        EXPECT_EQ(p.scaling, 0);
    }

    {
        Eigen::Matrix2d A =
            (0.9 * theta5) * Eigen::Matrix2d::Identity();

        const auto p = expm_pade_params(A);

        EXPECT_EQ(p.order, 5);
        EXPECT_EQ(p.scaling, 0);
    }

    {
        Eigen::Matrix2d A =
            (0.9 * theta7) * Eigen::Matrix2d::Identity();

        const auto p = expm_pade_params(A);

        EXPECT_EQ(p.order, 7);
        EXPECT_EQ(p.scaling, 0);
    }

    {
        Eigen::Matrix2d A =
            (0.9 * theta9) * Eigen::Matrix2d::Identity();

        const auto p = expm_pade_params(A);

        EXPECT_EQ(p.order, 9);
        EXPECT_EQ(p.scaling, 0);
    }

    {
        Eigen::Matrix2d A =
            (0.9 * theta13) * Eigen::Matrix2d::Identity();

        const auto p = expm_pade_params(A);

        EXPECT_EQ(p.order, 13);
        EXPECT_EQ(p.scaling, 0);
    }
}

TEST(ExpmTest, LargeCancellation){
    using namespace controlpp;
    Eigen::Matrix2d A;
    A << -147.0, 72.0,
         -192.0, 93.0;

    Eigen::Matrix2d expected;
    expected <<
        -0.09957413673572789,  0.07468060255179592,
        -0.19914827347145578,  0.14936120510359184;

    const Eigen::Matrix2d actual = expm(A);

    expect_matrix_near(actual, expected, 1e-13);
}

TEST(ExpmPadeTest, ScaledVsUnscaled){
    using namespace controlpp;
    Eigen::Matrix3d A;
    A << 1.0,  2.0,  0.0,
         0.0, -1.0,  3.0,
         2.0,  0.0,  1.0;

    const Eigen::Matrix3d reference =
        expm_pade_scaled(A, 13, 4);

    const Eigen::Matrix3d actual =
        expm(A);

    expect_matrix_near(actual, reference, 1e-13);
}