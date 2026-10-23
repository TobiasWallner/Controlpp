

// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/dare.hpp>

#include <iostream>

TEST(dare, qz_1x1_swap){

    const Eigen::Matrix<double, 2, 2> S({
        {1, 3},
        {0, 2}
    });
    const Eigen::Matrix<double, 2, 2> T({
        {2, 4},
        {0, 3}
    });

    const double lambda1 = S(0,0) / T(0,0);
    const double lambda2 = S(1,1) / T(1,1);

    const controlpp::QZSwap qz_swap = controlpp::create_local_qz_1x1_swap(S, T);

    Eigen::Matrix<double, 2, 2> S_ = qz_swap.U.transpose() * S * qz_swap.V;
    Eigen::Matrix<double, 2, 2> T_ = qz_swap.U.transpose() * T * qz_swap.V;

    // check if U and V are orthogonal
    EXPECT_TRUE((qz_swap.U.transpose() * qz_swap.U).isApprox(Eigen::Matrix2d::Identity(), 1e-12));
    EXPECT_TRUE((qz_swap.V.transpose() * qz_swap.V).isApprox(Eigen::Matrix2d::Identity(), 1e-12));

    // check if it is upper triangular:
    EXPECT_NEAR(S_(1, 0), 0, 1e-12);
    EXPECT_NEAR(T_(1, 0), 0, 1e-12);

    // check if the generalised eigenvalues swapped
    const double lambda1_ = S_(0,0) / T_(0,0);
    const double lambda2_ = S_(1,1) / T_(1,1);
    EXPECT_NEAR(lambda1, lambda2_, 1e-12);
    EXPECT_NEAR(lambda2, lambda1_, 1e-12);
}