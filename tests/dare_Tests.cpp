

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

    const controlpp::QZSwap qz_swap = controlpp::qz_1x1_swap_matrix(S, T);

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

TEST(dare, qz_1x2_swap) {
    const Eigen::Matrix3d S({
        {1, 2, 3},
        {0, 4, 1},
        {0, 2, 5}
    });

    const Eigen::Matrix3d T({
        {2, 1, 2},
        {0, 3, 1},
        {0, 0, 4}
    });

    // Eigenvalue of the original 1x1 block.
    const double lambda = S(0, 0) / T(0, 0);

    // Original trailing 2x2 generalized eigenvalue invariants.
    const Eigen::Matrix2d S_original = S.block<2, 2>(1, 1);
    const Eigen::Matrix2d T_original = T.block<2, 2>(1, 1);

    const Eigen::Vector2d original_invariants = controlpp::generalized_charpoly_invariants(S_original, T_original);

    const controlpp::QZSwap<double, 3> qz_swap = controlpp::qz_1x2_swap_matrix(S, T);
    const Eigen::Matrix3d S_ = qz_swap.U.transpose() * S * qz_swap.V;
    const Eigen::Matrix3d T_ = qz_swap.U.transpose() * T * qz_swap.V;

    // Check orthogonality of U and V.
    EXPECT_TRUE(qz_swap.U.transpose().isApprox(qz_swap.U.inverse(), 1e-12));
    EXPECT_TRUE(qz_swap.V.transpose().isApprox(qz_swap.V.inverse(), 1e-12));

    // Check generalized Schur structure:
    // S may have a 2x2 block in the leading diagonal,
    // while T must remain upper triangular.
    EXPECT_NEAR(S_(2, 0), 0.0, 1e-12);
    EXPECT_NEAR(T_(1, 0), 0.0, 1e-12);
    EXPECT_NEAR(T_(2, 0), 0.0, 1e-12);
    EXPECT_NEAR(T_(2, 1), 0.0, 1e-12);

    // The original 1x1 eigenvalue should now be in position (2,2).
    EXPECT_NEAR(S_(2, 2) / T_(2, 2), lambda, 1e-12);

    // The leading 2x2 block should retain the original generalized
    // eigenvalues, compared through its characteristic polynomial.
    const Eigen::Matrix2d S_after = S_.block<2, 2>(0, 0);
    const Eigen::Matrix2d T_after = T_.block<2, 2>(0, 0);

    const Eigen::Vector2d after_invariants = controlpp::generalized_charpoly_invariants(S_after, T_after);
    EXPECT_TRUE(original_invariants.isApprox(after_invariants, 1e-12));
}

TEST(dare, qz_2x1_swap) {
    const Eigen::Matrix3d S({
        {4, 1, 3},
        {2, 5, 1},
        {0, 0, 1}
    });

    const Eigen::Matrix3d T({
        {3, 1, 2},
        {0, 4, 1},
        {0, 0, 2}
    });

    // Eigenvalue of the original trailing 1x1 block.
    const double lambda = S(2, 2) / T(2, 2);

    // Generalized eigenvalue invariants of the original leading 2x2 block.
    const Eigen::Matrix2d S_original = S.block<2, 2>(0, 0);
    const Eigen::Matrix2d T_original = T.block<2, 2>(0, 0);

    const Eigen::Vector2d original_invariants = controlpp::generalized_charpoly_invariants(S_original, T_original);

    const controlpp::QZSwap<double, 3> qz_swap = controlpp::qz_2x1_swap_matrix(S, T);

    const Eigen::Matrix3d S_ = qz_swap.U.transpose() * S * qz_swap.V;

    const Eigen::Matrix3d T_ = qz_swap.U.transpose() * T * qz_swap.V;

    // Check orthogonality of U and V.
    EXPECT_TRUE(qz_swap.U.transpose().isApprox(qz_swap.U.inverse(), 1e-12));

    EXPECT_TRUE(qz_swap.V.transpose().isApprox(qz_swap.V.inverse(), 1e-12));

    // Check generalized Schur structure:
    // S may have a 2x2 block in the trailing diagonal,
    // while T must remain upper triangular.
    EXPECT_NEAR(S_(2, 0), 0.0, 1e-12);
    EXPECT_NEAR(T_(1, 0), 0.0, 1e-12);
    EXPECT_NEAR(T_(2, 0), 0.0, 1e-12);
    EXPECT_NEAR(T_(2, 1), 0.0, 1e-12);

    // The original 1x1 eigenvalue should now be in position (0,0).
    EXPECT_NEAR(S_(0, 0) / T_(0, 0), lambda, 1e-12);

    // The trailing 2x2 block should retain the original generalized
    // eigenvalues, compared through its characteristic polynomial.
    const Eigen::Matrix2d S_after = S_.block<2, 2>(1, 1);
    const Eigen::Matrix2d T_after = T_.block<2, 2>(1, 1);

    const Eigen::Vector2d after_invariants = controlpp::generalized_charpoly_invariants(S_after, T_after);

    EXPECT_TRUE(original_invariants.isApprox(after_invariants, 1e-12));
}