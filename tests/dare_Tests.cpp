

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

TEST(dare, qz_1x1_swap_v2){

    const Eigen::Matrix<double, 2, 2> S({
        {1, 1},
        {0, 2}
    });
    const Eigen::Matrix<double, 2, 2> T({
        {1, 1},
        {0, 0}
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


TEST(dare, qz_2x2_swap){
    using namespace controlpp;

    using Mat4 = Eigen::Matrix4d;
    using Mat2 = Eigen::Matrix2d;

    Mat4 S = Mat4::Zero();
    Mat4 T = Mat4::Zero();

    // S = [ A  C ]
    //     [ 0  D ]
    //
    // T = [ B  F ]
    //     [ 0  E ]

    // First generalized Schur block.
    S.block<2, 2>(0, 0) <<
        1.0, 2.0,
        0.0, 3.0;

    T.block<2, 2>(0, 0) <<
        1.0, 0.0,
        0.0, 1.0;

    // Second generalized Schur block.
    S.block<2, 2>(2, 2) <<
        8.0, 1.0,
        0.0, 10.0;

    T.block<2, 2>(2, 2) <<
        2.0, 0.0,
        0.0, 2.0;

    // Coupling blocks.
    S.block<2, 2>(0, 2) <<
        0.5, 1.0,
        1.5, 0.25;

    T.block<2, 2>(0, 2) <<
        0.1, 0.2,
        0.3, 0.15;

    const Mat4 Sb = S;
    const Mat4 Tb = T;

    // Eigenvalues before swapping:
    // Leading block: 1, 3
    // Trailing block: 4, 5
    const Mat2 S_first = S.block<2, 2>(0, 0);
    const Mat2 T_first = T.block<2, 2>(0, 0);
    const Mat2 S_second = S.block<2, 2>(2, 2);
    const Mat2 T_second = T.block<2, 2>(2, 2);

    // Perform the swap.
    const auto swap = qz_2x2_swap_matrix(Sb, Tb);

    const Mat4 U = swap.U;
    const Mat4 V = swap.V;

    const Mat4 Snew = U.transpose() * Sb * V;
    const Mat4 Tnew = U.transpose() * Tb * V;

    constexpr double tol = 1e-10;

    // 1. Orthogonality.
    EXPECT_NEAR((U.transpose() * U - Mat4::Identity()).norm(), 0.0, tol);
    EXPECT_NEAR((V.transpose() * V - Mat4::Identity()).norm(), 0.0, tol);

    // 2. Equivalence transformation reconstructs the original pair.
    EXPECT_NEAR((U * Snew * V.transpose() - Sb).norm(), 0.0, tol);
    EXPECT_NEAR((U * Tnew * V.transpose() - Tb).norm(),0.0, tol);

    // 3. Generalized real Schur structure:
    //    S is quasi-upper-triangular; T is upper-triangular.
    EXPECT_NEAR(Snew(2, 0), 0.0, tol);
    EXPECT_NEAR(Snew(2, 1), 0.0, tol);
    EXPECT_NEAR(Snew(3, 0), 0.0, tol);
    EXPECT_NEAR(Snew(3, 1), 0.0, tol);

    EXPECT_NEAR(Tnew(1, 0), 0.0, tol);
    EXPECT_NEAR(Tnew(2, 0), 0.0, tol);
    EXPECT_NEAR(Tnew(2, 1), 0.0, tol);
    EXPECT_NEAR(Tnew(3, 0), 0.0, tol);
    EXPECT_NEAR(Tnew(3, 1), 0.0, tol);
    EXPECT_NEAR(Tnew(3, 2), 0.0, tol);

    // 4. The trailing block after the swap should have eigenvalues 1, 3.
    const Mat2 S_trailing = Snew.block<2, 2>(2, 2);
    const Mat2 T_trailing = Tnew.block<2, 2>(2, 2);

    EXPECT_NEAR(
        generalized_charpoly_invariants(S_trailing, T_trailing)(0),
        generalized_charpoly_invariants(S_first, T_first)(0),
        tol);

    EXPECT_NEAR(
        generalized_charpoly_invariants(S_trailing, T_trailing)(1),
        generalized_charpoly_invariants(S_first, T_first)(1),
        tol);

    // 5. The leading block after the swap should have eigenvalues 4, 5.
    const Mat2 S_leading = Snew.block<2, 2>(0, 0);
    const Mat2 T_leading = Tnew.block<2, 2>(0, 0);

    EXPECT_NEAR(
        generalized_charpoly_invariants(S_leading, T_leading)(0),
        generalized_charpoly_invariants(S_second, T_second)(0),
        tol);

    EXPECT_NEAR(
        generalized_charpoly_invariants(S_leading, T_leading)(1),
        generalized_charpoly_invariants(S_second, T_second)(1),
        tol);
}