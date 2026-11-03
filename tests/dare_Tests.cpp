

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


struct Pt1Params {
	double K;    // static gain
	double Tau;  // time constant [s]
	double Ts;   // sample time [s]
	double q;    // state weight
	double r;    // input weight
};
 
// Closed form of the scalar DARE:
//   b^2 x^2 + (r(1-a^2) - q b^2) x - q r = 0, positive root
double scalar_dare(double a, double b, double q, double r){
	const double p = r * (1.0 - a * a) - q * b * b;
	return (-p + std::sqrt(p * p + 4.0 * b * b * q * r)) / (2.0 * b * b);
}

// For the scalar case the stable eigenvector of (H, J) is known in closed form:
//   w = [x; lambda; u] = [1; X; -K],  eigenvalue mu = a - b K.
// This test needs no QZ at all, so it isolates the builders.
TEST(dare, scalar_known_eigenvector){
	const double a = std::exp(-0.2), b = 2.0 * (1.0 - a), q = 1.0, r = 1.0;
 
	Eigen::Matrix<double, 1, 1> A, B, Q, R;
	A << a; B << b; Q << q; R << r;
 
	const double X  = scalar_dare(a, b, q, r);
	const double K  = a * b * X / (r + b * b * X);
	const double mu = a - b * K;
 
	const auto H = controlpp::build_dare_symplectic_pencil_H(A, B, Q, R);
	const auto J = controlpp::build_dare_symplectic_pencil_J(A, B);
 
	Eigen::Vector3d w(1.0, X, -K);
	EXPECT_NEAR((H * w - mu * (J * w)).norm(), 0.0, 1e-12);
}

// The pencil is symplectic: the finite spectrum is closed under lambda -> 1/lambda.
TEST(dare, spectrum_is_reciprocal){
	Eigen::Matrix2d A;  A << 1.0, 0.1, -0.2, 0.8;
	Eigen::Matrix<double, 2, 1> B;  B << 0.0, 0.1;
	Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
	Eigen::Matrix<double, 1, 1> R;  R << 1.0;
 
	const auto H = controlpp::build_dare_symplectic_pencil_H(A, B, Q, R);
	const auto J = controlpp::build_dare_symplectic_pencil_J(A, B);
 
	Eigen::GeneralizedEigenSolver<Eigen::Matrix<double, 5, 5>> ges;
	ges.compute(H, J);
 
	std::vector<std::complex<double>> finite;
	for (int i = 0; i < 5; ++i) {
		if (std::abs(ges.betas()(i)) > 1e-12) {
			finite.push_back(ges.alphas()(i) / ges.betas()(i));
		}
	}
	ASSERT_EQ(finite.size(), 4u);   // 2n finite + m infinite eigenvalues
 
	for (const auto& lam : finite) {
		bool found = false;
		for (const auto& mu : finite) {
			if (std::abs(mu - 1.0 / lam) < 1e-8) { found = true; }
		}
		EXPECT_TRUE(found) << "no reciprocal partner for " << lam;
	}
}

TEST(dare, is_stable_2x2){
	Eigen::Matrix2d T;  T << 1.0, 0.2, 0.0, 1.5;
 
	Eigen::Matrix2d S;  S << 0.3, -0.8, 0.9, 0.4;       // complex pair, |lambda|^2 = 0.84/1.5
	EXPECT_TRUE (controlpp::is_stable_pencil_block(S, T));
	EXPECT_FALSE(controlpp::is_stable_pencil_block(Eigen::Matrix2d(3.0 * S), T));
 
	Eigen::Matrix2d Ti = Eigen::Matrix2d::Identity();
	Eigen::Matrix2d Sr;  Sr << 0.2, 0.0, 0.1, 0.3;      // real eigenvalues 0.2, 0.3
	EXPECT_TRUE (controlpp::is_stable_pencil_block(Sr, Ti));
	Sr << 2.0, 0.0, 0.1, 3.0;                           // real eigenvalues 2, 3
	EXPECT_FALSE(controlpp::is_stable_pencil_block(Sr, Ti));
 
	Eigen::Matrix2d Tsing;  Tsing << 1.0, 0.0, 0.0, 0.0; // det(T) = 0 -> infinite eigenvalue
	EXPECT_FALSE(controlpp::is_stable_pencil_block(S, Tsing));
}


 
class DarePt1Test : public ::testing::TestWithParam<Pt1Params> {};
 
TEST_P(DarePt1Test, MatchesClosedFormAndStabilizes){
	using T = double;
	const Pt1Params p = GetParam();
 
	// ZOH discretisation
	const T a = std::exp(-p.Ts / p.Tau);
	const T b = p.K * (1.0 - a);
 
	Eigen::Matrix<T, 1, 1> A, B, Q, R;
	A << a;
	B << b;
	Q << p.q;
	R << p.r;
 
	Eigen::Matrix<T, 1, 1> X;
	ASSERT_NO_THROW((X = controlpp::dare_qz_solver<T, 1, 1>(A, B, Q, R)));
 
	const T x = X(0, 0);
	const T x_ref = scalar_dare(a, b, p.q, p.r);
 
	// solution matches the closed form
	EXPECT_GT(x, 0.0);
	EXPECT_NEAR(x, x_ref, 1e-9 * std::abs(x_ref));
 
	// DARE residual
	const T residual = a * x * a - (a * x * b) * (a * x * b) / (p.r + b * x * b) + p.q - x;
	EXPECT_NEAR(residual, 0.0, 1e-9 * std::max(T(1), std::abs(x)));
 
	// closed loop is stable: u = -k x
	const T k = (b * x * a) / (p.r + b * x * b);
	EXPECT_LT(std::abs(a - b * k), 1.0);
}
 
INSTANTIATE_TEST_SUITE_P(
	Pt1,
	DarePt1Test,
	::testing::Values(
		Pt1Params{2.0, 0.5,  0.1,  1.0,  1.0},   // baseline
		Pt1Params{1.0, 1.0,  0.1,  1.0,  1.0},   // unit gain
		Pt1Params{1.0, 1.0,  0.001, 1.0,  1.0},   // fast sampling, a close to 1
		Pt1Params{1.0, 0.1,  1.0,  1.0,  1.0},   // slow sampling, a close to 0
		Pt1Params{5.0, 2.0,  0.1,  10.0, 0.1},   // cheap control
		Pt1Params{0.5, 0.3,  0.05, 0.1,  10.0},  // expensive control
		Pt1Params{1.0, 1.0,  0.1,  1.0,  1e-3}   // very cheap control
	)
);
 

// ---- generated by dare_pt2_peak_reference.py ----

/*
#!/usr/bin/env python3
"""
Reference values for the C++ test of dare_qz_solver().

Plant (continuous time):   G(s) = PT2(s) * Peak(s)

    PT2(s)  = K * w0^2 / (s^2 + 2*zeta0*w0*s + w0^2)
    Peak(s) = (s^2 + 2*zeta_num*wn*s + wn^2) / (s^2 + 2*zeta_den*wn*s + wn^2)
              (inverse of a notch: numerator damping > denominator damping,
               gain at wn is zeta_num / zeta_den)

The plant is discretised with ZOH, then the DARE

    A' X A - A' X B (R + B' X B)^-1 B' X A + Q = X

is solved with scipy and cross-checked with an independent Newton (Kleinman)
iteration.  The script prints the matrices as C++ code for the GoogleTest.

Usage:  python3 dare_pt2_peak_reference.py
"""

import numpy as np
import scipy.linalg as sl
import scipy.signal as sig

# ---------------------------------------------------------------- parameters
K = 1.0
W0 = 10.0         # PT2 corner frequency [rad/s]
ZETA0 = 0.5       # PT2 damping
WN = 30.0         # peak frequency [rad/s]
ZETA_NUM = 0.5    # peak: numerator damping
ZETA_DEN = 0.05   # peak: denominator damping (lightly damped -> sharp peak)
TS = 0.01         # sample time [s]

Q_OUTPUT_WEIGHT = 1.0    # Q = w * C'C + q_reg * I
Q_REG = 0.01
R_WEIGHT = 1.0


def build_plant():
    pt2_num = [K * W0**2]
    pt2_den = [1.0, 2 * ZETA0 * W0, W0**2]
    peak_num = [1.0, 2 * ZETA_NUM * WN, WN**2]
    peak_den = [1.0, 2 * ZETA_DEN * WN, WN**2]

    num = np.polymul(pt2_num, peak_num)
    den = np.polymul(pt2_den, peak_den)

    Ac, Bc, Cc, Dc = sig.tf2ss(num, den)

    # The controllable canonical form has entries spanning ~10 orders of
    # magnitude (w^3 ... 1), which makes the DARE needlessly ill-conditioned.
    # A diagonal similarity transform (balancing) does not change the
    # transfer function but gives a well scaled problem.
    Ac, Tbal = sl.matrix_balance(Ac)
    Bc = np.linalg.solve(Tbal, Bc)
    Cc = Cc @ Tbal

    Ad, Bd, Cd, Dd, _ = sig.cont2discrete((Ac, Bc, Cc, Dc), TS, method="zoh")
    return Ad, Bd, Cd, Dd, (num, den)


def dare_residual(A, B, Q, R, X):
    G = R + B.T @ X @ B
    K_ = np.linalg.solve(G, B.T @ X @ A)
    return A.T @ X @ A - A.T @ X @ B @ K_ + Q - X


def newton_kleinman(A, B, Q, R, K0, iters=60):
    """Independent solver: Newton iteration, needs a stabilising start gain."""
    Kc = K0
    X = None
    for _ in range(iters):
        Acl = A - B @ Kc
        X = sl.solve_discrete_lyapunov(Acl.T, Q + Kc.T @ R @ Kc)
        Kc = np.linalg.solve(R + B.T @ X @ B, B.T @ X @ A)
    return X


def cpp_matrix(name, M, rows, cols):
    lines = [f"\tEigen::Matrix<double, {rows}, {cols}> {name};"]
    lines.append(f"\t{name} <<")
    for i in range(rows):
        row = ", ".join(repr(float(M[i, j])) for j in range(cols))
        lines.append("\t\t" + row + (",") if i < rows - 1 else "\t\t" + row + ";")
    return "\n".join(lines)


def main():
    A, B, C, D, (num, den) = build_plant()
    n, m = B.shape

    Q = Q_OUTPUT_WEIGHT * (C.T @ C) + Q_REG * np.eye(n)
    R = R_WEIGHT * np.eye(m)

    # ---- reference solution (scipy: QZ based) -------------------------------
    X = sl.solve_discrete_are(A, B, Q, R)
    X = 0.5 * (X + X.T)

    # ---- independent cross check (Newton, starts from K = 0: plant is stable)
    assert np.max(np.abs(np.linalg.eigvals(A))) < 1.0, "plant must be stable for K0 = 0"
    X_newton = newton_kleinman(A, B, Q, R, np.zeros((m, n)))

    Kgain = np.linalg.solve(R + B.T @ X @ B, B.T @ X @ A)
    cl_poles = np.linalg.eigvals(A - B @ Kgain)

    print("# plant poles (discrete):", np.round(np.linalg.eigvals(A), 6))
    print("# closed loop poles    :", np.round(cl_poles, 6))
    print("# |closed loop poles|  :", np.round(np.abs(cl_poles), 6))
    print("# residual (scipy)     :", np.linalg.norm(dare_residual(A, B, Q, R, X)))
    print("# residual (newton)    :", np.linalg.norm(dare_residual(A, B, Q, R, X_newton)))
    print("# |X_scipy - X_newton| / |X| :", np.linalg.norm(X - X_newton) / np.linalg.norm(X))
    print("# |X| :", np.linalg.norm(X), " cond(X):", np.linalg.cond(X))

    # eigenvalue structure of the symplectic pencil used by the C++ solver
    N = 2 * n + m
    H = np.zeros((N, N)); J = np.zeros((N, N))
    H[:n, :n] = A; H[:n, 2 * n:] = B
    H[n:2 * n, :n] = -Q; H[n:2 * n, n:2 * n] = np.eye(n)
    H[2 * n:, 2 * n:] = R
    J[:n, :n] = np.eye(n); J[n:2 * n, n:2 * n] = A.T; J[2 * n:, n:2 * n] = -B.T
    ev = sl.eigvals(H, J)
    finite = ev[np.isfinite(ev)]
    n_complex = int(np.sum(np.abs(finite.imag) > 1e-9))
    print("# pencil size:", N, " finite eigenvalues:", len(finite),
          " complex ones:", n_complex, " stable:", int(np.sum(np.abs(finite) < 1)))
    print()

    print("// ---- generated by dare_pt2_peak_reference.py ----")
    print(f"// PT2 (K={K}, w0={W0}, zeta={ZETA0}) * peak (wn={WN}, zeta_num={ZETA_NUM}, zeta_den={ZETA_DEN}), Ts={TS}")
    print(cpp_matrix("A", A, n, n))
    print(cpp_matrix("B", B, n, m))
    print(cpp_matrix("Q", Q, n, n))
    print(cpp_matrix("R", R, m, m))
    print(cpp_matrix("X_ref", X, n, n))
    print(cpp_matrix("K_ref", Kgain, m, n))


if __name__ == "__main__":
    main()
*/

struct Pt2PeakReference {
	Eigen::Matrix<double, 4, 4> A;
	Eigen::Matrix<double, 4, 1> B;
	Eigen::Matrix<double, 4, 4> Q;
	Eigen::Matrix<double, 1, 1> R;
	Eigen::Matrix<double, 4, 4> X_ref;
	Eigen::Matrix<double, 1, 4> K_ref;
 
	Pt2PeakReference(){
		A <<
			0.8297827133740366, -0.31081056892029296, -0.17570112054531462, -0.10122191762658658,
			0.29483020415893224, 0.9495574838136029, -0.028505688067733785, -0.01669510013871633,
			0.02431400183757585, 0.15729266532598135, 0.9984704171977887, -0.0009031225050561331,
			0.0013152674775857493, 0.012691328331557136, 0.1599386135719058, 0.9999635797394815;
		B <<
			0.0001439600606244787,
			2.3744142419507663e-05,
			1.2844408960798335e-06,
			5.179770384855041e-08;
		Q <<
			0.01, 0.0, 0.0, 0.0,
			0.0, 40000.01, 75000.0, 140625.0,
			0.0, 75000.0, 140625.01, 263671.875,
			0.0, 140625.0, 263671.875, 494384.775625;
		R <<
			1.0;
		X_ref <<
			2204418.0238277377, 1079839.2464192873, 3165232.5892979205, 1860975.0061927277,
			1079839.2464192873, 1272517.8264168578, 2171922.676930911, 1894647.1449130983,
			3165232.5892979205, 2171922.676930911, 5886815.830944924, 4046256.7640126324,
			1860975.0061927277, 1894647.1449130983, 4046256.7640126324, 3790948.7135033035;
		K_ref <<
			337.95970853191204, 148.03263222304633, 472.6971557154504, 264.91741107849253;
	}
};
// ---- end generated ----
 
double spectral_radius(const Eigen::Matrix4d& M){
	Eigen::EigenSolver<Eigen::Matrix4d> es(M);
	return es.eigenvalues().cwiseAbs().maxCoeff();
}
 
TEST(DarePt2Peak, MatchesPythonReference){
	const Pt2PeakReference ref;
 
	Eigen::Matrix4d X;
	ASSERT_NO_THROW((X = controlpp::dare_qz_solver<double, 4, 1>(ref.A, ref.B, ref.Q, ref.R)));
 
	// 1. same solution as scipy
	const double rel_err = (X - ref.X_ref).norm() / ref.X_ref.norm();
	EXPECT_LT(rel_err, 1e-7) << "relative error of X against scipy: " << rel_err;
 
	// 2. X is symmetric and positive semidefinite
	EXPECT_LT((X - X.transpose()).norm() / X.norm(), 1e-12);
	EXPECT_GT(X.selfadjointView<Eigen::Lower>().eigenvalues().minCoeff(), 0.0);
 
	// 3. residual of the DARE itself (independent of the reference)
	const Eigen::Matrix<double, 1, 1> G = ref.R + ref.B.transpose() * X * ref.B;
	const Eigen::Matrix<double, 1, 4> K = G.ldlt().solve(ref.B.transpose() * X * ref.A);
	const Eigen::Matrix4d residual = ref.A.transpose() * X * ref.A - ref.A.transpose() * X * ref.B * K + ref.Q - X;
	EXPECT_LT(residual.norm() / X.norm(), 1e-9) << "DARE residual";
 
	// 4. same feedback gain as scipy
	EXPECT_LT((K - ref.K_ref).norm() / ref.K_ref.norm(), 1e-7) << "gain K";
 
	// 5. same closed loop, and it is stable
	const double rho     = spectral_radius(ref.A - ref.B * K);
	const double rho_ref = spectral_radius(ref.A - ref.B * ref.K_ref);
	EXPECT_LT(rho, 1.0);
	EXPECT_NEAR(rho, rho_ref, 1e-7);
}