#pragma once

namespace controlpp{


/**
 * @brief Computes a diagonal similarity scaling that "balances" a square matrix A.
 *
 * We look for a diagonal matrix D = diag(d) such that
 *
 *     A' = D^-1 * A * D        <=>   A'(i,j) = A(i,j) * d(j) / d(i)
 *
 * has rows and columns of comparable size. A' has exactly the same eigenvalues
 * as A (it is a similarity transform), but badly scaled entries, for example
 * 1e3 next to 1e-3, are brought closer together. That makes QZ / eigenvalue
 * algorithms noticeably more accurate.
 *
 * Details:
 *  - Only the OFF-diagonal entries are balanced. A diagonal similarity never
 *    changes the diagonal (A'(i,i) = A(i,i) * d(i) / d(i)), so it is ignored.
 *  - Every d(i) is a power of two. Multiplying a floating point number by a
 *    power of two is exact, so the balancing itself adds no rounding error.
 *  - The size N must be known at compile time (fixed-size Eigen matrix).
 *
 * How to use the result (state coordinates x = D x'):
 *
 *     A' = D^-1 A D,   B' = D^-1 B,   Q' = D Q D,   X = D^-1 X' D^-1
 *
 * This is a simplified Parlett-Reinsch balancing (the idea behind LAPACK's
 * dgebal), scaling only, without the permutation step.
 *
 * @param A Square matrix to balance.
 * @return The scaling vector d (the diagonal of D), all entries powers of two.
 */
template<class Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>
balance_scaling(const Eigen::MatrixBase<Derived>& A){
	using T = typename Derived::Scalar;
	constexpr int N = Derived::RowsAtCompileTime;
 
	// d accumulates the scaling: after the loop, A'(i,j) = A(i,j) * d(j) / d(i).
	// Start with D = I (no scaling).
	Eigen::Matrix<T, N, 1> d = Eigen::Matrix<T, N, 1>::Ones();
 
	// M always holds the absolute values of the CURRENT balanced matrix |D^-1 A D|.
	// Using absolute values is enough: d is positive, so scaling never changes a sign,
	// and only the magnitudes matter for the balancing criterion.
	Eigen::Matrix<T, N, N> M = A.cwiseAbs();
 
	// Repeat until a full sweep over all indices changes nothing.
	// Scaling index i changes the row/column sums of the other indices too,
	// so several sweeps may be needed. 20 is only a safety cap; in practice it
	// converges within a few sweeps.
	for (int sweep = 0; sweep < 20; ++sweep) {
		bool changed = false;
 
		for (int i = 0; i < N; ++i) {
			// "Size" of column i and row i, excluding the diagonal entry
			// (the diagonal is not affected by the scaling, see above).
			// Here the size is the sum of absolute values (1-norm).
			const T c = M.col(i).sum() - M(i, i);   // how big column i is
			const T r = M.row(i).sum() - M(i, i);   // how big row i is
 
			// A row or column with only zeros cannot be balanced
			// (and r / c would divide by zero), so skip it.
			if (!(c > T(0)) || !(r > T(0))) { continue; }
 
			// Scaling index i by a factor f means d(i) *= f. This has two effects:
			//   column i is multiplied by f   (A'(j,i) = A(j,i) * d(i) / d(j))
			//   row    i is divided    by f   (A'(i,j) = A(i,j) * d(j) / d(i))
			// so afterwards:   column size = c * f,   row size = r / f.
			//
			// Balanced means c * f == r / f   <=>   f^2 = r / c   <=>   f = sqrt(r / c).
			// We want f to be a power of two, f = 2^e, so
			//   e = round( log2(sqrt(r / c)) ) = round( 0.5 * log2(r / c) ).
			const int e = static_cast<int>(std::lround(T(0.5) * std::log2(r / c)));
 
			// e == 0 means f == 1: this index is already balanced as well as
			// a power of two allows.
			if (e == 0) { continue; }
 
			// f = 1 * 2^e, computed exactly (ldexp only changes the exponent).
			const T f = std::ldexp(T(1), e);
 
			// Only apply the scaling if it clearly helps: the new size
			// c * f + r / f must be at least 5% smaller than the old size c + r.
			// Without this margin the algorithm could keep flipping indices back
			// and forth by tiny amounts and never terminate.
			if (c * f + r / f < T(0.95) * (c + r)) {
				d(i) *= f;          // remember the scaling (this is the result)
 
				// Update the working copy of |D^-1 A D| to match.
				// The diagonal entry M(i,i) is hit by both lines: it is multiplied
				// by f and then divided by f again, so it correctly stays unchanged.
				M.col(i) *= f;
				M.row(i) /= f;
 
				changed = true;
			}
		}
 
		// A whole sweep without any change: the matrix is balanced.
		if (!changed) { break; }
	}
 
	return d;
}

}