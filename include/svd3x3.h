
//              Copyright Manolis Lourakis 2024.
// Distributed under the Boost Software License, Version 1.0.
//    (See copy at https://www.boost.org/LICENSE_1_0.txt)

/* Fast SVD for 3x3 matrices based on the polar and eigen decompositions:
 * if A=Q*H and H=V*S*V', then an SVD is (Q*V)*S*V' \equiv U*S*V'
 *
 * See https://nhigham.com/2020/07/28/what-is-the-polar-decomposition/comment-page-1/
 * https://github.com/martinbis11/polar-decomposition-3x3/tree/master
 * and https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
 */

#ifndef __SVD_3X3_H__
#define __SVD_3X3_H__

#include <iostream>
#include <cmath>

#include "polar_decomposition_3x3.h"
#include "SymmetricEigensolver3x3.h"

namespace svd3x3
{

template <typename TReal>
class decomp
{
public:

// Compute the SVD of matrix A.
//
// This operator decomposes the input matrix into two orthogonal matrices U & V.
//
// Matrices are represented in row-major order.
//
// A [in]  : Matrix to decompose.
// U [out] : Orthogonal matrix of left singular vectors.
// s [out] : singular values vector.
// V [out] : Orthogonal matrix of right singular vectors.
void operator()(const TReal A[9], TReal U[9], TReal s[3], TReal V[9]) const
{
  TReal At[9], Q[9], H[9];
  TReal scfac = 1.0;
  const bool scale = true;

  // polar decomposition uses column-major order
  At[0] = A[0]; At[1] = A[3]; At[2] = A[6];
  At[3] = A[1]; At[4] = A[4]; At[5] = A[7];
  At[6] = A[2]; At[7] = A[5]; At[8] = A[8];
  if (scale) {
    // scale with the maximum absolute element of A
    scfac = std::fabs(At[0]);
    for (int i = 1; i < 9; ++i)
      if (std::fabs(At[i]) > scfac) scfac = std::fabs(At[i]);
    scfac = (scfac > 0) ? 1.0 / scfac : 1.0;

    At[0] *= scfac; At[1] *= scfac; At[2] *= scfac;
    At[3] *= scfac; At[4] *= scfac; At[5] *= scfac;
    At[6] *= scfac; At[7] *= scfac; At[8] *= scfac;
  }
  polar::polar_decomposition(Q, H, At); // A'=Q'*H'
  // back to row-major
  mattransp3x3(Q);
  // H is symmetric, no need to transpose

  std::array<TReal, 3> eval;
  std::array<std::array<TReal, 3>, 3> evec;
  gte::SymmetricEigensolver3x3<TReal> eig;

  eig(H[0], H[1], H[2], H[4], H[5], H[8], false, -1, eval, evec); // decreasing order

  bool flip = false;
  if (eval[0] * eval[1] * eval[2] >= 0.0) { // det(H)
    s[0] = eval[0];
    s[1] = eval[1];
    s[2] = eval[2];

    V[0] = evec[0][0]; V[1] = evec[1][0]; V[2] = evec[2][0];
    V[3] = evec[0][1]; V[4] = evec[1][1]; V[5] = evec[2][1];
    V[6] = evec[0][2]; V[7] = evec[1][2]; V[8] = evec[2][2];
  } else {
    // H is positive semidefinite, thus should be negated to make det(H) non-negative;
    // instead, the decomposition is adjusted by negating the evals and reordering evecs
    s[0] = -eval[2];
    s[1] = -eval[1];
    s[2] = -eval[0];

    V[0] = evec[2][0]; V[1] = evec[1][0]; V[2] = evec[0][0];
    V[3] = evec[2][1]; V[4] = evec[1][1]; V[5] = evec[0][1];
    V[6] = evec[2][2]; V[7] = evec[1][2]; V[8] = evec[0][2];
    flip = true;
  }

  matmul3x3(Q, V, U);

  if (flip) { // negate V
    V[0] = -V[0]; V[1] = -V[1]; V[2] = -V[2];
    V[3] = -V[3]; V[4] = -V[4]; V[5] = -V[5];
    V[6] = -V[6]; V[7] = -V[7]; V[8] = -V[8];
  }

  if (scale) {
    s[0] /= scfac;
    s[1] /= scfac;
    s[2] /= scfac;
  }
}

// compute U*s*V' for verification
inline void compose(const TReal U[9], const TReal s[3], const TReal V[9], TReal UsVt[9]) const
{
  TReal sVt[9];
  sVt[0] = V[0]*s[0]; sVt[1] = V[3]*s[0]; sVt[2] = V[6]*s[0];
  sVt[3] = V[1]*s[1]; sVt[4] = V[4]*s[1]; sVt[5] = V[7]*s[1];
  sVt[6] = V[2]*s[2]; sVt[7] = V[5]*s[2]; sVt[8] = V[8]*s[2];

  matmul3x3(U, sVt, UsVt);
}

private:

// A*B
static inline void matmul3x3(const TReal *A, const TReal *B, TReal *prod)
{
  prod[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
  prod[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
  prod[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];

  prod[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
  prod[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
  prod[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];

  prod[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
  prod[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
  prod[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
}

// transpose in place
static inline void mattransp3x3(TReal *A)
{
  TReal tmp;

  tmp = A[1]; A[1] = A[3]; A[3] = tmp;
  tmp = A[5]; A[5] = A[7]; A[7] = tmp;
  tmp = A[6]; A[6] = A[2]; A[2] = tmp;
}

};

} // namespace

#endif // __SVD_3X3_H__
