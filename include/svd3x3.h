
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
    for (int i = 1; i < 9; ++i) {
      TReal x = std::fabs(At[i]);
      scfac = (x > scfac) ? x : scfac;
    }
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
#if 1
  gte::SymmetricEigensolver3x3<TReal> eig;
  eig(H[0], H[1], H[2], H[4], H[5], H[8], false, -1, eval, evec); // decreasing order

#else
  // alternative eigensolver not relying on Geometric Tools
  diagonalizer3x3_desc(H, &evec[0][0], &eval[0]);
#endif

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
inline void compose(const TReal U[9], const TReal s[3], const TReal V[9], TReal UsVt[9]) const noexcept
{
  TReal sVt[9];
  sVt[0] = V[0]*s[0]; sVt[1] = V[3]*s[0]; sVt[2] = V[6]*s[0];
  sVt[3] = V[1]*s[1]; sVt[4] = V[4]*s[1]; sVt[5] = V[7]*s[1];
  sVt[6] = V[2]*s[2]; sVt[7] = V[5]*s[2]; sVt[8] = V[8]*s[2];

  matmul3x3(U, sVt, UsVt);
}

private:

// A*B
static inline void matmul3x3(const TReal *A, const TReal *B, TReal *prod) noexcept
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
static inline void mattransp3x3(TReal *A) noexcept
{
  TReal tmp;

  tmp = A[1]; A[1] = A[3]; A[3] = tmp;
  tmp = A[5]; A[5] = A[7]; A[7] = tmp;
  tmp = A[6]; A[6] = A[2]; A[2] = tmp;
}


/**
 * Diagonalization of a symmetric 3x3 matrix.
 * Based on code by S Melax modified by ML to avoid constructing a quaternion
 * See http://www.melax.com/diag
 *
 * Original logic preserved:
 * - Iterates to zero out largest off-diagonal elements.
 * - Returns orthogonal matrix Q such that its columns are A's eigenvectors: A = Q * D * Q^T.
 * - Eigenvalues are on the diagonal of matrix D = Q^T * A * Q and also returned in ev.
 */
static void diagonalizer3x3(const TReal *A, TReal *Q, TReal *ev)
{
    const int maxsteps = 24;

    // Q represents the accumulated rotations (eigenvectors in columns)
    TReal tmp[9], D[9];
    // start with no rotation: Q = I3
    Q[0] = Q[4] = Q[8] = 1.0;
    Q[1] = Q[2] = Q[3] = Q[5] = Q[6] = Q[7] = 0.0;

    // look-up tables to avoid modulo 3
    static const int k1s[3] = {1, 2, 0}, k2s[3] = {2, 0, 1};

    for(int i = 0; i < maxsteps; i++)
    {
        // calculate current matrix D = Q^T * A * Q
        // as Q holds eigenvectors in columns, we use the similarity transform Q^T * A * Q

        // tmp=Q^T*A
        tmp[0]=Q[0]*A[0] + Q[3]*A[3] + Q[6]*A[6];
        tmp[1]=Q[0]*A[1] + Q[3]*A[4] + Q[6]*A[7];
        tmp[2]=Q[0]*A[2] + Q[3]*A[5] + Q[6]*A[8];

        tmp[3]=Q[1]*A[0] + Q[4]*A[3] + Q[7]*A[6];
        tmp[4]=Q[1]*A[1] + Q[4]*A[4] + Q[7]*A[7];
        tmp[5]=Q[1]*A[2] + Q[4]*A[5] + Q[7]*A[8];

        tmp[6]=Q[2]*A[0] + Q[5]*A[3] + Q[8]*A[6];
        tmp[7]=Q[2]*A[1] + Q[5]*A[4] + Q[8]*A[7];
        tmp[8]=Q[2]*A[2] + Q[5]*A[5] + Q[8]*A[8];

        // D = tmp * Q  (where tmp is Q^T * A); upper triangle only!
        D[0]=tmp[0]*Q[0] + tmp[1]*Q[3] + tmp[2]*Q[6];
        D[1]=tmp[0]*Q[1] + tmp[1]*Q[4] + tmp[2]*Q[7];
        D[2]=tmp[0]*Q[2] + tmp[1]*Q[5] + tmp[2]*Q[8];
        //D[3]=tmp[3]*Q[0] + tmp[4]*Q[3] + tmp[5]*Q[6];
        D[4]=tmp[3]*Q[1] + tmp[4]*Q[4] + tmp[5]*Q[7];
        D[5]=tmp[3]*Q[2] + tmp[4]*Q[5] + tmp[5]*Q[8];
        //D[6]=tmp[6]*Q[0] + tmp[7]*Q[3] + tmp[8]*Q[6];
        //D[7]=tmp[6]*Q[1] + tmp[7]*Q[4] + tmp[8]*Q[7];
        D[8]=tmp[6]*Q[2] + tmp[7]*Q[5] + tmp[8]*Q[8];

        // off-diagonal elements and their magnitudes
        const TReal offdiag[3] = { D[5], D[2], D[1] };
        const TReal om[3] = { std::fabs(offdiag[0]), std::fabs(offdiag[1]), std::fabs(offdiag[2]) };

        // index of largest absolute off-diagonal element
        const int k = (om[0] > om[1] && om[0] > om[2])? 0 : (om[1] > om[2])? 1 : 2;

        if(om[k] < 1.E-10) break; // already diagonal

        const int k1 = k1s[k];  // (k + 1) % 3;
        const int k2 = k2s[k];  // (k + 2) % 3;

        // calculate rotation angle parameters (t, c, s)
        TReal theta = 0.5 * (D[k2*4]-D[k1*4])/offdiag[k]; // D[k*4] == D[k][k]

        const TReal sgn = (theta > 0.0) ? 1.0 : -1.0;
        theta *= sgn; // make positive

        // t = sgn / ( |theta| + sqrt(theta^2 + 1) )
        const TReal t = sgn / (theta + ((theta < 1.E6) ? std::sqrt(theta * theta + 1.0) : theta));
        const TReal c = 1.0 / std::sqrt(t * t + 1.0); //  c=1/(t^2+1) , t=s/c

        if(1.0 - c < 1.E-10) break; // reached machine precision

        const TReal s = t * c;
        // accumulate the rotation: Q_new = Q_old * J
        // for a Jacobi rotation matrix J

        // only cols k1, k2 of Q are updated; J is not constructed
        // the following loop is unrolled below:
        /***
        for (int i = 0; i < 3; ++i) {
            const TReal v1 = Q[i*3 + k1];
            const TReal v2 = Q[i*3 + k2];
            Q[i*3 + k1] =  c * v1 - s * v2;
            Q[i*3 + k2] =  s * v1 + c * v2;
        }
        ***/

        TReal& q0k1 = Q[0*3 + k1];
        TReal& q0k2 = Q[0*3 + k2];
        const TReal v1_0 = q0k1;
        const TReal v2_0 = q0k2;
        q0k1 = c * v1_0 - s * v2_0;
        q0k2 = s * v1_0 + c * v2_0;

        TReal& q1k1 = Q[1*3 + k1];
        TReal& q1k2 = Q[1*3 + k2];
        const TReal v1_1 = q1k1;
        const TReal v2_1 = q1k2;
        q1k1 = c * v1_1 - s * v2_1;
        q1k2 = s * v1_1 + c * v2_1;

        TReal& q2k1 = Q[2*3 + k1];
        TReal& q2k2 = Q[2*3 + k2];
        const TReal v1_2 = q2k1;
        const TReal v2_2 = q2k2;
        q2k1 = c * v1_2 - s * v2_2;
        q2k2 = s * v1_2 + c * v2_2;
    }

    ev[0] = D[0]; ev[1] = D[4]; ev[2] = D[8];
}

// as above with the eigenvalues sorted in descending order
// and the eigenvectors returned as rows of Qr
static inline void diagonalizer3x3_desc(const TReal *A, TReal *Qr, TReal *evs) noexcept
{
    TReal Q[9], ev[3];
    diagonalizer3x3(A, Q, ev);

    // sort in descending order
    int i0 = 0, i1 = 1, i2 = 2;
    // ensure ev[i0] >= ev[i1]
    if (ev[i0] < ev[i1]) std::swap(i0, i1);
    // ensure ev[i1] >= ev[i2]
    if (ev[i1] < ev[i2]) std::swap(i1, i2);
    // ensure ev[i0] >= ev[i1]
    if (ev[i0] < ev[i1]) std::swap(i0, i1);

    evs[0] = ev[i0]; evs[1] = ev[i1]; evs[2] = ev[i2];

    Qr[0] = Q[0 + i0];   Qr[3] = Q[0 + i1];   Qr[6] = Q[0 + i2];
    Qr[1] = Q[3 + i0];   Qr[4] = Q[3 + i1];   Qr[7] = Q[3 + i2];
    Qr[2] = Q[6 + i0];   Qr[5] = Q[6 + i1];   Qr[8] = Q[6 + i2];
}

};

} // namespace

#endif // __SVD_3X3_H__
